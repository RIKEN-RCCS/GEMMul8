#pragma once
#include "../../common/common.hpp"
#include "../general/helper_triangular.hpp"
#include "../general/helper_hermitian.hpp"

namespace gemmul8::scaling::fast {

inline constexpr float UNIT_ROUNDOFF = 0x1.0000000000000p-24F;

template <typename T>
__device__ __forceinline__ T reduction_nrm(T sum, T *ssum) {
    sum = common::inner_warp_sum(sum);

    if ((threadIdx.x & 31) == 0) {
        ssum[threadIdx.x >> 5] = sum;
    }
    __syncthreads();

    sum = common::Tconst<T>::zero();
    if (threadIdx.x < 32) {
        if (threadIdx.x < (blockDim.x >> 5)) {
            sum = ssum[threadIdx.x];
        }
        sum = common::inner_warp_sum(sum);
        if (threadIdx.x == 0) {
            ssum[0] = sum;
        }
    }
    __syncthreads();

    return ssum[0];
}

template <typename T,
          cublasFillMode_t UPLO = CUBLAS_FILL_MODE_FULL,
          cublasDiagType_t DIAG = CUBLAS_DIAG_NON_UNIT>
__device__ __forceinline__ common::underlying_t<T> calc_nrm(
    const unsigned length,
    const T *const __restrict__ ptr,
    common::underlying_t<T> *ssum //
) {
    using U = common::underlying_t<T>;
    U sum   = common::Tconst<U>::zero();

    const auto [begin, end] = general::column_load_range<UPLO, DIAG>(length);
    for (unsigned i = begin; i < end; i += blockDim.x) {
        const T tmp = ptr[i];
        sum         = common::Tsqr_add_ru<T>(tmp, sum);
    }

    return reduction_nrm(sum, ssum);
}

template <typename T,
          cublasFillMode_t UPLO = CUBLAS_FILL_MODE_FULL,
          cublasDiagType_t DIAG = CUBLAS_DIAG_NON_UNIT,
          bool HERM             = false>
__device__ __forceinline__ common::underlying_t<T> calc_nrm_tile(
    const unsigned m,
    const unsigned n,
    const T *const __restrict__ A,
    const size_t lda,
    common::underlying_t<T> ssum[][common::TILE_DIM + 1] //
) {
    using U          = common::underlying_t<T>;
    unsigned row_idx = blockIdx.x * common::TILE_DIM + threadIdx.x;
    U sum            = common::Tconst<U>::zero();

    if (row_idx < m) {
        const T *const __restrict__ row_ptr = A + row_idx;

        if constexpr (UPLO == CUBLAS_FILL_MODE_UPPER) {

            const unsigned row_base   = row_idx - threadIdx.x;
            constexpr unsigned unit   = (DIAG == CUBLAS_DIAG_UNIT) ? 1U : 0U;
            const unsigned full_begin = row_base + common::TILE_DIM - 1U + unit;
            const unsigned diag_end   = min(n, full_begin);

            for (unsigned col = row_base + threadIdx.y; col < diag_end; col += blockDim.y) {
                if (row_idx + unit <= col) {
                    T val = row_ptr[col * lda];
                    val   = general::hemm_diag_value<HERM, T>(val, row_idx, col);
                    sum   = common::Tsqr_add_ru<T>(val, sum);
                }
            }
            for (unsigned col = full_begin + threadIdx.y; col < n; col += blockDim.y) {
                const T tmp = row_ptr[col * lda];
                sum         = common::Tsqr_add_ru<T>(tmp, sum);
            }

        } else if constexpr (UPLO == CUBLAS_FILL_MODE_LOWER) {

            const unsigned row_base      = row_idx - threadIdx.x;
            constexpr unsigned keep_diag = (DIAG == CUBLAS_DIAG_UNIT) ? 0U : 1U;
            constexpr unsigned unit      = (DIAG == CUBLAS_DIAG_UNIT) ? 1U : 0U;
            const unsigned full_end      = min(n, row_base + keep_diag);
            const unsigned diag_end      = min(n, row_base + common::TILE_DIM);

            for (unsigned col = threadIdx.y; col < full_end; col += blockDim.y) {
                const T tmp = row_ptr[col * lda];
                sum         = common::Tsqr_add_ru<T>(tmp, sum);
            }
            for (unsigned col = full_end + threadIdx.y; col < diag_end; col += blockDim.y) {
                if (col + unit <= row_idx) {
                    T val = row_ptr[col * lda];
                    val   = general::hemm_diag_value<HERM, T>(val, row_idx, col);
                    sum   = common::Tsqr_add_ru<T>(val, sum);
                }
            }

        } else {

            for (unsigned col = threadIdx.y; col < n; col += blockDim.y) {
                const T tmp = row_ptr[col * lda];
                sum         = common::Tsqr_add_ru<T>(tmp, sum);
            }
        }
    }

    ssum[threadIdx.y][threadIdx.x] = sum;
    __syncthreads();

    sum = ssum[threadIdx.x][threadIdx.y];

    return common::inner_warp_sum<U, common::TILE_DIM>(sum);
}

template <typename T,
          cublasFillMode_t UPLO = CUBLAS_FILL_MODE_FULL,
          cublasDiagType_t DIAG = CUBLAS_DIAG_NON_UNIT,
          bool HERM             = false>
__device__ __forceinline__ common::underlying_t<T> calc_nrm_tile_range(
    const unsigned m,
    const unsigned n,
    const T *const __restrict__ A,
    const size_t lda,
    common::underlying_t<T> ssum[][common::TILE_DIM + 1],
    const unsigned col_begin,
    const unsigned col_end_in //
) {
    using U          = common::underlying_t<T>;
    unsigned row_idx = blockIdx.x * common::TILE_DIM + threadIdx.x;
    U sum            = common::Tconst<U>::zero();

    const unsigned col_end = min(n, col_end_in);

    if (row_idx < m && col_begin < col_end) {
        const T *const __restrict__ row_ptr = A + row_idx;

        if constexpr (UPLO == CUBLAS_FILL_MODE_UPPER) {

            const unsigned row_base   = row_idx - threadIdx.x;
            constexpr unsigned unit   = (DIAG == CUBLAS_DIAG_UNIT) ? 1U : 0U;
            const unsigned full_begin = row_base + common::TILE_DIM - 1U + unit;

            const unsigned diag_begin = max(col_begin, row_base);
            const unsigned diag_end   = min(col_end, min(n, full_begin));

            for (unsigned col = diag_begin + threadIdx.y; col < diag_end; col += blockDim.y) {
                if (row_idx + unit <= col) {
                    T val       = row_ptr[col * lda];
                    val         = general::hemm_diag_value<HERM, T>(val, row_idx, col);
                    sum         = common::Tsqr_add_ru<T>(val, sum);
                }
            }

            const unsigned dense_begin = max(col_begin, full_begin);
            for (unsigned col = dense_begin + threadIdx.y; col < col_end; col += blockDim.y) {
                const T tmp = row_ptr[col * lda];
                sum         = common::Tsqr_add_ru<T>(tmp, sum);
            }

        } else if constexpr (UPLO == CUBLAS_FILL_MODE_LOWER) {

            const unsigned row_base      = row_idx - threadIdx.x;
            constexpr unsigned keep_diag = (DIAG == CUBLAS_DIAG_UNIT) ? 0U : 1U;
            constexpr unsigned unit      = (DIAG == CUBLAS_DIAG_UNIT) ? 1U : 0U;
            const unsigned full_end      = min(n, row_base + keep_diag);
            const unsigned diag_end0     = min(n, row_base + common::TILE_DIM);

            const unsigned dense_begin = col_begin;
            const unsigned dense_end   = min(col_end, full_end);
            for (unsigned col = dense_begin + threadIdx.y; col < dense_end; col += blockDim.y) {
                const T tmp = row_ptr[col * lda];
                sum         = common::Tsqr_add_ru<T>(tmp, sum);
            }

            const unsigned diag_begin = max(col_begin, full_end);
            const unsigned diag_end   = min(col_end, diag_end0);
            for (unsigned col = diag_begin + threadIdx.y; col < diag_end; col += blockDim.y) {
                if (col + unit <= row_idx) {
                    T val       = row_ptr[col * lda];
                    val         = general::hemm_diag_value<HERM, T>(val, row_idx, col);
                    sum         = common::Tsqr_add_ru<T>(val, sum);
                }
            }

        } else {

            for (unsigned col = col_begin + threadIdx.y; col < col_end; col += blockDim.y) {
                const T tmp = row_ptr[col * lda];
                sum         = common::Tsqr_add_ru<T>(tmp, sum);
            }
        }
    }

    ssum[threadIdx.y][threadIdx.x] = sum;
    __syncthreads();

    sum    = ssum[threadIdx.x][threadIdx.y];
    return common::inner_warp_sum<U, common::TILE_DIM>(sum);
}

} // namespace gemmul8::scaling::fast
