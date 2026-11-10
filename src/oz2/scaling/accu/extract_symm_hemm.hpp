#pragma once
#include "../../common/launch_grid.hpp"
#include "config.hpp"
#include "calc_sft.hpp"
#include "store.hpp"

#include "../general/helper_temp.hpp"
#include "../general/helper_triangular.hpp"
#include "../general/helper_hermitian.hpp"
#include "../general/store.hpp"

namespace gemmul8::scaling::accu {

template <bool UPPER, bool HERM, bool STORE_TRANSPOSE, typename T, Backend BACKEND, bool EXTENDED = false>
__global__ void extract_symm_hemm_offdiag_kernel(
    const unsigned n,
    const T *const __restrict__ A, const size_t lda,
    common::matptr_t<common::low_t<BACKEND>, common::isComplex<T>> A_lo,
    const size_t lda_lo,
    const int16_t *const __restrict__ sftA //
) {
    __shared__ T tile[common::TILE_DIM][common::TILE_DIM + 1];
    __shared__ int16_t sft_row[common::TILE_DIM];
    __shared__ int16_t sft_col[common::TILE_DIM];

    const unsigned tile_i = blockIdx.x;
    const unsigned tile_j = common::block_y<EXTENDED>();

    if constexpr (UPPER) {
        if (tile_i >= tile_j) return;
    } else {
        if (tile_i <= tile_j) return;
    }

    const unsigned rowBase = tile_i * common::TILE_DIM;
    const unsigned colBase = tile_j * common::TILE_DIM;

    if (threadIdx.y == 0) {
        const unsigned r = rowBase + threadIdx.x;
        const unsigned c = colBase + threadIdx.x;

        sft_row[threadIdx.x] = (r < n) ? sftA[r] : 0;
        sft_col[threadIdx.x] = (c < n) ? sftA[c] : 0;
    }

    const unsigned row = rowBase + threadIdx.x;

#pragma unroll
    for (unsigned j = 0; j < common::TILE_DIM; j += threads_y_extract_symm_hemm<T>) {
        const unsigned yy = threadIdx.y + j;
        if (yy >= common::TILE_DIM) continue;

        const unsigned col    = colBase + yy;
        const T a             = (row < n && col < n) ? A[col * lda + row] : common::Tconst<T>::zero();
        tile[yy][threadIdx.x] = a;
    }
    __syncthreads();

#pragma unroll
    for (unsigned j = 0; j < common::TILE_DIM; j += threads_y_extract_symm_hemm<T>) {
        const unsigned yy = threadIdx.y + j;
        if (yy >= common::TILE_DIM) continue;

        const unsigned col = colBase + yy;
        if (row < n && col < n) {
            const size_t idx = col * lda_lo + row;
            const T a        = general::hemm_active_store_value<HERM, STORE_TRANSPOSE, T>(tile[yy][threadIdx.x], row, col);
            extract_store_one<T, BACKEND>(A_lo, idx, a, int32_t(sft_col[yy]));
        }

        const unsigned mrow = colBase + threadIdx.x;
        const unsigned mcol = rowBase + yy;
        if (mrow < n && mcol < n) {
            const size_t midx = mcol * lda_lo + mrow;
            const T b         = general::hemm_mirror_store_value<HERM, STORE_TRANSPOSE, T>(tile[threadIdx.x][yy]);
            extract_store_one<T, BACKEND>(A_lo, midx, b, int32_t(sft_row[yy]));
        }
    }
}

template <bool UPPER, bool HERM, bool STORE_TRANSPOSE, typename T, Backend BACKEND>
__global__ void extract_symm_hemm_diag_kernel(
    const unsigned n,
    const T *const __restrict__ A, const size_t lda,
    common::matptr_t<common::low_t<BACKEND>, common::isComplex<T>> A_lo,
    const size_t lda_lo,
    const int16_t *const __restrict__ sftA //
) {
    __shared__ T tile[common::TILE_DIM][common::TILE_DIM + 1];
    __shared__ int16_t sft_tile[common::TILE_DIM];

    const unsigned tile_id = blockIdx.x;
    const unsigned base    = tile_id * common::TILE_DIM;

    if (threadIdx.y == 0) {
        const unsigned idx    = base + threadIdx.x;
        sft_tile[threadIdx.x] = (idx < n) ? sftA[idx] : 0;
    }

#pragma unroll
    for (unsigned j = 0; j < common::TILE_DIM; j += threads_y_extract_symm_hemm<T>) {
        const unsigned yy = threadIdx.y + j;
        if (yy >= common::TILE_DIM) continue;

        const unsigned row  = base + threadIdx.x;
        const unsigned col  = base + yy;
        const bool in_range = (row < n && col < n);
        const bool active   = in_range && general::tri_elem_active<UPPER>(row, col);

        if (active) {
            T a                   = A[col * lda + row];
            a                     = general::hemm_diag_value<HERM, T>(a, row, col);
            tile[yy][threadIdx.x] = a;
        }
    }
    __syncthreads();

#pragma unroll
    for (unsigned j = 0; j < common::TILE_DIM; j += threads_y_extract_symm_hemm<T>) {
        const unsigned yy = threadIdx.y + j;
        if (yy >= common::TILE_DIM) continue;

        const unsigned row = base + threadIdx.x;
        const unsigned col = base + yy;
        if (row >= n || col >= n) continue;

        const bool active = general::tri_elem_active<UPPER>(row, col);
        const size_t idx  = col * lda_lo + row;
        if (active) {
            const T a = general::hemm_active_store_value<HERM, STORE_TRANSPOSE, T>(tile[yy][threadIdx.x], row, col);
            extract_store_one<T, BACKEND>(A_lo, idx, a, int32_t(sft_tile[yy]));
        } else {
            const T b = general::hemm_mirror_store_value<HERM, STORE_TRANSPOSE, T>(tile[threadIdx.x][yy]);
            extract_store_one<T, BACKEND>(A_lo, idx, b, int32_t(sft_tile[yy]));
        }
    }
}

template <typename T, Backend BACKEND,
          cublasFillMode_t UPLO, bool HERM,
          bool STORE_TRANSPOSE = false>
inline void extract_symm_hemm_launch(
    const cudaStream_t stream,
    const unsigned n,
    const T *const A, const size_t lda,
    common::matptr_t<common::low_t<BACKEND>, common::isComplex<T>> &A_lo,
    const size_t lda_lo, const size_t incA_lo,
    int16_t *const sftA //
) {
    static_assert(UPLO == CUBLAS_FILL_MODE_UPPER || UPLO == CUBLAS_FILL_MODE_LOWER,
                  "extract_symm_hemm_launch requires UPLO = UPPER or LOWER.");

    using U                       = common::underlying_t<T>;
    const unsigned num_col_blocks = (n < rowwise_sft_split_threshold) ? 1U : (n + rowwise_sft_col_tile - 1U) / rowwise_sft_col_tile;
    U *const partial_amax         = general::temporary_memory<common::low_t<BACKEND>, U>(A_lo.ptr0);
    U *const amaxA                = partial_amax + size_t(n) * num_col_blocks;

    constexpr dim3 threads_extract(threads_x_extract_rowwise, threads_y_extract_symm_hemm<T>);
    const unsigned nt      = (n + common::TILE_DIM - 1) / common::TILE_DIM;
    constexpr bool isUPPER = UPLO == CUBLAS_FILL_MODE_UPPER;

    calc_sft_before_sym_rowwise_launch<T, UPLO, HERM>(
        stream, n, A, lda, partial_amax, amaxA);

    calc_sft_before_sym_colwise<T, BACKEND, UPLO>
        <<<n, threads_accu, 0, stream>>>(
            n, A, lda, sftA, amaxA);

    common::launch_grid_y(dim3(nt, nt), [&]<bool EXTENDED>(dim3 part) {
        extract_symm_hemm_offdiag_kernel<isUPPER, HERM, STORE_TRANSPOSE, T, BACKEND, EXTENDED>
            <<<part, threads_extract, 0, stream>>>(
                n, A, lda, A_lo, lda_lo, sftA);
    });

    extract_symm_hemm_diag_kernel<isUPPER, HERM, STORE_TRANSPOSE, T, BACKEND>
        <<<nt, threads_extract, 0, stream>>>(
            n, A, lda, A_lo, lda_lo, sftA);

    general::memset_padding_low_mats_2d_async<T, BACKEND, 0U>(
        stream, A_lo, n, lda_lo, incA_lo / lda_lo);
}

template <typename T, Backend BACKEND, cublasFillMode_t UPLO>
void extract_symm(
    const cudaStream_t stream,
    const unsigned n,
    const T *const A, const size_t lda,
    common::matptr_t<common::low_t<BACKEND>, common::isComplex<T>> &A_lo,
    const size_t lda_lo, const size_t incA_lo,
    int16_t *const sftA //
) {
    extract_symm_hemm_launch<T, BACKEND, UPLO, false>(
        stream, n, A, lda, A_lo, lda_lo, incA_lo, sftA);
}

template <typename T, Backend BACKEND, cublasFillMode_t UPLO,
          bool STORE_TRANSPOSE>
void extract_hemm(
    const cudaStream_t stream,
    const unsigned n,
    const T *const A, const size_t lda,
    common::matptr_t<common::low_t<BACKEND>, common::isComplex<T>> &A_lo,
    const size_t lda_lo, const size_t incA_lo,
    int16_t *const sftA //
) {
    extract_symm_hemm_launch<T, BACKEND, UPLO, true, STORE_TRANSPOSE>(
        stream, n, A, lda, A_lo, lda_lo, incA_lo, sftA);
}

} // namespace gemmul8::scaling::accu
