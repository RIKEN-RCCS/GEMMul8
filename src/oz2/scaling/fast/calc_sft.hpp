#pragma once
#include "../../common/common.hpp"
#include "../../common/table.hpp"
#include "find_max.hpp"
#include "config.hpp"

namespace gemmul8::scaling::fast {

inline constexpr float u2 = 0x1.0000000000000p-23F; // 2^-23

template <Backend BACKEND, unsigned NUM_MODULI, typename T>
__device__ __forceinline__ int32_t calc_sft(T vecnrm) {
    if (vecnrm == T(0.0)) return 0;
    const int32_t exponent = common::Tilogb<T>(vecnrm);
    const float vecnrmf    = common::float_ru<T>(common::Tscalbn<T>(vecnrm, -exponent));
    const float exponent2  = exponent * 0.5f;
    constexpr float log2P  = common::table::log2P<BACKEND, NUM_MODULI>;
    const float mu0        = __fsub_rd(log2P, exponent2);
    const float mu1        = __fadd_ru(__log2f(vecnrmf) * 0.5f, u2);
    return __float2int_rd(__fsub_rd(mu0, mu1));
}

template <typename T, Backend BACKEND, unsigned NUM_MODULI,
          cublasFillMode_t UPLO, cublasDiagType_t DIAG>
__global__ void calc_sft_rowwise(
    const unsigned rows_A, const unsigned cols_A,
    const T *const __restrict__ A, const size_t lda,
    int16_t *const __restrict__ sftA //
) {
    using U = common::underlying_t<T>;
    __shared__ U ssum[common::TILE_DIM][common::TILE_DIM + 1];

    U sum = calc_nrm_tile<T, UPLO, DIAG>(rows_A, cols_A, A, lda, ssum);
    if constexpr (DIAG == CUBLAS_DIAG_UNIT) {
        constexpr U Uone = common::Tconst<U>::one();
        sum              = common::__Tadd_ru<U>(sum, Uone);
    }

    const unsigned row_idx = blockIdx.x * common::TILE_DIM + threadIdx.y;
    if (row_idx < rows_A && threadIdx.x == 0) {
        const int32_t sft = calc_sft<BACKEND, NUM_MODULI, U>(sum);
        sftA[row_idx]     = int16_t(-sft);
    }
}

template <typename T, Backend BACKEND, unsigned NUM_MODULI,
          cublasFillMode_t UPLO, cublasDiagType_t DIAG>
__device__ __forceinline__ int32_t calc_sft_colwise(
    const unsigned rows_A,
    const T *const __restrict__ in,
    int16_t *const __restrict__ sftA,
    common::underlying_t<T> *ssum //
) {
    using U = common::underlying_t<T>;

    U sum = calc_nrm<T, UPLO, DIAG>(rows_A, in, ssum);
    if constexpr (DIAG == CUBLAS_DIAG_UNIT) {
        constexpr U Uone = common::Tconst<U>::one();
        sum              = common::__Tadd_ru<U>(sum, Uone);
    }

    const int32_t sft = calc_sft<BACKEND, NUM_MODULI, U>(sum);
    if (threadIdx.x == 0) {
        const unsigned col_idx = blockIdx.x;
        sftA[col_idx]          = int16_t(-sft);
    }

    return sft;
}

template <typename T, cublasFillMode_t UPLO, bool HERM>
__global__ void calc_stat_sym_rowwise(
    const unsigned n,
    const T *const __restrict__ A, const size_t lda,
    common::underlying_t<T> *const __restrict__ sumA //
) {
    using U = common::underlying_t<T>;

    __shared__ U ssum[common::TILE_DIM][common::TILE_DIM + 1];

    const U sum = calc_nrm_tile<T, UPLO, CUBLAS_DIAG_NON_UNIT, HERM>(
        n, n, A, lda, ssum);

    const unsigned row_idx = blockIdx.x * common::TILE_DIM + threadIdx.y;
    if (row_idx < n && threadIdx.x == 0) {
        sumA[row_idx] = sum;
    }
}

template <typename T, Backend BACKEND, unsigned NUM_MODULI,
          cublasFillMode_t UPLO>
__global__ void calc_sft_sym_colwise(
    const unsigned n,
    const T *const __restrict__ A, const size_t lda,
    int16_t *const __restrict__ sftA,
    common::underlying_t<T> *const __restrict__ sumA //
) {
    using U = common::underlying_t<T>;

    __shared__ U ssum[32];

    const unsigned col_idx = blockIdx.x;
    if (col_idx >= n) return;

    const T *const __restrict__ in = A + col_idx * lda;

    const U sum_col = calc_nrm<T, UPLO, CUBLAS_DIAG_UNIT>(n, in, ssum);

    if (threadIdx.x == 0) {
        const U sum_all = common::__Tadd_ru<U>(sumA[col_idx], sum_col);

        const int32_t sft = calc_sft<BACKEND, NUM_MODULI, U>(sum_all);
        sftA[col_idx]     = int16_t(-sft);
    }
}

template <typename T, cublasFillMode_t UPLO, bool HERM>
__global__ void calc_stat_sym_rowwise_partial(
    const unsigned n,
    const T *const __restrict__ A, const size_t lda,
    common::underlying_t<T> *const __restrict__ partial_sum //
) {
    using U = common::underlying_t<T>;

    __shared__ U ssum[common::TILE_DIM][common::TILE_DIM + 1];

    const unsigned col_begin = blockIdx.y * rowwise_sft_col_tile;
    const unsigned col_end   = min(n, col_begin + rowwise_sft_col_tile);

    const U sum = calc_nrm_tile_range<T, UPLO, CUBLAS_DIAG_NON_UNIT, HERM>(
        n, n, A, lda, ssum, col_begin, col_end);

    const unsigned row_idx = blockIdx.x * common::TILE_DIM + threadIdx.y;

    if (row_idx < n && threadIdx.x == 0) {
        const size_t idx = size_t(row_idx) * gridDim.y + blockIdx.y;
        partial_sum[idx] = sum;
    }
}

template <typename T>
__global__ void calc_stat_sym_rowwise_reduce(
    const unsigned n,
    const unsigned num_col_blocks,
    const common::underlying_t<T> *const __restrict__ partial_sum,
    common::underlying_t<T> *const __restrict__ sumA //
) {
    using U = common::underlying_t<T>;

    const unsigned row_idx = blockIdx.x * blockDim.y + threadIdx.y;

    U sum = common::Tconst<U>::zero();

    if (row_idx < n) {
        for (unsigned cb = threadIdx.x; cb < num_col_blocks; cb += blockDim.x) {
            const size_t idx = size_t(row_idx) * num_col_blocks + cb;
            sum              = common::__Tadd_ru<U>(sum, partial_sum[idx]);
        }

        sum = common::inner_warp_sum<U>(sum);

        if (threadIdx.x == 0) {
            sumA[row_idx] = sum;
        }
    }
}

template <typename T, cublasFillMode_t UPLO, bool HERM>
inline void calc_stat_sym_rowwise_launch(
    const cudaStream_t stream,
    const unsigned n,
    const T *const A, const size_t lda,
    common::underlying_t<T> *const partial_sum,
    common::underlying_t<T> *const sumA //
) {
    constexpr dim3 threads_findmax(threads_x_findmax_tile,
                                   threads_y_findmax_tile);

    if (n < rowwise_sft_split_threshold) {
        const unsigned grid_findmax =
            (n + common::TILE_DIM - 1U) / common::TILE_DIM;

        calc_stat_sym_rowwise<T, UPLO, HERM>
            <<<grid_findmax, threads_findmax, 0, stream>>>(
                n, A, lda, sumA);
        return;
    }

    const unsigned num_col_blocks =
        (n + rowwise_sft_col_tile - 1U) / rowwise_sft_col_tile;

    const dim3 grid_partial(
        (n + common::TILE_DIM - 1U) / common::TILE_DIM,
        num_col_blocks);

    calc_stat_sym_rowwise_partial<T, UPLO, HERM>
        <<<grid_partial, threads_findmax, 0, stream>>>(
            n, A, lda, partial_sum);

    constexpr dim3 threads_reduce(threads_x_rowwise_sft_reduce,
                                  threads_y_rowwise_sft_reduce);

    const unsigned grid_reduce =
        (n + threads_y_rowwise_sft_reduce - 1U) / threads_y_rowwise_sft_reduce;

    calc_stat_sym_rowwise_reduce<T>
        <<<grid_reduce, threads_reduce, 0, stream>>>(
            n, num_col_blocks, partial_sum, sumA);
}

template <typename T, cublasFillMode_t UPLO, cublasDiagType_t DIAG>
__global__ void calc_sft_rowwise_partial(
    const unsigned rows_A, const unsigned cols_A,
    const T *const __restrict__ A, const size_t lda,
    common::underlying_t<T> *const __restrict__ partial_sum //
) {
    using U = common::underlying_t<T>;

    __shared__ U ssum[common::TILE_DIM][common::TILE_DIM + 1];

    const unsigned col_begin = blockIdx.y * rowwise_sft_col_tile;
    const unsigned col_end   = min(cols_A, col_begin + rowwise_sft_col_tile);

    const U sum = calc_nrm_tile_range<T, UPLO, DIAG>(
        rows_A, cols_A, A, lda, ssum, col_begin, col_end);

    const unsigned row_idx = blockIdx.x * common::TILE_DIM + threadIdx.y;
    if (row_idx < rows_A && threadIdx.x == 0) {
        const size_t out = size_t(row_idx) * gridDim.y + blockIdx.y;
        partial_sum[out] = sum;
    }
}

template <typename T, Backend BACKEND, unsigned NUM_MODULI, cublasDiagType_t DIAG>
__global__ void calc_sft_rowwise_reduce(
    const unsigned rows_A,
    const unsigned num_col_blocks,
    const common::underlying_t<T> *const __restrict__ partial_sum,
    int16_t *const __restrict__ sftA //
) {
    using U = common::underlying_t<T>;

    const unsigned row_idx = blockIdx.x * blockDim.y + threadIdx.y;

    U sum = common::Tconst<U>::zero();

    if (row_idx < rows_A) {
        for (unsigned cb = threadIdx.x; cb < num_col_blocks; cb += blockDim.x) {
            const size_t idx = size_t(row_idx) * num_col_blocks + cb;
            sum              = common::__Tadd_ru<U>(sum, partial_sum[idx]);
        }

        sum = common::inner_warp_sum<U>(sum);

        if constexpr (DIAG == CUBLAS_DIAG_UNIT) {
            constexpr U Uone = common::Tconst<U>::one();
            sum              = common::__Tadd_ru<U>(sum, Uone);
        }

        if (threadIdx.x == 0) {
            const int32_t sft = calc_sft<BACKEND, NUM_MODULI, U>(sum);
            sftA[row_idx]     = int16_t(-sft);
        }
    }
}

template <typename T, Backend BACKEND, unsigned NUM_MODULI,
          cublasFillMode_t UPLO, cublasDiagType_t DIAG>
inline void calc_sft_rowwise_launch(
    const cudaStream_t stream,
    const unsigned rows_A, const unsigned cols_A,
    const T *const A, const size_t lda,
    int16_t *const sftA,
    common::underlying_t<T> *const partial_sum //
) {
    constexpr dim3 threads_findmax(threads_x_findmax_tile,
                                   threads_y_findmax_tile);

    if (cols_A < rowwise_sft_split_threshold) {
        const unsigned grid_findmax = (rows_A + common::TILE_DIM - 1U) / common::TILE_DIM;
        calc_sft_rowwise<T, BACKEND, NUM_MODULI, UPLO, DIAG>
            <<<grid_findmax, threads_findmax, 0, stream>>>(
                rows_A, cols_A, A, lda, sftA);
        return;
    }

    const unsigned num_col_blocks =
        (cols_A + rowwise_sft_col_tile - 1U) / rowwise_sft_col_tile;
    dim3 grid_partial((rows_A + common::TILE_DIM - 1U) / common::TILE_DIM,
                      num_col_blocks);

    calc_sft_rowwise_partial<T, UPLO, DIAG>
        <<<grid_partial, threads_findmax, 0, stream>>>(
            rows_A, cols_A, A, lda, partial_sum);

    constexpr dim3 threads_reduce(threads_x_rowwise_sft_reduce,
                                  threads_y_rowwise_sft_reduce);
    const unsigned grid_reduce =
        (rows_A + threads_y_rowwise_sft_reduce - 1U) / threads_y_rowwise_sft_reduce;

    calc_sft_rowwise_reduce<T, BACKEND, NUM_MODULI, DIAG>
        <<<grid_reduce, threads_reduce, 0, stream>>>(
            rows_A, num_col_blocks, partial_sum, sftA);
}

} // namespace gemmul8::scaling::fast
