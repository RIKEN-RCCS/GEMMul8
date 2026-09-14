#pragma once
#include "../../common/common.hpp"
#include "../../common/table.hpp"
#include "find_max.hpp"
#include "config.hpp"

namespace gemmul8::scaling::accu {

//==========
// Number of bits for extraction
//==========
#if FP8_FNUZ
template <Backend b> inline constexpr int maxUFP = (b == Backend::FP8) ? 6 : 5;
#else
template <Backend b> inline constexpr int maxUFP = (b == Backend::FP8) ? 7 : 5;
#endif

template <Backend BACKEND, bool COMPLEX>
__host__ __device__ __forceinline__ int extraction_max_ufp(const unsigned k) {
    if constexpr (BACKEND == Backend::FP8) {
        return maxUFP<BACKEND>;
    } else {
#if defined(__CUDA_ARCH__) || (defined(__HIP_DEVICE_COMPILE__) && __HIP_DEVICE_COMPILE__)
        const int ufp = max(0, __clz(k) - 1 - int(COMPLEX)) / 2 - 1;
        return min(maxUFP<BACKEND>, ufp);
#else
        const int ufp = std::max(0, std::countl_zero(k) - 1 - int(COMPLEX)) / 2 - 1;
        return std::min(maxUFP<BACKEND>, ufp);
#endif
    }
}

template <Backend BACKEND, bool COMPLEX>
__host__ __device__ __forceinline__ bool norm_requires_fast_scaling(const size_t k) {
    return BACKEND == Backend::INT8 && k > size_t(0x7fffffffU >> unsigned(COMPLEX));
}

inline constexpr float mh4u_ru = -0x1.0000060000000p-1F; // -1 * round_up(0.5 / (1 - 4 * 2^-24))

template <unsigned NUM_MODULI, Backend BACKEND = Backend::INT8, bool COMPLEX = false>
__device__ __forceinline__ int32_t calc_sft(int32_t amax) {
    if (amax == 0) return 0;
    constexpr float log2P = common::table::log2P<Backend::INT8, NUM_MODULI, COMPLEX>;
    const float log2amax  = __log2f(__int2float_ru(amax));
    return __float2int_rd(__fmaf_rd(mh4u_ru, log2amax, log2P));
}

template <unsigned NUM_MODULI, Backend BACKEND = Backend::FP8, bool COMPLEX = false>
__device__ __forceinline__ int32_t calc_sft(float amax) {
    if (amax == 0.0f) return 0;
    constexpr float log2P = common::table::log2P<BACKEND, NUM_MODULI, COMPLEX>;
    const float log2amax  = __log2f(amax);
    return __float2int_rd(__fmaf_rd(mh4u_ru, log2amax, log2P));
}

template <Backend BACKEND, unsigned NUM_MODULI, bool COMPLEX = false>
__global__ void calc_sft_delta(
    const unsigned len,
    const int16_t *const __restrict__ sft_final,
    int16_t *const __restrict__ sft_saved_or_delta //
) {
    constexpr float log2P = common::table::log2P<BACKEND, NUM_MODULI, COMPLEX>;

    const unsigned idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= len) return;

    const int32_t sft1 = static_cast<int32_t>(sft_saved_or_delta[idx]);
    const int32_t sft2 = -static_cast<int32_t>(sft_final[idx]);

    const int32_t diff = sft2 - sft1;

    const float x       = static_cast<float>(diff) - log2P;
    const int32_t delta = __float2int_ru(x);

    sft_saved_or_delta[idx] = static_cast<int16_t>(delta);
}

template <typename T, Backend BACKEND, cublasFillMode_t UPLO, cublasDiagType_t DIAG>
__global__ void calc_sft_before_rowwise(
    const unsigned rows_A, const unsigned cols_A,
    const T *const __restrict__ A, const size_t lda,
    int16_t *const __restrict__ sftA //
) {
    using U = common::underlying_t<T>;
    __shared__ U samax[common::TILE_DIM][common::TILE_DIM + 1];

    U amax = find_amax_tile<T, UPLO, DIAG>(rows_A, cols_A, A, lda, samax);
    if constexpr (DIAG == CUBLAS_DIAG_UNIT) {
        constexpr U Uone = common::Tconst<U>::one();
        amax             = max(amax, Uone);
    }

    const unsigned row_idx = blockIdx.x * common::TILE_DIM + threadIdx.y;
    if (row_idx < rows_A && threadIdx.x == 0) {
        sftA[row_idx] = extraction_max_ufp<BACKEND, common::isComplex<T>>(cols_A) - common::Tilogb<U>(amax);
    }
}

template <typename T, cublasFillMode_t UPLO, bool HERM>
__global__ void calc_sft_before_sym_rowwise(
    const unsigned n,
    const T *const __restrict__ A, const size_t lda,
    common::underlying_t<T> *const __restrict__ amaxA //
) {
    using U = common::underlying_t<T>;
    __shared__ U samax[common::TILE_DIM][common::TILE_DIM + 1];

    const U amax = find_amax_tile<T, UPLO, CUBLAS_DIAG_NON_UNIT, HERM>(n, n, A, lda, samax);

    const unsigned row_idx = blockIdx.x * common::TILE_DIM + threadIdx.y;
    if (row_idx < n && threadIdx.x == 0) {
        amaxA[row_idx] = amax;
    }
}

template <typename T, Backend BACKEND, cublasFillMode_t UPLO, cublasDiagType_t DIAG>
__device__ __forceinline__ int32_t calc_sft_before_colwise(
    const unsigned rows_A,
    const T *const __restrict__ in,
    int16_t *const __restrict__ sftA,
    common::underlying_t<T> *samax //
) {
    using U = common::underlying_t<T>;

    U amax = find_amax<T, UPLO, DIAG>(rows_A, in, samax);
    if constexpr (DIAG == CUBLAS_DIAG_UNIT) {
        constexpr U Uone = common::Tconst<U>::one();
        amax             = max(amax, Uone);
    }

    const int32_t sft = extraction_max_ufp<BACKEND, common::isComplex<T>>(rows_A) - common::Tilogb<U>(amax);
    if (threadIdx.x == 0) {
        const unsigned col_idx = blockIdx.x;
        sftA[col_idx]          = int16_t(sft);
    }

    return sft;
}

template <typename T, Backend BACKEND, cublasFillMode_t UPLO>
__global__ void calc_sft_before_sym_colwise(
    const unsigned n,
    const T *const __restrict__ A, const size_t lda,
    int16_t *const __restrict__ sftA,
    common::underlying_t<T> *const __restrict__ amaxA //
) {
    using U = common::underlying_t<T>;
    __shared__ U samax[32];

    const unsigned col_idx = blockIdx.x;
    if (col_idx >= n) return;

    const T *const __restrict__ in = A + col_idx * lda;

    const U amax_col = find_amax<T, UPLO, CUBLAS_DIAG_UNIT>(n, in, samax);

    if (threadIdx.x == 0) {
        const U amax_all  = common::Tmax<U>(amaxA[col_idx], amax_col);
        const int32_t sft = extraction_max_ufp<BACKEND, common::isComplex<T>>(n) - common::Tilogb<U>(amax_all);
        sftA[col_idx]     = int16_t(sft);
    }
}

template <typename T, cublasFillMode_t UPLO, bool HERM>
__global__ void calc_sft_before_sym_rowwise_partial(
    const unsigned n,
    const T *const __restrict__ A, const size_t lda,
    common::underlying_t<T> *const __restrict__ partial_amax //
) {
    using U = common::underlying_t<T>;

    __shared__ U samax[common::TILE_DIM][common::TILE_DIM + 1];

    const unsigned col_begin = blockIdx.y * rowwise_sft_col_tile;
    const unsigned col_end   = min(n, col_begin + rowwise_sft_col_tile);

    const U amax = find_amax_tile_range<T, UPLO, CUBLAS_DIAG_NON_UNIT, HERM>(
        n, n, A, lda, samax, col_begin, col_end);

    const unsigned row_idx = blockIdx.x * common::TILE_DIM + threadIdx.y;

    if (row_idx < n && threadIdx.x == 0) {
        partial_amax[size_t(row_idx) * gridDim.y + blockIdx.y] = amax;
    }
}

template <typename T>
__global__ void calc_sft_before_sym_rowwise_reduce(
    const unsigned n,
    const unsigned num_col_blocks,
    const common::underlying_t<T> *const __restrict__ partial_amax,
    common::underlying_t<T> *const __restrict__ amaxA //
) {
    using U = common::underlying_t<T>;

    const unsigned row_idx = blockIdx.x * blockDim.y + threadIdx.y;

    U amax = common::Tconst<U>::zero();

    if (row_idx < n) {
        for (unsigned cb = threadIdx.x; cb < num_col_blocks; cb += blockDim.x) {
            amax = max(amax, partial_amax[size_t(row_idx) * num_col_blocks + cb]);
        }

        amax = common::inner_warp_max<U>(amax);

        if (threadIdx.x == 0) {
            amaxA[row_idx] = amax;
        }
    }
}

template <typename T, cublasFillMode_t UPLO, bool HERM>
inline void calc_sft_before_sym_rowwise_launch(
    const cudaStream_t stream,
    const unsigned n,
    const T *const A, const size_t lda,
    common::underlying_t<T> *const partial_amax,
    common::underlying_t<T> *const amaxA //
) {
    constexpr dim3 threads_findmax(threads_x_findmax_tile,
                                   threads_y_findmax_tile);

    if (n < rowwise_sft_split_threshold) {
        const unsigned grid_sft =
            (n + common::TILE_DIM - 1U) / common::TILE_DIM;

        calc_sft_before_sym_rowwise<T, UPLO, HERM>
            <<<grid_sft, threads_findmax, 0, stream>>>(
                n, A, lda, amaxA);
        return;
    }

    const unsigned num_col_blocks =
        (n + rowwise_sft_col_tile - 1U) / rowwise_sft_col_tile;

    const dim3 grid_partial((n + common::TILE_DIM - 1U) / common::TILE_DIM,
                            num_col_blocks);

    calc_sft_before_sym_rowwise_partial<T, UPLO, HERM>
        <<<grid_partial, threads_findmax, 0, stream>>>(
            n, A, lda, partial_amax);

    constexpr dim3 threads_reduce(threads_x_rowwise_sft_reduce,
                                  threads_y_rowwise_sft_reduce);

    const unsigned grid_reduce =
        (n + threads_y_rowwise_sft_reduce - 1U) / threads_y_rowwise_sft_reduce;

    calc_sft_before_sym_rowwise_reduce<T>
        <<<grid_reduce, threads_reduce, 0, stream>>>(
            n, num_col_blocks, partial_amax, amaxA);
}

template <typename T, Backend BACKEND, unsigned NUM_MODULI, cublasFillMode_t UPLO>
__global__ void calc_sft_after_rowwise(
    const unsigned m, const unsigned n, const unsigned k,
    common::matptr_t<common::hi_t<BACKEND>, common::isComplex<T>> C_hi,
    const size_t ldc_hi,
    int16_t *const __restrict__ sftA //
) {
    __shared__ common::hi_t<BACKEND> samax[common::TILE_DIM][common::TILE_DIM + 1];

    const common::hi_t<BACKEND> amax = find_max_tile<T, BACKEND, UPLO>(m, n, k, C_hi, ldc_hi, samax);

    const unsigned row_idx = blockIdx.x * common::TILE_DIM + threadIdx.y;
    if (row_idx < m && threadIdx.x == 0) {
        const int32_t sft_pre = sftA[row_idx];
        const int32_t sft_cor = calc_sft<NUM_MODULI, BACKEND, common::isComplex<T>>(amax);
        const int32_t sft_new = sft_pre + sft_cor;
        sftA[row_idx]         = static_cast<int16_t>(-sft_new);
    }
}

template <typename T, Backend BACKEND, unsigned NUM_MODULI, cublasFillMode_t UPLO>
__global__ void calc_sft_after_rowwise_with_delta(
    const unsigned m, const unsigned n, const unsigned k,
    common::matptr_t<common::hi_t<BACKEND>, common::isComplex<T>> C_hi,
    const size_t ldc_hi,
    const int16_t *const __restrict__ col_delta,
    int16_t *const __restrict__ sftA //
) {
    __shared__ float samax[common::TILE_DIM][common::TILE_DIM + 1];

    const float amax = find_max_tile_with_delta<T, BACKEND, UPLO>(m, n, k, C_hi, ldc_hi, col_delta, samax);

    const unsigned row_idx = blockIdx.x * common::TILE_DIM + threadIdx.y;
    if (row_idx < m && threadIdx.x == 0) {
        const int32_t sft_pre = sftA[row_idx];
        const int32_t sft_cor = calc_sft<NUM_MODULI, BACKEND, common::isComplex<T>>(amax);
        const int32_t sft_new = sft_pre + sft_cor;
        sftA[row_idx]         = static_cast<int16_t>(-sft_new);
    }
}

template <typename T, Backend BACKEND, unsigned NUM_MODULI, cublasFillMode_t UPLO>
__global__ void calc_sft_after_colwise(
    const unsigned m, const unsigned k,
    common::matptr_t<common::hi_t<BACKEND>, common::isComplex<T>> C_hi,
    const size_t ldc_hi,
    int16_t *const __restrict__ sftA //
) {
    __shared__ common::hi_t<BACKEND> samax[32];
    const common::hi_t<BACKEND> amax = find_max<T, BACKEND, UPLO>(m, k, C_hi, ldc_hi, samax);

    if (threadIdx.x == 0) {
        const unsigned col_idx = blockIdx.x;
        const int32_t sft_pre  = sftA[col_idx];
        const int32_t sft_cor  = calc_sft<NUM_MODULI, BACKEND, common::isComplex<T>>(amax);
        const int32_t sft_new  = sft_pre + sft_cor;
        sftA[col_idx]          = static_cast<int16_t>(-sft_new);
    }
}

template <typename T, Backend BACKEND, unsigned NUM_MODULI, cublasFillMode_t UPLO>
__global__ void calc_sft_after_colwise_with_delta(
    const unsigned m, const unsigned k,
    common::matptr_t<common::hi_t<BACKEND>, common::isComplex<T>> C_hi,
    const size_t ldc_hi,
    const int16_t *const __restrict__ row_delta,
    int16_t *const __restrict__ sftA //
) {
    __shared__ float samax[32];

    const float amax = find_max_with_delta<T, BACKEND, UPLO>(m, k, C_hi, ldc_hi, row_delta, samax);

    if (threadIdx.x == 0) {
        const unsigned col_idx = blockIdx.x;
        const int32_t sft_pre  = sftA[col_idx];
        const int32_t sft_cor  = calc_sft<NUM_MODULI, BACKEND, common::isComplex<T>>(amax);
        const int32_t sft_new  = sft_pre + sft_cor;
        sftA[col_idx]          = static_cast<int16_t>(-sft_new);
    }
}

template <typename T, Backend BACKEND, unsigned NUM_MODULI, cublasFillMode_t UPLO>
__global__ void calc_sft_after_rowwise_sym(
    const unsigned m, const unsigned n, const unsigned k,
    common::matptr_t<common::hi_t<BACKEND>, common::isComplex<T>> C_hi,
    const size_t ldc_hi,
    common::hi_t<BACKEND> *const __restrict__ maxC //
) {
    __shared__ common::hi_t<BACKEND> samax[common::TILE_DIM][common::TILE_DIM + 1];

    const common::hi_t<BACKEND> amax = find_max_tile<T, BACKEND, UPLO>(m, n, k, C_hi, ldc_hi, samax);

    const unsigned row_idx = blockIdx.x * common::TILE_DIM + threadIdx.y;
    if (row_idx < m && threadIdx.x == 0) {
        maxC[row_idx] = amax;
    }
}

template <typename T, Backend BACKEND, unsigned NUM_MODULI, cublasFillMode_t UPLO>
__global__ void calc_sft_after_colwise_sym(
    const unsigned m, const unsigned k,
    common::matptr_t<common::hi_t<BACKEND>, common::isComplex<T>> C_hi,
    const size_t ldc_hi,
    common::hi_t<BACKEND> *const __restrict__ maxC,
    int16_t *const __restrict__ sftA //
) {
    __shared__ common::hi_t<BACKEND> samax[32];
    const common::hi_t<BACKEND> amax = find_max<T, BACKEND, UPLO, true>(m, k, C_hi, ldc_hi, samax);

    if (threadIdx.x == 0) {
        const unsigned col_idx               = blockIdx.x;
        const common::hi_t<BACKEND> amax_pre = maxC[col_idx];
        const common::hi_t<BACKEND> amax_new = max(amax_pre, amax);
        const int32_t sft_pre                = sftA[col_idx];
        const int32_t sft_cor                = calc_sft<NUM_MODULI, BACKEND, common::isComplex<T>>(amax_new);
        const int32_t sft_new                = sft_pre + sft_cor;
        sftA[col_idx]                        = static_cast<int16_t>(-sft_new);
    }
}

template <typename T, Backend BACKEND, cublasFillMode_t UPLO>
__global__ void calc_sft_after_rowwise_sym_partial(
    const unsigned m, const unsigned n, const unsigned k,
    common::matptr_t<common::hi_t<BACKEND>, common::isComplex<T>> C_hi,
    const size_t ldc_hi,
    common::hi_t<BACKEND> *const __restrict__ partial_amax //
) {
    using HiT = common::hi_t<BACKEND>;

    __shared__ HiT samax[common::TILE_DIM][common::TILE_DIM + 1];

    const unsigned col_begin = blockIdx.y * rowwise_sft_col_tile;
    const unsigned col_end   = min(n, col_begin + rowwise_sft_col_tile);

    const HiT amax = find_max_tile_range<T, BACKEND, UPLO>(
        m, n, k, C_hi, ldc_hi, samax, col_begin, col_end);

    const unsigned row_idx = blockIdx.x * common::TILE_DIM + threadIdx.y;

    if (row_idx < m && threadIdx.x == 0) {
        partial_amax[size_t(row_idx) * gridDim.y + blockIdx.y] = amax;
    }
}

template <Backend BACKEND>
__global__ void calc_sft_after_rowwise_sym_reduce(
    const unsigned m,
    const unsigned num_col_blocks,
    const common::hi_t<BACKEND> *const __restrict__ partial_amax,
    common::hi_t<BACKEND> *const __restrict__ maxC //
) {
    using HiT = common::hi_t<BACKEND>;

    const unsigned row_idx = blockIdx.x * blockDim.y + threadIdx.y;

    HiT amax = common::Tconst<HiT>::zero();

    if (row_idx < m) {
        for (unsigned cb = threadIdx.x; cb < num_col_blocks; cb += blockDim.x) {
            amax = max(amax, partial_amax[size_t(row_idx) * num_col_blocks + cb]);
        }

        amax = common::inner_warp_max<HiT>(amax);

        if (threadIdx.x == 0) {
            maxC[row_idx] = amax;
        }
    }
}

template <typename T, Backend BACKEND, unsigned NUM_MODULI, cublasFillMode_t UPLO>
inline void calc_sft_after_rowwise_sym_launch(
    const cudaStream_t stream,
    const unsigned m, const unsigned n, const unsigned k,
    common::matptr_t<common::hi_t<BACKEND>, common::isComplex<T>> C_hi,
    const size_t ldc_hi,
    common::hi_t<BACKEND> *const partial_amax,
    common::hi_t<BACKEND> *const maxC //
) {
    constexpr dim3 threads_findmax(threads_x_findmax_tile,
                                   threads_y_findmax_tile);

    if (n < rowwise_sft_split_threshold) {
        const unsigned grid_findmax =
            (m + common::TILE_DIM - 1U) / common::TILE_DIM;

        calc_sft_after_rowwise_sym<T, BACKEND, NUM_MODULI, UPLO>
            <<<grid_findmax, threads_findmax, 0, stream>>>(
                m, n, k, C_hi, ldc_hi, maxC);
        return;
    }

    const unsigned num_col_blocks =
        (n + rowwise_sft_col_tile - 1U) / rowwise_sft_col_tile;

    const dim3 grid_partial((m + common::TILE_DIM - 1U) / common::TILE_DIM,
                            num_col_blocks);

    calc_sft_after_rowwise_sym_partial<T, BACKEND, UPLO>
        <<<grid_partial, threads_findmax, 0, stream>>>(
            m, n, k, C_hi, ldc_hi, partial_amax);

    constexpr dim3 threads_reduce(threads_x_rowwise_sft_reduce,
                                  threads_y_rowwise_sft_reduce);

    const unsigned grid_reduce =
        (m + threads_y_rowwise_sft_reduce - 1U) / threads_y_rowwise_sft_reduce;

    calc_sft_after_rowwise_sym_reduce<BACKEND>
        <<<grid_reduce, threads_reduce, 0, stream>>>(
            m, num_col_blocks, partial_amax, maxC);
}

template <typename T, cublasFillMode_t UPLO, cublasDiagType_t DIAG>
__global__ void calc_sft_before_rowwise_partial(
    const unsigned rows_A, const unsigned cols_A,
    const T *const __restrict__ A, const size_t lda,
    common::underlying_t<T> *const __restrict__ partial_amax //
) {
    using U = common::underlying_t<T>;
    __shared__ U samax[common::TILE_DIM][common::TILE_DIM + 1];

    const unsigned col_begin = blockIdx.y * rowwise_sft_col_tile;
    const unsigned col_end   = min(cols_A, col_begin + rowwise_sft_col_tile);

    const U amax = find_amax_tile_range<T, UPLO, DIAG>(
        rows_A, cols_A, A, lda, samax, col_begin, col_end);

    const unsigned row_idx = blockIdx.x * common::TILE_DIM + threadIdx.y;
    if (row_idx < rows_A && threadIdx.x == 0) {
        partial_amax[size_t(row_idx) * gridDim.y + blockIdx.y] = amax;
    }
}

template <typename T, Backend BACKEND, cublasDiagType_t DIAG>
__global__ void calc_sft_before_rowwise_reduce(
    const unsigned rows_A, const unsigned cols_A,
    const unsigned num_col_blocks,
    const common::underlying_t<T> *const __restrict__ partial_amax,
    int16_t *const __restrict__ sftA //
) {
    using U = common::underlying_t<T>;

    const unsigned row_idx = blockIdx.x * blockDim.y + threadIdx.y;

    U amax = common::Tconst<U>::zero();

    if (row_idx < rows_A) {
        for (unsigned cb = threadIdx.x; cb < num_col_blocks; cb += blockDim.x) {
            amax = max(amax, partial_amax[size_t(row_idx) * num_col_blocks + cb]);
        }

        amax = common::inner_warp_max<U>(amax);

        if constexpr (DIAG == CUBLAS_DIAG_UNIT) {
            constexpr U Uone = common::Tconst<U>::one();
            amax             = max(amax, Uone);
        }

        if (threadIdx.x == 0) {
            sftA[row_idx] = int16_t(extraction_max_ufp<BACKEND, common::isComplex<T>>(cols_A) - common::Tilogb<U>(amax));
        }
    }
}

template <typename T, Backend BACKEND, cublasFillMode_t UPLO, cublasDiagType_t DIAG>
inline void calc_sft_before_rowwise_launch(
    const cudaStream_t stream,
    const unsigned rows_A, const unsigned cols_A,
    const T *const A, const size_t lda,
    int16_t *const sftA,
    common::underlying_t<T> *const partial_amax //
) {
    constexpr dim3 threads_findmax(threads_x_findmax_tile, threads_y_findmax_tile);

    if (cols_A < rowwise_sft_split_threshold) {
        const unsigned grid_sft = (rows_A + common::TILE_DIM - 1U) / common::TILE_DIM;
        calc_sft_before_rowwise<T, BACKEND, UPLO, DIAG>
            <<<grid_sft, threads_findmax, 0, stream>>>(
                rows_A, cols_A, A, lda, sftA);
        return;
    }

    const unsigned num_col_blocks =
        (cols_A + rowwise_sft_col_tile - 1U) / rowwise_sft_col_tile;

    dim3 grid_partial((rows_A + common::TILE_DIM - 1U) / common::TILE_DIM,
                      num_col_blocks);

    calc_sft_before_rowwise_partial<T, UPLO, DIAG>
        <<<grid_partial, threads_findmax, 0, stream>>>(
            rows_A, cols_A, A, lda, partial_amax);

    constexpr dim3 threads_reduce(threads_x_rowwise_sft_reduce,
                                  threads_y_rowwise_sft_reduce);

    const unsigned grid_reduce =
        (rows_A + threads_y_rowwise_sft_reduce - 1U) / threads_y_rowwise_sft_reduce;

    calc_sft_before_rowwise_reduce<T, BACKEND, DIAG>
        <<<grid_reduce, threads_reduce, 0, stream>>>(
            rows_A, cols_A, num_col_blocks, partial_amax, sftA);
}

template <typename T, Backend BACKEND, cublasFillMode_t UPLO>
__global__ void calc_sft_after_rowwise_partial(
    const unsigned m, const unsigned n, const unsigned k,
    common::matptr_t<common::hi_t<BACKEND>, common::isComplex<T>> C_hi,
    const size_t ldc_hi,
    common::hi_t<BACKEND> *const __restrict__ partial_amax //
) {
    using HiT = common::hi_t<BACKEND>;
    __shared__ HiT samax[common::TILE_DIM][common::TILE_DIM + 1];

    const unsigned col_begin = blockIdx.y * rowwise_sft_col_tile;
    const unsigned col_end   = min(n, col_begin + rowwise_sft_col_tile);

    const HiT amax = find_max_tile_range<T, BACKEND, UPLO>(
        m, n, k, C_hi, ldc_hi, samax, col_begin, col_end);

    const unsigned row_idx = blockIdx.x * common::TILE_DIM + threadIdx.y;
    if (row_idx < m && threadIdx.x == 0) {
        partial_amax[size_t(row_idx) * gridDim.y + blockIdx.y] = amax;
    }
}

template <typename T, Backend BACKEND, unsigned NUM_MODULI>
__global__ void calc_sft_after_rowwise_reduce(
    const unsigned m,
    const unsigned num_col_blocks,
    const common::hi_t<BACKEND> *const __restrict__ partial_amax,
    int16_t *const __restrict__ sftA //
) {
    using HiT = common::hi_t<BACKEND>;

    const unsigned row_idx = blockIdx.x * blockDim.y + threadIdx.y;

    HiT amax = common::Tconst<HiT>::zero();

    if (row_idx < m) {
        for (unsigned cb = threadIdx.x; cb < num_col_blocks; cb += blockDim.x) {
            amax = max(amax, partial_amax[size_t(row_idx) * num_col_blocks + cb]);
        }

        amax = common::inner_warp_max<HiT>(amax);

        if (threadIdx.x == 0) {
            const int32_t sft_pre = sftA[row_idx];
            const int32_t sft_cor = calc_sft<NUM_MODULI, BACKEND, common::isComplex<T>>(amax);
            const int32_t sft_new = sft_pre + sft_cor;
            sftA[row_idx]         = static_cast<int16_t>(-sft_new);
        }
    }
}

template <typename T, Backend BACKEND, unsigned NUM_MODULI, cublasFillMode_t UPLO>
inline void calc_sft_after_rowwise_launch(
    const cudaStream_t stream,
    const unsigned m, const unsigned n, const unsigned k,
    common::matptr_t<common::hi_t<BACKEND>, common::isComplex<T>> C_hi,
    const size_t ldc_hi,
    int16_t *const sftA,
    common::hi_t<BACKEND> *const partial_amax //
) {
    constexpr dim3 threads_findmax(threads_x_findmax_tile, threads_y_findmax_tile);

    if (n < rowwise_sft_split_threshold) {
        const unsigned grid_findmax = (m + common::TILE_DIM - 1U) / common::TILE_DIM;
        calc_sft_after_rowwise<T, BACKEND, NUM_MODULI, UPLO>
            <<<grid_findmax, threads_findmax, 0, stream>>>(
                m, n, k, C_hi, ldc_hi, sftA);
        return;
    }

    const unsigned num_col_blocks =
        (n + rowwise_sft_col_tile - 1U) / rowwise_sft_col_tile;

    dim3 grid_partial((m + common::TILE_DIM - 1U) / common::TILE_DIM,
                      num_col_blocks);

    calc_sft_after_rowwise_partial<T, BACKEND, UPLO>
        <<<grid_partial, threads_findmax, 0, stream>>>(
            m, n, k, C_hi, ldc_hi, partial_amax);

    constexpr dim3 threads_reduce(threads_x_rowwise_sft_reduce,
                                  threads_y_rowwise_sft_reduce);

    const unsigned grid_reduce =
        (m + threads_y_rowwise_sft_reduce - 1U) / threads_y_rowwise_sft_reduce;

    calc_sft_after_rowwise_reduce<T, BACKEND, NUM_MODULI>
        <<<grid_reduce, threads_reduce, 0, stream>>>(
            m, num_col_blocks, partial_amax, sftA);
}

template <typename T, Backend BACKEND, cublasFillMode_t UPLO>
__global__ void calc_sft_after_rowwise_with_delta_partial(
    const unsigned m, const unsigned n, const unsigned k,
    common::matptr_t<common::hi_t<BACKEND>, common::isComplex<T>> C_hi,
    const size_t ldc_hi,
    const int16_t *const __restrict__ col_delta,
    float *const __restrict__ partial_amax //
) {
    __shared__ float samax[common::TILE_DIM][common::TILE_DIM + 1];

    const unsigned col_begin = blockIdx.y * rowwise_sft_col_tile;
    const unsigned col_end   = min(n, col_begin + rowwise_sft_col_tile);

    const float amax = find_max_tile_range_with_delta<T, BACKEND, UPLO>(
        m, n, k, C_hi, ldc_hi, samax, col_delta, col_begin, col_end);

    const unsigned row_idx = blockIdx.x * common::TILE_DIM + threadIdx.y;
    if (row_idx < m && threadIdx.x == 0) {
        partial_amax[size_t(row_idx) * gridDim.y + blockIdx.y] = amax;
    }
}

template <Backend BACKEND, unsigned NUM_MODULI, bool COMPLEX>
__global__ void calc_sft_after_rowwise_with_delta_reduce(
    const unsigned m,
    const unsigned num_col_blocks,
    const float *const __restrict__ partial_amax,
    int16_t *const __restrict__ sftA //
) {
    const unsigned row_idx = blockIdx.x * blockDim.y + threadIdx.y;

    float amax = 0.0F;

    if (row_idx < m) {
        for (unsigned cb = threadIdx.x; cb < num_col_blocks; cb += blockDim.x) {
            amax = max(amax, partial_amax[size_t(row_idx) * num_col_blocks + cb]);
        }

        amax = common::inner_warp_max<float>(amax);

        if (threadIdx.x == 0) {
            const int32_t sft_pre = sftA[row_idx];
            const int32_t sft_cor = calc_sft<NUM_MODULI, BACKEND, COMPLEX>(amax);
            const int32_t sft_new = sft_pre + sft_cor;
            sftA[row_idx]         = static_cast<int16_t>(-sft_new);
        }
    }
}

template <typename T, Backend BACKEND, unsigned NUM_MODULI, cublasFillMode_t UPLO>
inline void calc_sft_after_rowwise_with_delta_launch(
    const cudaStream_t stream,
    const unsigned m, const unsigned n, const unsigned k,
    common::matptr_t<common::hi_t<BACKEND>, common::isComplex<T>> C_hi,
    const size_t ldc_hi,
    int16_t *const sftA,
    const int16_t *const col_delta,
    float *const partial_amax //
) {
    constexpr dim3 threads_findmax(threads_x_findmax_tile, threads_y_findmax_tile);

    if (n < rowwise_sft_split_threshold) {
        const unsigned grid_findmax = (m + common::TILE_DIM - 1U) / common::TILE_DIM;
        calc_sft_after_rowwise_with_delta<T, BACKEND, NUM_MODULI, UPLO>
            <<<grid_findmax, threads_findmax, 0, stream>>>(
                m, n, k, C_hi, ldc_hi, col_delta, sftA);
        return;
    }

    const unsigned num_col_blocks =
        (n + rowwise_sft_col_tile - 1U) / rowwise_sft_col_tile;

    const dim3 grid_partial((m + common::TILE_DIM - 1U) / common::TILE_DIM,
                            num_col_blocks);

    calc_sft_after_rowwise_with_delta_partial<T, BACKEND, UPLO>
        <<<grid_partial, threads_findmax, 0, stream>>>(
            m, n, k, C_hi, ldc_hi, col_delta, partial_amax);

    constexpr dim3 threads_reduce(threads_x_rowwise_sft_reduce,
                                  threads_y_rowwise_sft_reduce);

    const unsigned grid_reduce =
        (m + threads_y_rowwise_sft_reduce - 1U) / threads_y_rowwise_sft_reduce;

    calc_sft_after_rowwise_with_delta_reduce<BACKEND, NUM_MODULI, common::isComplex<T>>
        <<<grid_reduce, threads_reduce, 0, stream>>>(
            m, num_col_blocks, partial_amax, sftA);
}

} // namespace gemmul8::scaling::accu
