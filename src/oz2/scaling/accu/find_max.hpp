#pragma once
#include "../../common/common.hpp"
#include "../general/helper_triangular.hpp"
#include "../general/helper_hermitian.hpp"

namespace gemmul8::scaling::accu {

inline constexpr float UNIT_ROUNDOFF = 0x1.0000000000000p-24F;

template <typename T>
__device__ __forceinline__ T reduction_max(T amax, T *smem) {
    amax = common::inner_warp_max<T>(amax);

    if ((threadIdx.x & 31) == 0) smem[threadIdx.x >> 5] = amax;
    __syncthreads();

    amax = common::Tconst<T>::zero();
    if (threadIdx.x < 32) {
        if (threadIdx.x < (blockDim.x >> 5)) amax = smem[threadIdx.x];
        amax = common::inner_warp_max<T>(amax);
        if (threadIdx.x == 0) smem[0] = amax;
    }
    __syncthreads();

    return smem[0];
}

__device__ __forceinline__ float fp8_upper_bound(const float nu, const float x) {
    return __fmaf_ru(nu, x, x);
}

__device__ __forceinline__ int32_t complex_int_bound(
    const int32_t x1, // (Re(A)-Im(A)) * (Re(B)-Im(B))
    const int32_t x2  // Re(A)*Im(B) + Im(A)*Re(B)
) {
    const int32_t x3 = x1 + x2; // Re(A)*Re(B) + Im(A)*Im(B)
    return max(x3, x2);
}

__device__ __forceinline__ float complex_fp8_bound(
    const float nu, // correction
    const float x1, // abs(Re(A)-Im(A)) * abs(Re(B)-Im(B))
    const float x2, // Re(A)*Im(B)
    const float x3  // Im(A)*Re(B)
) {
    const float x2_up = fp8_upper_bound(nu, x2); // Re(A)*Im(B)
    const float x3_up = fp8_upper_bound(nu, x3); // Im(A)*Re(B)
    const float im_up = __fadd_ru(x2_up, x3_up); // Re(A)*Im(B) + Im(A)*Re(B)
    const float x1_up = fp8_upper_bound(nu, x1); // abs(Re(A)-Im(A)) * abs(Re(B)-Im(B))
    const float re_up = __fadd_ru(x1_up, im_up); // upper bound of Re(A)*Re(B) + Im(A)*Im(B)
    return re_up;
}

//------------------------------
// accurate mode
//------------------------------
template <typename T,
          cublasFillMode_t UPLO = CUBLAS_FILL_MODE_FULL,
          cublasDiagType_t DIAG = CUBLAS_DIAG_NON_UNIT>
__device__ __forceinline__ common::underlying_t<T> find_amax(
    const unsigned length,
    const T *const __restrict__ ptr,
    common::underlying_t<T> *samax //
) {
    using U = common::underlying_t<T>;
    U amax  = common::Tconst<U>::zero();

    const auto [begin, end] = general::column_load_range<UPLO, DIAG>(length);
    for (unsigned i = begin; i < end; i += blockDim.x) {
        const T tmp = common::Tabs<T>(ptr[i]);
        amax        = common::Tmax<T>(tmp, amax);
    }

    return reduction_max<U>(amax, samax);
}

template <typename T,
          cublasFillMode_t UPLO = CUBLAS_FILL_MODE_FULL,
          cublasDiagType_t DIAG = CUBLAS_DIAG_NON_UNIT,
          bool HERM             = false>
__device__ __forceinline__ common::underlying_t<T> find_amax_tile(
    const unsigned m,
    const unsigned n,
    const T *const __restrict__ A,
    const size_t lda,
    common::underlying_t<T> samax[][common::TILE_DIM + 1] //
) {
    using U          = common::underlying_t<T>;
    unsigned row_idx = blockIdx.x * common::TILE_DIM + threadIdx.x;
    U amax           = common::Tconst<U>::zero();

    if (row_idx < m) {
        const T *const __restrict__ row_ptr = A + row_idx;

        if constexpr (UPLO == CUBLAS_FILL_MODE_UPPER) {

            const unsigned row_base   = row_idx - threadIdx.x;
            constexpr unsigned unit   = (DIAG == CUBLAS_DIAG_UNIT) ? 1U : 0U;
            const unsigned full_begin = row_base + common::TILE_DIM - 1U + unit;
            const unsigned diag_end   = min(n, full_begin);

            for (unsigned col = row_base + threadIdx.y; col < diag_end; col += blockDim.y) {
                if (row_idx + unit <= col) {
                    T val       = row_ptr[col * lda];
                    val         = general::hemm_diag_value<HERM, T>(val, row_idx, col);
                    const T tmp = common::Tabs<T>(val);
                    amax        = common::Tmax<T>(tmp, amax);
                }
            }
            for (unsigned col = full_begin + threadIdx.y; col < n; col += blockDim.y) {
                const T tmp = common::Tabs<T>(row_ptr[col * lda]);
                amax        = common::Tmax<T>(tmp, amax);
            }

        } else if constexpr (UPLO == CUBLAS_FILL_MODE_LOWER) {

            const unsigned row_base      = row_idx - threadIdx.x;
            constexpr unsigned keep_diag = (DIAG == CUBLAS_DIAG_UNIT) ? 0U : 1U;
            constexpr unsigned unit      = (DIAG == CUBLAS_DIAG_UNIT) ? 1U : 0U;
            const unsigned full_end      = min(n, row_base + keep_diag);
            const unsigned diag_end      = min(n, row_base + common::TILE_DIM);

            for (unsigned col = threadIdx.y; col < full_end; col += blockDim.y) {
                const T tmp = common::Tabs<T>(row_ptr[col * lda]);
                amax        = common::Tmax<T>(tmp, amax);
            }
            for (unsigned col = full_end + threadIdx.y; col < diag_end; col += blockDim.y) {
                if (col + unit <= row_idx) {
                    T val       = row_ptr[col * lda];
                    val         = general::hemm_diag_value<HERM, T>(val, row_idx, col);
                    const T tmp = common::Tabs<T>(val);
                    amax        = common::Tmax<T>(tmp, amax);
                }
            }

        } else {

            for (unsigned col = threadIdx.y; col < n; col += blockDim.y) {
                const T tmp = common::Tabs<T>(row_ptr[col * lda]);
                amax        = common::Tmax<T>(tmp, amax);
            }
        }
    }

    samax[threadIdx.y][threadIdx.x] = amax;
    __syncthreads();

    amax = samax[threadIdx.x][threadIdx.y];
    return common::inner_warp_max<U, common::TILE_DIM>(amax);
}

template <cublasFillMode_t UPLO = CUBLAS_FILL_MODE_FULL,
          cublasDiagType_t DIAG = CUBLAS_DIAG_NON_UNIT>
__device__ __forceinline__ int32_t find_max_real_8i(
    const unsigned length,
    const int32_t *const __restrict__ C_hi,
    int32_t *samax //
) {
    int32_t amax = 0;

    const auto [begin, end] = general::column_load_range<UPLO, DIAG>(length);
    for (unsigned i = begin; i < end; i += blockDim.x) {
        amax = max(C_hi[i], amax);
    }

    return reduction_max<int32_t>(amax, samax);
}

template <cublasFillMode_t UPLO = CUBLAS_FILL_MODE_FULL,
          cublasDiagType_t DIAG = CUBLAS_DIAG_NON_UNIT>
__device__ __forceinline__ float find_max_real_8f(
    const unsigned length,
    const unsigned k,
    const float *const __restrict__ C_hi,
    float *samax //
) {
    const float scale = __uint2float_ru(k) * UNIT_ROUNDOFF;
    float amax        = 0.0F;

    const auto [begin, end] = general::column_load_range<UPLO, DIAG>(length);
    for (unsigned i = begin; i < end; i += blockDim.x) {
        amax = max(fp8_upper_bound(scale, C_hi[i]), amax);
    }

    return reduction_max<float>(amax, samax);
}

template <cublasFillMode_t UPLO = CUBLAS_FILL_MODE_FULL,
          cublasDiagType_t DIAG = CUBLAS_DIAG_NON_UNIT>
__device__ __forceinline__ int32_t find_max_complex_8i(
    const unsigned length,
    const int32_t *const __restrict__ C_hi_1,
    const int32_t *const __restrict__ C_hi_2,
    int32_t *samax //
) {
    int32_t amax = 0;

    const auto [begin, end] = general::column_load_range<UPLO, DIAG>(length);
    for (unsigned i = begin; i < end; i += blockDim.x) {
        amax = max(complex_int_bound(C_hi_1[i], C_hi_2[i]), amax);
    }

    return reduction_max<int32_t>(amax, samax);
}

template <cublasFillMode_t UPLO = CUBLAS_FILL_MODE_FULL,
          cublasDiagType_t DIAG = CUBLAS_DIAG_NON_UNIT>
__device__ __forceinline__ float find_max_complex_8f(
    const unsigned length,
    const unsigned k,
    const float *const __restrict__ C_hi_1,
    const float *const __restrict__ C_hi_2,
    const float *const __restrict__ C_hi_3,
    float *samax //
) {
    const float scale = __uint2float_ru(k) * UNIT_ROUNDOFF;
    float amax        = 0.0F;

    const auto [begin, end] = general::column_load_range<UPLO, DIAG>(length);
    for (unsigned i = begin; i < end; i += blockDim.x) {
        amax = max(complex_fp8_bound(scale, C_hi_1[i], C_hi_2[i], C_hi_3[i]), amax);
    }

    return reduction_max<float>(amax, samax);
}

template <typename T, Backend BACKEND,
          cublasFillMode_t UPLO = CUBLAS_FILL_MODE_FULL,
          bool EXCLUDE_DIAG     = false>
__device__ __forceinline__ common::hi_t<BACKEND> find_max(
    const unsigned length,
    const unsigned k,
    common::matptr_t<common::hi_t<BACKEND>, common::isComplex<T>> C_hi,
    const size_t ldc,
    common::hi_t<BACKEND> *samax //
) {
    constexpr cublasDiagType_t DIAG = EXCLUDE_DIAG ? CUBLAS_DIAG_UNIT : CUBLAS_DIAG_NON_UNIT;
    const size_t idx                = blockIdx.x * ldc;
    if constexpr (common::isComplex<T>) {
        if constexpr (BACKEND == Backend::INT8) {
            return find_max_complex_8i<UPLO, DIAG>(length, C_hi.ptr0 + idx, C_hi.ptr1 + idx, samax);
        } else {
            return find_max_complex_8f<UPLO, DIAG>(length, k, C_hi.ptr0 + idx, C_hi.ptr1 + idx, C_hi.ptr2 + idx, samax);
        }
    } else {
        if constexpr (BACKEND == Backend::INT8) {
            return find_max_real_8i<UPLO, DIAG>(length, C_hi.ptr0 + idx, samax);
        } else {
            return find_max_real_8f<UPLO, DIAG>(length, k, C_hi.ptr0 + idx, samax);
        }
    }
}

template <cublasFillMode_t UPLO = CUBLAS_FILL_MODE_FULL>
__device__ __forceinline__ int32_t find_max_tile_real_8i(
    const unsigned m, const unsigned n,
    const int32_t *const __restrict__ C_hi, const size_t ldc,
    int32_t samax[][common::TILE_DIM + 1] //
) {
    unsigned row_idx = blockIdx.x * common::TILE_DIM + threadIdx.x;
    int32_t amax     = 0;

    if (row_idx < m) {
        const int32_t *const __restrict__ row_ptr = C_hi + row_idx;

        if constexpr (UPLO == CUBLAS_FILL_MODE_UPPER) {

            const unsigned row_base   = row_idx - threadIdx.x;
            const unsigned full_begin = row_base + common::TILE_DIM - 1U;
            const unsigned diag_end   = min(n, full_begin);

            for (unsigned col = row_base + threadIdx.y; col < diag_end; col += blockDim.y) {
                if (row_idx <= col) amax = max(row_ptr[col * ldc], amax);
            }
            for (unsigned col = full_begin + threadIdx.y; col < n; col += blockDim.y) {
                amax = max(row_ptr[col * ldc], amax);
            }

        } else if constexpr (UPLO == CUBLAS_FILL_MODE_LOWER) {

            const unsigned row_base = row_idx - threadIdx.x;
            const unsigned full_end = min(n, row_base + 1U);
            const unsigned diag_end = min(n, row_base + common::TILE_DIM);

            for (unsigned col = threadIdx.y; col < full_end; col += blockDim.y) {
                amax = max(row_ptr[col * ldc], amax);
            }
            for (unsigned col = full_end + threadIdx.y; col < diag_end; col += blockDim.y) {
                if (col <= row_idx) amax = max(row_ptr[col * ldc], amax);
            }

        } else {

            for (unsigned col = threadIdx.y; col < n; col += blockDim.y) {
                amax = max(row_ptr[col * ldc], amax);
            }
        }
    }

    samax[threadIdx.y][threadIdx.x] = amax;
    __syncthreads();

    amax = samax[threadIdx.x][threadIdx.y];
    return common::inner_warp_max<int32_t, common::TILE_DIM>(amax);
}

template <cublasFillMode_t UPLO = CUBLAS_FILL_MODE_FULL>
__device__ __forceinline__ float find_max_tile_real_8f(
    const unsigned m, const unsigned n, const unsigned k,
    const float *const __restrict__ C_hi, const size_t ldc,
    float samax[][common::TILE_DIM + 1] //
) {
    unsigned row_idx  = blockIdx.x * common::TILE_DIM + threadIdx.x;
    const float scale = __uint2float_ru(k) * UNIT_ROUNDOFF;
    float amax        = 0.0F;

    if (row_idx < m) {
        const float *const __restrict__ row_ptr = C_hi + row_idx;

        if constexpr (UPLO == CUBLAS_FILL_MODE_UPPER) {

            const unsigned row_base   = row_idx - threadIdx.x;
            const unsigned full_begin = row_base + common::TILE_DIM - 1U;
            const unsigned diag_end   = min(n, full_begin);

            for (unsigned col = row_base + threadIdx.y; col < diag_end; col += blockDim.y) {
                if (row_idx <= col) amax = max(fp8_upper_bound(scale, row_ptr[col * ldc]), amax);
            }
            for (unsigned col = full_begin + threadIdx.y; col < n; col += blockDim.y) {
                amax = max(fp8_upper_bound(scale, row_ptr[col * ldc]), amax);
            }

        } else if constexpr (UPLO == CUBLAS_FILL_MODE_LOWER) {

            const unsigned row_base = row_idx - threadIdx.x;
            const unsigned full_end = min(n, row_base + 1U);
            const unsigned diag_end = min(n, row_base + common::TILE_DIM);

            for (unsigned col = threadIdx.y; col < full_end; col += blockDim.y) {
                amax = max(fp8_upper_bound(scale, row_ptr[col * ldc]), amax);
            }
            for (unsigned col = full_end + threadIdx.y; col < diag_end; col += blockDim.y) {
                if (col <= row_idx) amax = max(fp8_upper_bound(scale, row_ptr[col * ldc]), amax);
            }

        } else {

            for (unsigned col = threadIdx.y; col < n; col += blockDim.y) {
                amax = max(fp8_upper_bound(scale, row_ptr[col * ldc]), amax);
            }
        }
    }

    samax[threadIdx.y][threadIdx.x] = amax;
    __syncthreads();

    amax = samax[threadIdx.x][threadIdx.y];
    return common::inner_warp_max<float, common::TILE_DIM>(amax);
}

template <cublasFillMode_t UPLO = CUBLAS_FILL_MODE_FULL>
__device__ __forceinline__ int32_t find_max_tile_complex_8i(
    const unsigned m, const unsigned n,
    const int32_t *const __restrict__ C_hi_1,
    const int32_t *const __restrict__ C_hi_2,
    const size_t ldc,
    int32_t samax[][common::TILE_DIM + 1] //
) {
    unsigned row_idx = blockIdx.x * common::TILE_DIM + threadIdx.x;
    int32_t amax     = 0;

    if (row_idx < m) {
        if constexpr (UPLO == CUBLAS_FILL_MODE_UPPER) {

            const unsigned row_base   = row_idx - threadIdx.x;
            const unsigned full_begin = row_base + common::TILE_DIM - 1U;
            const unsigned diag_end   = min(n, full_begin);

            for (unsigned col = row_base + threadIdx.y; col < diag_end; col += blockDim.y) {
                if (row_idx <= col) {
                    const size_t i = col * ldc + row_idx;
                    amax           = max(complex_int_bound(C_hi_1[i], C_hi_2[i]), amax);
                }
            }
            for (unsigned col = full_begin + threadIdx.y; col < n; col += blockDim.y) {
                const size_t i = col * ldc + row_idx;
                amax           = max(complex_int_bound(C_hi_1[i], C_hi_2[i]), amax);
            }

        } else if constexpr (UPLO == CUBLAS_FILL_MODE_LOWER) {

            const unsigned row_base = row_idx - threadIdx.x;
            const unsigned full_end = min(n, row_base + 1U);
            const unsigned diag_end = min(n, row_base + common::TILE_DIM);

            for (unsigned col = threadIdx.y; col < full_end; col += blockDim.y) {
                const size_t i = col * ldc + row_idx;
                amax           = max(complex_int_bound(C_hi_1[i], C_hi_2[i]), amax);
            }
            for (unsigned col = full_end + threadIdx.y; col < diag_end; col += blockDim.y) {
                if (col <= row_idx) {
                    const size_t i = col * ldc + row_idx;
                    amax           = max(complex_int_bound(C_hi_1[i], C_hi_2[i]), amax);
                }
            }

        } else {

            for (unsigned col = threadIdx.y; col < n; col += blockDim.y) {
                const size_t i = col * ldc + row_idx;
                amax           = max(complex_int_bound(C_hi_1[i], C_hi_2[i]), amax);
            }
        }
    }

    samax[threadIdx.y][threadIdx.x] = amax;
    __syncthreads();

    amax = samax[threadIdx.x][threadIdx.y];
    return common::inner_warp_max<int32_t, common::TILE_DIM>(amax);
}

template <cublasFillMode_t UPLO = CUBLAS_FILL_MODE_FULL>
__device__ __forceinline__ float find_max_tile_complex_8f(
    const unsigned m, const unsigned n, const unsigned k,
    const float *const __restrict__ C_hi_1,
    const float *const __restrict__ C_hi_2,
    const float *const __restrict__ C_hi_3,
    const size_t ldc,
    float samax[][common::TILE_DIM + 1] //
) {
    unsigned row_idx  = blockIdx.x * common::TILE_DIM + threadIdx.x;
    const float scale = __uint2float_ru(k) * UNIT_ROUNDOFF;
    float amax        = 0.0F;

    if (row_idx < m) {
        if constexpr (UPLO == CUBLAS_FILL_MODE_UPPER) {

            const unsigned row_base   = row_idx - threadIdx.x;
            const unsigned full_begin = row_base + common::TILE_DIM - 1U;
            const unsigned diag_end   = min(n, full_begin);

            for (unsigned col = row_base + threadIdx.y; col < diag_end; col += blockDim.y) {
                if (row_idx <= col) {
                    const size_t i = col * ldc + row_idx;
                    amax           = max(complex_fp8_bound(scale, C_hi_1[i], C_hi_2[i], C_hi_3[i]), amax);
                }
            }
            for (unsigned col = full_begin + threadIdx.y; col < n; col += blockDim.y) {
                const size_t i = col * ldc + row_idx;
                amax           = max(complex_fp8_bound(scale, C_hi_1[i], C_hi_2[i], C_hi_3[i]), amax);
            }

        } else if constexpr (UPLO == CUBLAS_FILL_MODE_LOWER) {

            const unsigned row_base = row_idx - threadIdx.x;
            const unsigned full_end = min(n, row_base + 1U);
            const unsigned diag_end = min(n, row_base + common::TILE_DIM);

            for (unsigned col = threadIdx.y; col < full_end; col += blockDim.y) {
                const size_t i = col * ldc + row_idx;
                amax           = max(complex_fp8_bound(scale, C_hi_1[i], C_hi_2[i], C_hi_3[i]), amax);
            }
            for (unsigned col = full_end + threadIdx.y; col < diag_end; col += blockDim.y) {
                if (col <= row_idx) {
                    const size_t i = col * ldc + row_idx;
                    amax           = max(complex_fp8_bound(scale, C_hi_1[i], C_hi_2[i], C_hi_3[i]), amax);
                }
            }

        } else {

            for (unsigned col = threadIdx.y; col < n; col += blockDim.y) {
                const size_t i = col * ldc + row_idx;
                amax           = max(complex_fp8_bound(scale, C_hi_1[i], C_hi_2[i], C_hi_3[i]), amax);
            }
        }
    }

    samax[threadIdx.y][threadIdx.x] = amax;
    __syncthreads();

    amax = samax[threadIdx.x][threadIdx.y];
    return common::inner_warp_max<float, common::TILE_DIM>(amax);
}

template <typename T, Backend BACKEND, cublasFillMode_t UPLO = CUBLAS_FILL_MODE_FULL>
__device__ __forceinline__ common::hi_t<BACKEND> find_max_tile(
    const unsigned m, const unsigned n, const unsigned k,
    common::matptr_t<common::hi_t<BACKEND>, common::isComplex<T>> C_hi, const size_t ldc,
    common::hi_t<BACKEND> samax[][common::TILE_DIM + 1] //
) {
    if constexpr (common::isComplex<T>) {
        if constexpr (BACKEND == Backend::INT8) {
            return find_max_tile_complex_8i<UPLO>(m, n, C_hi.ptr0, C_hi.ptr1, ldc, samax);
        } else {
            return find_max_tile_complex_8f<UPLO>(m, n, k, C_hi.ptr0, C_hi.ptr1, C_hi.ptr2, ldc, samax);
        }
    } else {
        if constexpr (BACKEND == Backend::INT8) {
            return find_max_tile_real_8i<UPLO>(m, n, C_hi.ptr0, ldc, samax);
        } else {
            return find_max_tile_real_8f<UPLO>(m, n, k, C_hi.ptr0, ldc, samax);
        }
    }
}

template <typename T,
          cublasFillMode_t UPLO = CUBLAS_FILL_MODE_FULL,
          cublasDiagType_t DIAG = CUBLAS_DIAG_NON_UNIT,
          bool HERM             = false>
__device__ __forceinline__ common::underlying_t<T> find_amax_tile_range(
    const unsigned m,
    const unsigned n,
    const T *const __restrict__ A,
    const size_t lda,
    common::underlying_t<T> samax[][common::TILE_DIM + 1],
    const unsigned col_begin,
    const unsigned col_end_in //
) {
    using U          = common::underlying_t<T>;
    unsigned row_idx = blockIdx.x * common::TILE_DIM + threadIdx.x;
    U amax           = common::Tconst<U>::zero();

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
                    const T tmp = common::Tabs<T>(val);
                    amax        = common::Tmax<T>(tmp, amax);
                }
            }

            const unsigned dense_begin = max(col_begin, full_begin);
            for (unsigned col = dense_begin + threadIdx.y; col < col_end; col += blockDim.y) {
                const T tmp = common::Tabs<T>(row_ptr[col * lda]);
                amax        = common::Tmax<T>(tmp, amax);
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
                const T tmp = common::Tabs<T>(row_ptr[col * lda]);
                amax        = common::Tmax<T>(tmp, amax);
            }

            const unsigned diag_begin = max(col_begin, full_end);
            const unsigned diag_end   = min(col_end, diag_end0);
            for (unsigned col = diag_begin + threadIdx.y; col < diag_end; col += blockDim.y) {
                if (col + unit <= row_idx) {
                    T val       = row_ptr[col * lda];
                    val         = general::hemm_diag_value<HERM, T>(val, row_idx, col);
                    const T tmp = common::Tabs<T>(val);
                    amax        = common::Tmax<T>(tmp, amax);
                }
            }

        } else {

            for (unsigned col = col_begin + threadIdx.y; col < col_end; col += blockDim.y) {
                const T tmp = common::Tabs<T>(row_ptr[col * lda]);
                amax        = common::Tmax<T>(tmp, amax);
            }
        }
    }

    samax[threadIdx.y][threadIdx.x] = amax;
    __syncthreads();

    amax = samax[threadIdx.x][threadIdx.y];
    return common::inner_warp_max<U, common::TILE_DIM>(amax);
}

template <cublasFillMode_t UPLO = CUBLAS_FILL_MODE_FULL>
__device__ __forceinline__ int32_t find_max_tile_range_real_8i(
    const unsigned m, const unsigned n,
    const int32_t *const __restrict__ C_hi, const size_t ldc,
    int32_t samax[][common::TILE_DIM + 1],
    const unsigned col_begin,
    const unsigned col_end_in //
) {
    unsigned row_idx = blockIdx.x * common::TILE_DIM + threadIdx.x;
    int32_t amax     = 0;

    const unsigned col_end = min(n, col_end_in);

    if (row_idx < m && col_begin < col_end) {
        const int32_t *const __restrict__ row_ptr = C_hi + row_idx;

        if constexpr (UPLO == CUBLAS_FILL_MODE_UPPER) {

            const unsigned row_base   = row_idx - threadIdx.x;
            const unsigned full_begin = row_base + common::TILE_DIM - 1U;

            const unsigned diag_begin = max(col_begin, row_base);
            const unsigned diag_end   = min(col_end, full_begin);

            for (unsigned col = diag_begin + threadIdx.y; col < diag_end; col += blockDim.y) {
                if (row_idx <= col) amax = max(row_ptr[col * ldc], amax);
            }

            const unsigned dense_begin = max(col_begin, full_begin);
            for (unsigned col = dense_begin + threadIdx.y; col < col_end; col += blockDim.y) {
                amax = max(row_ptr[col * ldc], amax);
            }

        } else if constexpr (UPLO == CUBLAS_FILL_MODE_LOWER) {

            const unsigned row_base  = row_idx - threadIdx.x;
            const unsigned full_end  = min(n, row_base + 1U);
            const unsigned diag_end0 = min(n, row_base + common::TILE_DIM);

            const unsigned dense_begin = col_begin;
            const unsigned dense_end   = min(col_end, full_end);

            for (unsigned col = dense_begin + threadIdx.y; col < dense_end; col += blockDim.y) {
                amax = max(row_ptr[col * ldc], amax);
            }

            const unsigned diag_begin = max(col_begin, full_end);
            const unsigned diag_end   = min(col_end, diag_end0);

            for (unsigned col = diag_begin + threadIdx.y; col < diag_end; col += blockDim.y) {
                if (col <= row_idx) amax = max(row_ptr[col * ldc], amax);
            }

        } else {

            for (unsigned col = col_begin + threadIdx.y; col < col_end; col += blockDim.y) {
                amax = max(row_ptr[col * ldc], amax);
            }
        }
    }

    samax[threadIdx.y][threadIdx.x] = amax;
    __syncthreads();

    amax = samax[threadIdx.x][threadIdx.y];
    return common::inner_warp_max<int32_t, common::TILE_DIM>(amax);
}

template <cublasFillMode_t UPLO = CUBLAS_FILL_MODE_FULL>
__device__ __forceinline__ float find_max_tile_range_real_8f(
    const unsigned m, const unsigned n, const unsigned k,
    const float *const __restrict__ C_hi, const size_t ldc,
    float samax[][common::TILE_DIM + 1],
    const unsigned col_begin,
    const unsigned col_end_in //
) {
    unsigned row_idx  = blockIdx.x * common::TILE_DIM + threadIdx.x;
    const float scale = __uint2float_ru(k) * UNIT_ROUNDOFF;
    float amax        = 0.0F;

    const unsigned col_end = min(n, col_end_in);

    if (row_idx < m && col_begin < col_end) {
        const float *const __restrict__ row_ptr = C_hi + row_idx;

        if constexpr (UPLO == CUBLAS_FILL_MODE_UPPER) {

            const unsigned row_base   = row_idx - threadIdx.x;
            const unsigned full_begin = row_base + common::TILE_DIM - 1U;

            const unsigned diag_begin = max(col_begin, row_base);
            const unsigned diag_end   = min(col_end, full_begin);

            for (unsigned col = diag_begin + threadIdx.y; col < diag_end; col += blockDim.y) {
                if (row_idx <= col) {
                    amax = max(fp8_upper_bound(scale, row_ptr[col * ldc]), amax);
                }
            }

            const unsigned dense_begin = max(col_begin, full_begin);
            for (unsigned col = dense_begin + threadIdx.y; col < col_end; col += blockDim.y) {
                amax = max(fp8_upper_bound(scale, row_ptr[col * ldc]), amax);
            }

        } else if constexpr (UPLO == CUBLAS_FILL_MODE_LOWER) {

            const unsigned row_base  = row_idx - threadIdx.x;
            const unsigned full_end  = min(n, row_base + 1U);
            const unsigned diag_end0 = min(n, row_base + common::TILE_DIM);

            const unsigned dense_begin = col_begin;
            const unsigned dense_end   = min(col_end, full_end);

            for (unsigned col = dense_begin + threadIdx.y; col < dense_end; col += blockDim.y) {
                amax = max(fp8_upper_bound(scale, row_ptr[col * ldc]), amax);
            }

            const unsigned diag_begin = max(col_begin, full_end);
            const unsigned diag_end   = min(col_end, diag_end0);

            for (unsigned col = diag_begin + threadIdx.y; col < diag_end; col += blockDim.y) {
                if (col <= row_idx) {
                    amax = max(fp8_upper_bound(scale, row_ptr[col * ldc]), amax);
                }
            }

        } else {

            for (unsigned col = col_begin + threadIdx.y; col < col_end; col += blockDim.y) {
                amax = max(fp8_upper_bound(scale, row_ptr[col * ldc]), amax);
            }
        }
    }

    samax[threadIdx.y][threadIdx.x] = amax;
    __syncthreads();

    amax = samax[threadIdx.x][threadIdx.y];
    return common::inner_warp_max<float, common::TILE_DIM>(amax);
}

template <cublasFillMode_t UPLO = CUBLAS_FILL_MODE_FULL>
__device__ __forceinline__ int32_t find_max_tile_range_complex_8i(
    const unsigned m, const unsigned n,
    const int32_t *const __restrict__ C_hi_1,
    const int32_t *const __restrict__ C_hi_2,
    const size_t ldc,
    int32_t samax[][common::TILE_DIM + 1],
    const unsigned col_begin,
    const unsigned col_end_in //
) {
    unsigned row_idx = blockIdx.x * common::TILE_DIM + threadIdx.x;
    int32_t amax     = 0;

    const unsigned col_end = min(n, col_end_in);

    if (row_idx < m && col_begin < col_end) {
        if constexpr (UPLO == CUBLAS_FILL_MODE_UPPER) {

            const unsigned row_base   = row_idx - threadIdx.x;
            const unsigned full_begin = row_base + common::TILE_DIM - 1U;

            const unsigned diag_begin = max(col_begin, row_base);
            const unsigned diag_end   = min(col_end, full_begin);

            for (unsigned col = diag_begin + threadIdx.y; col < diag_end; col += blockDim.y) {
                if (row_idx <= col) {
                    const size_t i = col * ldc + row_idx;
                    amax           = max(complex_int_bound(C_hi_1[i], C_hi_2[i]), amax);
                }
            }

            const unsigned dense_begin = max(col_begin, full_begin);
            for (unsigned col = dense_begin + threadIdx.y; col < col_end; col += blockDim.y) {
                const size_t i = col * ldc + row_idx;
                amax           = max(complex_int_bound(C_hi_1[i], C_hi_2[i]), amax);
            }

        } else if constexpr (UPLO == CUBLAS_FILL_MODE_LOWER) {

            const unsigned row_base  = row_idx - threadIdx.x;
            const unsigned full_end  = min(n, row_base + 1U);
            const unsigned diag_end0 = min(n, row_base + common::TILE_DIM);

            const unsigned dense_begin = col_begin;
            const unsigned dense_end   = min(col_end, full_end);

            for (unsigned col = dense_begin + threadIdx.y; col < dense_end; col += blockDim.y) {
                const size_t i = col * ldc + row_idx;
                amax           = max(complex_int_bound(C_hi_1[i], C_hi_2[i]), amax);
            }

            const unsigned diag_begin = max(col_begin, full_end);
            const unsigned diag_end   = min(col_end, diag_end0);

            for (unsigned col = diag_begin + threadIdx.y; col < diag_end; col += blockDim.y) {
                if (col <= row_idx) {
                    const size_t i = col * ldc + row_idx;
                    amax           = max(complex_int_bound(C_hi_1[i], C_hi_2[i]), amax);
                }
            }

        } else {

            for (unsigned col = col_begin + threadIdx.y; col < col_end; col += blockDim.y) {
                const size_t i = col * ldc + row_idx;
                amax           = max(complex_int_bound(C_hi_1[i], C_hi_2[i]), amax);
            }
        }
    }

    samax[threadIdx.y][threadIdx.x] = amax;
    __syncthreads();

    amax = samax[threadIdx.x][threadIdx.y];
    return common::inner_warp_max<int32_t, common::TILE_DIM>(amax);
}

template <cublasFillMode_t UPLO = CUBLAS_FILL_MODE_FULL>
__device__ __forceinline__ float find_max_tile_range_complex_8f(
    const unsigned m, const unsigned n, const unsigned k,
    const float *const __restrict__ C_hi_1,
    const float *const __restrict__ C_hi_2,
    const float *const __restrict__ C_hi_3,
    const size_t ldc,
    float samax[][common::TILE_DIM + 1],
    const unsigned col_begin,
    const unsigned col_end_in //
) {
    unsigned row_idx  = blockIdx.x * common::TILE_DIM + threadIdx.x;
    const float scale = __uint2float_ru(k) * UNIT_ROUNDOFF;
    float amax        = 0.0F;

    const unsigned col_end = min(n, col_end_in);

    if (row_idx < m && col_begin < col_end) {
        if constexpr (UPLO == CUBLAS_FILL_MODE_UPPER) {

            const unsigned row_base   = row_idx - threadIdx.x;
            const unsigned full_begin = row_base + common::TILE_DIM - 1U;

            const unsigned diag_begin = max(col_begin, row_base);
            const unsigned diag_end   = min(col_end, full_begin);

            for (unsigned col = diag_begin + threadIdx.y; col < diag_end; col += blockDim.y) {
                if (row_idx <= col) {
                    const size_t i = col * ldc + row_idx;
                    amax           = max(complex_fp8_bound(scale, C_hi_1[i], C_hi_2[i], C_hi_3[i]), amax);
                }
            }

            const unsigned dense_begin = max(col_begin, full_begin);
            for (unsigned col = dense_begin + threadIdx.y; col < col_end; col += blockDim.y) {
                const size_t i = col * ldc + row_idx;
                amax           = max(complex_fp8_bound(scale, C_hi_1[i], C_hi_2[i], C_hi_3[i]), amax);
            }

        } else if constexpr (UPLO == CUBLAS_FILL_MODE_LOWER) {

            const unsigned row_base  = row_idx - threadIdx.x;
            const unsigned full_end  = min(n, row_base + 1U);
            const unsigned diag_end0 = min(n, row_base + common::TILE_DIM);

            const unsigned dense_begin = col_begin;
            const unsigned dense_end   = min(col_end, full_end);

            for (unsigned col = dense_begin + threadIdx.y; col < dense_end; col += blockDim.y) {
                const size_t i = col * ldc + row_idx;
                amax           = max(complex_fp8_bound(scale, C_hi_1[i], C_hi_2[i], C_hi_3[i]), amax);
            }

            const unsigned diag_begin = max(col_begin, full_end);
            const unsigned diag_end   = min(col_end, diag_end0);

            for (unsigned col = diag_begin + threadIdx.y; col < diag_end; col += blockDim.y) {
                if (col <= row_idx) {
                    const size_t i = col * ldc + row_idx;
                    amax           = max(complex_fp8_bound(scale, C_hi_1[i], C_hi_2[i], C_hi_3[i]), amax);
                }
            }

        } else {

            for (unsigned col = col_begin + threadIdx.y; col < col_end; col += blockDim.y) {
                const size_t i = col * ldc + row_idx;
                amax           = max(complex_fp8_bound(scale, C_hi_1[i], C_hi_2[i], C_hi_3[i]), amax);
            }
        }
    }

    samax[threadIdx.y][threadIdx.x] = amax;
    __syncthreads();

    amax = samax[threadIdx.x][threadIdx.y];
    return common::inner_warp_max<float, common::TILE_DIM>(amax);
}

template <typename T, Backend BACKEND, cublasFillMode_t UPLO = CUBLAS_FILL_MODE_FULL>
__device__ __forceinline__ common::hi_t<BACKEND> find_max_tile_range(
    const unsigned m, const unsigned n, const unsigned k,
    common::matptr_t<common::hi_t<BACKEND>, common::isComplex<T>> C_hi, const size_t ldc,
    common::hi_t<BACKEND> samax[][common::TILE_DIM + 1],
    const unsigned col_begin,
    const unsigned col_end //
) {
    if constexpr (common::isComplex<T>) {
        if constexpr (BACKEND == Backend::INT8) {
            return find_max_tile_range_complex_8i<UPLO>(
                m, n, C_hi.ptr0, C_hi.ptr1, ldc, samax, col_begin, col_end);
        } else {
            return find_max_tile_range_complex_8f<UPLO>(
                m, n, k, C_hi.ptr0, C_hi.ptr1, C_hi.ptr2, ldc, samax, col_begin, col_end);
        }
    } else {
        if constexpr (BACKEND == Backend::INT8) {
            return find_max_tile_range_real_8i<UPLO>(
                m, n, C_hi.ptr0, ldc, samax, col_begin, col_end);
        } else {
            return find_max_tile_range_real_8f<UPLO>(
                m, n, k, C_hi.ptr0, ldc, samax, col_begin, col_end);
        }
    }
}

//=====
// for skip-scaling
//=====
template <cublasFillMode_t UPLO = CUBLAS_FILL_MODE_FULL,
          cublasDiagType_t DIAG = CUBLAS_DIAG_NON_UNIT>
__device__ __forceinline__ float find_max_real_8i_with_delta(
    const unsigned length,
    const int32_t *const __restrict__ C_hi,
    const int16_t *const __restrict__ delta,
    float *samax //
) {
    float amax = 0.0F;

    const auto [begin, end] = general::column_load_range<UPLO, DIAG>(length);
    for (unsigned i = begin; i < end; i += blockDim.x) {
        const float v = scalbnf(__int2float_ru(C_hi[i]), -delta[i]);
        amax          = max(v, amax);
    }

    return reduction_max<float>(amax, samax);
}

template <cublasFillMode_t UPLO = CUBLAS_FILL_MODE_FULL,
          cublasDiagType_t DIAG = CUBLAS_DIAG_NON_UNIT>
__device__ __forceinline__ float find_max_real_8f_with_delta(
    const unsigned length,
    const unsigned k,
    const float *const __restrict__ C_hi,
    const int16_t *const __restrict__ delta,
    float *samax //
) {
    const float scale = __uint2float_ru(k) * UNIT_ROUNDOFF;
    float amax        = 0.0F;

    const auto [begin, end] = general::column_load_range<UPLO, DIAG>(length);
    for (unsigned i = begin; i < end; i += blockDim.x) {
        const float v = scalbnf(fp8_upper_bound(scale, C_hi[i]), -delta[i]);
        amax          = max(v, amax);
    }

    return reduction_max<float>(amax, samax);
}

template <cublasFillMode_t UPLO = CUBLAS_FILL_MODE_FULL,
          cublasDiagType_t DIAG = CUBLAS_DIAG_NON_UNIT>
__device__ __forceinline__ float find_max_complex_8i_with_delta(
    const unsigned length,
    const int32_t *const __restrict__ C_hi_1,
    const int32_t *const __restrict__ C_hi_2,
    const int16_t *const __restrict__ delta,
    float *samax //
) {
    float amax = 0.0F;

    const auto [begin, end] = general::column_load_range<UPLO, DIAG>(length);
    for (unsigned i = begin; i < end; i += blockDim.x) {
        const float v = scalbnf(__int2float_ru(complex_int_bound(C_hi_1[i], C_hi_2[i])), -delta[i]);
        amax          = max(v, amax);
    }

    return reduction_max<float>(amax, samax);
}

template <cublasFillMode_t UPLO = CUBLAS_FILL_MODE_FULL,
          cublasDiagType_t DIAG = CUBLAS_DIAG_NON_UNIT>
__device__ __forceinline__ float find_max_complex_8f_with_delta(
    const unsigned length,
    const unsigned k,
    const float *const __restrict__ C_hi_1,
    const float *const __restrict__ C_hi_2,
    const float *const __restrict__ C_hi_3,
    const int16_t *const __restrict__ delta,
    float *samax //
) {
    const float scale = __uint2float_ru(k) * UNIT_ROUNDOFF;
    float amax        = 0.0F;

    const auto [begin, end] = general::column_load_range<UPLO, DIAG>(length);
    for (unsigned i = begin; i < end; i += blockDim.x) {
        const float v = scalbnf(complex_fp8_bound(scale, C_hi_1[i], C_hi_2[i], C_hi_3[i]), -delta[i]);
        amax          = max(v, amax);
    }

    return reduction_max<float>(amax, samax);
}

template <typename T, Backend BACKEND,
          cublasFillMode_t UPLO = CUBLAS_FILL_MODE_FULL,
          bool EXCLUDE_DIAG     = false>
__device__ __forceinline__ float find_max_with_delta(
    const unsigned length,
    const unsigned k,
    common::matptr_t<common::hi_t<BACKEND>, common::isComplex<T>> C_hi,
    const size_t ldc,
    const int16_t *const __restrict__ delta,
    float *samax //
) {
    constexpr cublasDiagType_t DIAG = EXCLUDE_DIAG ? CUBLAS_DIAG_UNIT : CUBLAS_DIAG_NON_UNIT;
    const size_t idx                = blockIdx.x * ldc;
    if constexpr (common::isComplex<T>) {
        if constexpr (BACKEND == Backend::INT8) {
            return find_max_complex_8i_with_delta<UPLO, DIAG>(
                length, C_hi.ptr0 + idx, C_hi.ptr1 + idx, delta, samax);
        } else {
            return find_max_complex_8f_with_delta<UPLO, DIAG>(
                length, k, C_hi.ptr0 + idx, C_hi.ptr1 + idx, C_hi.ptr2 + idx, delta, samax);
        }
    } else {
        if constexpr (BACKEND == Backend::INT8) {
            return find_max_real_8i_with_delta<UPLO, DIAG>(
                length, C_hi.ptr0 + idx, delta, samax);
        } else {
            return find_max_real_8f_with_delta<UPLO, DIAG>(
                length, k, C_hi.ptr0 + idx, delta, samax);
        }
    }
}

template <cublasFillMode_t UPLO = CUBLAS_FILL_MODE_FULL>
__device__ __forceinline__ float find_max_tile_real_8i_with_delta(
    const unsigned m, const unsigned n,
    const int32_t *const __restrict__ C_hi, const size_t ldc,
    const int16_t *const __restrict__ delta,
    float samax[][common::TILE_DIM + 1] //
) {
    unsigned row_idx = blockIdx.x * common::TILE_DIM + threadIdx.x;
    float amax       = 0.0F;

    if (row_idx < m) {
        const int32_t *const __restrict__ row_ptr = C_hi + row_idx;

        if constexpr (UPLO == CUBLAS_FILL_MODE_UPPER) {

            const unsigned row_base   = row_idx - threadIdx.x;
            const unsigned full_begin = row_base + common::TILE_DIM - 1U;
            const unsigned diag_end   = min(n, full_begin);

            for (unsigned col = row_base + threadIdx.y; col < diag_end; col += blockDim.y) {
                if (row_idx <= col) {
                    const float v = scalbnf(__int2float_ru(row_ptr[col * ldc]), -delta[col]);
                    amax          = max(v, amax);
                }
            }
            for (unsigned col = full_begin + threadIdx.y; col < n; col += blockDim.y) {
                const float v = scalbnf(__int2float_ru(row_ptr[col * ldc]), -delta[col]);
                amax          = max(v, amax);
            }

        } else if constexpr (UPLO == CUBLAS_FILL_MODE_LOWER) {

            const unsigned row_base = row_idx - threadIdx.x;
            const unsigned full_end = min(n, row_base + 1U);
            const unsigned diag_end = min(n, row_base + common::TILE_DIM);

            for (unsigned col = threadIdx.y; col < full_end; col += blockDim.y) {
                const float v = scalbnf(__int2float_ru(row_ptr[col * ldc]), -delta[col]);
                amax          = max(v, amax);
            }
            for (unsigned col = full_end + threadIdx.y; col < diag_end; col += blockDim.y) {
                if (col <= row_idx) {
                    const float v = scalbnf(__int2float_ru(row_ptr[col * ldc]), -delta[col]);
                    amax          = max(v, amax);
                }
            }

        } else {

            for (unsigned col = threadIdx.y; col < n; col += blockDim.y) {
                const float v = scalbnf(__int2float_ru(row_ptr[col * ldc]), -delta[col]);
                amax          = max(v, amax);
            }
        }
    }

    samax[threadIdx.y][threadIdx.x] = amax;
    __syncthreads();

    amax = samax[threadIdx.x][threadIdx.y];
    return common::inner_warp_max<float, common::TILE_DIM>(amax);
}

template <cublasFillMode_t UPLO = CUBLAS_FILL_MODE_FULL>
__device__ __forceinline__ float find_max_tile_real_8f_with_delta(
    const unsigned m, const unsigned n, const unsigned k,
    const float *const __restrict__ C_hi, const size_t ldc,
    const int16_t *const __restrict__ delta,
    float samax[][common::TILE_DIM + 1] //
) {
    unsigned row_idx  = blockIdx.x * common::TILE_DIM + threadIdx.x;
    const float scale = __uint2float_ru(k) * UNIT_ROUNDOFF;
    float amax        = 0.0F;

    if (row_idx < m) {
        const float *const __restrict__ row_ptr = C_hi + row_idx;

        if constexpr (UPLO == CUBLAS_FILL_MODE_UPPER) {

            const unsigned row_base   = row_idx - threadIdx.x;
            const unsigned full_begin = row_base + common::TILE_DIM - 1U;
            const unsigned diag_end   = min(n, full_begin);

            for (unsigned col = row_base + threadIdx.y; col < diag_end; col += blockDim.y) {
                if (row_idx <= col) {
                    const float v = scalbnf(fp8_upper_bound(scale, row_ptr[col * ldc]), -delta[col]);
                    amax          = max(v, amax);
                }
            }
            for (unsigned col = full_begin + threadIdx.y; col < n; col += blockDim.y) {
                const float v = scalbnf(fp8_upper_bound(scale, row_ptr[col * ldc]), -delta[col]);
                amax          = max(v, amax);
            }

        } else if constexpr (UPLO == CUBLAS_FILL_MODE_LOWER) {

            const unsigned row_base = row_idx - threadIdx.x;
            const unsigned full_end = min(n, row_base + 1U);
            const unsigned diag_end = min(n, row_base + common::TILE_DIM);

            for (unsigned col = threadIdx.y; col < full_end; col += blockDim.y) {
                const float v = scalbnf(fp8_upper_bound(scale, row_ptr[col * ldc]), -delta[col]);
                amax          = max(v, amax);
            }
            for (unsigned col = full_end + threadIdx.y; col < diag_end; col += blockDim.y) {
                if (col <= row_idx) {
                    const float v = scalbnf(fp8_upper_bound(scale, row_ptr[col * ldc]), -delta[col]);
                    amax          = max(v, amax);
                }
            }

        } else {

            for (unsigned col = threadIdx.y; col < n; col += blockDim.y) {
                const float v = scalbnf(fp8_upper_bound(scale, row_ptr[col * ldc]), -delta[col]);
                amax          = max(v, amax);
            }
        }
    }

    samax[threadIdx.y][threadIdx.x] = amax;
    __syncthreads();

    amax = samax[threadIdx.x][threadIdx.y];
    return common::inner_warp_max<float, common::TILE_DIM>(amax);
}

template <cublasFillMode_t UPLO = CUBLAS_FILL_MODE_FULL>
__device__ __forceinline__ float find_max_tile_complex_8i_with_delta(
    const unsigned m, const unsigned n,
    const int32_t *const __restrict__ C_hi_1,
    const int32_t *const __restrict__ C_hi_2,
    const size_t ldc,
    const int16_t *const __restrict__ delta,
    float samax[][common::TILE_DIM + 1] //
) {
    unsigned row_idx = blockIdx.x * common::TILE_DIM + threadIdx.x;
    float amax       = 0.0F;

    if (row_idx < m) {
        if constexpr (UPLO == CUBLAS_FILL_MODE_UPPER) {

            const unsigned row_base   = row_idx - threadIdx.x;
            const unsigned full_begin = row_base + common::TILE_DIM - 1U;
            const unsigned diag_end   = min(n, full_begin);

            for (unsigned col = row_base + threadIdx.y; col < diag_end; col += blockDim.y) {
                if (row_idx <= col) {
                    const size_t i = col * ldc + row_idx;
                    const float v  = scalbnf(__int2float_ru(complex_int_bound(C_hi_1[i], C_hi_2[i])), -delta[col]);
                    amax           = max(v, amax);
                }
            }
            for (unsigned col = full_begin + threadIdx.y; col < n; col += blockDim.y) {
                const size_t i = col * ldc + row_idx;
                const float v  = scalbnf(__int2float_ru(complex_int_bound(C_hi_1[i], C_hi_2[i])), -delta[col]);
                amax           = max(v, amax);
            }

        } else if constexpr (UPLO == CUBLAS_FILL_MODE_LOWER) {

            const unsigned row_base = row_idx - threadIdx.x;
            const unsigned full_end = min(n, row_base + 1U);
            const unsigned diag_end = min(n, row_base + common::TILE_DIM);

            for (unsigned col = threadIdx.y; col < full_end; col += blockDim.y) {
                const size_t i = col * ldc + row_idx;
                const float v  = scalbnf(__int2float_ru(complex_int_bound(C_hi_1[i], C_hi_2[i])), -delta[col]);
                amax           = max(v, amax);
            }
            for (unsigned col = full_end + threadIdx.y; col < diag_end; col += blockDim.y) {
                if (col <= row_idx) {
                    const size_t i = col * ldc + row_idx;
                    const float v  = scalbnf(__int2float_ru(complex_int_bound(C_hi_1[i], C_hi_2[i])), -delta[col]);
                    amax           = max(v, amax);
                }
            }

        } else {

            for (unsigned col = threadIdx.y; col < n; col += blockDim.y) {
                const size_t i = col * ldc + row_idx;
                const float v  = scalbnf(__int2float_ru(complex_int_bound(C_hi_1[i], C_hi_2[i])), -delta[col]);
                amax           = max(v, amax);
            }
        }
    }

    samax[threadIdx.y][threadIdx.x] = amax;
    __syncthreads();

    amax = samax[threadIdx.x][threadIdx.y];
    return common::inner_warp_max<float, common::TILE_DIM>(amax);
}

template <cublasFillMode_t UPLO = CUBLAS_FILL_MODE_FULL>
__device__ __forceinline__ float find_max_tile_complex_8f_with_delta(
    const unsigned m, const unsigned n, const unsigned k,
    const float *const __restrict__ C_hi_1,
    const float *const __restrict__ C_hi_2,
    const float *const __restrict__ C_hi_3,
    const size_t ldc,
    const int16_t *const __restrict__ delta,
    float samax[][common::TILE_DIM + 1] //
) {
    unsigned row_idx  = blockIdx.x * common::TILE_DIM + threadIdx.x;
    const float scale = __uint2float_ru(k) * UNIT_ROUNDOFF;
    float amax        = 0.0F;

    if (row_idx < m) {
        if constexpr (UPLO == CUBLAS_FILL_MODE_UPPER) {

            const unsigned row_base   = row_idx - threadIdx.x;
            const unsigned full_begin = row_base + common::TILE_DIM - 1U;
            const unsigned diag_end   = min(n, full_begin);

            for (unsigned col = row_base + threadIdx.y; col < diag_end; col += blockDim.y) {
                if (row_idx <= col) {
                    const size_t i = col * ldc + row_idx;
                    const float v  = scalbnf(complex_fp8_bound(scale, C_hi_1[i], C_hi_2[i], C_hi_3[i]), -delta[col]);
                    amax           = max(v, amax);
                }
            }
            for (unsigned col = full_begin + threadIdx.y; col < n; col += blockDim.y) {
                const size_t i = col * ldc + row_idx;
                const float v  = scalbnf(complex_fp8_bound(scale, C_hi_1[i], C_hi_2[i], C_hi_3[i]), -delta[col]);
                amax           = max(v, amax);
            }

        } else if constexpr (UPLO == CUBLAS_FILL_MODE_LOWER) {

            const unsigned row_base = row_idx - threadIdx.x;
            const unsigned full_end = min(n, row_base + 1U);
            const unsigned diag_end = min(n, row_base + common::TILE_DIM);

            for (unsigned col = threadIdx.y; col < full_end; col += blockDim.y) {
                const size_t i = col * ldc + row_idx;
                const float v  = scalbnf(complex_fp8_bound(scale, C_hi_1[i], C_hi_2[i], C_hi_3[i]), -delta[col]);
                amax           = max(v, amax);
            }
            for (unsigned col = full_end + threadIdx.y; col < diag_end; col += blockDim.y) {
                if (col <= row_idx) {
                    const size_t i = col * ldc + row_idx;
                    const float v  = scalbnf(complex_fp8_bound(scale, C_hi_1[i], C_hi_2[i], C_hi_3[i]), -delta[col]);
                    amax           = max(v, amax);
                }
            }

        } else {

            for (unsigned col = threadIdx.y; col < n; col += blockDim.y) {
                const size_t i = col * ldc + row_idx;
                const float v  = scalbnf(complex_fp8_bound(scale, C_hi_1[i], C_hi_2[i], C_hi_3[i]), -delta[col]);
                amax           = max(v, amax);
            }
        }
    }

    samax[threadIdx.y][threadIdx.x] = amax;
    __syncthreads();

    amax = samax[threadIdx.x][threadIdx.y];
    return common::inner_warp_max<float, common::TILE_DIM>(amax);
}

template <typename T, Backend BACKEND, cublasFillMode_t UPLO = CUBLAS_FILL_MODE_FULL>
__device__ __forceinline__ float find_max_tile_with_delta(
    const unsigned m, const unsigned n, const unsigned k,
    common::matptr_t<common::hi_t<BACKEND>, common::isComplex<T>> C_hi, const size_t ldc,
    const int16_t *const __restrict__ delta,
    float samax[][common::TILE_DIM + 1] //
) {
    if constexpr (common::isComplex<T>) {
        if constexpr (BACKEND == Backend::INT8) {
            return find_max_tile_complex_8i_with_delta<UPLO>(
                m, n, C_hi.ptr0, C_hi.ptr1, ldc, delta, samax);
        } else {
            return find_max_tile_complex_8f_with_delta<UPLO>(
                m, n, k, C_hi.ptr0, C_hi.ptr1, C_hi.ptr2, ldc, delta, samax);
        }
    } else {
        if constexpr (BACKEND == Backend::INT8) {
            return find_max_tile_real_8i_with_delta<UPLO>(
                m, n, C_hi.ptr0, ldc, delta, samax);
        } else {
            return find_max_tile_real_8f_with_delta<UPLO>(
                m, n, k, C_hi.ptr0, ldc, delta, samax);
        }
    }
}

template <cublasFillMode_t UPLO = CUBLAS_FILL_MODE_FULL>
__device__ __forceinline__ float find_max_tile_range_real_8i_with_delta(
    const unsigned m, const unsigned n,
    const int32_t *const __restrict__ C_hi, const size_t ldc,
    float samax[][common::TILE_DIM + 1],
    const int16_t *const __restrict__ delta,
    const unsigned col_begin,
    const unsigned col_end_in //
) {
    unsigned row_idx = blockIdx.x * common::TILE_DIM + threadIdx.x;
    float amax       = 0.0F;

    const unsigned col_end = min(n, col_end_in);

    if (row_idx < m && col_begin < col_end) {
        const int32_t *const __restrict__ row_ptr = C_hi + row_idx;

        if constexpr (UPLO == CUBLAS_FILL_MODE_UPPER) {

            const unsigned row_base   = row_idx - threadIdx.x;
            const unsigned full_begin = row_base + common::TILE_DIM - 1U;

            const unsigned diag_begin = max(col_begin, row_base);
            const unsigned diag_end   = min(col_end, full_begin);

            for (unsigned col = diag_begin + threadIdx.y; col < diag_end; col += blockDim.y) {
                if (row_idx <= col) {
                    const float v = scalbnf(__int2float_ru(row_ptr[col * ldc]), -delta[col]);
                    amax          = max(v, amax);
                }
            }

            const unsigned dense_begin = max(col_begin, full_begin);
            for (unsigned col = dense_begin + threadIdx.y; col < col_end; col += blockDim.y) {
                const float v = scalbnf(__int2float_ru(row_ptr[col * ldc]), -delta[col]);
                amax          = max(v, amax);
            }

        } else if constexpr (UPLO == CUBLAS_FILL_MODE_LOWER) {

            const unsigned row_base  = row_idx - threadIdx.x;
            const unsigned full_end  = min(n, row_base + 1U);
            const unsigned diag_end0 = min(n, row_base + common::TILE_DIM);

            const unsigned dense_begin = col_begin;
            const unsigned dense_end   = min(col_end, full_end);

            for (unsigned col = dense_begin + threadIdx.y; col < dense_end; col += blockDim.y) {
                const float v = scalbnf(__int2float_ru(row_ptr[col * ldc]), -delta[col]);
                amax          = max(v, amax);
            }

            const unsigned diag_begin = max(col_begin, full_end);
            const unsigned diag_end   = min(col_end, diag_end0);

            for (unsigned col = diag_begin + threadIdx.y; col < diag_end; col += blockDim.y) {
                if (col <= row_idx) {
                    const float v = scalbnf(__int2float_ru(row_ptr[col * ldc]), -delta[col]);
                    amax          = max(v, amax);
                }
            }

        } else {

            for (unsigned col = col_begin + threadIdx.y; col < col_end; col += blockDim.y) {
                const float v = scalbnf(__int2float_ru(row_ptr[col * ldc]), -delta[col]);
                amax          = max(v, amax);
            }
        }
    }

    samax[threadIdx.y][threadIdx.x] = amax;
    __syncthreads();

    amax = samax[threadIdx.x][threadIdx.y];
    return common::inner_warp_max<float, common::TILE_DIM>(amax);
}

template <cublasFillMode_t UPLO = CUBLAS_FILL_MODE_FULL>
__device__ __forceinline__ float find_max_tile_range_real_8f_with_delta(
    const unsigned m, const unsigned n, const unsigned k,
    const float *const __restrict__ C_hi, const size_t ldc,
    float samax[][common::TILE_DIM + 1],
    const int16_t *const __restrict__ delta,
    const unsigned col_begin,
    const unsigned col_end_in //
) {
    unsigned row_idx  = blockIdx.x * common::TILE_DIM + threadIdx.x;
    const float scale = __uint2float_ru(k) * UNIT_ROUNDOFF;
    float amax        = 0.0F;

    const unsigned col_end = min(n, col_end_in);

    if (row_idx < m && col_begin < col_end) {
        const float *const __restrict__ row_ptr = C_hi + row_idx;

        if constexpr (UPLO == CUBLAS_FILL_MODE_UPPER) {

            const unsigned row_base   = row_idx - threadIdx.x;
            const unsigned full_begin = row_base + common::TILE_DIM - 1U;

            const unsigned diag_begin = max(col_begin, row_base);
            const unsigned diag_end   = min(col_end, full_begin);

            for (unsigned col = diag_begin + threadIdx.y; col < diag_end; col += blockDim.y) {
                if (row_idx <= col) {
                    const float v = scalbnf(fp8_upper_bound(scale, row_ptr[col * ldc]), -delta[col]);
                    amax          = max(v, amax);
                }
            }

            const unsigned dense_begin = max(col_begin, full_begin);
            for (unsigned col = dense_begin + threadIdx.y; col < col_end; col += blockDim.y) {
                const float v = scalbnf(fp8_upper_bound(scale, row_ptr[col * ldc]), -delta[col]);
                amax          = max(v, amax);
            }

        } else if constexpr (UPLO == CUBLAS_FILL_MODE_LOWER) {

            const unsigned row_base  = row_idx - threadIdx.x;
            const unsigned full_end  = min(n, row_base + 1U);
            const unsigned diag_end0 = min(n, row_base + common::TILE_DIM);

            const unsigned dense_begin = col_begin;
            const unsigned dense_end   = min(col_end, full_end);

            for (unsigned col = dense_begin + threadIdx.y; col < dense_end; col += blockDim.y) {
                const float v = scalbnf(fp8_upper_bound(scale, row_ptr[col * ldc]), -delta[col]);
                amax          = max(v, amax);
            }

            const unsigned diag_begin = max(col_begin, full_end);
            const unsigned diag_end   = min(col_end, diag_end0);

            for (unsigned col = diag_begin + threadIdx.y; col < diag_end; col += blockDim.y) {
                if (col <= row_idx) {
                    const float v = scalbnf(fp8_upper_bound(scale, row_ptr[col * ldc]), -delta[col]);
                    amax          = max(v, amax);
                }
            }

        } else {

            for (unsigned col = col_begin + threadIdx.y; col < col_end; col += blockDim.y) {
                const float v = scalbnf(fp8_upper_bound(scale, row_ptr[col * ldc]), -delta[col]);
                amax          = max(v, amax);
            }
        }
    }

    samax[threadIdx.y][threadIdx.x] = amax;
    __syncthreads();

    amax = samax[threadIdx.x][threadIdx.y];
    return common::inner_warp_max<float, common::TILE_DIM>(amax);
}

template <cublasFillMode_t UPLO = CUBLAS_FILL_MODE_FULL>
__device__ __forceinline__ float find_max_tile_range_complex_8i_with_delta(
    const unsigned m, const unsigned n,
    const int32_t *const __restrict__ C_hi_1,
    const int32_t *const __restrict__ C_hi_2,
    const size_t ldc,
    float samax[][common::TILE_DIM + 1],
    const int16_t *const __restrict__ delta,
    const unsigned col_begin,
    const unsigned col_end_in //
) {
    unsigned row_idx = blockIdx.x * common::TILE_DIM + threadIdx.x;
    float amax       = 0.0F;

    const unsigned col_end = min(n, col_end_in);

    if (row_idx < m && col_begin < col_end) {
        if constexpr (UPLO == CUBLAS_FILL_MODE_UPPER) {

            const unsigned row_base   = row_idx - threadIdx.x;
            const unsigned full_begin = row_base + common::TILE_DIM - 1U;

            const unsigned diag_begin = max(col_begin, row_base);
            const unsigned diag_end   = min(col_end, full_begin);

            for (unsigned col = diag_begin + threadIdx.y; col < diag_end; col += blockDim.y) {
                if (row_idx <= col) {
                    const size_t i = col * ldc + row_idx;
                    const float v  = scalbnf(__int2float_ru(complex_int_bound(C_hi_1[i], C_hi_2[i])), -delta[col]);
                    amax           = max(v, amax);
                }
            }

            const unsigned dense_begin = max(col_begin, full_begin);
            for (unsigned col = dense_begin + threadIdx.y; col < col_end; col += blockDim.y) {
                const size_t i = col * ldc + row_idx;
                const float v  = scalbnf(__int2float_ru(complex_int_bound(C_hi_1[i], C_hi_2[i])), -delta[col]);
                amax           = max(v, amax);
            }

        } else if constexpr (UPLO == CUBLAS_FILL_MODE_LOWER) {

            const unsigned row_base  = row_idx - threadIdx.x;
            const unsigned full_end  = min(n, row_base + 1U);
            const unsigned diag_end0 = min(n, row_base + common::TILE_DIM);

            const unsigned dense_begin = col_begin;
            const unsigned dense_end   = min(col_end, full_end);

            for (unsigned col = dense_begin + threadIdx.y; col < dense_end; col += blockDim.y) {
                const size_t i = col * ldc + row_idx;
                const float v  = scalbnf(__int2float_ru(complex_int_bound(C_hi_1[i], C_hi_2[i])), -delta[col]);
                amax           = max(v, amax);
            }

            const unsigned diag_begin = max(col_begin, full_end);
            const unsigned diag_end   = min(col_end, diag_end0);

            for (unsigned col = diag_begin + threadIdx.y; col < diag_end; col += blockDim.y) {
                if (col <= row_idx) {
                    const size_t i = col * ldc + row_idx;
                    const float v  = scalbnf(__int2float_ru(complex_int_bound(C_hi_1[i], C_hi_2[i])), -delta[col]);
                    amax           = max(v, amax);
                }
            }

        } else {

            for (unsigned col = col_begin + threadIdx.y; col < col_end; col += blockDim.y) {
                const size_t i = col * ldc + row_idx;
                const float v  = scalbnf(__int2float_ru(complex_int_bound(C_hi_1[i], C_hi_2[i])), -delta[col]);
                amax           = max(v, amax);
            }
        }
    }

    samax[threadIdx.y][threadIdx.x] = amax;
    __syncthreads();

    amax = samax[threadIdx.x][threadIdx.y];
    return common::inner_warp_max<float, common::TILE_DIM>(amax);
}

template <cublasFillMode_t UPLO = CUBLAS_FILL_MODE_FULL>
__device__ __forceinline__ float find_max_tile_range_complex_8f_with_delta(
    const unsigned m, const unsigned n, const unsigned k,
    const float *const __restrict__ C_hi_1,
    const float *const __restrict__ C_hi_2,
    const float *const __restrict__ C_hi_3,
    const size_t ldc,
    float samax[][common::TILE_DIM + 1],
    const int16_t *const __restrict__ delta,
    const unsigned col_begin,
    const unsigned col_end_in //
) {
    unsigned row_idx  = blockIdx.x * common::TILE_DIM + threadIdx.x;
    const float scale = __uint2float_ru(k) * UNIT_ROUNDOFF;
    float amax        = 0.0F;

    const unsigned col_end = min(n, col_end_in);

    if (row_idx < m && col_begin < col_end) {
        if constexpr (UPLO == CUBLAS_FILL_MODE_UPPER) {

            const unsigned row_base   = row_idx - threadIdx.x;
            const unsigned full_begin = row_base + common::TILE_DIM - 1U;

            const unsigned diag_begin = max(col_begin, row_base);
            const unsigned diag_end   = min(col_end, full_begin);

            for (unsigned col = diag_begin + threadIdx.y; col < diag_end; col += blockDim.y) {
                if (row_idx <= col) {
                    const size_t i = col * ldc + row_idx;
                    const float v  = scalbnf(complex_fp8_bound(scale, C_hi_1[i], C_hi_2[i], C_hi_3[i]), -delta[col]);
                    amax           = max(v, amax);
                }
            }

            const unsigned dense_begin = max(col_begin, full_begin);
            for (unsigned col = dense_begin + threadIdx.y; col < col_end; col += blockDim.y) {
                const size_t i = col * ldc + row_idx;
                const float v  = scalbnf(complex_fp8_bound(scale, C_hi_1[i], C_hi_2[i], C_hi_3[i]), -delta[col]);
                amax           = max(v, amax);
            }

        } else if constexpr (UPLO == CUBLAS_FILL_MODE_LOWER) {

            const unsigned row_base  = row_idx - threadIdx.x;
            const unsigned full_end  = min(n, row_base + 1U);
            const unsigned diag_end0 = min(n, row_base + common::TILE_DIM);

            const unsigned dense_begin = col_begin;
            const unsigned dense_end   = min(col_end, full_end);

            for (unsigned col = dense_begin + threadIdx.y; col < dense_end; col += blockDim.y) {
                const size_t i = col * ldc + row_idx;
                const float v  = scalbnf(complex_fp8_bound(scale, C_hi_1[i], C_hi_2[i], C_hi_3[i]), -delta[col]);
                amax           = max(v, amax);
            }

            const unsigned diag_begin = max(col_begin, full_end);
            const unsigned diag_end   = min(col_end, diag_end0);

            for (unsigned col = diag_begin + threadIdx.y; col < diag_end; col += blockDim.y) {
                if (col <= row_idx) {
                    const size_t i = col * ldc + row_idx;
                    const float v  = scalbnf(complex_fp8_bound(scale, C_hi_1[i], C_hi_2[i], C_hi_3[i]), -delta[col]);
                    amax           = max(v, amax);
                }
            }

        } else {

            for (unsigned col = col_begin + threadIdx.y; col < col_end; col += blockDim.y) {
                const size_t i = col * ldc + row_idx;
                const float v  = scalbnf(complex_fp8_bound(scale, C_hi_1[i], C_hi_2[i], C_hi_3[i]), -delta[col]);
                amax           = max(v, amax);
            }
        }
    }

    samax[threadIdx.y][threadIdx.x] = amax;
    __syncthreads();

    amax = samax[threadIdx.x][threadIdx.y];
    return common::inner_warp_max<float, common::TILE_DIM>(amax);
}

template <typename T, Backend BACKEND, cublasFillMode_t UPLO = CUBLAS_FILL_MODE_FULL>
__device__ __forceinline__ float find_max_tile_range_with_delta(
    const unsigned m, const unsigned n, const unsigned k,
    common::matptr_t<common::hi_t<BACKEND>, common::isComplex<T>> C_hi, const size_t ldc,
    float samax[][common::TILE_DIM + 1],
    const int16_t *const __restrict__ delta,
    const unsigned col_begin,
    const unsigned col_end //
) {
    if constexpr (common::isComplex<T>) {
        if constexpr (BACKEND == Backend::INT8) {
            return find_max_tile_range_complex_8i_with_delta<UPLO>(
                m, n, C_hi.ptr0, C_hi.ptr1, ldc, samax, delta, col_begin, col_end);
        } else {
            return find_max_tile_range_complex_8f_with_delta<UPLO>(
                m, n, k, C_hi.ptr0, C_hi.ptr1, C_hi.ptr2, ldc, samax, delta, col_begin, col_end);
        }
    } else {
        if constexpr (BACKEND == Backend::INT8) {
            return find_max_tile_range_real_8i_with_delta<UPLO>(
                m, n, C_hi.ptr0, ldc, samax, delta, col_begin, col_end);
        } else {
            return find_max_tile_range_real_8f_with_delta<UPLO>(
                m, n, k, C_hi.ptr0, ldc, samax, delta, col_begin, col_end);
        }
    }
}

} // namespace gemmul8::scaling::accu
