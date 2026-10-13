#pragma once
#include "../common/common.hpp"
#include "mod_core.hpp"
#include "mod_reduce_matprod_fp8.hpp"
#include <cassert>

#define GEMMUL8_FP8_FOR_EACH_REDUCTION_PERIOD(M) \
    M(13)                                        \
    M(17)                                        \
    M(19)                                        \
    M(23)                                        \
    M(27)                                        \
    M(29)                                        \
    M(31)                                        \
    M(32)                                        \
    M(37)                                        \
    M(39)                                        \
    M(41)                                        \
    M(43)                                        \
    M(45)                                        \
    M(47)                                        \
    M(81)                                        \
    M(135)

namespace gemmul8::mod {

template <Backend BACKEND, unsigned IDX, bool COMPLEX>
__device__ __forceinline__
    common::hi_t<BACKEND>
    mod_reduce_matprod_core(const common::hi_t<BACKEND> c) {
    if constexpr (BACKEND == Backend::INT8) {
        return mod_small<BACKEND, IDX, COMPLEX>(c);
    } else if constexpr (BACKEND == Backend::FP8) {
        constexpr unsigned p = common::table::moduli<Backend::FP8, IDX, COMPLEX>;
        return fp8_matprod_reduction<p>::apply(c);
    } else {
        return common::hi_t<BACKEND>(0);
    }
}

template <Backend BACKEND, unsigned IDX, bool COMPLEX>
__global__ void mod_reduce_matprod_kernel(
    common::hi_t<BACKEND> *__restrict__ C,
    const int64_t m,
    const int64_t n,
    const int64_t ldc //
) {
    const int64_t i = int64_t(blockIdx.x * blockDim.x + threadIdx.x);
    const int64_t j = int64_t(blockIdx.y * blockDim.y + threadIdx.y);

    if (i >= m || j >= n) return;

    common::hi_t<BACKEND> *ptr = C + j * ldc + i;

    *ptr = mod_reduce_matprod_core<BACKEND, IDX, COMPLEX>(*ptr);
}

template <Backend BACKEND, unsigned IDX, bool COMPLEX>
inline void mod_reduce_matprod_launch(
    const cudaStream_t stream,
    common::hi_t<BACKEND> *C,
    const int m,
    const int n,
    const size_t ldc //
) {
    if (m == 0 || n == 0) return;

    if constexpr (BACKEND == Backend::FP8) {
        if ((m % 2) == 0 &&
            (ldc % 2U) == 0U &&
            (reinterpret_cast<uintptr_t>(C) % alignof(float2)) == 0U) {
            constexpr unsigned p = common::table::moduli<BACKEND, IDX, COMPLEX>;
            if (mod_reduce_matprod_fp8_is_flat(m, ldc)) {
                return mod_reduce_matprod_fp8_flat_launch<p, false>(
                    stream, C, ldc * n / 2U, 0U, 1U);
            }
            return mod_reduce_matprod_fp8_tile_launch<p, false>(
                stream, C, m, n, ldc, 0U, 1U);
        }
    }

    constexpr int TX = 32;
    constexpr int TY = 8;

    constexpr int max_columns = 65535 * TY;
    if (n > max_columns) {
        for (int first = 0; first < n;) {
            const int count = std::min(max_columns, n - first);
            mod_reduce_matprod_launch<BACKEND, IDX, COMPLEX>(
                stream, C + first * ldc, m, count, ldc);
            first += count;
        }
        return;
    }

    const dim3 block(TX, TY);
    const dim3 grid((m + TX - 1) / TX, (n + TY - 1) / TY);

    mod_reduce_matprod_kernel<BACKEND, IDX, COMPLEX><<<grid, block, 0, stream>>>(
        C, int64_t(m), int64_t(n), int64_t(ldc));
}

#define GEMMUL8_FP8_REDUCE_CASE(P)                                                         \
    case P + 20U: {                                                                        \
        return mod_reduce_matprod_launch<BACKEND, P + 20U, COMPLEX>(stream, C, m, n, ldc); \
    }

template <Backend BACKEND, bool COMPLEX>
void mod_reduce_matprod_dispatch(
    const cudaStream_t stream,
    common::hi_t<BACKEND> *C,
    const int m,
    const int n,
    const size_t ldc,
    const unsigned modulus_idx //
) {
    if constexpr (BACKEND == Backend::FP8) {
        switch (modulus_idx) {
            GEMMUL8_FP8_FOR_EACH_MODULUS(GEMMUL8_FP8_REDUCE_CASE)
            GEMMUL8_FP8_FOR_EACH_REDUCTION_PERIOD(GEMMUL8_FP8_REDUCE_CASE)
        }
    } else {
        switch (modulus_idx) {
        case 0U: return mod_reduce_matprod_launch<BACKEND, 0U, COMPLEX>(stream, C, m, n, ldc);
        case 1U: return mod_reduce_matprod_launch<BACKEND, 1U, COMPLEX>(stream, C, m, n, ldc);
        case 2U: return mod_reduce_matprod_launch<BACKEND, 2U, COMPLEX>(stream, C, m, n, ldc);
        case 3U: return mod_reduce_matprod_launch<BACKEND, 3U, COMPLEX>(stream, C, m, n, ldc);
        case 4U: return mod_reduce_matprod_launch<BACKEND, 4U, COMPLEX>(stream, C, m, n, ldc);
        case 5U: return mod_reduce_matprod_launch<BACKEND, 5U, COMPLEX>(stream, C, m, n, ldc);
        case 6U: return mod_reduce_matprod_launch<BACKEND, 6U, COMPLEX>(stream, C, m, n, ldc);
        case 7U: return mod_reduce_matprod_launch<BACKEND, 7U, COMPLEX>(stream, C, m, n, ldc);
        case 8U: return mod_reduce_matprod_launch<BACKEND, 8U, COMPLEX>(stream, C, m, n, ldc);
        case 9U: return mod_reduce_matprod_launch<BACKEND, 9U, COMPLEX>(stream, C, m, n, ldc);
        case 10U: return mod_reduce_matprod_launch<BACKEND, 10U, COMPLEX>(stream, C, m, n, ldc);
        case 11U: return mod_reduce_matprod_launch<BACKEND, 11U, COMPLEX>(stream, C, m, n, ldc);
        case 12U: return mod_reduce_matprod_launch<BACKEND, 12U, COMPLEX>(stream, C, m, n, ldc);
        case 13U: return mod_reduce_matprod_launch<BACKEND, 13U, COMPLEX>(stream, C, m, n, ldc);
        case 14U: return mod_reduce_matprod_launch<BACKEND, 14U, COMPLEX>(stream, C, m, n, ldc);
        case 15U: return mod_reduce_matprod_launch<BACKEND, 15U, COMPLEX>(stream, C, m, n, ldc);
        case 16U: return mod_reduce_matprod_launch<BACKEND, 16U, COMPLEX>(stream, C, m, n, ldc);
        case 17U: return mod_reduce_matprod_launch<BACKEND, 17U, COMPLEX>(stream, C, m, n, ldc);
        case 18U: return mod_reduce_matprod_launch<BACKEND, 18U, COMPLEX>(stream, C, m, n, ldc);
        case 19U: return mod_reduce_matprod_launch<BACKEND, 19U, COMPLEX>(stream, C, m, n, ldc);
        }
    }
}

#undef GEMMUL8_FP8_REDUCE_CASE

template <Backend BACKEND, unsigned IDX, bool COMPLEX>
__global__ void mod_reduce_matprod_strided_kernel(
    common::hi_t<BACKEND> *__restrict__ C,
    const int64_t m,
    const int64_t n,
    const int64_t ldc,
    const int64_t strideC,
    const int batchCount //
) {
    const int64_t i = int64_t(blockIdx.x * blockDim.x + threadIdx.x);
    const int64_t j = int64_t(blockIdx.y * blockDim.y + threadIdx.y);
    const int b     = int(blockIdx.z);

    if (i >= m || j >= n || b >= batchCount) return;

    common::hi_t<BACKEND> *ptr = C + b * strideC + j * ldc + i;

    *ptr = mod_reduce_matprod_core<BACKEND, IDX, COMPLEX>(*ptr);
}

template <Backend BACKEND, unsigned IDX, bool COMPLEX>
inline void mod_reduce_matprod_strided_launch(
    const cudaStream_t stream,
    common::hi_t<BACKEND> *C,
    const int m,
    const int n,
    const size_t ldc,
    const int64_t strideC,
    const int batchCount //
) {
    if (m == 0 || n == 0 || batchCount == 0) return;

    if constexpr (BACKEND == Backend::FP8) {
        if ((m % 2) == 0 &&
            (ldc % 2U) == 0U &&
            strideC >= 0 &&
            (strideC % 2) == 0 &&
            (reinterpret_cast<uintptr_t>(C) % alignof(float2)) == 0U) {
            constexpr unsigned p = common::table::moduli<BACKEND, IDX, COMPLEX>;
            if (!mod_reduce_matprod_fp8_is_flat(m, ldc)) {
                return mod_reduce_matprod_fp8_tile_launch<p, true>(
                    stream, C, m, n, ldc, size_t(strideC), unsigned(batchCount));
            }

            const size_t pairs = ldc * n / 2U;

            if (size_t(strideC) == 2U * pairs) {
                return mod_reduce_matprod_fp8_flat_launch<p, false>(
                    stream, C, pairs * batchCount, 0U, 1U);
            }

            return mod_reduce_matprod_fp8_flat_launch<p, true>(
                stream, C, pairs, size_t(strideC) / 2U, unsigned(batchCount));
        }
    }

    constexpr int TX = 32;
    constexpr int TY = 8;

    constexpr int max_columns = 65535 * TY;
    if (n > max_columns) {
        for (int first = 0; first < n;) {
            const int count = std::min(max_columns, n - first);
            mod_reduce_matprod_strided_launch<BACKEND, IDX, COMPLEX>(
                stream, C + first * ldc, m, count, ldc, strideC, batchCount);
            first += count;
        }
        return;
    }

    const dim3 block(TX, TY);
    const dim3 grid((m + TX - 1) / TX, (n + TY - 1) / TY, batchCount);

    mod_reduce_matprod_strided_kernel<BACKEND, IDX, COMPLEX><<<grid, block, 0, stream>>>(
        C, int64_t(m), int64_t(n), int64_t(ldc),
        strideC, batchCount);
}

#define GEMMUL8_FP8_REDUCE_CASE(P)                                           \
    case P + 20U: {                                                          \
        return mod_reduce_matprod_strided_launch<BACKEND, P + 20U, COMPLEX>( \
            stream, C, m, n, ldc, strideC, batchCount);                      \
    }

template <Backend BACKEND, bool COMPLEX>
void mod_reduce_matprod_strided_dispatch(
    const cudaStream_t stream,
    common::hi_t<BACKEND> *C,
    const int m,
    const int n,
    const size_t ldc,
    const int64_t strideC,
    const int batchCount,
    const unsigned modulus_idx //
) {
    if constexpr (BACKEND == Backend::FP8) {
        switch (modulus_idx) {
            GEMMUL8_FP8_FOR_EACH_MODULUS(GEMMUL8_FP8_REDUCE_CASE)
            GEMMUL8_FP8_FOR_EACH_REDUCTION_PERIOD(GEMMUL8_FP8_REDUCE_CASE)
        }
    } else {
        switch (modulus_idx) {
        case 0U: return mod_reduce_matprod_strided_launch<BACKEND, 0U, COMPLEX>(stream, C, m, n, ldc, strideC, batchCount);
        case 1U: return mod_reduce_matprod_strided_launch<BACKEND, 1U, COMPLEX>(stream, C, m, n, ldc, strideC, batchCount);
        case 2U: return mod_reduce_matprod_strided_launch<BACKEND, 2U, COMPLEX>(stream, C, m, n, ldc, strideC, batchCount);
        case 3U: return mod_reduce_matprod_strided_launch<BACKEND, 3U, COMPLEX>(stream, C, m, n, ldc, strideC, batchCount);
        case 4U: return mod_reduce_matprod_strided_launch<BACKEND, 4U, COMPLEX>(stream, C, m, n, ldc, strideC, batchCount);
        case 5U: return mod_reduce_matprod_strided_launch<BACKEND, 5U, COMPLEX>(stream, C, m, n, ldc, strideC, batchCount);
        case 6U: return mod_reduce_matprod_strided_launch<BACKEND, 6U, COMPLEX>(stream, C, m, n, ldc, strideC, batchCount);
        case 7U: return mod_reduce_matprod_strided_launch<BACKEND, 7U, COMPLEX>(stream, C, m, n, ldc, strideC, batchCount);
        case 8U: return mod_reduce_matprod_strided_launch<BACKEND, 8U, COMPLEX>(stream, C, m, n, ldc, strideC, batchCount);
        case 9U: return mod_reduce_matprod_strided_launch<BACKEND, 9U, COMPLEX>(stream, C, m, n, ldc, strideC, batchCount);
        case 10U: return mod_reduce_matprod_strided_launch<BACKEND, 10U, COMPLEX>(stream, C, m, n, ldc, strideC, batchCount);
        case 11U: return mod_reduce_matprod_strided_launch<BACKEND, 11U, COMPLEX>(stream, C, m, n, ldc, strideC, batchCount);
        case 12U: return mod_reduce_matprod_strided_launch<BACKEND, 12U, COMPLEX>(stream, C, m, n, ldc, strideC, batchCount);
        case 13U: return mod_reduce_matprod_strided_launch<BACKEND, 13U, COMPLEX>(stream, C, m, n, ldc, strideC, batchCount);
        case 14U: return mod_reduce_matprod_strided_launch<BACKEND, 14U, COMPLEX>(stream, C, m, n, ldc, strideC, batchCount);
        case 15U: return mod_reduce_matprod_strided_launch<BACKEND, 15U, COMPLEX>(stream, C, m, n, ldc, strideC, batchCount);
        case 16U: return mod_reduce_matprod_strided_launch<BACKEND, 16U, COMPLEX>(stream, C, m, n, ldc, strideC, batchCount);
        case 17U: return mod_reduce_matprod_strided_launch<BACKEND, 17U, COMPLEX>(stream, C, m, n, ldc, strideC, batchCount);
        case 18U: return mod_reduce_matprod_strided_launch<BACKEND, 18U, COMPLEX>(stream, C, m, n, ldc, strideC, batchCount);
        case 19U: return mod_reduce_matprod_strided_launch<BACKEND, 19U, COMPLEX>(stream, C, m, n, ldc, strideC, batchCount);
        }
    }
}

#undef GEMMUL8_FP8_REDUCE_CASE

template <Backend BACKEND, unsigned IDX, bool COMPLEX>
__global__ void mod_reduce_matprod_pointer_and_advance_kernel(
    void **Aarray,
    void **Barray,
    void **Carray,
    const int64_t m,
    const int64_t n,
    const int64_t ldc,
    const int batchCount,
    const int k_advance //
) {
    const int64_t i = int64_t(blockIdx.x * blockDim.x + threadIdx.x);
    const int64_t j = int64_t(blockIdx.y * blockDim.y + threadIdx.y);
    const int b     = int(blockIdx.z);

    if (b >= batchCount) return;

    if (blockIdx.x == 0 && blockIdx.y == 0 && threadIdx.x == 0 && threadIdx.y == 0) {
        using LowT = common::low_t<BACKEND>;
        Aarray[b]  = reinterpret_cast<LowT *>(Aarray[b]) + k_advance;
        Barray[b]  = reinterpret_cast<LowT *>(Barray[b]) + k_advance;
    }

    if (i >= m || j >= n) return;

    common::hi_t<BACKEND> *C   = reinterpret_cast<common::hi_t<BACKEND> *>(Carray[b]);
    common::hi_t<BACKEND> *ptr = C + j * ldc + i;

    *ptr = mod_reduce_matprod_core<BACKEND, IDX, COMPLEX>(*ptr);
}

template <Backend BACKEND, unsigned IDX, bool COMPLEX>
__global__ void mod_reduce_matprod_pointer_and_advance_large_n_kernel(
    void **Aarray, void **Barray, void **Carray,
    const int64_t m, const int64_t n, const int64_t ldc,
    const int batchCount, const int k_advance //
) {
    const int64_t i  = int64_t(blockIdx.x * blockDim.x + threadIdx.x);
    const int64_t j0 = int64_t(blockIdx.y * blockDim.y + threadIdx.y);
    const int b      = int(blockIdx.z);
    if (b >= batchCount) return;

    // Advance once per batch, not once per column iteration.
    if (blockIdx.x == 0 && blockIdx.y == 0 && threadIdx.x == 0 && threadIdx.y == 0) {
        using LowT = common::low_t<BACKEND>;
        Aarray[b]  = reinterpret_cast<LowT *>(Aarray[b]) + k_advance;
        Barray[b]  = reinterpret_cast<LowT *>(Barray[b]) + k_advance;
    }
    if (i >= m) return;
    auto *C = reinterpret_cast<common::hi_t<BACKEND> *>(Carray[b]);
    for (int64_t j = j0; j < n; j += int64_t(gridDim.y) * blockDim.y) {
        auto *ptr = C + j * ldc + i;
        *ptr      = mod_reduce_matprod_core<BACKEND, IDX, COMPLEX>(*ptr);
    }
}

template <Backend BACKEND, unsigned IDX, bool COMPLEX>
inline void mod_reduce_matprod_pointer_and_advance_launch(
    const cudaStream_t stream,
    void **Aarray,
    void **Barray,
    void **Carray,
    const int m,
    const int n,
    const size_t ldc,
    const int batchCount,
    const int k_advance //
) {
    if (m == 0 || n == 0 || batchCount == 0) return;

    if constexpr (BACKEND == Backend::FP8) {
        if (mod_reduce_matprod_fp8_is_flat(m, ldc)) {
            constexpr unsigned p = common::table::moduli<BACKEND, IDX, COMPLEX>;
            const size_t pairs   = ldc * n / 2U;
            const dim3 grid(unsigned(pairs / FP8_REDUCE_THREADS), unsigned(batchCount));
            if (pairs <= size_t(UINT32_MAX)) {
                mod_reduce_matprod_fp8_pointer_kernel<p, uint32_t>
                    <<<grid, FP8_REDUCE_THREADS, 0, stream>>>(
                        Aarray, Barray, Carray, unsigned(k_advance));
            } else {
                mod_reduce_matprod_fp8_pointer_kernel<p, size_t>
                    <<<grid, FP8_REDUCE_THREADS, 0, stream>>>(
                        Aarray, Barray, Carray, unsigned(k_advance));
            }
            return;
        }
    }

    constexpr int TX = 32;
    constexpr int TY = 8;

    if (n > 65535 * TY) {
        const dim3 grid((m + TX - 1) / TX, 65535U, batchCount);
        mod_reduce_matprod_pointer_and_advance_large_n_kernel<BACKEND, IDX, COMPLEX>
            <<<grid, dim3(TX, TY), 0, stream>>>(
                Aarray, Barray, Carray, int64_t(m), int64_t(n), int64_t(ldc), batchCount, k_advance);
        return;
    }

    const dim3 block(TX, TY);
    const dim3 grid((m + TX - 1) / TX, (n + TY - 1) / TY, batchCount);

    mod_reduce_matprod_pointer_and_advance_kernel<BACKEND, IDX, COMPLEX><<<
        grid, block, 0, stream>>>(
        Aarray, Barray, Carray,
        int64_t(m), int64_t(n), int64_t(ldc),
        batchCount, k_advance);
}

#define GEMMUL8_FP8_REDUCE_CASE(P)                                                       \
    case P + 20U: {                                                                      \
        return mod_reduce_matprod_pointer_and_advance_launch<BACKEND, P + 20U, COMPLEX>( \
            stream, Aarray, Barray, Carray, m, n, ldc, batchCount, k_advance);           \
    }

template <Backend BACKEND, bool COMPLEX>
void mod_reduce_matprod_pointer_and_advance_dispatch(
    const cudaStream_t stream,
    void **Aarray,
    void **Barray,
    void **Carray,
    const int m,
    const int n,
    const size_t ldc,
    const int batchCount,
    const int k_advance,
    const unsigned modulus_idx //
) {
    if constexpr (BACKEND == Backend::FP8) {
        switch (modulus_idx) {
            GEMMUL8_FP8_FOR_EACH_MODULUS(GEMMUL8_FP8_REDUCE_CASE)
            GEMMUL8_FP8_FOR_EACH_REDUCTION_PERIOD(GEMMUL8_FP8_REDUCE_CASE)
        }
    } else {
        switch (modulus_idx) {
        case 0U: return mod_reduce_matprod_pointer_and_advance_launch<BACKEND, 0U, COMPLEX>(stream, Aarray, Barray, Carray, m, n, ldc, batchCount, k_advance);
        case 1U: return mod_reduce_matprod_pointer_and_advance_launch<BACKEND, 1U, COMPLEX>(stream, Aarray, Barray, Carray, m, n, ldc, batchCount, k_advance);
        case 2U: return mod_reduce_matprod_pointer_and_advance_launch<BACKEND, 2U, COMPLEX>(stream, Aarray, Barray, Carray, m, n, ldc, batchCount, k_advance);
        case 3U: return mod_reduce_matprod_pointer_and_advance_launch<BACKEND, 3U, COMPLEX>(stream, Aarray, Barray, Carray, m, n, ldc, batchCount, k_advance);
        case 4U: return mod_reduce_matprod_pointer_and_advance_launch<BACKEND, 4U, COMPLEX>(stream, Aarray, Barray, Carray, m, n, ldc, batchCount, k_advance);
        case 5U: return mod_reduce_matprod_pointer_and_advance_launch<BACKEND, 5U, COMPLEX>(stream, Aarray, Barray, Carray, m, n, ldc, batchCount, k_advance);
        case 6U: return mod_reduce_matprod_pointer_and_advance_launch<BACKEND, 6U, COMPLEX>(stream, Aarray, Barray, Carray, m, n, ldc, batchCount, k_advance);
        case 7U: return mod_reduce_matprod_pointer_and_advance_launch<BACKEND, 7U, COMPLEX>(stream, Aarray, Barray, Carray, m, n, ldc, batchCount, k_advance);
        case 8U: return mod_reduce_matprod_pointer_and_advance_launch<BACKEND, 8U, COMPLEX>(stream, Aarray, Barray, Carray, m, n, ldc, batchCount, k_advance);
        case 9U: return mod_reduce_matprod_pointer_and_advance_launch<BACKEND, 9U, COMPLEX>(stream, Aarray, Barray, Carray, m, n, ldc, batchCount, k_advance);
        case 10U: return mod_reduce_matprod_pointer_and_advance_launch<BACKEND, 10U, COMPLEX>(stream, Aarray, Barray, Carray, m, n, ldc, batchCount, k_advance);
        case 11U: return mod_reduce_matprod_pointer_and_advance_launch<BACKEND, 11U, COMPLEX>(stream, Aarray, Barray, Carray, m, n, ldc, batchCount, k_advance);
        case 12U: return mod_reduce_matprod_pointer_and_advance_launch<BACKEND, 12U, COMPLEX>(stream, Aarray, Barray, Carray, m, n, ldc, batchCount, k_advance);
        case 13U: return mod_reduce_matprod_pointer_and_advance_launch<BACKEND, 13U, COMPLEX>(stream, Aarray, Barray, Carray, m, n, ldc, batchCount, k_advance);
        case 14U: return mod_reduce_matprod_pointer_and_advance_launch<BACKEND, 14U, COMPLEX>(stream, Aarray, Barray, Carray, m, n, ldc, batchCount, k_advance);
        case 15U: return mod_reduce_matprod_pointer_and_advance_launch<BACKEND, 15U, COMPLEX>(stream, Aarray, Barray, Carray, m, n, ldc, batchCount, k_advance);
        case 16U: return mod_reduce_matprod_pointer_and_advance_launch<BACKEND, 16U, COMPLEX>(stream, Aarray, Barray, Carray, m, n, ldc, batchCount, k_advance);
        case 17U: return mod_reduce_matprod_pointer_and_advance_launch<BACKEND, 17U, COMPLEX>(stream, Aarray, Barray, Carray, m, n, ldc, batchCount, k_advance);
        case 18U: return mod_reduce_matprod_pointer_and_advance_launch<BACKEND, 18U, COMPLEX>(stream, Aarray, Barray, Carray, m, n, ldc, batchCount, k_advance);
        case 19U: return mod_reduce_matprod_pointer_and_advance_launch<BACKEND, 19U, COMPLEX>(stream, Aarray, Barray, Carray, m, n, ldc, batchCount, k_advance);
        }
    }
}

#undef GEMMUL8_FP8_REDUCE_CASE

template <Backend BACKEND>
void mod_reduce_matprod(
    const cudaStream_t stream,
    common::hi_t<BACKEND> *C,
    const int m,
    const int n,
    const size_t ldc,
    const unsigned modulus_idx,
    const bool complex_moduli //
) {
    if constexpr (BACKEND == Backend::FP8) {
        mod_reduce_matprod_dispatch<BACKEND, false>(
            stream, C, m, n, ldc, modulus_idx);
    } else {
        if (complex_moduli) {
            mod_reduce_matprod_dispatch<BACKEND, true>(
                stream, C, m, n, ldc, modulus_idx);
        } else {
            mod_reduce_matprod_dispatch<BACKEND, false>(
                stream, C, m, n, ldc, modulus_idx);
        }
    }
}

template <Backend BACKEND>
void mod_reduce_matprod_strided(
    const cudaStream_t stream,
    common::hi_t<BACKEND> *C,
    const int m,
    const int n,
    const size_t ldc,
    const int64_t strideC,
    const int batchCount,
    const unsigned modulus_idx,
    const bool complex_moduli //
) {
    if constexpr (BACKEND == Backend::FP8) {
        mod_reduce_matprod_strided_dispatch<BACKEND, false>(
            stream, C, m, n, ldc, strideC, batchCount, modulus_idx);
    } else {
        if (complex_moduli) {
            mod_reduce_matprod_strided_dispatch<BACKEND, true>(
                stream, C, m, n, ldc, strideC, batchCount, modulus_idx);
        } else {
            mod_reduce_matprod_strided_dispatch<BACKEND, false>(
                stream, C, m, n, ldc, strideC, batchCount, modulus_idx);
        }
    }
}

template <Backend BACKEND>
void mod_reduce_matprod_pointer_and_advance(
    const cudaStream_t stream,
    void **Aarray,
    void **Barray,
    void **Carray,
    const int m,
    const int n,
    const size_t ldc,
    const int batchCount,
    const int k_advance,
    const unsigned modulus_idx,
    const bool complex_moduli //
) {
    if constexpr (BACKEND == Backend::FP8) {
        mod_reduce_matprod_pointer_and_advance_dispatch<BACKEND, false>(
            stream, Aarray, Barray, Carray, m, n, ldc, batchCount, k_advance, modulus_idx);
    } else {
        if (complex_moduli) {
            mod_reduce_matprod_pointer_and_advance_dispatch<BACKEND, true>(
                stream, Aarray, Barray, Carray, m, n, ldc, batchCount, k_advance, modulus_idx);
        } else {
            mod_reduce_matprod_pointer_and_advance_dispatch<BACKEND, false>(
                stream, Aarray, Barray, Carray, m, n, ldc, batchCount, k_advance, modulus_idx);
        }
    }
}

template <unsigned P>
inline void mod_reduce_matprod_fp8_batch_launch(
    cudaStream_t stream,
    float *C,
    int m, int n,
    size_t ldc,
    int64_t strideC,
    unsigned batches,
    unsigned mask //
) {
    constexpr auto s        = common::fp8_plan::get(P);
    constexpr unsigned full = (1U << s.products) - 1U;

    if constexpr (s.q[0] == P &&
                  (s.products < 2U || s.q[1] == P) &&
                  (s.products < 3U || s.q[2] == P)) {
        if ((mask & full) == full) {
            return mod_reduce_matprod_strided_launch<Backend::FP8, P + 20U, false>(
                stream, C, m, n, ldc, strideC, int(batches));
        }
    }

    const auto selected = make_fp8_reduction_batch<P>(batches, mask);
    if (selected.count == 0) return;

    if (mod_reduce_matprod_fp8_is_flat(m, ldc) &&
        strideC >= 0 && (strideC % 2) == 0 &&
        (reinterpret_cast<uintptr_t>(C) % alignof(float2)) == 0U) {
        return mod_reduce_matprod_fp8_batch_flat_launch<false>(
            stream, C, ldc * n / 2U, size_t(strideC) / 2U,
            nullptr, nullptr, nullptr, batches, 0U, selected);
    }

    for (unsigned i = 0; i < selected.count; ++i) {
        const auto plane = selected.planes[i];
        mod_reduce_matprod<Backend::FP8>(
            stream,
            C + int64_t(plane.index) * strideC,
            m, n, ldc,
            common::fp8_plan::modulus_index(unsigned(plane.p)),
            false);
    }
}

template <unsigned P>
inline void mod_reduce_matprod_fp8_pointer_batch_launch(
    cudaStream_t stream,
    void **Aarray,
    void **Barray,
    void **Carray,
    int m, int n,
    size_t ldc,
    unsigned batches,
    unsigned k_advance,
    unsigned mask //
) {
    constexpr auto s        = common::fp8_plan::get(P);
    constexpr unsigned full = (1U << s.products) - 1U;

    if constexpr (s.q[0] == P &&
                  (s.products < 2U || s.q[1] == P) &&
                  (s.products < 3U || s.q[2] == P)) {
        if ((mask & full) == full) {
            return mod_reduce_matprod_pointer_and_advance_launch<Backend::FP8, P + 20U, false>(
                stream, Aarray, Barray, Carray, m, n, ldc, int(batches), int(k_advance));
        }
    }

    const auto selected = make_fp8_reduction_batch<P>(batches, mask);
    if (selected.count == 0) {
        mod_reduce_matprod_fp8_advance_kernel<32><<<1, 32, 0, stream>>>(
            Aarray, Barray, batches, k_advance);
    } else if (mod_reduce_matprod_fp8_is_flat(m, ldc)) {
        mod_reduce_matprod_fp8_batch_flat_launch<true>(
            stream, nullptr, ldc * n / 2U, 0U,
            Aarray, Barray, Carray, batches, k_advance, selected);
    } else {
        const dim3 grid((unsigned(m) + 31U) / 32U, (unsigned(n) + 7U) / 8U, selected.count);
        mod_reduce_matprod_fp8_pointer_tile_kernel<32, 8><<<grid, dim3(32, 8), 0, stream>>>(
            Aarray, Barray, Carray,
            unsigned(m), unsigned(n),
            ldc, batches, k_advance, selected);
    }
}

#define GEMMUL8_FP8_REDUCE_CASE(P)                                      \
    case P: {                                                           \
        return mod_reduce_matprod_fp8_batch_launch<P>(                  \
            stream, C, m, n, ldc, strideC, unsigned(batchCount), mask); \
    }

template <Backend BACKEND>
void mod_reduce_matprod_fp8_batch(
    cudaStream_t stream,
    float *C,
    int m, int n,
    size_t ldc,
    int64_t strideC,
    int batchCount,
    unsigned p,
    unsigned mask //
) {
    if constexpr (BACKEND == Backend::FP8) {
        if (m <= 0 || n <= 0 || batchCount <= 0 || mask == 0) return;
        assert(batchCount <= 9);
        switch (p) {
            GEMMUL8_FP8_FOR_EACH_MODULUS(GEMMUL8_FP8_REDUCE_CASE)
        }
    }
}

#undef GEMMUL8_FP8_REDUCE_CASE

#define GEMMUL8_FP8_REDUCE_CASE(P)                                                                       \
    case P: {                                                                                            \
        return mod_reduce_matprod_fp8_pointer_batch_launch<P>(                                           \
            stream, Aarray, Barray, Carray, m, n, ldc, unsigned(batchCount), unsigned(k_advance), mask); \
    }

template <Backend BACKEND>
void mod_reduce_matprod_fp8_pointer_batch(
    cudaStream_t stream,
    void **Aarray,
    void **Barray,
    void **Carray,
    int m, int n,
    size_t ldc,
    int batchCount,
    int k_advance,
    unsigned p,
    unsigned mask //
) {
    if constexpr (BACKEND == Backend::FP8) {
        if (m <= 0 || n <= 0 || batchCount <= 0) return;
        assert(batchCount <= 9);
        switch (p) {
            GEMMUL8_FP8_FOR_EACH_MODULUS(GEMMUL8_FP8_REDUCE_CASE)
        }
    }
}

#undef GEMMUL8_FP8_REDUCE_CASE

} // namespace gemmul8::mod

#undef GEMMUL8_FP8_FOR_EACH_REDUCTION_PERIOD
