#pragma once
#include "../common/common.hpp"
#include "mod_core.hpp"

namespace gemmul8::mod {

template <Backend BACKEND, unsigned IDX, bool COMPLEX>
__device__ __forceinline__
    common::hi_t<BACKEND>
    mod_reduce_matprod_core(const common::hi_t<BACKEND> c) {
    if constexpr (BACKEND == Backend::INT8) {
        return mod_small<BACKEND, IDX, COMPLEX>(c);
    } else if constexpr (BACKEND == Backend::FP8) {
        constexpr float p         = common::table::moduli<Backend::FP8, IDX, COMPLEX>;
        constexpr float neg_inv_p = float(-1.0 / double(p));
        constexpr float half_p    = 0.5f * p;

        float r = fmaf(rintf(c * neg_inv_p), p, c);
        return (r > half_p) ? (r - p) : ((r < -half_p) ? (r + p) : r);
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

    if (i >= m || j >= n) {
        return;
    }

    common::hi_t<BACKEND> *ptr = C + j * ldc + i;
    *ptr                       = mod_reduce_matprod_core<BACKEND, IDX, COMPLEX>(*ptr);
}

template <Backend BACKEND, unsigned IDX, bool COMPLEX>
inline void mod_reduce_matprod_launch(
    const cudaStream_t stream,
    common::hi_t<BACKEND> *C,
    const int m,
    const int n,
    const size_t ldc //
) {
    if (m == 0 || n == 0) {
        return;
    }

    constexpr int TX = 32;
    constexpr int TY = 8;

    const dim3 block(TX, TY);
    const dim3 grid((m + TX - 1) / TX, (n + TY - 1) / TY);

    mod_reduce_matprod_kernel<BACKEND, IDX, COMPLEX><<<grid, block, 0, stream>>>(
        C, int64_t(m), int64_t(n), int64_t(ldc));
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
    default: break;
    }
}

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

    if (i >= m || j >= n || b >= batchCount) {
        return;
    }

    common::hi_t<BACKEND> *ptr = C + b * strideC + j * ldc + i;
    *ptr                       = mod_reduce_matprod_core<BACKEND, IDX, COMPLEX>(*ptr);
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
    if (m == 0 || n == 0 || batchCount == 0) {
        return;
    }

    constexpr int TX = 32;
    constexpr int TY = 8;

    const dim3 block(TX, TY);
    const dim3 grid((m + TX - 1) / TX, (n + TY - 1) / TY, batchCount);

    mod_reduce_matprod_strided_kernel<BACKEND, IDX, COMPLEX><<<grid, block, 0, stream>>>(
        C, int64_t(m), int64_t(n), int64_t(ldc),
        strideC, batchCount);
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
    default: break;
    }
}

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

    if (b >= batchCount) {
        return;
    }

    if (blockIdx.x == 0 && blockIdx.y == 0 && threadIdx.x == 0 && threadIdx.y == 0) {
        using LowT = common::low_t<BACKEND>;
        Aarray[b]  = reinterpret_cast<LowT *>(Aarray[b]) + k_advance;
        Barray[b]  = reinterpret_cast<LowT *>(Barray[b]) + k_advance;
    }

    if (i >= m || j >= n) {
        return;
    }

    common::hi_t<BACKEND> *C   = reinterpret_cast<common::hi_t<BACKEND> *>(Carray[b]);
    common::hi_t<BACKEND> *ptr = C + j * ldc + i;
    *ptr                       = mod_reduce_matprod_core<BACKEND, IDX, COMPLEX>(*ptr);
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
    if (m == 0 || n == 0 || batchCount == 0) {
        return;
    }

    constexpr int TX = 32;
    constexpr int TY = 8;

    const dim3 block(TX, TY);
    const dim3 grid((m + TX - 1) / TX, (n + TY - 1) / TY, batchCount);

    mod_reduce_matprod_pointer_and_advance_kernel<BACKEND, IDX, COMPLEX><<<grid, block, 0, stream>>>(
        Aarray, Barray, Carray,
        int64_t(m), int64_t(n), int64_t(ldc),
        batchCount, k_advance);
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
    default: break;
    }
}

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
    if (complex_moduli) {
        mod_reduce_matprod_dispatch<BACKEND, true>(stream, C, m, n, ldc, modulus_idx);
    } else {
        mod_reduce_matprod_dispatch<BACKEND, false>(stream, C, m, n, ldc, modulus_idx);
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
    if (complex_moduli) {
        mod_reduce_matprod_strided_dispatch<BACKEND, true>(stream, C, m, n, ldc, strideC, batchCount, modulus_idx);
    } else {
        mod_reduce_matprod_strided_dispatch<BACKEND, false>(stream, C, m, n, ldc, strideC, batchCount, modulus_idx);
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
    if (complex_moduli) {
        mod_reduce_matprod_pointer_and_advance_dispatch<BACKEND, true>(stream, Aarray, Barray, Carray, m, n, ldc, batchCount, k_advance, modulus_idx);
    } else {
        mod_reduce_matprod_pointer_and_advance_dispatch<BACKEND, false>(stream, Aarray, Barray, Carray, m, n, ldc, batchCount, k_advance, modulus_idx);
    }
}

} // namespace gemmul8::mod
