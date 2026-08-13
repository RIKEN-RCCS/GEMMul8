#pragma once
#include "../common/common.hpp"

namespace gemmul8::common {

__device__ __forceinline__ float mod_f32_core(
    const float c,
    const float p,
    const float neg_inv_p,
    const float half_p //
) {
    float r = fmaf(rintf(c * neg_inv_p), p, c);
    return (r > half_p) ? (r - p) : ((r < -half_p) ? (r + p) : r);
}

static __global__ void mod_f32_kernel(
    float *__restrict__ C,
    const int64_t m,
    const int64_t n,
    const int64_t ldc,
    const float p,
    const float neg_inv_p,
    const float half_p //
) {
    const int64_t i = int64_t(blockIdx.x * blockDim.x + threadIdx.x);
    const int64_t j = int64_t(blockIdx.y * blockDim.y + threadIdx.y);

    if (i >= m || j >= n) {
        return;
    }

    float *__restrict__ ptr = C + j * ldc + i;
    *ptr                    = mod_f32_core(*ptr, p, neg_inv_p, half_p);
}

inline void mod_f32(
    const cudaStream_t stream,
    float *C,
    const int m,
    const int n,
    const size_t ldc,
    const Handle_t &h //
) {
    if (m == 0 || n == 0) {
        return;
    }

    constexpr int TX = 32;
    constexpr int TY = 8;

    const dim3 block(TX, TY);
    const dim3 grid((m + TX - 1) / TX, (n + TY - 1) / TY);

    mod_f32_kernel<<<grid, block, 0, stream>>>(
        C, int64_t(m), int64_t(n), int64_t(ldc),
        h.fp8_modulus,
        h.fp8_neg_inv_modulus,
        h.fp8_half_modulus);
}

static __global__ void mod_f32_strided_kernel(
    float *__restrict__ C,
    const int64_t m,
    const int64_t n,
    const int64_t ldc,
    const int64_t strideC,
    const int batchCount,
    const float p,
    const float neg_inv_p,
    const float half_p //
) {
    const int64_t i = int64_t(blockIdx.x * blockDim.x + threadIdx.x);
    const int64_t j = int64_t(blockIdx.y * blockDim.y + threadIdx.y);
    const int b     = int(blockIdx.z);

    if (i >= m || j >= n || b >= batchCount) {
        return;
    }

    float *__restrict__ ptr = C + b * strideC + j * ldc + i;
    *ptr                    = mod_f32_core(*ptr, p, neg_inv_p, half_p);
}

inline void mod_f32_strided(
    const cudaStream_t stream,
    float *C,
    const int m,
    const int n,
    const size_t ldc,
    const int64_t strideC,
    const int batchCount,
    const Handle_t &h //
) {
    if (m == 0 || n == 0 || batchCount == 0) {
        return;
    }

    constexpr int TX = 32;
    constexpr int TY = 8;

    const dim3 block(TX, TY);
    const dim3 grid((m + TX - 1) / TX, (n + TY - 1) / TY, batchCount);

    mod_f32_strided_kernel<<<grid, block, 0, stream>>>(
        C, int64_t(m), int64_t(n), int64_t(ldc),
        strideC,
        batchCount,
        h.fp8_modulus,
        h.fp8_neg_inv_modulus,
        h.fp8_half_modulus);
}

static __global__ void mod_f32_pointer_and_advance_kernel(
    void **Aarray,
    void **Barray,
    void **Carray,
    const int64_t m,
    const int64_t n,
    const int64_t ldc,
    const int batchCount,
    const int k_advance,
    const float p,
    const float neg_inv_p,
    const float half_p //
) {
    const int64_t i = int64_t(blockIdx.x * blockDim.x + threadIdx.x);
    const int64_t j = int64_t(blockIdx.y * blockDim.y + threadIdx.y);
    const int b     = int(blockIdx.z);

    if (b >= batchCount) {
        return;
    }

    if (blockIdx.x == 0 && blockIdx.y == 0 && threadIdx.x == 0 && threadIdx.y == 0) {
        using LowT = low_t<Backend::FP8>;
        Aarray[b]  = reinterpret_cast<LowT *>(Aarray[b]) + k_advance;
        Barray[b]  = reinterpret_cast<LowT *>(Barray[b]) + k_advance;
    }

    if (i >= m || j >= n) {
        return;
    }

    float *__restrict__ C   = reinterpret_cast<float *>(Carray[b]);
    float *__restrict__ ptr = C + j * ldc + i;
    *ptr                    = mod_f32_core(*ptr, p, neg_inv_p, half_p);
}

inline void mod_f32_pointer_and_advance(
    const cudaStream_t stream,
    void **Aarray,
    void **Barray,
    void **Carray,
    const int m,
    const int n,
    const size_t ldc,
    const int batchCount,
    const int k_advance,
    const Handle_t &h //
) {
    if (m == 0 || n == 0 || batchCount == 0) {
        return;
    }

    constexpr int TX = 32;
    constexpr int TY = 8;

    const dim3 block(TX, TY);
    const dim3 grid((m + TX - 1) / TX, (n + TY - 1) / TY, batchCount);

    mod_f32_pointer_and_advance_kernel<<<grid, block, 0, stream>>>(
        Aarray, Barray, Carray,
        int64_t(m), int64_t(n), int64_t(ldc),
        batchCount,
        k_advance,
        h.fp8_modulus,
        h.fp8_neg_inv_modulus,
        h.fp8_half_modulus);
}

} // namespace gemmul8::common
