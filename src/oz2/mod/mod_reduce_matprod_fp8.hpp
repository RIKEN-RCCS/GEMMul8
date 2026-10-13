#pragma once
#include "../common/common.hpp"
#include "../common/fp8_plan.hpp"

namespace gemmul8::mod {

template <unsigned P>
struct fp8_matprod_reduction {
    static_assert(P >= 5U && P < (1U << 24));
    static constexpr float neg_inv_p = float(-1.0 / double(P));
    static constexpr float bias      = 0x1.8p23f; // 1.5*2^23

    __device__ __forceinline__ static float apply(float c) {
        const float q = __fadd_rn(__fmaf_rn(c, neg_inv_p, bias), -bias);
        return __fmaf_rn(q, float(P), c); // <= floor(p/2) + 1
    }

#if defined(__CUDA_ARCH__) && __CUDA_ARCH__ >= 1000 && __CUDA_ARCH__ < 1100 && !defined(__HIPCC__)

    __device__ __forceinline__ static float2 apply(float2 c) {
        float2 q = __ffma2_rn(c, {neg_inv_p, neg_inv_p}, {bias, bias});
        q        = __fadd2_rn(q, {-bias, -bias});
        return __ffma2_rn(q, {float(P), float(P)}, c); // <= floor(p/2) + 1
    }

#else

    __device__ __forceinline__ static float2 apply(float2 c) {
        return {apply(c.x), apply(c.y)}; // <= floor(p/2) + 1
    }

#endif
};

__device__ __forceinline__ float fp8_matprod_reduce(float c, float p, float neg_inv_p) {
    constexpr float bias = 0x1.8p23f; // 1.5*2^23
    const float q        = __fadd_rn(__fmaf_rn(c, neg_inv_p, bias), -bias);
    return __fmaf_rn(q, p, c);
}

#if defined(__CUDA_ARCH__) && __CUDA_ARCH__ >= 1000 && __CUDA_ARCH__ < 1100 && !defined(__HIPCC__)

__device__ __forceinline__ float2 fp8_matprod_reduce(float2 c, float p, float neg_inv_p) {
    constexpr float bias = 0x1.8p23f; // 1.5*2^23
    float2 q             = __ffma2_rn(c, {neg_inv_p, neg_inv_p}, {bias, bias});
    q                    = __fadd2_rn(q, {-bias, -bias});
    return __ffma2_rn(q, {p, p}, c);
}

#else

__device__ __forceinline__ float2 fp8_matprod_reduce(float2 c, float p, float neg_inv_p) {
    return {fp8_matprod_reduce(c.x, p, neg_inv_p),
            fp8_matprod_reduce(c.y, p, neg_inv_p)};
}

#endif

inline constexpr unsigned FP8_REDUCE_THREADS = 128U;

template <unsigned P, class Index, bool STRIDED>
__global__ void mod_reduce_matprod_fp8_flat_kernel(
    float2 *__restrict__ C,
    const size_t stride2 //
) {
    const Index i = Index(blockIdx.x) * FP8_REDUCE_THREADS + threadIdx.x;
    if constexpr (STRIDED) C += size_t(blockIdx.y) * stride2;
    C[i] = fp8_matprod_reduction<P>::apply(C[i]);
}

template <unsigned P, bool STRIDED>
__global__ void mod_reduce_matprod_fp8_tile_kernel(
    float2 *__restrict__ C,
    const unsigned m2,
    const unsigned n,
    const size_t ldc2,
    const size_t stride2 //
) {
    const unsigned i = blockIdx.x * 32U + threadIdx.x;
    const unsigned j = blockIdx.y * 8U + threadIdx.y;
    if (i >= m2 || j >= n) return;
    size_t offset = j * ldc2 + i;
    if constexpr (STRIDED) offset += blockIdx.z * stride2;
    C[offset] = fp8_matprod_reduction<P>::apply(C[offset]);
}

template <unsigned P, bool STRIDED>
inline void mod_reduce_matprod_fp8_tile_launch(
    cudaStream_t stream,
    float *C,
    int m, int n,
    size_t ldc,
    size_t stride,
    unsigned batches //
) {
    constexpr int max_columns = 65535 * 8;
    if (n > max_columns) {
        for (int first = 0; first < n;) {
            const int count = std::min(max_columns, n - first);
            mod_reduce_matprod_fp8_tile_launch<P, STRIDED>(
                stream, C + first * ldc, m, count, ldc, stride, batches);
            first += count;
        }
        return;
    }
    const unsigned m2 = unsigned(m) / 2U;
    const dim3 grid((m2 + 31U) / 32U, (unsigned(n) + 7U) / 8U, batches);
    mod_reduce_matprod_fp8_tile_kernel<P, STRIDED><<<grid, dim3(32, 8), 0, stream>>>(
        reinterpret_cast<float2 *>(C), m2, unsigned(n), ldc / 2U, stride / 2U);
}

template <unsigned P, class Index>
__global__ void mod_reduce_matprod_fp8_pointer_kernel(
    void **Aarray,
    void **Barray,
    void **Carray,
    const unsigned k_advance //
) {
    const unsigned b = blockIdx.y;
    const Index i    = Index(blockIdx.x) * FP8_REDUCE_THREADS + threadIdx.x;
    auto *C          = reinterpret_cast<float2 *>(Carray[b]);
    C[i]             = fp8_matprod_reduction<P>::apply(C[i]);
    if (blockIdx.x == 0 && threadIdx.x == 0) {
        using LowT = common::low_t<Backend::FP8>;
        Aarray[b]  = reinterpret_cast<LowT *>(Aarray[b]) + k_advance;
        Barray[b]  = reinterpret_cast<LowT *>(Barray[b]) + k_advance;
    }
}

template <unsigned P, bool STRIDED>
inline void mod_reduce_matprod_fp8_flat_launch(
    cudaStream_t stream,
    float *C,
    size_t pairs,
    size_t stride2,
    unsigned batches //
) {
    const dim3 grid(unsigned(pairs / FP8_REDUCE_THREADS), batches);
    if (pairs <= size_t(UINT32_MAX)) {
        mod_reduce_matprod_fp8_flat_kernel<P, uint32_t, STRIDED>
            <<<grid, FP8_REDUCE_THREADS, 0, stream>>>(
                reinterpret_cast<float2 *>(C), stride2);
    } else {
        mod_reduce_matprod_fp8_flat_kernel<P, size_t, STRIDED>
            <<<grid, FP8_REDUCE_THREADS, 0, stream>>>(
                reinterpret_cast<float2 *>(C), stride2);
    }
}

inline bool mod_reduce_matprod_fp8_is_flat(int m, size_t ldc) {
    return size_t(m) == ldc && (ldc % (2U * FP8_REDUCE_THREADS)) == 0U;
}

struct fp8_reduction_batch {
    struct plane {
        unsigned index;
        float p, neg_inv_p;
    } planes[9];
    unsigned count;
};

template <unsigned P>
inline fp8_reduction_batch make_fp8_reduction_batch(unsigned batches, unsigned mask) {
    constexpr auto s = common::fp8_plan::get(P);

    constexpr fp8_reduction_batch::plane periods[] = {
        {0U, float(s.q[0]), float(-1.0 / double(s.q[0]))},
        {0U, float(s.q[1]), float(-1.0 / double(s.q[1]))},
        {0U, float(s.q[2]), float(-1.0 / double(s.q[2]))}
    };

    fp8_reduction_batch result{};
    for (unsigned b = 0; b < batches; ++b) {
        const unsigned j = b % s.products;
        if (mask & (1U << j)) {
            auto plane                    = periods[j];
            plane.index                   = b;
            result.planes[result.count++] = plane;
        }
    }
    return result;
}

__device__ __forceinline__ void fp8_advance_pointer(
    void **Aarray, void **Barray,
    unsigned b, unsigned k_advance //
) {
    using LowT = common::low_t<Backend::FP8>;
    Aarray[b]  = reinterpret_cast<LowT *>(Aarray[b]) + k_advance;
    Barray[b]  = reinterpret_cast<LowT *>(Barray[b]) + k_advance;
}

#if defined(__CUDACC__) && !defined(__HIPCC__)
    #define GEMMUL8_FP8_GRID_CONSTANT const __grid_constant__
#else
    #define GEMMUL8_FP8_GRID_CONSTANT
#endif

template <bool POINTER, class Index>
__global__ void mod_reduce_matprod_fp8_batch_kernel(
    float2 *__restrict__ C,
    size_t stride2,
    void **Aarray,
    void **Barray,
    void **Carray,
    unsigned batches,
    unsigned k_advance,
    GEMMUL8_FP8_GRID_CONSTANT fp8_reduction_batch selected //
) {
    if constexpr (POINTER) {
        if (blockIdx.x == 0 && blockIdx.y == 0 && threadIdx.x < batches) {
            fp8_advance_pointer(Aarray, Barray, threadIdx.x, k_advance);
        }
    }

    const auto plane = selected.planes[blockIdx.y];
    const Index i    = Index(blockIdx.x) * FP8_REDUCE_THREADS + threadIdx.x;

    if constexpr (POINTER) {
        C = reinterpret_cast<float2 *>(Carray[plane.index]);
    } else {
        C += size_t(plane.index) * stride2;
    }
    C[i] = fp8_matprod_reduce(C[i], plane.p, plane.neg_inv_p);
}

template <unsigned TX, unsigned TY>
__global__ void mod_reduce_matprod_fp8_pointer_tile_kernel(
    void **Aarray,
    void **Barray,
    void **Carray,
    unsigned m, unsigned n,
    size_t ldc,
    unsigned batches,
    unsigned k_advance,
    GEMMUL8_FP8_GRID_CONSTANT fp8_reduction_batch selected //
) {
    if (blockIdx.x == 0 && blockIdx.y == 0 && blockIdx.z == 0 &&
        threadIdx.y == 0 && threadIdx.x < batches) {
        fp8_advance_pointer(Aarray, Barray, threadIdx.x, k_advance);
    }

    const unsigned i = blockIdx.x * TX + threadIdx.x;
    const unsigned j = blockIdx.y * TY + threadIdx.y;
    if (i >= m || j >= n) return;

    const auto plane = selected.planes[blockIdx.z];
    float *C         = static_cast<float *>(Carray[plane.index]) + size_t(j) * ldc + i;
    *C               = fp8_matprod_reduce(*C, plane.p, plane.neg_inv_p);
}

#undef GEMMUL8_FP8_GRID_CONSTANT

template <unsigned THREADS>
__global__ void mod_reduce_matprod_fp8_advance_kernel(
    void **Aarray,
    void **Barray,
    unsigned batches,
    unsigned k_advance //
) {
    const unsigned b = blockIdx.x * THREADS + threadIdx.x;
    if (b < batches) fp8_advance_pointer(Aarray, Barray, b, k_advance);
}

template <bool POINTER>
inline void mod_reduce_matprod_fp8_batch_flat_launch(
    cudaStream_t stream,
    float *C,
    size_t pairs, size_t stride2,
    void **Aarray,
    void **Barray,
    void **Carray,
    unsigned batches,
    unsigned k_advance,
    fp8_reduction_batch selected //
) {
    const dim3 grid(unsigned(pairs / FP8_REDUCE_THREADS), selected.count);
    if (pairs <= size_t(UINT32_MAX)) {
        mod_reduce_matprod_fp8_batch_kernel<POINTER, uint32_t>
            <<<grid, FP8_REDUCE_THREADS, 0, stream>>>(
                reinterpret_cast<float2 *>(C), stride2,
                Aarray, Barray, Carray, batches, k_advance, selected);
    } else {
        mod_reduce_matprod_fp8_batch_kernel<POINTER, size_t>
            <<<grid, FP8_REDUCE_THREADS, 0, stream>>>(
                reinterpret_cast<float2 *>(C), stride2,
                Aarray, Barray, Carray, batches, k_advance, selected);
    }
}

} // namespace gemmul8::mod
