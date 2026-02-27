#pragma once
#include "self_hipify.hpp"

namespace makemat {

template <typename Tin, typename Tout> __forceinline__ __device__ Tout casting(Tin in);

// --- Tin = double ---
template <> __forceinline__ __device__ float casting<double, float>(double in) { return __double2float_rn(in); }
template <> __forceinline__ __device__ double casting<double, double>(double in) { return in; }
template <> __forceinline__ __device__ float2 casting<double, float2>(double in) { return float2{__double2float_rn(in), 0.0f}; }
template <> __forceinline__ __device__ double2 casting<double, double2>(double in) { return double2{in, 0.0}; }

// --- Tin = float ---
template <> __forceinline__ __device__ float casting<float, float>(float in) { return in; }
template <> __forceinline__ __device__ double casting<float, double>(float in) { return static_cast<double>(in); }
template <> __forceinline__ __device__ float2 casting<float, float2>(float in) { return float2{in, 0.0f}; }
template <> __forceinline__ __device__ double2 casting<float, double2>(float in) { return double2{double(in), 0.0}; }

// --- Tin = float2 ---
template <> __forceinline__ __device__ float casting<float2, float>(float2 in) { return in.x; }
template <> __forceinline__ __device__ double casting<float2, double>(float2 in) { return static_cast<double>(in.x); }
template <> __forceinline__ __device__ float2 casting<float2, float2>(float2 in) { return in; }
template <> __forceinline__ __device__ double2 casting<float2, double2>(float2 in) { return double2{static_cast<double>(in.x), static_cast<double>(in.y)}; }

// --- Tin = double2 ---
template <> __forceinline__ __device__ float casting<double2, float>(double2 in) { return __double2float_rn(in.x); }
template <> __forceinline__ __device__ double casting<double2, double>(double2 in) { return in.x; }
template <> __forceinline__ __device__ float2 casting<double2, float2>(double2 in) { return float2{__double2float_rn(in.x), __double2float_rn(in.y)}; }
template <> __forceinline__ __device__ double2 casting<double2, double2>(double2 in) { return in; }

#pragma clang optimize off
template <typename T>
__global__ void randmat_kernel(size_t m,                      // rows of A
                               size_t n,                      // columns of A
                               T *const A,                    // output
                               double phi,                    // difficulty for matrix multiplication
                               const unsigned long long seed) // seed for random numbers
{
    const size_t idx = size_t(threadIdx.x) + size_t(blockIdx.x) * size_t(blockDim.x);
    if (idx >= m * n) return;
    curandState state;
    curand_init(seed, idx, 0, &state);

    T out;
    if constexpr (std::is_same_v<T, cuDoubleComplex>) {
        const double rand_r  = curand_uniform_double(&state);
        const double rand_i  = curand_uniform_double(&state);
        const double randn_r = curand_normal_double(&state);
        const double randn_i = curand_normal_double(&state);
        if (phi < 0) {
            out.x = randn_r;
            out.y = randn_i;
        } else {
            out.x = (rand_r - 0.5) * exp(randn_r * phi);
            out.y = (rand_i - 0.5) * exp(randn_i * phi);
        }

    } else if constexpr (std::is_same_v<T, cuFloatComplex>) {
        const double rand_r  = curand_uniform_double(&state);
        const double rand_i  = curand_uniform_double(&state);
        const double randn_r = curand_normal_double(&state);
        const double randn_i = curand_normal_double(&state);
        if (phi < 0) {
            out.x = static_cast<float>(randn_r);
            out.y = static_cast<float>(randn_i);
        } else {
            out.x = static_cast<float>((rand_r - 0.5) * exp(randn_r * phi));
            out.y = static_cast<float>((rand_i - 0.5) * exp(randn_i * phi));
        }

    } else {
        const double rand  = curand_uniform_double(&state);
        const double randn = curand_normal_double(&state);
        if (phi < 0) {
            out = static_cast<T>(randn);
        } else {
            out = static_cast<T>((rand - 0.5) * exp(randn * phi));
        }
    }
    A[idx] = out;
}
#pragma clang optimize on

template <typename T>
void randmat(size_t m,                     // rows of A
             size_t n,                     // columns of A
             T *const A,                   // output
             double phi,                   // difficulty for matrix multiplication
             const unsigned long long seed // seed for random numbers
) {
    constexpr size_t block_size = 256;
    const size_t grid_size      = (m * n + block_size - 1) / block_size;
    randmat_kernel<T><<<grid_size, block_size>>>(m, n, A, phi, seed);
    cudaDeviceSynchronize();
}

template <typename Tin, typename Tout>
__global__ void casting_kernel(size_t sizeA,
                               const Tin *const __restrict__ in,
                               Tout *const __restrict__ out //
) {
    const size_t idx = size_t(threadIdx.x) + size_t(blockIdx.x) * size_t(blockDim.x);
    if (idx >= sizeA) return;
    out[idx] = casting<Tin, Tout>(in[idx]);
}

template <typename Tin, typename Tout>
void casting(size_t m,            // rows of A
             size_t n,            // columns of A
             const Tin *const in, // input
             Tout *const out      // output
) {
    constexpr size_t block_size = 256;
    const size_t grid_size      = (m * n + block_size - 1) / block_size;
    casting_kernel<Tin, Tout><<<grid_size, block_size>>>(m * n, in, out);
}

} // namespace makemat
