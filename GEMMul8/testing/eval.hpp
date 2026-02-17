#pragma once
#include "self_hipify.hpp"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>
#include <type_traits>

namespace eval {

//------------------------------
// dd
//------------------------------
namespace dd {

struct double2_complex {
    double2 x, y;
};

__host__ __device__ __forceinline__ double FMA(double a, double b, double c) {
#ifdef __CUDA_ARCH__
    return fma(a, b, c);
#else
    return std::fma(a, b, c);
#endif
}

#pragma clang optimize off
__device__ __host__ __forceinline__ void two_sum(
    const double a, const double b,
    double &c, double &d //
) {
    c        = a + b;
    double s = c - a;
    double t = b - s;
    double u = c - s;
    d        = (a - u) + t;
}

__device__ __host__ __forceinline__ void two_sub(
    const double a, const double b,
    double &c, double &d //
) {
    c        = a - b;
    double s = c - a;
    double t = b + s;
    double u = c - s;
    d        = (a - u) - t;
}

__device__ __host__ __forceinline__ void fast_two_sum(
    const double a, const double b,
    double &c, double &d //
) {
    c = a + b;
    d = (a - c) + b;
}

__device__ __host__ __forceinline__ void two_prod(
    const double a, const double b,
    double &c, double &d //
) {
    c = a * b;
    d = FMA(a, b, -c);
}

__device__ __host__ __forceinline__ double2 add(const double2 a, const double2 b) {
    double2 c;
    two_sum(a.x, b.x, c.x, c.y);
    c.y += a.y;
    c.y += b.y;
    fast_two_sum(c.x, c.y, c.x, c.y);
    return c;
}

__device__ __host__ __forceinline__ double2 add(const double a, const double2 b) {
    double2 c;
    two_sum(a, b.x, c.x, c.y);
    c.y += b.y;
    fast_two_sum(c.x, c.y, c.x, c.y);
    return c;
}

__device__ __host__ __forceinline__ double2 add(const double2 a, const double b) {
    double2 c;
    two_sum(a.x, b, c.x, c.y);
    c.y += a.y;
    fast_two_sum(c.x, c.y, c.x, c.y);
    return c;
}

__device__ __host__ __forceinline__ double2 sub(const double2 a, const double2 b) {
    double2 c;
    two_sub(a.x, b.x, c.x, c.y);
    c.y += a.y;
    c.y -= b.y;
    fast_two_sum(c.x, c.y, c.x, c.y);
    return c;
}

__device__ __host__ __forceinline__ double2 sub(const double a, const double2 b) {
    double2 c;
    two_sub(a, b.x, c.x, c.y);
    c.y -= b.y;
    fast_two_sum(c.x, c.y, c.x, c.y);
    return c;
}

__device__ __host__ __forceinline__ double2 sub(const double a, const double b) {
    double2 c;
    two_sub(a, b, c.x, c.y);
    return c;
}

__device__ __host__ __forceinline__ double2 mul(const double2 a, const double2 b) {
    double2 c;
    two_prod(a.x, b.x, c.x, c.y);
    c.y = FMA(a.y, b.x, FMA(a.x, b.y, c.y));
    fast_two_sum(c.x, c.y, c.x, c.y);
    return c;
}

__device__ __host__ __forceinline__ double2 mul(const double a, const double2 b) {
    double2 c;
    two_prod(a, b.x, c.x, c.y);
    c.y = FMA(a, b.y, c.y);
    fast_two_sum(c.x, c.y, c.x, c.y);
    return c;
}

__device__ __host__ __forceinline__ double2 mul(const double a, const double b) {
    double2 c;
    two_prod(a, b, c.x, c.y);
    return c;
}

__device__ __host__ __forceinline__ double2 div(const double2 a, const double2 b) {
    double2 c;
    double s, t;
    c.x = a.x / b.x;
    two_prod(c.x, b.x, s, t);
    double u = a.x - s;
    u -= t;
    u += a.y;
    u = FMA(-c.x, b.y, u);
    u /= b.x;
    fast_two_sum(c.x, u, c.x, c.y);
    return c;
}

__device__ __host__ __forceinline__ double2 div(const double2 a, const double b) {
    double2 c;
    double s, t;
    c.x = a.x / b;
    two_prod(c.x, b, s, t);
    double u = a.x - s;
    u -= t;
    u += a.y;
    u /= b;
    fast_two_sum(c.x, u, c.x, c.y);
    return c;
}
#pragma clang optimize on

template <typename T, bool trans>
__device__ __forceinline__ T load_A_element(const T *A,
                                            size_t m,
                                            size_t k,
                                            size_t row,
                                            size_t t //
) {
    if constexpr (!trans) {
        return __ldg(A + row + t * m);
    } else {
        return __ldg(A + t + row * k);
    }
}

template <typename T, bool trans>
__device__ __forceinline__ T load_B_element(const T *B,
                                            size_t k,
                                            size_t n,
                                            size_t t,
                                            size_t col //
) {
    if constexpr (!trans) {
        return __ldg(B + t + col * k);
    } else {
        return __ldg(B + col + t * n);
    }
}

template <bool transA, bool transB>
__global__ void simple_gemm_device(size_t m,
                                   size_t n,
                                   size_t k,
                                   const float *__restrict__ A,
                                   const float *__restrict__ B,
                                   double *__restrict__ C //
) {
    const int TILE = 32;
    __shared__ float Asub[TILE][TILE + 1];
    __shared__ float Bsub[TILE][TILE + 1];

    size_t row = blockIdx.y * blockDim.y + threadIdx.y;
    size_t col = blockIdx.x * blockDim.x + threadIdx.x;

    double sum = 0.0;

    int numTiles = (int)((k + TILE - 1) / TILE);

    for (int t = 0; t < numTiles; ++t) {
        size_t a_col = t * TILE + threadIdx.x;
        if (row < m && a_col < k) {
            Asub[threadIdx.y][threadIdx.x] = load_A_element<float, transA>(A, m, k, row, a_col);
        } else {
            Asub[threadIdx.y][threadIdx.x] = 0.0f;
        }

        size_t b_row = t * TILE + threadIdx.y;
        if (col < n && b_row < k) {
            Bsub[threadIdx.y][threadIdx.x] = load_B_element<float, transB>(B, k, n, b_row, col);
        } else {
            Bsub[threadIdx.y][threadIdx.x] = 0.0f;
        }

        __syncthreads();

#pragma unroll
        for (int i = 0; i < TILE; ++i) {
            double a = static_cast<double>(Asub[threadIdx.y][i]);
            double b = static_cast<double>(Bsub[i][threadIdx.x]);
            sum      = fma(a, b, sum);
        }

        __syncthreads();
    }

    if (row < m && col < n) {
        C[col * m + row] = sum;
    }
}

template <bool transA, bool transB>
__global__ void simple_gemm_device(size_t m,
                                   size_t n,
                                   size_t k,
                                   const double *__restrict__ A,
                                   const double *__restrict__ B,
                                   double2 *__restrict__ C //
) {
    const int TILE = 32;
    __shared__ double Asub[TILE][TILE + 1];
    __shared__ double Bsub[TILE][TILE + 1];

    size_t row = blockIdx.y * blockDim.y + threadIdx.y;
    size_t col = blockIdx.x * blockDim.x + threadIdx.x;

    double2 sum{};

    int numTiles = (int)((k + TILE - 1) / TILE);

    for (int t = 0; t < numTiles; ++t) {
        size_t a_col = t * TILE + threadIdx.x;
        if (row < m && a_col < k) {
            Asub[threadIdx.y][threadIdx.x] = load_A_element<double, transA>(A, m, k, row, a_col);
        } else {
            Asub[threadIdx.y][threadIdx.x] = 0.0;
        }

        size_t b_row = t * TILE + threadIdx.y;
        if (col < n && b_row < k) {
            Bsub[threadIdx.y][threadIdx.x] = load_B_element<double, transB>(B, k, n, b_row, col);
        } else {
            Bsub[threadIdx.y][threadIdx.x] = 0.0;
        }

        __syncthreads();

#pragma unroll
        for (int i = 0; i < TILE; ++i) {
            double a = Asub[threadIdx.y][i];
            double b = Bsub[i][threadIdx.x];
            sum      = add(sum, mul(a, b));
        }

        __syncthreads();
    }

    if (row < m && col < n) {
        C[col * m + row] = sum;
    }
}

template <bool transA, bool transB>
__global__ void simple_gemm_device(size_t m,
                                   size_t n,
                                   size_t k,
                                   const cuFloatComplex *__restrict__ A,
                                   const cuFloatComplex *__restrict__ B,
                                   cuDoubleComplex *__restrict__ C //
) {
    const int TILE = 32;
    __shared__ cuFloatComplex Asub[TILE][TILE + 1];
    __shared__ cuFloatComplex Bsub[TILE][TILE + 1];

    size_t row = blockIdx.y * blockDim.y + threadIdx.y;
    size_t col = blockIdx.x * blockDim.x + threadIdx.x;

    double2 sum{};

    int numTiles = (int)((k + TILE - 1) / TILE);

    for (int t = 0; t < numTiles; ++t) {
        size_t a_col = t * TILE + threadIdx.x;
        if (row < m && a_col < k) {
            Asub[threadIdx.y][threadIdx.x] = load_A_element<cuFloatComplex, transA>(A, m, k, row, a_col);
        } else {
            Asub[threadIdx.y][threadIdx.x] = cuFloatComplex{0.0f, 0.0f};
        }

        size_t b_row = t * TILE + threadIdx.y;
        if (col < n && b_row < k) {
            Bsub[threadIdx.y][threadIdx.x] = load_B_element<cuFloatComplex, transB>(B, k, n, b_row, col);
        } else {
            Bsub[threadIdx.y][threadIdx.x] = cuFloatComplex{0.0f, 0.0f};
        }

        __syncthreads();

#pragma unroll
        for (int i = 0; i < TILE; ++i) {
            cuFloatComplex a = Asub[threadIdx.y][i];
            cuFloatComplex b = Bsub[i][threadIdx.x];

            // (ar + i ai)(br + i bi)
            double ar = double(a.x), ai = double(a.y);
            double br = double(b.x), bi = double(b.y);

            // --- real contribution: ar*br - ai*bi ---
            sum.x = fma(ar, br, fma(-ai, bi, sum.x));

            // --- imag contribution: ar*bi + ai*br ---
            sum.y = fma(ar, bi, fma(ai, br, sum.y));
        }

        __syncthreads();
    }

    if (row < m && col < n) {
        C[col * m + row] = sum;
    }
}

template <bool transA, bool transB>
__global__ void simple_gemm_device(size_t m,
                                   size_t n,
                                   size_t k,
                                   const cuDoubleComplex *__restrict__ A,
                                   const cuDoubleComplex *__restrict__ B,
                                   eval::dd::double2_complex *__restrict__ C //
) {
    const int TILE = 32;
    __shared__ cuDoubleComplex Asub[TILE][TILE + 1];
    __shared__ cuDoubleComplex Bsub[TILE][TILE + 1];

    size_t row = blockIdx.y * blockDim.y + threadIdx.y;
    size_t col = blockIdx.x * blockDim.x + threadIdx.x;

    double2 sum_r{}, sum_i{};

    int numTiles = (int)((k + TILE - 1) / TILE);

    for (int t = 0; t < numTiles; ++t) {
        size_t a_col = t * TILE + threadIdx.x;
        if (row < m && a_col < k) {
            Asub[threadIdx.y][threadIdx.x] = load_A_element<cuDoubleComplex, transA>(A, m, k, row, a_col);
        } else {
            Asub[threadIdx.y][threadIdx.x] = cuDoubleComplex{0.0, 0.0};
        }

        size_t b_row = t * TILE + threadIdx.y;
        if (col < n && b_row < k) {
            Bsub[threadIdx.y][threadIdx.x] = load_B_element<cuDoubleComplex, transB>(B, k, n, b_row, col);
        } else {
            Bsub[threadIdx.y][threadIdx.x] = cuDoubleComplex{0.0, 0.0};
        }

        __syncthreads();

#pragma unroll
        for (int i = 0; i < TILE; ++i) {
            cuDoubleComplex a = Asub[threadIdx.y][i];
            cuDoubleComplex b = Bsub[i][threadIdx.x];

            // (ar + i ai)(br + i bi)
            double ar = double(a.x), ai = double(a.y);
            double br = double(b.x), bi = double(b.y);

            // --- real contribution: ar*br - ai*bi ---
            double2 re = sub(mul(ar, br), mul(ai, bi));
            sum_r      = add(sum_r, re);

            // --- imag contribution: ar*bi + ai*br ---
            double2 im = add(mul(ar, bi), mul(ai, br));
            sum_i      = add(sum_i, im);
        }

        __syncthreads();
    }

    if (row < m && col < n) {
        C[col * m + row] = eval::dd::double2_complex{sum_r, sum_i};
    }
}

template <typename Tin, typename Tout>
void simple_gemm(size_t m,
                 size_t n,
                 size_t k,
                 const Tin *A,
                 const Tin *B,
                 Tout *C,
                 cublasOperation_t op_A = CUBLAS_OP_N, //
                 cublasOperation_t op_B = CUBLAS_OP_N  //
) {
    dim3 threadsPerBlock(32, 32);
    dim3 numBlocks((n + threadsPerBlock.x - 1) / threadsPerBlock.x,
                   (m + threadsPerBlock.y - 1) / threadsPerBlock.y);

    if (op_A == CUBLAS_OP_N) {
        if (op_B == CUBLAS_OP_N) {
            simple_gemm_device<false, false><<<numBlocks, threadsPerBlock>>>(m, n, k, A, B, C);
        } else {
            simple_gemm_device<false, true><<<numBlocks, threadsPerBlock>>>(m, n, k, A, B, C);
        }
    } else {
        if (op_B == CUBLAS_OP_N) {
            simple_gemm_device<true, false><<<numBlocks, threadsPerBlock>>>(m, n, k, A, B, C);
        } else {
            simple_gemm_device<true, true><<<numBlocks, threadsPerBlock>>>(m, n, k, A, B, C);
        }
    }
    cudaDeviceSynchronize();
}

} // namespace dd

//------------------------------
// evaluate error
//------------------------------
namespace err {

__global__ void gemm_err_kernel(const size_t m,
                                const size_t n,
                                float *const __restrict__ C,             // calculated value
                                const double *const __restrict__ C_exact // true value
) {
    const size_t idx = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (idx >= m * n) return;
    double2 gap = dd::sub(double(C[idx]), C_exact[idx]);
    double2 err = dd::div(gap, C_exact[idx]);
    C[idx]      = static_cast<float>(fabs(err.x));
}

__global__ void gemm_err_kernel(const size_t m,
                                const size_t n,
                                double *const __restrict__ C,             // calculated value
                                const double2 *const __restrict__ C_exact // true value
) {
    const size_t idx = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (idx >= m * n) return;
    double2 gap = dd::sub(C[idx], C_exact[idx]);
    double2 err = dd::div(gap, C_exact[idx]);
    C[idx]      = fabs(err.x);
}

__global__ void gemm_err_kernel(const size_t m,
                                const size_t n,
                                cuFloatComplex *const __restrict__ C,             // calculated value
                                const cuDoubleComplex *const __restrict__ C_exact // true value
) {
    const size_t idx = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (idx >= m * n) return;
    double2 gap_r = dd::sub(double(C[idx].x), C_exact[idx].x);
    double2 gap_i = dd::sub(double(C[idx].y), C_exact[idx].y);
    double2 err_r = dd::div(gap_r, C_exact[idx].x);
    double2 err_i = dd::div(gap_i, C_exact[idx].y);
    C[idx]        = {static_cast<float>(fabs(err_r.x)), static_cast<float>(fabs(err_i.x))};
}

__global__ void gemm_err_kernel(const size_t m,
                                const size_t n,
                                cuDoubleComplex *const __restrict__ C,                      // calculated value
                                const eval::dd::double2_complex *const __restrict__ C_exact // true value
) {
    const size_t idx = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (idx >= m * n) return;
    double2 gap_r = dd::sub(double(C[idx].x), C_exact[idx].x);
    double2 gap_i = dd::sub(double(C[idx].y), C_exact[idx].y);
    double2 err_r = dd::div(gap_r, C_exact[idx].x);
    double2 err_i = dd::div(gap_i, C_exact[idx].y);
    C[idx]        = {static_cast<double>(fabs(err_r.x)), static_cast<double>(fabs(err_i.x))};
}

template <typename T_lo, typename T_hi>
double2 gemm_err(const size_t m,
                 const size_t n,
                 T_lo *const C,            // calculated value
                 const T_hi *const C_exact // true value
) {
    cudaDeviceSynchronize();
    gemm_err_kernel<<<(m * n + 255) / 256, 256>>>(m, n, C, C_exact);
    if constexpr (std::is_same_v<T_lo, cuDoubleComplex> || std::is_same_v<T_lo, cuFloatComplex>) {
        size_t sizeC = m * n * 2;
        using U      = std::conditional_t<(std::is_same_v<T_lo, cuDoubleComplex>), double, float>;
        std::vector<U> hC(sizeC);
        cudaMemcpy(hC.data(), C, sizeC * sizeof(U), cudaMemcpyDeviceToHost);
        std::sort(hC.begin(), hC.end());
        double err1 = double(hC[sizeC - 1]);
        double err2 = (sizeC & 1) ? double(hC[sizeC / 2]) : ((double(hC[sizeC / 2]) + double(hC[sizeC / 2 - 1])) * 0.5);
        return {err1, err2};
    } else {
        size_t sizeC = m * n;
        std::vector<T_lo> hC(sizeC);
        cudaMemcpy(hC.data(), C, sizeC * sizeof(T_lo), cudaMemcpyDeviceToHost);
        std::sort(hC.begin(), hC.end());
        double err1 = double(hC[sizeC - 1]);
        double err2 = (sizeC & 1) ? double(hC[sizeC / 2]) : ((double(hC[sizeC / 2]) + double(hC[sizeC / 2 - 1])) * 0.5);
        return {err1, err2};
    }
}

} // namespace err

} // namespace eval
