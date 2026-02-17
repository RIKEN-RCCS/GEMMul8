#include "gemmul8.hpp"
#include <cuComplex.h>
#include <cublasLt.h>
#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <iomanip>
#include <iostream>
#include <vector>
#include <random>

//
//
//==========
int imax        = 5;
int total_tests = 0;
int idx_tests   = 0;
int nmin        = 32;
int nmax        = 47;

#if !defined(SEED)
    #define SEED 9999
#endif

inline constexpr uint32_t seed = SEED;

std::vector<cublasOperation_t> opsA = {CUBLAS_OP_N, CUBLAS_OP_T, CUBLAS_OP_C};
std::vector<cublasOperation_t> opsB = {CUBLAS_OP_N, CUBLAS_OP_T, CUBLAS_OP_C};

#if defined(UseFP8)
constexpr gemmul8::Backend backend = gemmul8::Backend::FP8;
#else
constexpr gemmul8::Backend backend = gemmul8::Backend::INT8;
#endif

#if defined(UseLT)
    #define HANDLE_t      cublasLtHandle_t
    #define HANDLECreate  cublasLtCreate
    #define HANDLEDestroy cublasLtDestroy
#else
    #define HANDLE_t      cublasHandle_t
    #define HANDLECreate  cublasCreate
    #define HANDLEDestroy cublasDestroy
#endif
//==========
//
//

__global__ void f2d_kernel(size_t sizeA, const float *const __restrict__ in, double *const __restrict__ out) {
    const auto idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx >= sizeA) return;
    out[idx] = static_cast<double>(in[idx]);
}

void f2d(size_t m,              // rows of A
         size_t n,              // columns of A
         const float *const in, // input
         double *const out)     // output
{
    constexpr size_t block_size = 256;
    const size_t grid_size      = (m * n + block_size - 1) / block_size;
    f2d_kernel<<<grid_size, block_size>>>(m * n, in, out);
    cudaDeviceSynchronize();
}

__global__ void f2d_kernel(size_t sizeA, const cuComplex *const __restrict__ in, cuDoubleComplex *const __restrict__ out) {
    const auto idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx >= sizeA) return;
    cuComplex intmp = in[idx];
    cuDoubleComplex outtmp;
    outtmp.x = static_cast<double>(intmp.x);
    outtmp.y = static_cast<double>(intmp.y);
    out[idx] = outtmp;
}

void f2d(size_t m,                   // rows of A
         size_t n,                   // columns of A
         const cuComplex *const in,  // input
         cuDoubleComplex *const out) // output
{
    constexpr size_t block_size = 256;
    const size_t grid_size      = (m * n + block_size - 1) / block_size;
    f2d_kernel<<<grid_size, block_size>>>(m * n, in, out);
    cudaDeviceSynchronize();
}

template <typename T> void fill_random(T *data, int n) {
    std::mt19937 mt(seed);
    for (int i = 0; i < n; ++i)
        data[i] = static_cast<T>(mt()) / RAND_MAX;
}

template <> void fill_random(cuComplex *data, int n) {
    std::mt19937 mt(seed);
    for (int i = 0; i < n; ++i)
        data[i] = make_cuComplex(static_cast<float>(mt()) / RAND_MAX,
                                 static_cast<float>(mt()) / RAND_MAX);
}

template <> void fill_random(cuDoubleComplex *data, int n) {
    std::mt19937 mt(seed);
    for (int i = 0; i < n; ++i)
        data[i] = make_cuDoubleComplex(static_cast<double>(mt()) / RAND_MAX,
                                       static_cast<double>(mt()) / RAND_MAX);
}

template <typename T> struct Tconst {
    static constexpr T alpha[5] = {T(1), T(1), T(-1), T(-1), T(-1.5)};
    static constexpr T beta[5]  = {T(0), T(1), T(0), T(1), T(1.5)};
};
template <> struct Tconst<cuDoubleComplex> {
    static constexpr cuDoubleComplex alpha[5] = {
        { 1.0, 0.0},
        { 1.0, 0.0},
        {-1.0, 0.0},
        {-1.0, 0.0},
        {-1.5, 1.2}
    };
    static constexpr cuDoubleComplex beta[5] = {
        {0.0, 0.0},
        {1.0, 0.0},
        {0.0, 0.0},
        {1.0, 0.0},
        {1.5, 1.2}
    };
};
template <> struct Tconst<cuFloatComplex> {
    static constexpr cuFloatComplex alpha[5] = {
        { 1.0f, 0.0f},
        { 1.0f, 0.0f},
        {-1.0f, 0.0f},
        {-1.0f, 0.0f},
        {-1.5f, 1.2f}
    };
    static constexpr cuFloatComplex beta[5] = {
        {0.0f, 0.0f},
        {1.0f, 0.0f},
        {0.0f, 0.0f},
        {1.0f, 0.0f},
        {1.5f, 1.2f}
    };
};

template <typename T>
__global__ void ones_kernel(int64_t n, T *__restrict__ x) {
    const auto idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;
    x[idx] = Tconst<T>::alpha[4];
}

template <typename T>
__forceinline__ void ones(int64_t m, int64_t n, T *x) {
    constexpr size_t block_size = 256;
    const size_t grid_size      = (m * n + block_size - 1) / block_size;
    ones_kernel<T><<<grid_size, block_size>>>(m * n, x);
}

template <typename T> inline constexpr bool isComplex = (std::is_same_v<T, cuComplex> || std::is_same_v<T, cuDoubleComplex>);
template <typename T> inline constexpr bool isFloat   = (std::is_same_v<T, float> || std::is_same_v<T, cuComplex>);

template <typename T> double err_check(const int sizeC, const T *hC1, const T *hC2) {
    double err = 0.;
    for (int i = 0; i < sizeC; i++) {
        if constexpr (isComplex<T>) {
            double refx  = static_cast<double>(hC1[i].x);
            double refy  = static_cast<double>(hC1[i].y);
            double diffx = (refx == 0.0) ? fabs(static_cast<double>(hC1[i].x) - static_cast<double>(hC2[i].x))
                                         : fabs((static_cast<double>(hC1[i].x) - static_cast<double>(hC2[i].x)) / refx);
            double diffy = (refy == 0.0) ? fabs(static_cast<double>(hC1[i].y) - static_cast<double>(hC2[i].y))
                                         : fabs((static_cast<double>(hC1[i].y) - static_cast<double>(hC2[i].y)) / refy);
            err          = std::max(err, std::max(diffx, diffy));
        } else {
            double ref  = static_cast<double>(hC1[i]);
            double diff = (ref == 0.0) ? fabs(static_cast<double>(hC1[i]) - static_cast<double>(hC2[i]))
                                       : fabs((static_cast<double>(hC1[i]) - static_cast<double>(hC2[i])) / ref);
            err         = std::max(err, diff);
        }
    }
    return err;
}

const char *op_to_char(cublasOperation_t op) {
    switch (op) {
    case CUBLAS_OP_N:
        return "N";
    case CUBLAS_OP_T:
        return "T";
    case CUBLAS_OP_C:
        return "C";
    default:
        return "?";
    }
}

template <typename T> void run_gemm(HANDLE_t handle,
                                    cublasHandle_t handle_gemm,
                                    int m,
                                    int n,
                                    int k,
                                    cublasOperation_t transa,
                                    cublasOperation_t transb,
                                    T alpha,
                                    T beta,
                                    const T *dA,
                                    const T *dB,
                                    T *dC1,
                                    T *dC2,
                                    T *hC1,
                                    void *work) {

    float progress = float(idx_tests) / total_tests;
    int barWidth   = 50;

    std::cout << "\r[";
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0)
              << "% (" << idx_tests << "/" << total_tests
              << " NOW: "
              << (isFloat<T> ? (isComplex<T> ? "CGEMM" : "SGEMM") : (isComplex<T> ? "ZGEMM" : "DGEMM"))
              << ")"
              << std::flush;

    const int sizeC = m * n;

    cudaDeviceSynchronize();
    if constexpr (std::is_same_v<T, float>) {
        cudaMemcpy(dC1, hC1, sizeof(T) * sizeC, cudaMemcpyHostToDevice);
        cublasSgemm(handle_gemm, transa, transb, m, n, k, &alpha, dA, (transa == CUBLAS_OP_N ? m : k), dB, (transb == CUBLAS_OP_N ? k : n), &beta, dC1, m);
    } else if constexpr (std::is_same_v<T, double>) {
        cudaMemcpy(dC1, hC1, sizeof(T) * sizeC, cudaMemcpyHostToDevice);
        cublasDgemm(handle_gemm, transa, transb, m, n, k, &alpha, dA, (transa == CUBLAS_OP_N ? m : k), dB, (transb == CUBLAS_OP_N ? k : n), &beta, dC1, m);
    } else if constexpr (std::is_same_v<T, cuComplex>) {
        cudaMemcpy(dC1, hC1, sizeof(T) * sizeC, cudaMemcpyHostToDevice);
        cublasCgemm(handle_gemm, transa, transb, m, n, k, &alpha, dA, (transa == CUBLAS_OP_N ? m : k), dB, (transb == CUBLAS_OP_N ? k : n), &beta, dC1, m);
    } else {
        cudaMemcpy(dC1, hC1, sizeof(T) * sizeC, cudaMemcpyHostToDevice);
        cublasZgemm(handle_gemm, transa, transb, m, n, k, &alpha, dA, (transa == CUBLAS_OP_N ? m : k), dB, (transb == CUBLAS_OP_N ? k : n), &beta, dC1, m);
    }

    std::vector<T> hC3(sizeC), hC4(sizeC);
    cudaMemcpy(hC3.data(), dC1, sizeof(T) * sizeC, cudaMemcpyDeviceToHost);

    const unsigned num_moduli_start = (std::is_same_v<T, float> || std::is_same_v<T, cuComplex>) ? 6u : 7u;
    const unsigned num_moduli_end   = (std::is_same_v<T, float> || std::is_same_v<T, cuComplex>) ? 15u : 20u;

    bool breakflag = true;
    for (int pass = 0; pass < 2; ++pass) {
        const bool fastmode = (pass == 0);
        for (unsigned num_moduli = num_moduli_start; num_moduli < num_moduli_end; ++num_moduli) {
            cudaMemcpy(dC2, hC1, sizeof(T) * sizeC, cudaMemcpyHostToDevice);
            cudaDeviceSynchronize();
            gemmul8::gemm<T, backend>(handle, transa, transb, m, n, k, &alpha, dA, (transa == CUBLAS_OP_N ? m : k), dB, (transb == CUBLAS_OP_N ? k : n), &beta, dC2, m, num_moduli, fastmode, work);
            cudaDeviceSynchronize();
            cudaMemcpy(hC4.data(), dC2, sizeof(T) * sizeC, cudaMemcpyDeviceToHost);

            double err = err_check(sizeC, hC3.data(), hC4.data());

            if (err > 1) {

                if (breakflag) {
                    std::cout << std::endl;
                    breakflag = false;
                }
                if constexpr (std::is_same_v<T, float>) {
                    std::cout << std::left << std::setw(11) << (fastmode ? "SGEMM fast" : "SGEMM accu");
                } else if (std::is_same_v<T, double>) {
                    std::cout << std::left << std::setw(11) << (fastmode ? "DGEMM fast" : "DGEMM accu");
                } else if (std::is_same_v<T, cuComplex>) {
                    std::cout << std::left << std::setw(11) << (fastmode ? "CGEMM fast" : "CGEMM accu");
                } else {
                    std::cout << std::left << std::setw(11) << (fastmode ? "ZGEMM fast" : "ZGEMM accu");
                }
                double a, b;
                if constexpr (isComplex<T>) {
                    a = alpha.x;
                    b = beta.x;
                } else {
                    a = alpha;
                    b = beta;
                }
                std::cout << std::left
                          << " ( m=" << std::setw(2) << m
                          << ", n=" << std::setw(2) << n
                          << ", k=" << std::setw(2) << k
                          << ", transa=" << op_to_char(transa)
                          << ", transb=" << op_to_char(transb)
                          << ", alpha=" << std::scientific << std::setprecision(1) << a
                          << ", beta=" << std::scientific << std::setprecision(1) << b
                          << ", num_mod=" << std::setw(2) << num_moduli
                          << ", err=" << std::scientific << std::setprecision(1) << err
                          << " )"
                          << std::endl;
            }
        }
    }

    // std::cout << std::right << std::setw(3) << idx_tests;
    // std::cout << "/";
    // std::cout << std::left << std::setw(3) << total_tests << std::endl;
}

int main() {
#if defined(SGEMM)
    total_tests += (nmax - nmin + 1) * 4;
#endif
#if defined(DGEMM)
    total_tests += (nmax - nmin + 1) * 4;
#endif
#if defined(CGEMM)
    total_tests += (nmax - nmin + 1) * 9;
#endif
#if defined(ZGEMM)
    total_tests += (nmax - nmin + 1) * 9;
#endif
    total_tests *= imax;

    srand(0);
    cudaSetDevice(0);
    HANDLE_t handle;
    HANDLECreate(&handle);
    cublasHandle_t handle_gemm;
    cublasCreate(&handle_gemm);

    const int maxsize = nmax * nmax;

    std::cout << std::endl;
#if defined(SGEMM)
    {
        using Type = float;
        std::vector<Type> hA(maxsize), hB(maxsize), hC1(maxsize);
        fill_random(hA.data(), maxsize);
        fill_random(hB.data(), maxsize);
        fill_random(hC1.data(), maxsize);
        Type *dA, *dB, *dC, *dD;
        cudaMalloc(&dA, sizeof(Type) * maxsize);
        cudaMalloc(&dB, sizeof(Type) * maxsize);
        cudaMalloc(&dC, sizeof(Type) * maxsize);
        cudaMalloc(&dD, sizeof(Type) * maxsize);
        cudaMemcpy(dA, hA.data(), sizeof(Type) * maxsize, cudaMemcpyHostToDevice);
        cudaMemcpy(dB, hB.data(), sizeof(Type) * maxsize, cudaMemcpyHostToDevice);

        const size_t worksize = gemmul8::workSize<isComplex<Type>, backend>(maxsize, maxsize, maxsize, 20u);
        void *work;
        cudaMalloc(&work, worksize);

        for (auto transa : opsA)
            for (auto transb : opsB)
                for (int i = 0; i < imax; ++i)
                    for (int m = nmin; m <= nmax; ++m) {
                        Type alpha = Tconst<Type>::alpha[i];
                        Type beta  = Tconst<Type>::beta[i];
                        int n      = m;
                        int k      = m;
                        if (transa != CUBLAS_OP_C && transb != CUBLAS_OP_C) {
                            idx_tests++;
                            run_gemm<Type>(handle, handle_gemm, m, n, k, transa, transb, alpha, beta, dA, dB, dC, dD, hC1.data(), work);
                        } else if (isComplex<Type>) {
                            idx_tests++;
                            run_gemm<Type>(handle, handle_gemm, m, n, k, transa, transb, alpha, beta, dA, dB, dC, dD, hC1.data(), work);
                        }
                    }

        cudaFree(dA);
        cudaFree(dB);
        cudaFree(dC);
        cudaFree(dD);
        cudaFree(work);
    }
#endif

#if defined(DGEMM)
    {
        using Type = double;
        std::vector<Type> hA(maxsize), hB(maxsize), hC1(maxsize);
        fill_random(hA.data(), maxsize);
        fill_random(hB.data(), maxsize);
        fill_random(hC1.data(), maxsize);
        Type *dA, *dB, *dC, *dD;
        cudaMalloc(&dA, sizeof(Type) * maxsize);
        cudaMalloc(&dB, sizeof(Type) * maxsize);
        cudaMalloc(&dC, sizeof(Type) * maxsize);
        cudaMalloc(&dD, sizeof(Type) * maxsize);
        cudaMemcpy(dA, hA.data(), sizeof(Type) * maxsize, cudaMemcpyHostToDevice);
        cudaMemcpy(dB, hB.data(), sizeof(Type) * maxsize, cudaMemcpyHostToDevice);

        const size_t worksize = gemmul8::workSize<isComplex<Type>, backend>(maxsize, maxsize, maxsize, 20u);
        void *work;
        cudaMalloc(&work, worksize);

        for (auto transa : opsA)
            for (auto transb : opsB)
                for (int i = 0; i < imax; ++i)
                    for (int m = nmin; m <= nmax; ++m) {
                        Type alpha = Tconst<Type>::alpha[i];
                        Type beta  = Tconst<Type>::beta[i];
                        int n      = m;
                        int k      = m;
                        if (transa != CUBLAS_OP_C && transb != CUBLAS_OP_C) {
                            idx_tests++;
                            run_gemm<Type>(handle, handle_gemm, m, n, k, transa, transb, alpha, beta, dA, dB, dC, dD, hC1.data(), work);
                        } else if (isComplex<Type>) {
                            idx_tests++;
                            run_gemm<Type>(handle, handle_gemm, m, n, k, transa, transb, alpha, beta, dA, dB, dC, dD, hC1.data(), work);
                        }
                    }

        cudaFree(dA);
        cudaFree(dB);
        cudaFree(dC);
        cudaFree(dD);
        cudaFree(work);
    }
#endif

#if defined(CGEMM)
    {
        using Type = cuComplex;
        std::vector<Type> hA(maxsize), hB(maxsize), hC1(maxsize);
        fill_random(hA.data(), maxsize);
        fill_random(hB.data(), maxsize);
        fill_random(hC1.data(), maxsize);
        Type *dA, *dB, *dC, *dD;
        cudaMalloc(&dA, sizeof(Type) * maxsize);
        cudaMalloc(&dB, sizeof(Type) * maxsize);
        cudaMalloc(&dC, sizeof(Type) * maxsize);
        cudaMalloc(&dD, sizeof(Type) * maxsize);
        cudaMemcpy(dA, hA.data(), sizeof(Type) * maxsize, cudaMemcpyHostToDevice);
        cudaMemcpy(dB, hB.data(), sizeof(Type) * maxsize, cudaMemcpyHostToDevice);

        const size_t worksize = gemmul8::workSize<isComplex<Type>, backend>(maxsize, maxsize, maxsize, 20u);
        void *work;
        cudaMalloc(&work, worksize);

        for (auto transa : opsA)
            for (auto transb : opsB)
                for (int i = 0; i < imax; ++i)
                    for (int m = nmin; m <= nmax; ++m) {
                        Type alpha = Tconst<Type>::alpha[i];
                        Type beta  = Tconst<Type>::beta[i];
                        int n      = m;
                        int k      = m;
                        if (transa != CUBLAS_OP_C && transb != CUBLAS_OP_C) {
                            idx_tests++;
                            run_gemm<Type>(handle, handle_gemm, m, n, k, transa, transb, alpha, beta, dA, dB, dC, dD, hC1.data(), work);
                        } else if (isComplex<Type>) {
                            idx_tests++;
                            run_gemm<Type>(handle, handle_gemm, m, n, k, transa, transb, alpha, beta, dA, dB, dC, dD, hC1.data(), work);
                        }
                    }

        cudaFree(dA);
        cudaFree(dB);
        cudaFree(dC);
        cudaFree(dD);
        cudaFree(work);
    }
#endif

#if defined(ZGEMM)
    {
        using Type = cuDoubleComplex;
        std::vector<Type> hA(maxsize), hB(maxsize), hC1(maxsize);
        fill_random(hA.data(), maxsize);
        fill_random(hB.data(), maxsize);
        fill_random(hC1.data(), maxsize);
        Type *dA, *dB, *dC, *dD;
        cudaMalloc(&dA, sizeof(Type) * maxsize);
        cudaMalloc(&dB, sizeof(Type) * maxsize);
        cudaMalloc(&dC, sizeof(Type) * maxsize);
        cudaMalloc(&dD, sizeof(Type) * maxsize);
        cudaMemcpy(dA, hA.data(), sizeof(Type) * maxsize, cudaMemcpyHostToDevice);
        cudaMemcpy(dB, hB.data(), sizeof(Type) * maxsize, cudaMemcpyHostToDevice);

        const size_t worksize = gemmul8::workSize<isComplex<Type>, backend>(maxsize, maxsize, maxsize, 20u);
        void *work;
        cudaMalloc(&work, worksize);

        for (auto transa : opsA)
            for (auto transb : opsB)
                for (int i = 0; i < imax; ++i)
                    for (int m = nmin; m <= nmax; ++m) {
                        Type alpha = Tconst<Type>::alpha[i];
                        Type beta  = Tconst<Type>::beta[i];
                        int n      = m;
                        int k      = m;
                        if (transa != CUBLAS_OP_C && transb != CUBLAS_OP_C) {
                            idx_tests++;
                            run_gemm<Type>(handle, handle_gemm, m, n, k, transa, transb, alpha, beta, dA, dB, dC, dD, hC1.data(), work);
                        } else if (isComplex<Type>) {
                            idx_tests++;
                            run_gemm<Type>(handle, handle_gemm, m, n, k, transa, transb, alpha, beta, dA, dB, dC, dD, hC1.data(), work);
                        }
                    }

        cudaFree(dA);
        cudaFree(dB);
        cudaFree(dC);
        cudaFree(dD);
        cudaFree(work);
    }
#endif

    std::cout << std::endl;
    std::cout << std::endl;
    HANDLEDestroy(handle);
    cublasDestroy(handle_gemm);
    return 0;
}
