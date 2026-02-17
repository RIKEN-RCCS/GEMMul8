#include "gemmul8.hpp"
#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <dlfcn.h>
#include <iostream>

#define CHECK_CUDA(call)                                               \
    do {                                                               \
        cudaError_t err = call;                                        \
        if (err != cudaSuccess) {                                      \
            std::cerr << "CUDA error: " << cudaGetErrorString(err)     \
                      << " (" << __FILE__ << ":" << __LINE__ << ")\n"; \
            exit(EXIT_FAILURE);                                        \
        }                                                              \
    } while (0)

#define CHECK_CUBLAS(call)                                                          \
    do {                                                                            \
        cublasStatus_t status = call;                                               \
        if (status != CUBLAS_STATUS_SUCCESS) {                                      \
            std::cerr << "cuBLAS error at " << __FILE__ << ":" << __LINE__ << "\n"; \
            exit(EXIT_FAILURE);                                                     \
        }                                                                           \
    } while (0)

void fill_random(double *data, int n) {
    for (int i = 0; i < n; ++i) data[i] = static_cast<double>(rand()) / RAND_MAX;
}
void fill_random(float *data, int n) {
    for (int i = 0; i < n; ++i) data[i] = static_cast<float>(rand()) / RAND_MAX;
}

void run_dgemm(cublasHandle_t handle, const double *dA, const double *dB, double *dC, double *dC2, int m, int n, int k) {
    const double alpha = 1.0, beta = 0.0;

    std::cout << "HIJACK DGEMM" << std::endl;
    cudaDeviceSynchronize();
    CHECK_CUBLAS(cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, m, n, k, &alpha, dA, m, dB, k, &beta, dC, m));
    cudaDeviceSynchronize();

    const unsigned num_moduli = 18u;                                           // Accuracy knob: 2 <= num_moduli <= 20
    const bool fastmode       = false;                                         // true (fast mode) or false (accurate mode)
    const size_t worksize     = gemmul8::workSize<false>(m, n, k, num_moduli); // calculate required memory (Byte)
    void *work;
    cudaMalloc(&work, worksize);
    std::cout << "gemmul8 DGEMM" << std::endl;
    cudaDeviceSynchronize();
    gemmul8::gemm<double>(handle, CUBLAS_OP_N, CUBLAS_OP_N, m, n, k, &alpha, dA, m, dB, k, &beta, dC2, m, num_moduli, fastmode, work, nullptr, nullptr, false, false, false, false);
    cudaDeviceSynchronize();

    std::vector<double> hC(m * n), hC2(m * n);
    CHECK_CUDA(cudaMemcpy(hC.data(), dC, sizeof(double) * m * n, cudaMemcpyDeviceToHost));
    CHECK_CUDA(cudaMemcpy(hC2.data(), dC2, sizeof(double) * m * n, cudaMemcpyDeviceToHost));

    double err_chk = 0.;
    for (int i = 0; i < m * n; i++) {
        double diff = hC[i] - hC2[i];
        err_chk += diff * diff;
    }
    std::cout << "[DGEMM] L2 error: " << err_chk;
    if (err_chk > 1e-6) std::cout << "  L2 Error Too High!";
    std::cout << std::endl;

    cudaFree(work);
}

void run_sgemm(cublasHandle_t handle, const float *dA, const float *dB, float *dC, float *dC2, int m, int n, int k) {
    const float alpha = 1.0f, beta = 0.0f;

    std::cout << "HIJACK SGEMM" << std::endl;
    cudaDeviceSynchronize();
    CHECK_CUBLAS(cublasSgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, m, n, k, &alpha, dA, m, dB, k, &beta, dC, m));
    cudaDeviceSynchronize();

    const unsigned num_moduli = 13u;                                           // Accuracy knob: 2 <= num_moduli <= 20
    const bool fastmode       = false;                                         // true (fast mode) or false (accurate mode)
    const size_t worksize     = gemmul8::workSize<false>(m, n, k, num_moduli); // calculate required memory (Byte)
    void *work;
    cudaMalloc(&work, worksize);
    std::cout << "gemmul8 SGEMM" << std::endl;
    cudaDeviceSynchronize();
    gemmul8::gemm<float>(handle, CUBLAS_OP_N, CUBLAS_OP_N, m, n, k, &alpha, dA, m, dB, k, &beta, dC2, m, num_moduli, fastmode, work, nullptr, nullptr, false, false, false, false);
    cudaDeviceSynchronize();

    std::vector<float> hC(m * n), hC2(m * n);
    CHECK_CUDA(cudaMemcpy(hC.data(), dC, sizeof(float) * m * n, cudaMemcpyDeviceToHost));
    CHECK_CUDA(cudaMemcpy(hC2.data(), dC2, sizeof(float) * m * n, cudaMemcpyDeviceToHost));

    double err_chk = 0.;
    for (int i = 0; i < m * n; i++) {
        double diff = static_cast<double>(hC[i]) - static_cast<double>(hC2[i]);
        err_chk += diff * diff;
    }
    std::cout << "[SGEMM] L2 error: " << err_chk;
    if (err_chk > 1e-3) std::cout << "  L2 Error Too High!";
    std::cout << std::endl;

    cudaFree(work);
}

int main() {
    srand(0);
    CHECK_CUDA(cudaSetDevice(0));

    int m1 = 100, k1 = 80, n1 = 90;  // A1: m1×k1, B1: k1×n1
    int m2 = 120, k2 = 70, n2 = 110; // A2: m2×k2, B2: k2×n2

    double *A1, *A2, *B1, *B2;
    A1 = (double *)malloc(sizeof(double) * 120 * 120);
    A2 = (double *)malloc(sizeof(double) * 120 * 120);
    B1 = (double *)malloc(sizeof(double) * 120 * 120);
    B2 = (double *)malloc(sizeof(double) * 120 * 120);
    fill_random(A1, 120 * 120);
    fill_random(A2, 120 * 120);
    fill_random(B1, 120 * 120);
    fill_random(B2, 120 * 120);

    float *sA1, *sA2, *sB1, *sB2;
    sA1 = (float *)malloc(sizeof(float) * 120 * 120);
    sA2 = (float *)malloc(sizeof(float) * 120 * 120);
    sB1 = (float *)malloc(sizeof(float) * 120 * 120);
    sB2 = (float *)malloc(sizeof(float) * 120 * 120);
    fill_random(sA1, 120 * 120);
    fill_random(sA2, 120 * 120);
    fill_random(sB1, 120 * 120);
    fill_random(sB2, 120 * 120);

    double *dA1, *dA2, *dB1, *dB2, *dCd, *dCd2;
    float *dsA1, *dsA2, *dsB1, *dsB2, *dsC, *dsC2;

    CHECK_CUDA(cudaMalloc(&dA1, sizeof(double) * 120 * 120));
    CHECK_CUDA(cudaMalloc(&dA2, sizeof(double) * 120 * 120));
    CHECK_CUDA(cudaMalloc(&dB1, sizeof(double) * 120 * 120));
    CHECK_CUDA(cudaMalloc(&dB2, sizeof(double) * 120 * 120));

    CHECK_CUDA(cudaMalloc(&dsA1, sizeof(float) * 120 * 120));
    CHECK_CUDA(cudaMalloc(&dsA2, sizeof(float) * 120 * 120));
    CHECK_CUDA(cudaMalloc(&dsB1, sizeof(float) * 120 * 120));
    CHECK_CUDA(cudaMalloc(&dsB2, sizeof(float) * 120 * 120));

    int maxCd = 120 * 120;
    int maxCs = maxCd;
    CHECK_CUDA(cudaMalloc(&dCd, sizeof(double) * maxCd));
    CHECK_CUDA(cudaMalloc(&dsC, sizeof(float) * maxCs));
    CHECK_CUDA(cudaMalloc(&dCd2, sizeof(double) * maxCd));
    CHECK_CUDA(cudaMalloc(&dsC2, sizeof(float) * maxCs));

    CHECK_CUDA(cudaMemcpy(dA1, A1, sizeof(double) * 120 * 120, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(dA2, A2, sizeof(double) * 120 * 120, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(dB1, B1, sizeof(double) * 120 * 120, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(dB2, B2, sizeof(double) * 120 * 120, cudaMemcpyHostToDevice));

    CHECK_CUDA(cudaMemcpy(dsA1, sA1, sizeof(float) * 120 * 120, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(dsA2, sA2, sizeof(float) * 120 * 120, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(dsB1, sB1, sizeof(float) * 120 * 120, cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(dsB2, sB2, sizeof(float) * 120 * 120, cudaMemcpyHostToDevice));

    cublasHandle_t handle;
    CHECK_CUBLAS(cublasCreate(&handle));

    std::cout << "Running GEMMs...\n";

    // ---- 2〜7: DGEMM ----
    run_dgemm(handle, dA1, dB1, dCd, dCd2, m1, n1, k1); // 2
    run_dgemm(handle, dA1, dB2, dCd, dCd2, m1, n2, k2); // 3
    run_dgemm(handle, dA1, dB1, dCd, dCd2, m1, n1, k1); // 4
    run_dgemm(handle, dA2, dB1, dCd, dCd2, m2, n1, k2); // 5
    run_dgemm(handle, dA2, dB1, dCd, dCd2, m2, n1, k2); // 6
    run_dgemm(handle, dA2, dB2, dCd, dCd2, m2, n2, k2); // 7

    // ---- 9〜14: SGEMM ----
    run_sgemm(handle, dsA1, dsB1, dsC, dsC2, m1, n1, k1); // 9
    run_sgemm(handle, dsA1, dsB2, dsC, dsC2, m1, n2, k2); // 10
    run_dgemm(handle, dA1, dB2, dCd, dCd2, m1, n2, k2);   // 11
    run_sgemm(handle, dsA2, dsB1, dsC, dsC2, m2, n1, k2); // 12
    run_sgemm(handle, dsA2, dsB1, dsC, dsC2, m2, n1, k2); // 13
    run_sgemm(handle, dsA1, dsB1, dsC, dsC2, m1, n1, k1); // 14

    std::cout << "All GEMMs completed.\n";

    CHECK_CUBLAS(cublasDestroy(handle));

    cudaFree(dA1);
    cudaFree(dA2);
    cudaFree(dB1);
    cudaFree(dB2);
    cudaFree(dsA1);
    cudaFree(dsA2);
    cudaFree(dsB1);
    cudaFree(dsB2);
    cudaFree(dCd);
    cudaFree(dsC);
    cudaFree(dCd2);
    cudaFree(dsC2);

    free(A1);
    free(A2);
    free(B1);
    free(B2);
    free(sA1);
    free(sA2);
    free(sB1);
    free(sB2);
    
    CHECK_CUDA(cudaDeviceReset());

    return 0;
}
