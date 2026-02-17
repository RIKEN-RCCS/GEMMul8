#include "../include/gemmul8.hpp"
#include "../testing/self_hipify.hpp"
#include <iostream>

void disp_mat(int m, int n, double *Mat) {
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j)
            printf("%24.16e  ", Mat[j * m + i]);
        printf("\n");
    }
    printf("\n");
}

int main(int argc, char **argv) {
    const gemmul8::Backend backend = gemmul8::Backend::INT8;
    const unsigned num_moduli      = 13;
    const bool fastmode            = false;
    const bool is_complex          = false;

    cudaSetDevice(0);
    cublasLtHandle_t handle;
    cublasLtCreate(&handle);

    const int m = 4, n = 3, k = 5;
    const int lda = m, ldb = k, ldc = m;
    std::vector<double> hA = {0x1.13491b78f7ff1p-1, 0x1.d5797d024f750p+0, -0x1.2121e4d9576a2p+1, 0x1.b96ec80cedfb6p-1,
                              0x1.466a65212f053p-2, -0x1.4ec4a901fe3c4p+0, -0x1.bbff8c0e700a1p-2, 0x1.5ed8f2ba5f2dbp-2,
                              0x1.ca08e9321d439p+1, 0x1.627ce99fd7ed1p+1, -0x1.599230c5450f8p+0, 0x1.84785f44e10f1p+1,
                              0x1.73682ebd0c291p-1, -0x1.0245d3a33f7d8p-4, 0x1.6df2c829f659fp-1, -0x1.a3c53ea980203p-3,
                              -0x1.fc7ec8b9281f7p-4, 0x1.7d5cd28a5e35bp+0, 0x1.68b67bfca10cfp+0, 0x1.6acd1f3bd1cafp+0};

    std::vector<double> hB = {0x1.57ce78e868ad7p-1, -0x1.351ddceb47a8bp+0, 0x1.6f39e78dc4de4p-1, 0x1.a1571993bf63bp+0,
                              0x1.f4a0918ad43eep-2, 0x1.08e1a41eff3c4p+0, 0x1.742a49c7a8c1fp-1, -0x1.36b937c0e54f0p-2,
                              0x1.2ceca451a1789p-2, -0x1.9316bb4db16cfp-1, 0x1.c6dbcad09ddd8p-1, -0x1.25a662f3a6d75p+0,
                              -0x1.11a17e8d7e02fp+0, -0x1.9e769ce56b489p-1, -0x1.78de4dacf30d6p+1};

    std::vector<double> hC(lda * n, 0.0);
    std::vector<double> hC_exact = {0x1.d51136ef01e9dp+1, 0x1.5b07528da2db2p+2, -0x1.b7d034d197c42p-4, 0x1.59b0e0e988db5p+1,
                                    0x1.ad784e3b16dc5p-7, -0x1.15b1323003b06p+0, -0x1.922e5c1c4b38bp+1, -0x1.e95843f74c224p-1,
                                    -0x1.f79e85fefa19bp+1, -0x1.0a9fa599dc6d9p+2, -0x1.32cc3fa2fc921p+2, -0x1.b82c3fad3ab16p+2};

    double *A, *B, *C;
    cudaMalloc(reinterpret_cast<void **>(&A), lda * k * sizeof(double));
    cudaMalloc(reinterpret_cast<void **>(&B), ldb * n * sizeof(double));
    cudaMalloc(reinterpret_cast<void **>(&C), ldc * n * sizeof(double));

    cudaMemcpy(A, hA.data(), lda * k * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(B, hB.data(), ldb * n * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(C, hC.data(), ldc * n * sizeof(double), cudaMemcpyHostToDevice);

    void *work;
    const size_t lwork = gemmul8::workSize<is_complex, backend>(m, n, k, num_moduli);
    cudaMalloc(&work, lwork);

    const double alpha = 1.0;
    const double beta  = 0.0;
    gemmul8::gemm<double, backend>(handle, CUBLAS_OP_N, CUBLAS_OP_N,
                                   m, n, k, &alpha, A, lda, B, ldb, &beta, C, ldc,
                                   num_moduli, fastmode, work);
    cudaMemcpy(hC.data(), C, ldc * n * sizeof(double), cudaMemcpyDeviceToHost);

    double nrm = 0.0;
    for (int i = 0; i < ldc * n; ++i) {
        const double err = hC_exact[i] - hC[i];
        nrm              = std::fma(err, err, nrm);
    }
    nrm = std::sqrt(nrm);

    printf("===== A: %d * %d =====\n", m, k);
    disp_mat(m, k, hA.data());

    printf("===== B: %d * %d =====\n", k, n);
    disp_mat(k, n, hB.data());

    printf("===== C_exact: %d * %d =====\n", m, n);
    disp_mat(m, n, hC_exact.data());

    printf("===== C: %d * %d =====\n", m, n);
    disp_mat(m, n, hC.data());

    printf("error = %e\n", nrm);

    cudaFree(work);
    cudaFree(C);
    cudaFree(B);
    cudaFree(A);
    cublasLtDestroy(handle);
    cudaDeviceReset();
    return 0;
}
