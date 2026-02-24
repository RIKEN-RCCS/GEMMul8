#pragma once
#include "getWatt.hpp"

template <typename T, gemmul8::Backend backend = gemmul8::Backend::FP8>
__inline__ void watt_check(std::string &deviceName, std::string &dateTime) {
    std::string fileName = std::string("oz2_results_") + backendType<backend> + "_" + gemmTraits<T>::prefix + "_watt_" + deviceName + "_" + dateTime + ".csv";
    std::ofstream outFile(fileName);

    CHECK_CUDA(cudaSetDevice(0));
    cublasHandle_t handle;
    CHECK_CUBLAS(cublasCreate(&handle));
    cublasLtHandle_t handleLt;
    CHECK_CUBLASLT(cublasLtCreate(&handleLt));

    outFile << std::scientific;
    std::cout << std::scientific;
    outFile << "phi,m,n,k,function,err_max,err_med,watt,GFLOPS/watt," << std::endl;
    std::cout << "phi,m,n,k,function,err_max,err_med,watt,GFLOPS/watt," << std::endl;

    cudaEvent_t start, stop;
    CHECK_CUDA(cudaEventCreate(&start));
    CHECK_CUDA(cudaEventCreate(&stop));
    std::vector<float> times(mainloop, 0.0);

    using accu_t        = typename gemmTraits<T>::ACCU_TYPE;
    size_t total_memory = size_t(GPU_MEM_MB);

    const double phi = 0.5;
    const T alpha    = gemmTraits<T>::one();
    const T beta     = gemmTraits<T>::zero();

    for (size_t m = 1024; m <= size_max; m *= 2) {
        size_t n = m;
        for (size_t k = 1024; k <= size_max; k *= 2) {

            const size_t size_A           = m * k;
            const size_t size_B           = k * n;
            const size_t size_C           = m * n;
            const unsigned num_moduli_max = NUM_MODULI_MAX<T>;
            const size_t lwork_gemmul8    = gemmul8::workSize<gemmTraits<T>::is_complex, backend>(m, n, k, num_moduli_max);
            const size_t lwork            = std::max(size_C * sizeof(accu_t), lwork_gemmul8);
            const size_t total_work       = lwork + (size_A + size_B + size_C) * sizeof(T);
            if ((total_work + 256ULL * sizeof(accu_t)) * 1.e-9 > total_memory) {
                continue;
            }

            T *A, *B, *C;
            accu_t *C_hi;
            std::vector<accu_t> C_hi_h(size_C);
            void *work;

            CHECK_CUDA(cudaMalloc(reinterpret_cast<void **>(&A), size_A * sizeof(T)));
            CHECK_CUDA(cudaMalloc(reinterpret_cast<void **>(&B), size_B * sizeof(T)));
            CHECK_CUDA(cudaMalloc(reinterpret_cast<void **>(&C), size_C * sizeof(T)));
            CHECK_CUDA(cudaMalloc(&work, lwork_gemmul8));
            C_hi = reinterpret_cast<accu_t *>(work);

            const int64_t mi = static_cast<int64_t>(m);
            const int64_t ni = static_cast<int64_t>(n);
            const int64_t ki = static_cast<int64_t>(k);

            // test matrices
            makemat::randmat<T>(m, k, A, phi, seedA);
            CHECK_CUDA(cudaGetLastError());
            CHECK_CUDA(cudaDeviceSynchronize());
            makemat::randmat<T>(k, n, B, phi, seedB);
            CHECK_CUDA(cudaGetLastError());
            CHECK_CUDA(cudaDeviceSynchronize());

            // high-precision AB
            eval::dd::simple_gemm(m, n, k, A, B, C_hi);
            CHECK_CUDA(cudaGetLastError());
            CHECK_CUDA(cudaDeviceSynchronize());
            CHECK_CUDA(cudaMemcpy(C_hi_h.data(), C_hi, size_C * sizeof(accu_t), cudaMemcpyDeviceToHost));

            // native gemm
            {
                CHECK_CUBLAS(gemmTraits<T>::gemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, mi, ni, ki, &alpha, A, mi, B, ki, &beta, C, mi));
                auto [err_max, err_med] = eval::err::gemm_err(m, n, C, C_hi);
                CHECK_CUDA(cudaGetLastError());
                CHECK_CUDA(cudaDeviceSynchronize());

                for (int i = 1; i < warmup; ++i) {
                    CHECK_CUBLAS(gemmTraits<T>::gemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, mi, ni, ki, &alpha, A, mi, B, ki, &beta, C, mi));
                }

                std::vector<double> res = getWatt::getWatt(
                    [&]() { CHECK_CUBLAS(gemmTraits<T>::gemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, mi, ni, ki, &alpha, A, mi, B, ki, &beta, C, mi)); },
                    m, n, k);

                outFile << phi << "," << m << "," << n << "," << k << "," << gemmTraits<T>::prefix_upper() << "GEMM" << ",";
                outFile << err_max << "," << err_med << "," << res[0] << "," << ((gemmTraits<T>::is_complex) ? 4.0 : 1.0) * res[1] * 1.e-9 << "," << std::endl;
                std::cout << phi << "," << m << "," << n << "," << k << "," << gemmTraits<T>::prefix_upper() << "GEMM" << ",";
                std::cout << err_max << "," << err_med << "," << res[0] << "," << ((gemmTraits<T>::is_complex) ? 4.0 : 1.0) * res[1] * 1.e-9 << "," << std::endl;
            }

            // fast mode
            for (unsigned num_moduli = NUM_MODULI_MIN<T>; num_moduli <= NUM_MODULI_MAX<T>; ++num_moduli) {

                gemmul8::gemm<T, backend>(handleLt, CUBLAS_OP_N, CUBLAS_OP_N, m, n, k, &alpha, A, m, B, k, &beta, C, m, num_moduli, true, work);
                CHECK_CUDA(cudaGetLastError());
                CHECK_CUDA(cudaDeviceSynchronize());
                CHECK_CUDA(cudaMemcpy(C_hi, C_hi_h.data(), size_C * sizeof(accu_t), cudaMemcpyHostToDevice));
                auto [err_max, err_med] = eval::err::gemm_err(m, n, C, C_hi);
                CHECK_CUDA(cudaGetLastError());
                CHECK_CUDA(cudaDeviceSynchronize());

                for (int i = 1; i < warmup; ++i) {
                    gemmul8::gemm<T, backend>(handleLt, CUBLAS_OP_N, CUBLAS_OP_N, m, n, k, &alpha, A, m, B, k, &beta, C, m, num_moduli, true, work);
                }

                std::vector<double> res = getWatt::getWatt(
                    [&]() { gemmul8::gemm<T, backend>(handleLt, CUBLAS_OP_N, CUBLAS_OP_N, m, n, k, &alpha, A, m, B, k, &beta, C, m, num_moduli, true, work); },
                    m, n, k);

                outFile << phi << "," << m << "," << n << "," << k << "," << "OS2-fast-" << num_moduli << ",";
                outFile << err_max << "," << err_med << "," << res[0] << "," << ((gemmTraits<T>::is_complex) ? 4.0 : 1.0) * res[1] * 1.e-9 << "," << std::endl;
                std::cout << phi << "," << m << "," << n << "," << k << "," << "OS2-fast-" << num_moduli << ",";
                std::cout << err_max << "," << err_med << "," << res[0] << "," << ((gemmTraits<T>::is_complex) ? 4.0 : 1.0) * res[1] * 1.e-9 << "," << std::endl;
            }

            // accu mode
            for (unsigned num_moduli = NUM_MODULI_MIN<T>; num_moduli <= NUM_MODULI_MAX<T>; ++num_moduli) {

                gemmul8::gemm<T, backend>(handleLt, CUBLAS_OP_N, CUBLAS_OP_N, m, n, k, &alpha, A, m, B, k, &beta, C, m, num_moduli, false, work);
                CHECK_CUDA(cudaGetLastError());
                CHECK_CUDA(cudaDeviceSynchronize());
                CHECK_CUDA(cudaMemcpy(C_hi, C_hi_h.data(), size_C * sizeof(accu_t), cudaMemcpyHostToDevice));
                auto [err_max, err_med] = eval::err::gemm_err(m, n, C, C_hi);
                CHECK_CUDA(cudaGetLastError());
                CHECK_CUDA(cudaDeviceSynchronize());

                for (int i = 1; i < warmup; ++i) {
                    gemmul8::gemm<T, backend>(handleLt, CUBLAS_OP_N, CUBLAS_OP_N, m, n, k, &alpha, A, m, B, k, &beta, C, m, num_moduli, false, work);
                }

                std::vector<double> res = getWatt::getWatt(
                    [&]() { gemmul8::gemm<T, backend>(handleLt, CUBLAS_OP_N, CUBLAS_OP_N, m, n, k, &alpha, A, m, B, k, &beta, C, m, num_moduli, false, work); },
                    m, n, k);

                outFile << phi << "," << m << "," << n << "," << k << "," << "OS2-accu-" << num_moduli << ",";
                outFile << err_max << "," << err_med << "," << res[0] << "," << ((gemmTraits<T>::is_complex) ? 4.0 : 1.0) * res[1] * 1.e-9 << "," << std::endl;
                std::cout << phi << "," << m << "," << n << "," << k << "," << "OS2-accu-" << num_moduli << ",";
                std::cout << err_max << "," << err_med << "," << res[0] << "," << ((gemmTraits<T>::is_complex) ? 4.0 : 1.0) * res[1] * 1.e-9 << "," << std::endl;
            }

            CHECK_CUDA(cudaFree(work));
            CHECK_CUDA(cudaFree(C));
            CHECK_CUDA(cudaFree(B));
            CHECK_CUDA(cudaFree(A));
        }
    }

    std::cout << std::endl;
    CHECK_CUDA(cudaEventDestroy(stop));
    CHECK_CUDA(cudaEventDestroy(start));
    CHECK_CUBLASLT(cublasLtDestroy(handleLt));
    CHECK_CUBLAS(cublasDestroy(handle));
    outFile.close();
}
