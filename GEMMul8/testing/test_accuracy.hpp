#pragma once

template <typename T>
__inline__ void accuracy_check(std::string &deviceName, std::string &dateTime) {
    std::string fileName = std::string("oz2_results_") + gemmTraits<T>::prefix + "_accuracy_" + deviceName + "_" + dateTime + ".csv";
    std::ofstream outFile(fileName);

    CHECK_CUDA(cudaSetDevice(0));
    cublasHandle_t handle;
    CHECK_CUBLAS(cublasCreate(&handle));
    cublasLtHandle_t handleLt;
    CHECK_CUBLASLT(cublasLtCreate(&handleLt));

    outFile << std::scientific;
    std::cout << std::scientific;

    const size_t m = 128;
    const size_t n = 128;
    const size_t k = *max_element(begin(size_list), end(size_list));

    const int64_t mi = static_cast<int64_t>(m);
    const int64_t ni = static_cast<int64_t>(n);
    const int64_t ki = static_cast<int64_t>(k);

    const size_t size_A           = m * k;
    const size_t size_B           = k * n;
    const size_t size_C           = m * n;
    const unsigned num_moduli_max = NUM_MODULI_MAX<T>;
    const size_t lwork            = gemmul8::workSize<gemmTraits<T>::is_complex, backend>(m, n, k, num_moduli_max);

    using accu_t = typename gemmTraits<T>::ACCU_TYPE;
    T *A, *B, *C;
    accu_t *C_hi;
    std::vector<accu_t> C_hi_h(size_C);
    void *work;
    const T alpha = gemmTraits<T>::one();
    const T beta  = gemmTraits<T>::zero();

    CHECK_CUDA(cudaMalloc(reinterpret_cast<void **>(&A), size_A * sizeof(T)));
    CHECK_CUDA(cudaMalloc(reinterpret_cast<void **>(&B), size_B * sizeof(T)));
    CHECK_CUDA(cudaMalloc(reinterpret_cast<void **>(&C), size_C * sizeof(T)));
    CHECK_CUDA(cudaMalloc(&work, std::max(lwork, size_C * sizeof(accu_t))));
    C_hi = reinterpret_cast<accu_t *>(work);

    outFile << "phi,function,";
    std::cout << "phi,function,";
    for (unsigned num_moduli = NUM_MODULI_MIN<T>; num_moduli <= NUM_MODULI_MAX<T>; ++num_moduli) {
        outFile << num_moduli << ",";
        std::cout << num_moduli << ",";
    }
    outFile << std::endl;
    std::cout << std::endl;

    for (auto &phi : phi_list) {
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
            outFile << phi << "," << gemmTraits<T>::prefix_upper() << "GEMM (k=" << k << "),";
            std::cout << phi << "," << gemmTraits<T>::prefix_upper() << "GEMM (k=" << k << "),";
            for (unsigned num_moduli = NUM_MODULI_MIN<T>; num_moduli <= NUM_MODULI_MAX<T>; ++num_moduli) {
                outFile << err_max << ",";
                std::cout << err_max << ",";
            }
            outFile << std::endl;
            std::cout << std::endl;
        }

        // fast mode
        outFile << phi << ",OS2-fast (k=" << k << "),";
        std::cout << phi << ",OS2-fast (k=" << k << "),";
        for (unsigned num_moduli = NUM_MODULI_MIN<T>; num_moduli <= NUM_MODULI_MAX<T>; ++num_moduli) {
            gemmul8::gemm<T, backend>(handleLt, CUBLAS_OP_N, CUBLAS_OP_N, m, n, k, &alpha, A, m, B, k, &beta, C, m, num_moduli, true, work);
            CHECK_CUDA(cudaGetLastError());
            CHECK_CUDA(cudaDeviceSynchronize());
            CHECK_CUDA(cudaMemcpy(C_hi, C_hi_h.data(), size_C * sizeof(accu_t), cudaMemcpyHostToDevice));
            auto [err_max, err_med] = eval::err::gemm_err(m, n, C, C_hi);
            CHECK_CUDA(cudaGetLastError());
            CHECK_CUDA(cudaDeviceSynchronize());
            outFile << err_max << ",";
            std::cout << err_max << ",";
        }
        outFile << std::endl;
        std::cout << std::endl;

        // accu mode
        outFile << phi << ",OS2-accu (k=" << k << "),";
        std::cout << phi << ",OS2-accu (k=" << k << "),";
        for (unsigned num_moduli = NUM_MODULI_MIN<T>; num_moduli <= NUM_MODULI_MAX<T>; ++num_moduli) {
            gemmul8::gemm<T, backend>(handleLt, CUBLAS_OP_N, CUBLAS_OP_N, m, n, k, &alpha, A, m, B, k, &beta, C, m, num_moduli, false, work);
            CHECK_CUDA(cudaGetLastError());
            CHECK_CUDA(cudaDeviceSynchronize());
            CHECK_CUDA(cudaMemcpy(C_hi, C_hi_h.data(), size_C * sizeof(accu_t), cudaMemcpyHostToDevice));
            auto [err_max, err_med] = eval::err::gemm_err(m, n, C, C_hi);
            CHECK_CUDA(cudaGetLastError());
            CHECK_CUDA(cudaDeviceSynchronize());
            outFile << err_max << ",";
            std::cout << err_max << ",";
        }
        outFile << std::endl;
        std::cout << std::endl;
    }

    std::cout << std::endl;
    CHECK_CUDA(cudaFree(work));
    CHECK_CUDA(cudaFree(C));
    CHECK_CUDA(cudaFree(B));
    CHECK_CUDA(cudaFree(A));
    CHECK_CUBLASLT(cublasLtDestroy(handleLt));
    CHECK_CUBLAS(cublasDestroy(handle));
    outFile.close();
}
