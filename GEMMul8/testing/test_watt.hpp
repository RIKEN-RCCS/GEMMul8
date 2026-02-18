#pragma once
#include "getWatt.hpp"

template <typename T>
__inline__ void watt_check(std::string &deviceName, std::string &dateTime) {
    std::string fileName = std::string("oz2_results_") + gemmTraits<T>::prefix + "_watt_" + deviceName + "_" + dateTime + ".csv";
    std::ofstream outFile(fileName);

    CHECK_CUDA(cudaSetDevice(0));
    cublasHandle_t handle;
    CHECK_CUBLAS(cublasCreate(&handle));
    cublasLtHandle_t handleLt;
    CHECK_CUBLASLT(cublasLtCreate(&handleLt));

    const size_t n_max            = max_size<T>();
    const size_t size_A           = n_max * n_max;
    const size_t size_B           = n_max * n_max;
    const size_t size_C           = n_max * n_max;
    const unsigned num_moduli_max = NUM_MODULI_MAX<T>;
    const size_t lwork            = gemmul8::workSize<gemmTraits<T>::is_complex, backend>(n_max, n_max, n_max, num_moduli_max);
    std::vector<size_t> size_list_new;
    for (auto &size : size_list) {
        if (size <= n_max) size_list_new.push_back(size);
    }
    if (size_list_new.back() != n_max) size_list_new.push_back(n_max);
    
    std::cout << std::endl;
    std::cout << "n_max = " << n_max << std::endl;
    std::cout << std::endl;

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

    outFile << std::scientific;
    std::cout << std::scientific;
    outFile << "phi,m,n,k,function,err_max,err_med,watt,GFLOPS/watt," << std::endl;
    std::cout << "phi,m,n,k,function,err_max,err_med,watt,GFLOPS/watt," << std::endl;

    cudaEvent_t start, stop;
    CHECK_CUDA(cudaEventCreate(&start));
    CHECK_CUDA(cudaEventCreate(&stop));
    std::vector<float> times(mainloop, 0.0);

    const double phi = 0.5;
    for (auto &size : size_list_new) {
        if (size > n_max) continue;

        const size_t m   = size;
        const size_t n   = size;
        const size_t k   = size;
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
    }

    std::cout << std::endl;
    CHECK_CUDA(cudaEventDestroy(stop));
    CHECK_CUDA(cudaEventDestroy(start));
    CHECK_CUDA(cudaFree(work));
    CHECK_CUDA(cudaFree(C));
    CHECK_CUDA(cudaFree(B));
    CHECK_CUDA(cudaFree(A));
    CHECK_CUBLASLT(cublasLtDestroy(handleLt));
    CHECK_CUBLAS(cublasDestroy(handle));
    outFile.close();
}
