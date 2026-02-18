#pragma once

template <typename T>
__inline__ void time_check(std::string &deviceName, std::string &dateTime) {
    std::string fileName = std::string("oz2_results_") + gemmTraits<T>::prefix + "_time_" + deviceName + "_" + dateTime + ".csv";
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
    outFile << "phi,m,n,k,function,err_max,err_med,TFLOPS,total_time[sec],quantization,low_prec_gemm,requantization,dequantization," << std::endl;
    std::cout << "phi,m,n,k,function,err_max,err_med,TFLOPS,total_time[sec],quantization,low_prec_gemm,requantization,dequantization," << std::endl;

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

            for (int i = 0; i < mainloop; ++i) {
                CHECK_CUDA(cudaEventRecord(start));
                CHECK_CUBLAS(gemmTraits<T>::gemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, mi, ni, ki, &alpha, A, mi, B, ki, &beta, C, mi));
                CHECK_CUDA(cudaEventRecord(stop));
                CHECK_CUDA(cudaEventSynchronize(stop));
                CHECK_CUDA(cudaEventElapsedTime(&times[i], start, stop));
            }

            sort(times.begin(), times.end());
            double time_med = (mainloop & 1) ? double(times[mainloop / 2]) : ((double(times[mainloop / 2]) + double(times[mainloop / 2 - 1])) * 0.5);
            time_med *= 1.e-3;
            double TFLOPS = 2.0 * m * n * k * ((gemmTraits<T>::is_complex) ? 4.0 : 1.0) / time_med * 1.0e-12;

            outFile << phi << "," << m << "," << n << "," << k << "," << gemmTraits<T>::prefix_upper() << "GEMM" << ",";
            outFile << err_max << "," << err_med << "," << TFLOPS << "," << time_med << "," << "," << "," << "," << "," << std::endl;
            std::cout << phi << "," << m << "," << n << "," << k << "," << gemmTraits<T>::prefix_upper() << "GEMM" << ",";
            std::cout << err_max << "," << err_med << "," << TFLOPS << "," << time_med << "," << "," << "," << "," << "," << std::endl;
        }

        // fast mode
        for (unsigned num_moduli = NUM_MODULI_MIN<T>; num_moduli <= NUM_MODULI_MAX<T>; ++num_moduli) {

            std::vector<double> time0(mainloop, 0.0), time1(mainloop, 0.0), time2(mainloop, 0.0), time3(mainloop, 0.0), timestmp(4, 0.0);
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

            for (int i = 0; i < mainloop; ++i) {
                CHECK_CUDA(cudaEventRecord(start));
                timestmp = gemmul8::gemm<T, backend>(handleLt, CUBLAS_OP_N, CUBLAS_OP_N, m, n, k, &alpha, A, m, B, k, &beta, C, m, num_moduli, true, work);
                CHECK_CUDA(cudaEventRecord(stop));
                CHECK_CUDA(cudaEventSynchronize(stop));
                CHECK_CUDA(cudaEventElapsedTime(&times[i], start, stop));
                time0[i] = timestmp[0];
                time1[i] = timestmp[1];
                time2[i] = timestmp[2];
                time3[i] = timestmp[3];
            }

            sort(times.begin(), times.end());
            sort(time0.begin(), time0.end());
            sort(time1.begin(), time1.end());
            sort(time2.begin(), time2.end());
            sort(time3.begin(), time3.end());
            double time_med  = (mainloop & 1) ? double(times[mainloop / 2]) : ((double(times[mainloop / 2]) + double(times[mainloop / 2 - 1])) * 0.5);
            double time0_med = (mainloop & 1) ? double(time0[mainloop / 2]) : ((double(time0[mainloop / 2]) + double(time0[mainloop / 2 - 1])) * 0.5);
            double time1_med = (mainloop & 1) ? double(time1[mainloop / 2]) : ((double(time1[mainloop / 2]) + double(time1[mainloop / 2 - 1])) * 0.5);
            double time2_med = (mainloop & 1) ? double(time2[mainloop / 2]) : ((double(time2[mainloop / 2]) + double(time2[mainloop / 2 - 1])) * 0.5);
            double time3_med = (mainloop & 1) ? double(time3[mainloop / 2]) : ((double(time3[mainloop / 2]) + double(time3[mainloop / 2 - 1])) * 0.5);
            time_med *= 1.e-3;
            time0_med *= 1.e-9;
            time1_med *= 1.e-9;
            time2_med *= 1.e-9;
            time3_med *= 1.e-9;
            double TFLOPS = 2.0 * m * n * k * ((gemmTraits<T>::is_complex) ? 4.0 : 1.0) / time_med * 1.0e-12;

            outFile << phi << "," << m << "," << n << "," << k << "," << "OS2-fast-" << num_moduli << ",";
            outFile << err_max << "," << err_med << "," << TFLOPS << "," << time_med << ",";
            outFile << time0_med << "," << time1_med << "," << time2_med << "," << time3_med << "," << std::endl;
            std::cout << phi << "," << m << "," << n << "," << k << "," << "OS2-fast-" << num_moduli << ",";
            std::cout << err_max << "," << err_med << "," << TFLOPS << "," << time_med << ",";
            std::cout << time0_med << "," << time1_med << "," << time2_med << "," << time3_med << "," << std::endl;
        }

        // accu mode
        for (unsigned num_moduli = NUM_MODULI_MIN<T>; num_moduli <= NUM_MODULI_MAX<T>; ++num_moduli) {

            std::vector<double> time0(mainloop, 0.0), time1(mainloop, 0.0), time2(mainloop, 0.0), time3(mainloop, 0.0), timestmp(4, 0.0);
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

            for (int i = 0; i < mainloop; ++i) {
                CHECK_CUDA(cudaEventRecord(start));
                timestmp = gemmul8::gemm<T, backend>(handleLt, CUBLAS_OP_N, CUBLAS_OP_N, m, n, k, &alpha, A, m, B, k, &beta, C, m, num_moduli, false, work);
                CHECK_CUDA(cudaEventRecord(stop));
                CHECK_CUDA(cudaEventSynchronize(stop));
                CHECK_CUDA(cudaEventElapsedTime(&times[i], start, stop));
                time0[i] = timestmp[0];
                time1[i] = timestmp[1];
                time2[i] = timestmp[2];
                time3[i] = timestmp[3];
            }

            sort(times.begin(), times.end());
            sort(time0.begin(), time0.end());
            sort(time1.begin(), time1.end());
            sort(time2.begin(), time2.end());
            sort(time3.begin(), time3.end());
            double time_med  = (mainloop & 1) ? double(times[mainloop / 2]) : ((double(times[mainloop / 2]) + double(times[mainloop / 2 - 1])) * 0.5);
            double time0_med = (mainloop & 1) ? double(time0[mainloop / 2]) : ((double(time0[mainloop / 2]) + double(time0[mainloop / 2 - 1])) * 0.5);
            double time1_med = (mainloop & 1) ? double(time1[mainloop / 2]) : ((double(time1[mainloop / 2]) + double(time1[mainloop / 2 - 1])) * 0.5);
            double time2_med = (mainloop & 1) ? double(time2[mainloop / 2]) : ((double(time2[mainloop / 2]) + double(time2[mainloop / 2 - 1])) * 0.5);
            double time3_med = (mainloop & 1) ? double(time3[mainloop / 2]) : ((double(time3[mainloop / 2]) + double(time3[mainloop / 2 - 1])) * 0.5);
            time_med *= 1.e-3;
            time0_med *= 1.e-9;
            time1_med *= 1.e-9;
            time2_med *= 1.e-9;
            time3_med *= 1.e-9;
            double TFLOPS = 2.0 * m * n * k * ((gemmTraits<T>::is_complex) ? 4.0 : 1.0) / time_med * 1.0e-12;

            outFile << phi << "," << m << "," << n << "," << k << "," << "OS2-accu-" << num_moduli << ",";
            outFile << err_max << "," << err_med << "," << TFLOPS << "," << time_med << ",";
            outFile << time0_med << "," << time1_med << "," << time2_med << "," << time3_med << "," << std::endl;
            std::cout << phi << "," << m << "," << n << "," << k << "," << "OS2-accu-" << num_moduli << ",";
            std::cout << err_max << "," << err_med << "," << TFLOPS << "," << time_med << ",";
            std::cout << time0_med << "," << time1_med << "," << time2_med << "," << time3_med << "," << std::endl;
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
