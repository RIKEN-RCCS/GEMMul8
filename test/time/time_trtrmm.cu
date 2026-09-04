#include "../common/common.hpp"
#include "../common/ozaki1.hpp"
#include "time_trtrmm.hpp"

namespace bench::time::trtrmm {

inline constexpr gemmul8::Func func = gemmul8::Func::trtrmm;

inline cublasFillMode_t flip_uplo(const cublasFillMode_t uplo) {
    return (uplo == CUBLAS_FILL_MODE_UPPER) ? CUBLAS_FILL_MODE_LOWER : CUBLAS_FILL_MODE_UPPER;
}

inline cublasFillMode_t effective_uplo(
    const cublasFillMode_t uplo,
    const cublasOperation_t trans) {
    return (trans == CUBLAS_OP_N) ? uplo : flip_uplo(uplo);
}

inline cublasFillMode_t output_uplo(
    const cublasFillMode_t uplo_A,
    const cublasFillMode_t uplo_B,
    const cublasOperation_t trans_A,
    const cublasOperation_t trans_B) {
    const cublasFillMode_t eff_A = effective_uplo(uplo_A, trans_A);
    const cublasFillMode_t eff_B = effective_uplo(uplo_B, trans_B);
    return (eff_A == eff_B) ? eff_A : CUBLAS_FILL_MODE_FULL;
}

template <bool isOzaki2 = false, class F>
inline void evaluate_time(
    std::ofstream &outFile,
    const double computational_cost,
    const std::string &func_name,
    F &&f,
    const cudaStream_t stream //
) {
    cudaEvent_t start, stop;
    CHECK_CUDA(cudaEventCreate(&start));
    CHECK_CUDA(cudaEventCreate(&stop));
    std::vector<double> times_rep(repetitions, 0.0);
    std::vector<double> time0_rep(repetitions, 0.0);
    std::vector<double> time1_rep(repetitions, 0.0);
    std::vector<double> time2_rep(repetitions, 0.0);
    std::vector<double> time3_rep(repetitions, 0.0);

    for (unsigned rep = 0; rep < repetitions; ++rep) {
        std::vector<float> times;
        times.reserve(mainloop_max);
        float time_tmp;

        std::vector<double> time0, time1, time2, time3, timestmp(4, 0.0);
        time0.reserve(mainloop_max);
        time1.reserve(mainloop_max);
        time2.reserve(mainloop_max);
        time3.reserve(mainloop_max);

        double elapsed_time = 0.0;
        for (unsigned i = 0; i < warmup_max; ++i) {
            CHECK_CUDA(cudaStreamSynchronize(stream));
            CHECK_CUDA(cudaEventRecord(start, stream));
            f();
            CHECK_CUDA(cudaEventRecord(stop, stream));
            CHECK_CUDA(cudaEventSynchronize(stop));
            CHECK_CUDA(cudaEventElapsedTime(&time_tmp, start, stop));
            elapsed_time += double(time_tmp);
            if (i + 1 >= warmup_min && elapsed_time > warmup_ms_min) break;
        }

        elapsed_time = 0.0;
        for (unsigned i = 0; i < mainloop_max; ++i) {
            CHECK_CUDA(cudaStreamSynchronize(stream));
            CHECK_CUDA(cudaEventRecord(start, stream));
            if constexpr (isOzaki2) {
                timestmp = f();
            } else {
                f();
            }
            CHECK_CUDA(cudaEventRecord(stop, stream));
            CHECK_CUDA(cudaEventSynchronize(stop));
            CHECK_CUDA(cudaEventElapsedTime(&time_tmp, start, stop));
            times.push_back(time_tmp);
            if constexpr (isOzaki2) {
                time0.push_back(timestmp[0]);
                time1.push_back(timestmp[1]);
                time2.push_back(timestmp[2]);
                time3.push_back(timestmp[3]);
            }
            elapsed_time += double(time_tmp);
            if (i + 1 >= mainloop_min && elapsed_time > mainloop_ms_min) break;
        }

        double times_med = calc_median(times);
        times_med *= 1.e-3;
        times_rep[rep] = times_med;

        if constexpr (isOzaki2) {
            double time0_med = calc_median(time0);
            double time1_med = calc_median(time1);
            double time2_med = calc_median(time2);
            double time3_med = calc_median(time3);
            time0_rep[rep]   = time0_med;
            time1_rep[rep]   = time1_med;
            time2_rep[rep]   = time2_med;
            time3_rep[rep]   = time3_med;
        }
    }

    double times_rep_med = calc_median(times_rep);
    double TFLOPS        = computational_cost / times_rep_med * 1.0e-12;
    if constexpr (isOzaki2) {
        double time0_rep_med = calc_median(time0_rep);
        double time1_rep_med = calc_median(time1_rep);
        double time2_rep_med = calc_median(time2_rep);
        double time3_rep_med = calc_median(time3_rep);
        PRINT(outFile, func_name << "," << TFLOPS << "," << times_rep_med << "," << time0_rep_med << "," << time1_rep_med << "," << time2_rep_med << "," << time3_rep_med << ",");
    } else {
        PRINT(outFile, func_name << "," << TFLOPS << "," << times_rep_med << "," << "," << "," << "," << ",");
    }

    CHECK_CUDA(cudaEventDestroy(stop));
    CHECK_CUDA(cudaEventDestroy(start));
}

template <typename T>
void check_time(
    std::string &deviceName,
    std::string &dateTime,
    size_t memory_limit,
    cublasFillMode_t uplo_A,
    cublasFillMode_t uplo_B,
    cublasOperation_t trans_A,
    cublasOperation_t trans_B,
    cublasDiagType_t diag_A,
    cublasDiagType_t diag_B,
    const bool run_Ozaki2_I8,
    const bool run_Ozaki2_F8,
    const bool run_Ozaki1_I8 //
) {

    if (uplo_A != CUBLAS_FILL_MODE_UPPER && uplo_A != CUBLAS_FILL_MODE_LOWER) {
        assert(false && "TRTRMM requires uplo_A = UPPER or LOWER.");
        return;
    }
    if (uplo_B != CUBLAS_FILL_MODE_UPPER && uplo_B != CUBLAS_FILL_MODE_LOWER) {
        assert(false && "TRTRMM requires uplo_B = UPPER or LOWER.");
        return;
    }
    if (!testTraits<T>::is_complex &&
        (trans_A == CUBLAS_OP_C || trans_B == CUBLAS_OP_C)) {
        assert(false && "real TRTRMM accuracy test does not use trans = C.");
        return;
    }

    auto trans_tag = [](const cublasOperation_t trans) -> const char * {
        return (trans == CUBLAS_OP_N) ? "n_" : ((trans == CUBLAS_OP_T) ? "t_" : "c_");
    };

    std::string fileName = std::string("oz2_results_") +
                           testTraits<T>::prefix +
                           std::string(getFuncName<func>()) +
                           std::string("_time_") +
                           ((uplo_A == CUBLAS_FILL_MODE_UPPER) ? "upper_" : "lower_") +
                           ((uplo_B == CUBLAS_FILL_MODE_UPPER) ? "upper_" : "lower_") +
                           trans_tag(trans_A) +
                           trans_tag(trans_B) +
                           std::string("nonunit_") +
                           memTag(memory_limit) +
                           deviceName + "_" + dateTime + ".csv";

    std::ofstream outFile(fileName);

    PRINT(outFile,
          "phi,n,function,TFLOPS,total_time[sec],scaling,low_prec_gemm,mod_hi2mid,undo_scaling,");

    CHECK_CUDA(cudaSetDevice(0));

    cublasHandle_t handle;
    cublasHandle_t handle_emu;
    cublasLtHandle_t handleLt;
    cudaStream_t stream;

    CHECK_CUBLAS(cublasCreate(&handle));
    CHECK_CUBLAS(cublasCreate(&handle_emu));
    CHECK_CUBLAS(cublasLtCreate(&handleLt));
    CHECK_CUDA(cudaStreamCreate(&stream));
    CHECK_CUDA(cudaStreamSynchronize(stream));

    CHECK_CUBLAS(cublasSetStream(handle, stream));
    CHECK_CUBLAS(cublasSetStream(handle_emu, stream));
    CHECK_CUDA(cudaStreamSynchronize(stream));

    constexpr bool isDouble   = testTraits<T>::is_double;
    constexpr bool isSingle   = !testTraits<T>::is_double;
    constexpr bool use_ozaki1 = bool(avail_Ozaki1) && isDouble;
    constexpr bool use_bf16x9 = bool(avail_BF16x9) && isSingle;

    const double phi              = -1.0;
    const T alpha                 = testTraits<T>::one();
    const T beta                  = testTraits<T>::zero();
    const unsigned num_moduli_min = testTraits<T>::NUM_MODULI_MIN;
    const unsigned num_moduli_max = testTraits<T>::NUM_MODULI_MAX;

    auto free_async_if_needed = [&](auto *&ptr, const cudaStream_t stream) {
        if (ptr) {
            CHECK_CUDA(cudaFreeAsync(ptr, stream));
            ptr = nullptr;
        }
    };

    auto free_async_if_needed_void = [&](void *&ptr, const cudaStream_t stream) {
        if (ptr) {
            CHECK_CUDA(cudaFreeAsync(ptr, stream));
            ptr = nullptr;
        }
    };

    if (memory_limit > 0) {
        gemmul8::set_memory_saving(handle, true);
        gemmul8::set_max_worksize(handle, memory_limit);
        gemmul8::set_memory_savingLt(handleLt, true);
        gemmul8::set_max_worksizeLt(handleLt, memory_limit);
    }

    for (auto &n : N_list) {
        const size_t k = n;

        const int64_t ni = static_cast<int64_t>(n);
        const size_t lda = n;
        const size_t ldb = n;
        const size_t ldc = n;

        T *A            = nullptr;
        T *B            = nullptr;
        T *C            = nullptr;
        void *work_emu  = nullptr;
        void *work_blas = nullptr;

        bool run_oz2_i8 = run_Ozaki2_I8;
        bool run_oz2_f8 = run_Ozaki2_F8;
        bool run_oz1_i8 = run_Ozaki1_I8 && use_ozaki1;

        const cublasFillMode_t uplo_out = output_uplo(uplo_A, uplo_B, trans_A, trans_B);
        const double comp_cost_base     = getFuncCost<func>(n, n, n, testTraits<T>::is_complex);
        const double comp_cost          = (uplo_out == CUBLAS_FILL_MODE_FULL) ? (2.0 * comp_cost_base) : comp_cost_base;

        const size_t size_A               = n * n * sizeof(T);
        const size_t size_B               = n * n * sizeof(T);
        const size_t size_C               = n * n * sizeof(T);
        constexpr size_t lwork_blas       = size_t(32) << 20;  // 32 MiB
        constexpr size_t safety_margin    = size_t(256) << 20; // 256 MiB
        constexpr size_t alignment_margin = size_t(1) << 20;   // 1 MiB
        size_t lwork_ABC                  = lwork_blas + size_A + size_B + size_C;

        size_t free_bytes  = 0;
        size_t total_bytes = 0;
        CHECK_CUDA(cudaStreamSynchronize(stream));
        CHECK_CUDA(cudaMemGetInfo(&free_bytes, &total_bytes));

        const size_t lwork_total = lwork_ABC + safety_margin + alignment_margin;
        bool alloc_ok            = (lwork_total <= free_bytes);
        CHECK_CUDA(cudaStreamSynchronize(stream));
        if (alloc_ok) alloc_ok = alloc_ok && (cudaMallocAsync(reinterpret_cast<void **>(&A), size_A, stream) == cudaSuccess);
        CHECK_CUDA(cudaStreamSynchronize(stream));
        if (alloc_ok) alloc_ok = alloc_ok && (cudaMallocAsync(reinterpret_cast<void **>(&B), size_B, stream) == cudaSuccess);
        CHECK_CUDA(cudaStreamSynchronize(stream));
        if (alloc_ok) alloc_ok = alloc_ok && (cudaMallocAsync(reinterpret_cast<void **>(&C), size_C, stream) == cudaSuccess);
        CHECK_CUDA(cudaStreamSynchronize(stream));
        if (alloc_ok) alloc_ok = alloc_ok && (cudaMallocAsync(&work_blas, lwork_blas, stream) == cudaSuccess);
        CHECK_CUDA(cudaStreamSynchronize(stream));

        if (!alloc_ok) {
            free_async_if_needed(A, stream);
            free_async_if_needed(B, stream);
            free_async_if_needed(C, stream);
            free_async_if_needed_void(work_blas, stream);
            CHECK_CUDA(cudaStreamSynchronize(stream));
            continue;
        }

        CHECK_CUBLAS(cublasSetWorkspace(handle, work_blas, lwork_blas));
        CHECK_CUDA(cudaStreamSynchronize(stream));

        // test matrices
        makemat::randmat<T>(n, n, A, phi, seedA, stream, uplo_A, diag_A);
        makemat::randmat<T>(n, n, B, phi, seedB, stream, uplo_B, diag_B);

        eval::tri_2_full<T>(n, A, lda, stream, uplo_A, diag_A);
        eval::tri_2_full<T>(n, B, ldb, stream, uplo_B, diag_B);

        bool run_native = false;

        //-------------------------------
        // fast mode int8
        //-------------------------------
        if (run_oz2_i8) {
            for (unsigned num_moduli = num_moduli_min; num_moduli <= num_moduli_max; ++num_moduli) {

                const size_t lwork_gemmul8_i8 = (memory_limit > 0) ? memory_limit : gemmul8::workSize<testTraits<T>::is_complex, gemmul8::Backend::INT8, func>(n, n, n, num_moduli);
                bool alloc_ok                 = (lwork_total + lwork_gemmul8_i8 <= free_bytes);
                if (alloc_ok) alloc_ok = alloc_ok && (cudaMallocAsync(&work_emu, lwork_gemmul8_i8, stream) == cudaSuccess);
                CHECK_CUDA(cudaStreamSynchronize(stream));
                if (!alloc_ok) {
                    free_async_if_needed_void(work_emu, stream);
                    break;
                }

                run_native = true;

                const std::string funcname =
                    std::to_string(phi) + "," +
                    std::to_string(n) + "," +
                    std::string("OS2-i8-fast-") + std::to_string(num_moduli);

                constexpr bool isOzaki2            = true;
                constexpr bool fastmode            = true;
                constexpr gemmul8::Backend backend = gemmul8::Backend::INT8;

                evaluate_time<isOzaki2>(
                    outFile,
                    comp_cost,
                    funcname,
                    [&]() { return gemmul8::trtrmm<T, backend, T, T>(
                                handle,
                                uplo_A, uplo_B,
                                trans_A, trans_B,
                                diag_A, diag_B,
                                n,
                                &alpha,
                                A, lda,
                                B, ldb,
                                &beta,
                                C, ldc,
                                num_moduli, fastmode,
                                work_emu,
                                nullptr, nullptr,
                                false, false,
                                false, false); },
                    stream);

                CHECK_CUDA(cudaStreamSynchronize(stream));
                free_async_if_needed_void(work_emu, stream);
            }
        }

        //-------------------------------
        // accu mode int8
        //-------------------------------
        if (run_oz2_i8) {
            for (unsigned num_moduli = num_moduli_min; num_moduli <= num_moduli_max; ++num_moduli) {

                const size_t lwork_gemmul8_i8 = (memory_limit > 0) ? memory_limit : gemmul8::workSize<testTraits<T>::is_complex, gemmul8::Backend::INT8, func>(n, n, n, num_moduli);
                bool alloc_ok                 = (lwork_total + lwork_gemmul8_i8 <= free_bytes);
                if (alloc_ok) alloc_ok = alloc_ok && (cudaMallocAsync(&work_emu, lwork_gemmul8_i8, stream) == cudaSuccess);
                CHECK_CUDA(cudaStreamSynchronize(stream));
                if (!alloc_ok) {
                    free_async_if_needed_void(work_emu, stream);
                    break;
                }

                run_native = true;

                const std::string funcname =
                    std::to_string(phi) + "," +
                    std::to_string(n) + "," +
                    std::string("OS2-i8-accu-") + std::to_string(num_moduli);

                constexpr bool isOzaki2            = true;
                constexpr bool fastmode            = false;
                constexpr gemmul8::Backend backend = gemmul8::Backend::INT8;

                evaluate_time<isOzaki2>(
                    outFile,
                    comp_cost,
                    funcname,
                    [&]() { return gemmul8::trtrmm<T, backend, T, T>(
                                handle,
                                uplo_A, uplo_B,
                                trans_A, trans_B,
                                diag_A, diag_B,
                                n,
                                &alpha,
                                A, lda,
                                B, ldb,
                                &beta,
                                C, ldc,
                                num_moduli, fastmode,
                                work_emu,
                                nullptr, nullptr,
                                false, false,
                                false, false); },
                    stream);

                CHECK_CUDA(cudaStreamSynchronize(stream));
                free_async_if_needed_void(work_emu, stream);
            }
        }

        //-------------------------------
        // fast mode fp8
        //-------------------------------
        if (run_oz2_f8) {
            for (unsigned num_moduli = num_moduli_min; num_moduli <= num_moduli_max; ++num_moduli) {

                const size_t lwork_gemmul8_f8 = (memory_limit > 0) ? memory_limit : gemmul8::workSize<testTraits<T>::is_complex, gemmul8::Backend::FP8, func>(n, n, n, num_moduli);
                bool alloc_ok                 = (lwork_total + lwork_gemmul8_f8 <= free_bytes);
                if (alloc_ok) alloc_ok = alloc_ok && (cudaMallocAsync(&work_emu, lwork_gemmul8_f8, stream) == cudaSuccess);
                CHECK_CUDA(cudaStreamSynchronize(stream));
                if (!alloc_ok) {
                    free_async_if_needed_void(work_emu, stream);
                    break;
                }

                run_native = true;

                const std::string funcname =
                    std::to_string(phi) + "," +
                    std::to_string(n) + "," +
                    std::string("OS2-f8-fast-") + std::to_string(num_moduli);

                constexpr bool isOzaki2            = true;
                constexpr bool fastmode            = true;
                constexpr gemmul8::Backend backend = gemmul8::Backend::FP8;

                evaluate_time<isOzaki2>(
                    outFile,
                    comp_cost,
                    funcname,
                    [&]() { return gemmul8::trtrmmLt<T, backend, T, T>(
                                handleLt,
                                uplo_A, uplo_B,
                                trans_A, trans_B,
                                diag_A, diag_B,
                                n,
                                &alpha,
                                A, lda,
                                B, ldb,
                                &beta,
                                C, ldc,
                                num_moduli, fastmode,
                                work_emu,
                                nullptr, nullptr,
                                false, false,
                                false, false,
                                stream); },
                    stream);

                CHECK_CUDA(cudaStreamSynchronize(stream));
                free_async_if_needed_void(work_emu, stream);
            }
        }

        //-------------------------------
        // accu mode fp8
        //-------------------------------
        if (run_oz2_f8) {
            for (unsigned num_moduli = num_moduli_min; num_moduli <= num_moduli_max; ++num_moduli) {

                const size_t lwork_gemmul8_f8 = (memory_limit > 0) ? memory_limit : gemmul8::workSize<testTraits<T>::is_complex, gemmul8::Backend::FP8, func>(n, n, n, num_moduli);
                bool alloc_ok                 = (lwork_total + lwork_gemmul8_f8 <= free_bytes);
                if (alloc_ok) alloc_ok = alloc_ok && (cudaMallocAsync(&work_emu, lwork_gemmul8_f8, stream) == cudaSuccess);
                CHECK_CUDA(cudaStreamSynchronize(stream));
                if (!alloc_ok) {
                    free_async_if_needed_void(work_emu, stream);
                    break;
                }

                run_native = true;

                const std::string funcname =
                    std::to_string(phi) + "," +
                    std::to_string(n) + "," +
                    std::string("OS2-f8-accu-") + std::to_string(num_moduli);

                constexpr bool isOzaki2            = true;
                constexpr bool fastmode            = false;
                constexpr gemmul8::Backend backend = gemmul8::Backend::FP8;

                evaluate_time<isOzaki2>(
                    outFile,
                    comp_cost,
                    funcname,
                    [&]() { return gemmul8::trtrmmLt<T, backend, T, T>(
                                handleLt,
                                uplo_A, uplo_B,
                                trans_A, trans_B,
                                diag_A, diag_B,
                                n,
                                &alpha,
                                A, lda,
                                B, ldb,
                                &beta,
                                C, ldc,
                                num_moduli, fastmode,
                                work_emu,
                                nullptr, nullptr,
                                false, false,
                                false, false,
                                stream); },
                    stream);

                CHECK_CUDA(cudaStreamSynchronize(stream));
                free_async_if_needed_void(work_emu, stream);
            }
        }

        //-------------------------------
        // cuBLAS Ozaki-I
        //-------------------------------
#if avail_Ozaki1
        if (run_oz1_i8) {
            for (auto &num_slice : oz1_slice_list) {

                int mantissaBitCount      = num_slice * 8 - 1;
                const size_t lwork_ozaki1 = ozaki1::workSize(n, n, n, 1, testTraits<T>::is_complex, mantissaBitCount);
                bool alloc_ok             = (lwork_total + lwork_ozaki1 <= free_bytes);
                if (alloc_ok) alloc_ok = alloc_ok && (cudaMallocAsync(&work_emu, lwork_ozaki1, stream) == cudaSuccess);
                CHECK_CUDA(cudaStreamSynchronize(stream));
                if (!alloc_ok) {
                    free_async_if_needed_void(work_emu, stream);
                    break;
                }

                run_native = true;

                CHECK_CUBLAS(cublasSetMathMode(handle_emu, CUBLAS_FP64_EMULATED_FIXEDPOINT_MATH));
                CHECK_CUBLAS(cublasSetEmulationStrategy(handle_emu, CUBLAS_EMULATION_STRATEGY_EAGER));
                CHECK_CUBLAS(cublasSetFixedPointEmulationMantissaControl(handle_emu, CUDA_EMULATION_MANTISSA_CONTROL_FIXED));
                CHECK_CUBLAS(cublasSetFixedPointEmulationMaxMantissaBitCount(handle_emu, mantissaBitCount));
                CHECK_CUBLAS(cublasSetWorkspace(handle_emu, work_emu, lwork_ozaki1));

                const std::string funcname =
                    std::to_string(phi) + "," +
                    std::to_string(n) + "," +
                    std::string("OS1-") + std::to_string(num_slice);

                constexpr bool isOzaki2 = false;

                evaluate_time<isOzaki2>(
                    outFile,
                    comp_cost,
                    funcname,
                    [&]() { CHECK_CUBLAS(testTraits<T>::trtrmm(
                                handle_emu,
                                trans_A, trans_B,
                                ni, ni, ni,
                                &alpha,
                                A, static_cast<int64_t>(n),
                                B, static_cast<int64_t>(n),
                                &beta,
                                C, static_cast<int64_t>(ldc))); },
                    stream);

                CHECK_CUDA(cudaStreamSynchronize(stream));
                CHECK_CUBLAS(cublasSetWorkspace(handle_emu, nullptr, 0));
                CHECK_CUBLAS(cublasSetMathMode(handle_emu, CUBLAS_DEFAULT_MATH));
                free_async_if_needed_void(work_emu, stream);
            }
        }
#endif

        //-------------------------------
        // cuBLAS BF16x9
        //-------------------------------
#if avail_BF16x9
        if (use_bf16x9) {
            CHECK_CUBLAS(cublasSetMathMode(handle_emu, CUBLAS_FP32_EMULATED_BF16X9_MATH));
            CHECK_CUBLAS(cublasSetWorkspace(handle_emu, work_blas, lwork_blas));

            run_native = true;

            const std::string funcname =
                std::to_string(phi) + "," +
                std::to_string(n) + "," +
                std::string("BF16x9");

            constexpr bool isOzaki2 = false;

            evaluate_time<isOzaki2>(
                outFile,
                comp_cost,
                funcname,
                [&]() { CHECK_CUBLAS(testTraits<T>::trtrmm(
                            handle_emu,
                            trans_A, trans_B,
                            ni, ni, ni,
                            &alpha,
                            A, static_cast<int64_t>(n),
                            B, static_cast<int64_t>(n),
                            &beta,
                            C, static_cast<int64_t>(ldc))); },
                stream);

            CHECK_CUDA(cudaStreamSynchronize(stream));
            CHECK_CUBLAS(cublasSetWorkspace(handle_emu, nullptr, 0));
            CHECK_CUBLAS(cublasSetMathMode(handle_emu, CUBLAS_DEFAULT_MATH));
        }
#endif

        if (!run_native) {
            CHECK_CUDA(cudaStreamSynchronize(stream));
            free_async_if_needed(A, stream);
            free_async_if_needed(B, stream);
            free_async_if_needed(C, stream);
            free_async_if_needed_void(work_blas, stream);
            CHECK_CUDA(cudaStreamSynchronize(stream));
            continue;
        }

        //-------------------------------
        // native routine
        //-------------------------------
        if constexpr (isDouble || isSingle) {

            const std::string funcname =
                std::to_string(phi) + "," +
                std::to_string(n) + "," +
                std::string(1, testTraits<T>::prefix_upper) + std::string(getFuncName_upper<func>());

            constexpr bool isOzaki2 = false;

            evaluate_time<isOzaki2>(
                outFile,
                comp_cost,
                funcname,
                [&]() { CHECK_CUBLAS(testTraits<T>::trtrmm(
                            handle,
                            trans_A, trans_B,
                            ni, ni, ni,
                            &alpha,
                            A, static_cast<int64_t>(n),
                            B, static_cast<int64_t>(n),
                            &beta,
                            C, static_cast<int64_t>(ldc))); },
                stream);
            CHECK_CUDA(cudaStreamSynchronize(stream));
        }

        CHECK_CUDA(cudaStreamSynchronize(stream));
        free_async_if_needed(A, stream);
        free_async_if_needed(B, stream);
        free_async_if_needed(C, stream);
        free_async_if_needed_void(work_blas, stream);
        CHECK_CUDA(cudaStreamSynchronize(stream));
    }

    std::cout << std::endl;
    CHECK_CUDA(cudaStreamSynchronize(stream));
    if (memory_limit > 0) {
        gemmul8::clear_config(handle);
        gemmul8::clear_configLt(handleLt);
    }
    CHECK_CUBLAS(cublasLtDestroy(handleLt));
    CHECK_CUBLAS(cublasDestroy(handle_emu));
    CHECK_CUBLAS(cublasDestroy(handle));
    cudaStreamDestroy(stream);
    outFile.close();
}

template void check_time<float>(
    std::string &, std::string &, size_t, cublasFillMode_t, cublasFillMode_t,
    cublasOperation_t, cublasOperation_t, cublasDiagType_t, cublasDiagType_t,
    const bool, const bool, const bool);

template void check_time<double>(
    std::string &, std::string &, size_t, cublasFillMode_t, cublasFillMode_t,
    cublasOperation_t, cublasOperation_t, cublasDiagType_t, cublasDiagType_t,
    const bool, const bool, const bool);

template void check_time<cuFloatComplex>(
    std::string &, std::string &, size_t, cublasFillMode_t, cublasFillMode_t,
    cublasOperation_t, cublasOperation_t, cublasDiagType_t, cublasDiagType_t,
    const bool, const bool, const bool);

template void check_time<cuDoubleComplex>(
    std::string &, std::string &, size_t, cublasFillMode_t, cublasFillMode_t,
    cublasOperation_t, cublasOperation_t, cublasDiagType_t, cublasDiagType_t,
    const bool, const bool, const bool);

} // namespace bench::time::trtrmm
