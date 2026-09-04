#include "../common/common.hpp"
#include "../common/ozaki1.hpp"
#include "accuracy_trtrmm.hpp"

namespace bench::accuracy::trtrmm {

inline constexpr gemmul8::Func func = gemmul8::Func::trtrmm;

inline cublasFillMode_t flip_uplo(const cublasFillMode_t uplo) {
    return (uplo == CUBLAS_FILL_MODE_UPPER) ? CUBLAS_FILL_MODE_LOWER : CUBLAS_FILL_MODE_UPPER;
}

inline cublasFillMode_t effective_uplo(
    const cublasFillMode_t uplo,
    const cublasOperation_t trans //
) {
    return (trans == CUBLAS_OP_N) ? uplo : flip_uplo(uplo);
}

inline cublasFillMode_t output_uplo(
    const cublasFillMode_t uplo_A,
    const cublasFillMode_t uplo_B,
    const cublasOperation_t trans_A,
    const cublasOperation_t trans_B //
) {
    const cublasFillMode_t eff_A = effective_uplo(uplo_A, trans_A);
    const cublasFillMode_t eff_B = effective_uplo(uplo_B, trans_B);
    return (eff_A == eff_B) ? eff_A : CUBLAS_FILL_MODE_FULL;
}

template <bool isOzaki2 = false, typename T, typename accu_t, class F>
inline double evaluate_error(
    std::ofstream &outFile,
    F &&f,
    const cublasFillMode_t uplo_out,
    const size_t n,
    T *const C, const size_t ldc,
    const accu_t *const C_exact, const size_t ldc_exact,
    const cudaStream_t stream //
) {
    CHECK_CUDA(cudaStreamSynchronize(stream));
    CHECK_CUDA(cudaMemsetAsync(C, 0, ldc * n * sizeof(T), stream));
    CHECK_CUDA(cudaStreamSynchronize(stream));

    f();

    CHECK_CUDA(cudaStreamSynchronize(stream));

    double2 err;
    if (uplo_out == CUBLAS_FILL_MODE_UPPER) {
        err = eval::calc_err<T, CUBLAS_FILL_MODE_UPPER>(
            n, n, C, ldc, C_exact, ldc_exact, stream);
    } else if (uplo_out == CUBLAS_FILL_MODE_LOWER) {
        err = eval::calc_err<T, CUBLAS_FILL_MODE_LOWER>(
            n, n, C, ldc, C_exact, ldc_exact, stream);
    } else {
        err = eval::calc_err<T, CUBLAS_FILL_MODE_FULL>(
            n, n, C, ldc, C_exact, ldc_exact, stream);
    }

    PRINT_nobreak(outFile, err.x << ",");
    return err.x;
}

template <typename T>
void check_accuracy(
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

    const cublasFillMode_t uplo_out = output_uplo(uplo_A, uplo_B, trans_A, trans_B);

    auto trans_tag = [](const cublasOperation_t trans) -> const char * {
        return (trans == CUBLAS_OP_N) ? "n_" : ((trans == CUBLAS_OP_T) ? "t_" : "c_");
    };

    std::string fileName = std::string("oz2_results_") +
                           testTraits<T>::prefix +
                           std::string(getFuncName<func>()) +
                           std::string("_accuracy_") +
                           ((uplo_A == CUBLAS_FILL_MODE_UPPER) ? "upper_" : "lower_") +
                           ((uplo_B == CUBLAS_FILL_MODE_UPPER) ? "upper_" : "lower_") +
                           trans_tag(trans_A) +
                           trans_tag(trans_B) +
                           std::string("nonunit_") +
                           memTag(memory_limit) +
                           deviceName + "_" + dateTime + ".csv";

    std::ofstream outFile(fileName);

    const unsigned num_moduli_min = testTraits<T>::NUM_MODULI_MIN_accuracy;
    const unsigned num_moduli_max = testTraits<T>::NUM_MODULI_MAX_accuracy;

    std::string num_moduli_str = std::string("");
    for (unsigned num_moduli = num_moduli_min; num_moduli <= num_moduli_max; ++num_moduli) {
        num_moduli_str += std::to_string(num_moduli) + ",";
    }

    PRINT(outFile, "phi,n,function," + num_moduli_str);

    CHECK_CUDA(cudaSetDevice(0));

    cublasHandle_t handle;
    cublasHandle_t handle_emu;
    cublasLtHandle_t handleLt;
    cudaStream_t stream;

    CHECK_CUBLAS(cublasCreate(&handle));
    CHECK_CUBLAS(cublasCreate(&handle_emu));
    CHECK_CUBLAS(cublasLtCreate(&handleLt));
    CHECK_CUDA(cudaStreamCreate(&stream));

    CHECK_CUBLAS(cublasSetStream(handle, stream));
    CHECK_CUBLAS(cublasSetStream(handle_emu, stream));
    CHECK_CUDA(cudaStreamSynchronize(stream));

    constexpr bool isDouble   = testTraits<T>::is_double;
    constexpr bool isSingle   = !testTraits<T>::is_double;
    constexpr bool use_ozaki1 = bool(avail_Ozaki1) && isDouble;
    constexpr bool use_bf16x9 = bool(avail_BF16x9) && isSingle;

    const T alpha = testTraits<T>::one();
    const T beta  = testTraits<T>::zero();

    const int num_slice_max = *std::max_element(begin(oz1_slice_list), end(oz1_slice_list));

    using accu_t = typename testTraits<T>::accu_t;

    T *A         = nullptr;
    T *B         = nullptr;
    T *A_full    = nullptr;
    T *B_full    = nullptr;
    T *C         = nullptr;
    accu_t *C_hi = nullptr;

    void *work_emu  = nullptr;
    void *work_blas = nullptr;

    bool run_oz2_i8 = run_Ozaki2_I8;
    bool run_oz2_f8 = run_Ozaki2_F8;
    bool run_oz1_i8 = run_Ozaki1_I8 && use_ozaki1;

    std::vector<size_t> K_list;
    for (auto &n : N_list) {
        if (n <= 8192) {
            K_list.push_back(n);
        }
    }

    const size_t n_max  = *std::max_element(begin(K_list), end(K_list));
    const size_t size_A = n_max * n_max;
    const size_t size_B = n_max * n_max;
    const size_t size_C = n_max * n_max;

    const size_t lwork_blas = size_t(32) << 20;

    if (memory_limit > 0) {
        gemmul8::set_memory_saving(handle, true);
        gemmul8::set_max_worksize(handle, memory_limit);
        gemmul8::set_memory_savingLt(handleLt, true);
        gemmul8::set_max_worksizeLt(handleLt, memory_limit);
    }
    const size_t lwork_gemmul8_i8_ = gemmul8::workSize<testTraits<T>::is_complex, gemmul8::Backend::INT8, func>(n_max, n_max, n_max, num_moduli_max);
    const size_t lwork_gemmul8_f8_ = gemmul8::workSize<testTraits<T>::is_complex, gemmul8::Backend::FP8, func>(n_max, n_max, n_max, num_moduli_max);
    const size_t lwork_ozaki1_     = ozaki1::workSize(n_max, n_max, n_max, 1, testTraits<T>::is_complex, 8 * num_slice_max - 1);
    const size_t lwork_gemmul8_i8  = (run_oz2_i8) ? std::max(memory_limit, lwork_gemmul8_i8_) : 0;
    const size_t lwork_gemmul8_f8  = (run_oz2_f8) ? std::max(memory_limit, lwork_gemmul8_f8_) : 0;
    const size_t lwork_ozaki1      = (run_oz1_i8) ? lwork_ozaki1_ : 0;
    const size_t lwork_emu         = std::max(std::max(lwork_gemmul8_i8, lwork_gemmul8_f8), lwork_ozaki1);

    CHECK_CUDA(cudaMallocAsync(reinterpret_cast<void **>(&A), size_A * sizeof(T), stream));
    CHECK_CUDA(cudaMallocAsync(reinterpret_cast<void **>(&B), size_B * sizeof(T), stream));
    CHECK_CUDA(cudaMallocAsync(reinterpret_cast<void **>(&A_full), size_A * sizeof(T), stream));
    CHECK_CUDA(cudaMallocAsync(reinterpret_cast<void **>(&B_full), size_B * sizeof(T), stream));
    CHECK_CUDA(cudaMallocAsync(reinterpret_cast<void **>(&C), size_C * sizeof(T), stream));
    CHECK_CUDA(cudaMallocAsync(reinterpret_cast<void **>(&C_hi), size_C * sizeof(accu_t), stream));

    CHECK_CUDA(cudaMallocAsync(&work_emu, lwork_emu, stream));
    CHECK_CUDA(cudaMallocAsync(&work_blas, lwork_blas, stream));
    CHECK_CUDA(cudaStreamSynchronize(stream));

    CHECK_CUBLAS(cublasSetWorkspace(handle, work_blas, lwork_blas));
    CHECK_CUDA(cudaStreamSynchronize(stream));

    for (auto &phi : phi_list) {
        std::vector<size_t> N_test_list;
        N_test_list = K_list;
        for (auto &n : N_test_list) {

            const int64_t ni = static_cast<int64_t>(n);
            const size_t lda = n;
            const size_t ldb = n;
            const size_t ldc = n;

            makemat::randmat<T>(n, n, A, phi, seedA, stream, uplo_A, diag_A);
            makemat::randmat<T>(n, n, B, phi, seedB, stream, uplo_B, diag_B);

            eval::tri_2_full<T>(n, A, lda, A_full, n, stream, uplo_A, diag_A);
            eval::tri_2_full<T>(n, B, ldb, B_full, n, stream, uplo_B, diag_B);

            eval::DDtrtrmm<T, T, accu_t>(
                uplo_A, uplo_B,
                trans_A, trans_B,
                diag_A, diag_B,
                n,
                accu_t(alpha),
                A, lda,
                B, ldb,
                accu_t(beta),
                C_hi, ldc,
                stream);

            //-------------------------------
            // native routine via cuBLAS GEMM
            //-------------------------------
            if constexpr (isDouble || isSingle) {
                const std::string funcname =
                    std::to_string(phi) + "," +
                    std::to_string(n) + "," +
                    std::string(1, testTraits<T>::prefix_upper) +
                    "GEMM,";

                PRINT_nobreak(outFile, funcname);

                constexpr bool isOzaki2 = false;

                const double err = evaluate_error<isOzaki2, T, accu_t>(
                    outFile,
                    [&]() {
                        CHECK_CUBLAS(testTraits<T>::trtrmm(
                            handle,
                            trans_A, trans_B,
                            ni, ni, ni,
                            &alpha,
                            A_full, static_cast<int64_t>(n),
                            B_full, static_cast<int64_t>(n),
                            &beta,
                            C, static_cast<int64_t>(ldc)));
                    },
                    uplo_out,
                    n, C, ldc, C_hi, ldc, stream);

                for (unsigned num_moduli = num_moduli_min + 1; num_moduli <= num_moduli_max; ++num_moduli) {
                    PRINT_nobreak(outFile, err << ",");
                }
                PRINT(outFile, "");
                CHECK_CUDA(cudaStreamSynchronize(stream));
            }

            //-------------------------------
            // cuBLAS Ozaki-I via GEMM
            //-------------------------------
#if avail_Ozaki1
            if (run_oz1_i8) {
                CHECK_CUBLAS(cublasSetMathMode(handle_emu, CUBLAS_FP64_EMULATED_FIXEDPOINT_MATH));
                CHECK_CUBLAS(cublasSetEmulationStrategy(handle_emu, CUBLAS_EMULATION_STRATEGY_EAGER));
                CHECK_CUBLAS(cublasSetFixedPointEmulationMantissaControl(handle_emu, CUDA_EMULATION_MANTISSA_CONTROL_FIXED));
                CHECK_CUBLAS(cublasSetWorkspace(handle_emu, work_emu, lwork_emu));

                for (auto &num_slice : oz1_slice_list) {
                    int mantissaBitCount = num_slice * 8 - 1;
                    CHECK_CUBLAS(cublasSetFixedPointEmulationMaxMantissaBitCount(handle_emu, mantissaBitCount));

                    const std::string funcname =
                        std::to_string(phi) + "," +
                        std::to_string(n) + "," +
                        std::string("OS1-") + std::to_string(num_slice) + ",";

                    PRINT_nobreak(outFile, funcname);

                    constexpr bool isOzaki2 = false;

                    const double err = evaluate_error<isOzaki2, T, accu_t>(
                        outFile,
                        [&]() {
                            CHECK_CUBLAS(testTraits<T>::trtrmm(
                                handle_emu,
                                trans_A, trans_B,
                                ni, ni, ni,
                                &alpha,
                                A_full, static_cast<int64_t>(n),
                                B_full, static_cast<int64_t>(n),
                                &beta,
                                C, static_cast<int64_t>(ldc)));
                        },
                        uplo_out,
                        n, C, ldc, C_hi, ldc, stream);

                    for (unsigned num_moduli = num_moduli_min + 1; num_moduli <= num_moduli_max; ++num_moduli) {
                        PRINT_nobreak(outFile, err << ",");
                    }
                    PRINT(outFile, "");
                }

                CHECK_CUBLAS(cublasSetWorkspace(handle_emu, nullptr, 0));
                CHECK_CUBLAS(cublasSetMathMode(handle_emu, CUBLAS_DEFAULT_MATH));
                CHECK_CUDA(cudaStreamSynchronize(stream));
            }
#endif

            //-------------------------------
            // cuBLAS BF16x9 via GEMM
            //-------------------------------
#if avail_BF16x9
            if (use_bf16x9) {
                const std::string funcname =
                    std::to_string(phi) + "," +
                    std::to_string(n) + "," +
                    std::string("BF16x9") + ",";

                PRINT_nobreak(outFile, funcname);

                CHECK_CUBLAS(cublasSetMathMode(handle_emu, CUBLAS_FP32_EMULATED_BF16X9_MATH));

                constexpr bool isOzaki2 = false;

                const double err = evaluate_error<isOzaki2, T, accu_t>(
                    outFile,
                    [&]() {
                        CHECK_CUBLAS(testTraits<T>::trtrmm(
                            handle_emu,
                            trans_A, trans_B,
                            ni, ni, ni,
                            &alpha,
                            A_full, static_cast<int64_t>(n),
                            B_full, static_cast<int64_t>(n),
                            &beta,
                            C, static_cast<int64_t>(ldc)));
                    },
                    uplo_out,
                    n, C, ldc, C_hi, ldc, stream);

                for (unsigned num_moduli = num_moduli_min + 1; num_moduli <= num_moduli_max; ++num_moduli) {
                    PRINT_nobreak(outFile, err << ",");
                }
                PRINT(outFile, "");
                CHECK_CUBLAS(cublasSetMathMode(handle_emu, CUBLAS_DEFAULT_MATH));
                CHECK_CUDA(cudaStreamSynchronize(stream));
            }
#endif

            //-------------------------------
            // fast mode int8
            //-------------------------------
            if (run_oz2_i8) {
                const std::string funcname =
                    std::to_string(phi) + "," +
                    std::to_string(n) + "," +
                    std::string("OS2-i8-fast") + ",";

                PRINT_nobreak(outFile, funcname);

                for (unsigned num_moduli = num_moduli_min; num_moduli <= num_moduli_max; ++num_moduli) {
                    constexpr bool isOzaki2            = true;
                    constexpr bool fastmode            = true;
                    constexpr gemmul8::Backend backend = gemmul8::Backend::INT8;

                    evaluate_error<isOzaki2, T, accu_t>(
                        outFile,
                        [&]() {
                            gemmul8::trtrmm<T, backend, T, T>(
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
                                false, false);
                        },
                        uplo_out,
                        n, C, ldc, C_hi, ldc, stream);
                }

                PRINT(outFile, "");
                CHECK_CUDA(cudaStreamSynchronize(stream));
            }

            //-------------------------------
            // accu mode int8
            //-------------------------------
            if (run_oz2_i8) {
                const std::string funcname =
                    std::to_string(phi) + "," +
                    std::to_string(n) + "," +
                    std::string("OS2-i8-accu") + ",";

                PRINT_nobreak(outFile, funcname);

                for (unsigned num_moduli = num_moduli_min; num_moduli <= num_moduli_max; ++num_moduli) {
                    constexpr bool isOzaki2            = true;
                    constexpr bool fastmode            = false;
                    constexpr gemmul8::Backend backend = gemmul8::Backend::INT8;

                    evaluate_error<isOzaki2, T, accu_t>(
                        outFile,
                        [&]() {
                            gemmul8::trtrmm<T, backend, T, T>(
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
                                false, false);
                        },
                        uplo_out,
                        n, C, ldc, C_hi, ldc, stream);
                }

                PRINT(outFile, "");
                CHECK_CUDA(cudaStreamSynchronize(stream));
            }

            //-------------------------------
            // fast mode fp8
            //-------------------------------
            if (run_oz2_f8) {
                const std::string funcname =
                    std::to_string(phi) + "," +
                    std::to_string(n) + "," +
                    std::string("OS2-f8-fast") + ",";

                PRINT_nobreak(outFile, funcname);

                for (unsigned num_moduli = num_moduli_min; num_moduli <= num_moduli_max; ++num_moduli) {
                    constexpr bool isOzaki2            = true;
                    constexpr bool fastmode            = true;
                    constexpr gemmul8::Backend backend = gemmul8::Backend::FP8;

                    evaluate_error<isOzaki2, T, accu_t>(
                        outFile,
                        [&]() {
                            gemmul8::trtrmmLt<T, backend, T, T>(
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
                                stream);
                        },
                        uplo_out,
                        n, C, ldc, C_hi, ldc, stream);
                }

                PRINT(outFile, "");
                CHECK_CUDA(cudaStreamSynchronize(stream));
            }

            //-------------------------------
            // accu mode fp8
            //-------------------------------
            if (run_oz2_f8) {
                const std::string funcname =
                    std::to_string(phi) + "," +
                    std::to_string(n) + "," +
                    std::string("OS2-f8-accu") + ",";

                PRINT_nobreak(outFile, funcname);

                for (unsigned num_moduli = num_moduli_min; num_moduli <= num_moduli_max; ++num_moduli) {
                    constexpr bool isOzaki2            = true;
                    constexpr bool fastmode            = false;
                    constexpr gemmul8::Backend backend = gemmul8::Backend::FP8;

                    evaluate_error<isOzaki2, T, accu_t>(
                        outFile,
                        [&]() {
                            gemmul8::trtrmmLt<T, backend, T, T>(
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
                                stream);
                        },
                        uplo_out,
                        n, C, ldc, C_hi, ldc, stream);
                }

                PRINT(outFile, "");
                CHECK_CUDA(cudaStreamSynchronize(stream));
            }
        }
    }

    std::cout << std::endl;

    CHECK_CUDA(cudaStreamSynchronize(stream));

    if (memory_limit > 0) {
        gemmul8::clear_config(handle);
        gemmul8::clear_configLt(handleLt);
    }

    CHECK_CUDA(cudaFreeAsync(work_blas, stream));
    CHECK_CUDA(cudaFreeAsync(work_emu, stream));
    CHECK_CUDA(cudaFreeAsync(C_hi, stream));
    CHECK_CUDA(cudaFreeAsync(C, stream));
    CHECK_CUDA(cudaFreeAsync(B_full, stream));
    CHECK_CUDA(cudaFreeAsync(A_full, stream));
    CHECK_CUDA(cudaFreeAsync(B, stream));
    CHECK_CUDA(cudaFreeAsync(A, stream));

    CHECK_CUBLAS(cublasLtDestroy(handleLt));
    CHECK_CUBLAS(cublasDestroy(handle_emu));
    CHECK_CUBLAS(cublasDestroy(handle));
    cudaStreamDestroy(stream);

    outFile.close();
}

template void check_accuracy<float>(
    std::string &, std::string &, size_t, cublasFillMode_t, cublasFillMode_t,
    cublasOperation_t, cublasOperation_t, cublasDiagType_t, cublasDiagType_t,
    const bool, const bool, const bool);

template void check_accuracy<double>(
    std::string &, std::string &, size_t, cublasFillMode_t, cublasFillMode_t,
    cublasOperation_t, cublasOperation_t, cublasDiagType_t, cublasDiagType_t,
    const bool, const bool, const bool);

template void check_accuracy<cuFloatComplex>(
    std::string &, std::string &, size_t, cublasFillMode_t, cublasFillMode_t,
    cublasOperation_t, cublasOperation_t, cublasDiagType_t, cublasDiagType_t,
    const bool, const bool, const bool);

template void check_accuracy<cuDoubleComplex>(
    std::string &, std::string &, size_t, cublasFillMode_t, cublasFillMode_t,
    cublasOperation_t, cublasOperation_t, cublasDiagType_t, cublasDiagType_t,
    const bool, const bool, const bool);

} // namespace bench::accuracy::trtrmm
