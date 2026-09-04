#include "../common/common.hpp"
#include "../common/ozaki1.hpp"
#include "accuracy_trsm.hpp"

namespace bench::accuracy::trsm {

inline constexpr gemmul8::Func func = gemmul8::Func::trsm;

template <typename T>
__global__ void stabilize_trsm_triangular_kernel(
    const size_t n,
    T *const A,
    const size_t lda,
    const cublasFillMode_t uplo,
    const cublasDiagType_t diag //
) {
    const size_t idx = blockIdx.x * size_t(blockDim.x) + threadIdx.x;
    if (idx >= n * n) return;

    const size_t col = idx / n;
    const size_t row = idx - col * n;

    const bool active =
        (uplo == CUBLAS_FILL_MODE_UPPER) ? (row <= col)
                                         : ((uplo == CUBLAS_FILL_MODE_LOWER) ? (row >= col)
                                                                             : true);

    if (!active) return;

    const size_t pos = row + col * lda;

    if (row == col) {
        if (diag == CUBLAS_DIAG_UNIT) {
            if constexpr (std::is_same_v<T, cuFloatComplex>) {
                A[pos] = cuFloatComplex{1.0f, 0.0f};
            } else if constexpr (std::is_same_v<T, cuDoubleComplex>) {
                A[pos] = cuDoubleComplex{1.0, 0.0};
            } else {
                A[pos] = T(1);
            }
        } else {
            if constexpr (std::is_same_v<T, cuFloatComplex>) {
                A[pos] = cuFloatComplex{2.0f, 0.0f};
            } else if constexpr (std::is_same_v<T, cuDoubleComplex>) {
                A[pos] = cuDoubleComplex{2.0, 0.0};
            } else {
                A[pos] = T(2);
            }
        }
        return;
    }

    constexpr double scale = 1.0 / 64.0;

    if constexpr (std::is_same_v<T, cuFloatComplex>) {
        A[pos].x *= static_cast<float>(scale);
        A[pos].y *= static_cast<float>(scale);
    } else if constexpr (std::is_same_v<T, cuDoubleComplex>) {
        A[pos].x *= scale;
        A[pos].y *= scale;
    } else {
        A[pos] = T(double(A[pos]) * scale);
    }
}

template <typename T>
inline void stabilize_trsm_triangular(
    const size_t n,
    T *const A,
    const size_t lda,
    const cublasFillMode_t uplo,
    const cublasDiagType_t diag,
    const cudaStream_t stream //
) {
    constexpr size_t block_size = 256;
    const size_t grid_size      = (n * n + block_size - 1) / block_size;

    stabilize_trsm_triangular_kernel<T>
        <<<grid_size, block_size, 0, stream>>>(
            n, A, lda, uplo, diag);
}

template <bool isOzaki2 = false, typename T, typename accu_t, class F>
inline double evaluate_error(
    std::ofstream &outFile,
    F &&f,
    const size_t rows, const size_t cols,
    const T *const B0,
    T *const B, const size_t ldb,
    const accu_t *const B_exact, const size_t ldb_exact,
    const cudaStream_t stream //
) {
    CHECK_CUDA(cudaStreamSynchronize(stream));
    CHECK_CUDA(cudaMemcpyAsync(B, B0, ldb * cols * sizeof(T), cudaMemcpyDeviceToDevice, stream));
    CHECK_CUDA(cudaStreamSynchronize(stream));

    f();

    CHECK_CUDA(cudaStreamSynchronize(stream));
    auto [err_max, err_med] =
        eval::calc_err(rows, cols, B, ldb, B_exact, ldb_exact, stream);

    PRINT_nobreak(outFile, err_max << ",");
    return err_max;
}

inline const char *side_tag(const cublasSideMode_t side) {
    return (side == CUBLAS_SIDE_LEFT) ? "left_" : "right_";
}

inline const char *uplo_tag(const cublasFillMode_t uplo) {
    return (uplo == CUBLAS_FILL_MODE_UPPER) ? "upper_" : "lower_";
}

inline const char *trans_tag(const cublasOperation_t trans) {
    return (trans == CUBLAS_OP_N) ? "n_" : ((trans == CUBLAS_OP_T) ? "t_" : "c_");
}

inline const char *diag_tag(const cublasDiagType_t diag) {
    return (diag == CUBLAS_DIAG_UNIT) ? "unit_" : "nonunit_";
}

template <typename T>
void check_accuracy(
    std::string &deviceName,
    std::string &dateTime,
    size_t memory_limit,
    cublasSideMode_t side,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    cublasDiagType_t diag,
    const bool run_Ozaki2_I8,
    const bool run_Ozaki2_F8,
    const bool run_Ozaki1_I8,
    const bool is_square //
) {
    static_assert(std::is_same_v<T, float> || std::is_same_v<T, double> ||
                      std::is_same_v<T, cuFloatComplex> || std::is_same_v<T, cuDoubleComplex>,
                  "Unsupported type for TRSM accuracy test.");

    if (side != CUBLAS_SIDE_LEFT && side != CUBLAS_SIDE_RIGHT) {
        assert(false && "TRSM requires side = LEFT or RIGHT.");
        return;
    }
    if (uplo != CUBLAS_FILL_MODE_UPPER && uplo != CUBLAS_FILL_MODE_LOWER) {
        assert(false && "TRSM requires uplo = UPPER or LOWER.");
        return;
    }
    if (!testTraits<T>::is_complex && trans == CUBLAS_OP_C) {
        assert(false && "real TRSM accuracy test does not use trans = C.");
        return;
    }

    std::string square_tag = is_square ? std::string("square_") : std::string("");
    std::string fileName   = std::string("oz2_results_") +
                             testTraits<T>::prefix +
                             std::string(getFuncName<func>()) +
                             std::string("_accuracy_") +
                             square_tag +
                             std::string(side_tag(side)) +
                             std::string(uplo_tag(uplo)) +
                             std::string(trans_tag(trans)) +
                             std::string(diag_tag(diag)) +
                             memTag(memory_limit) +
                             deviceName + "_" + dateTime + ".csv";

    std::ofstream outFile(fileName);

    const unsigned num_moduli_min = testTraits<T>::NUM_MODULI_MIN_accuracy;
    const unsigned num_moduli_max = testTraits<T>::NUM_MODULI_MAX_accuracy;

    std::string num_moduli_str = std::string("");
    for (unsigned num_moduli = num_moduli_min; num_moduli <= num_moduli_max; ++num_moduli) {
        num_moduli_str += std::to_string(num_moduli) + ",";
    }

    PRINT(outFile, "phi,m,n,function," + num_moduli_str);

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

    const T alpha           = testTraits<T>::one();
    const int num_slice_max = *std::max_element(begin(oz1_slice_list), end(oz1_slice_list));

    using accu_t = typename testTraits<T>::accu_t;

    T *A         = nullptr;
    T *B0        = nullptr;
    T *B         = nullptr;
    accu_t *B_hi = nullptr;

    void *work_emu  = nullptr;
    void *work_blas = nullptr;

    bool run_oz2_i8 = run_Ozaki2_I8;
    bool run_oz2_f8 = run_Ozaki2_F8;
    bool run_oz1_i8 = run_Ozaki1_I8 && use_ozaki1;

    const size_t fixed_dim = is_square ? size_t(8192) : size_t(128);
    const size_t k_max     = is_square ? size_t(8192) : *std::max_element(begin(N_list), end(N_list));

    const size_t m_max = (side == CUBLAS_SIDE_LEFT) ? k_max : fixed_dim;
    const size_t n_max = (side == CUBLAS_SIDE_LEFT) ? fixed_dim : k_max;

    const size_t size_A = k_max * k_max;
    const size_t size_B = m_max * n_max;

    const size_t lwork_blas = size_t(32) << 20;

    if (memory_limit > 0) {
        gemmul8::set_memory_saving(handle, true);
        gemmul8::set_max_worksize(handle, memory_limit);
        gemmul8::set_memory_savingLt(handleLt, true);
        gemmul8::set_max_worksizeLt(handleLt, memory_limit);
    }
    const size_t lwork_gemmul8_i8_ = gemmul8::workSizeTrsm<T, gemmul8::Backend::INT8>(handle, side, m_max, n_max, num_moduli_max);
    const size_t lwork_gemmul8_f8_ = gemmul8::workSizeTrsmLt<T, gemmul8::Backend::FP8>(handleLt, side, m_max, n_max, num_moduli_max);
    const size_t lwork_ozaki1_     = std::max<size_t>(lwork_blas, ozaki1::workSize(m_max, n_max, k_max, 1, testTraits<T>::is_complex, 8 * num_slice_max - 1));
    const size_t lwork_gemmul8_i8  = (run_oz2_i8) ? std::max(memory_limit, lwork_gemmul8_i8_) : 0;
    const size_t lwork_gemmul8_f8  = (run_oz2_f8) ? std::max(memory_limit, lwork_gemmul8_f8_) : 0;
    const size_t lwork_ozaki1      = (run_oz1_i8) ? lwork_ozaki1_ : 0;
    const size_t lwork_emu         = std::max(std::max(lwork_gemmul8_i8, lwork_gemmul8_f8), lwork_ozaki1);

    CHECK_CUDA(cudaMallocAsync(reinterpret_cast<void **>(&A), size_A * sizeof(T), stream));
    CHECK_CUDA(cudaMallocAsync(reinterpret_cast<void **>(&B0), size_B * sizeof(T), stream));
    CHECK_CUDA(cudaMallocAsync(reinterpret_cast<void **>(&B), size_B * sizeof(T), stream));
    CHECK_CUDA(cudaMallocAsync(reinterpret_cast<void **>(&B_hi), size_B * sizeof(accu_t), stream));

    CHECK_CUDA(cudaMallocAsync(&work_emu, lwork_emu, stream));
    CHECK_CUDA(cudaMallocAsync(&work_blas, lwork_blas, stream));
    CHECK_CUDA(cudaStreamSynchronize(stream));

    CHECK_CUBLAS(cublasSetWorkspace(handle, work_blas, lwork_blas));
    CHECK_CUDA(cudaStreamSynchronize(stream));

    double phi = -1.0;
    std::vector<size_t> K_list;
    if (is_square) {
        K_list = std::vector<size_t>{size_t(8192)};
    } else {
        K_list = N_list;
    }

    for (auto &k : K_list) {

        const size_t m = (side == CUBLAS_SIDE_LEFT) ? k : fixed_dim;
        const size_t n = (side == CUBLAS_SIDE_LEFT) ? fixed_dim : k;

        const int64_t mi = static_cast<int64_t>(m);
        const int64_t ni = static_cast<int64_t>(n);

        const size_t lda = k;
        const size_t ldb = m;

        makemat::randmat<T>(k, k, A, phi, seedA, stream, uplo, diag);
        stabilize_trsm_triangular<T>(k, A, lda, uplo, diag, stream);

        makemat::set_ones<T>(m, n, B, ldb, stream);

        const T trmm_alpha = testTraits<T>::one();
        CHECK_CUBLAS(testTraits<T>::trmm(
            handle,
            side, uplo, trans, diag,
            mi, ni,
            &trmm_alpha,
            A, static_cast<int64_t>(lda),
            B, static_cast<int64_t>(ldb),
            B0, static_cast<int64_t>(ldb)));

        // Accurate result
        CHECK_CUDA(cudaMemsetAsync(B_hi, 0, size_B * sizeof(accu_t), stream));
        eval::addvec<T, accu_t>(m * n, B0, B_hi, stream);
        eval::DDtrsm<T, accu_t>(
            handle,
            side, uplo, trans, diag,
            m, n,
            accu_t(alpha),
            A, lda,
            B_hi, ldb,
            stream);

        //-------------------------------
        // native routine
        //-------------------------------
        if constexpr (isDouble || isSingle) {
            const std::string funcname =
                std::to_string(phi) + "," +
                std::to_string(m) + "," +
                std::to_string(n) + "," +
                std::string(1, testTraits<T>::prefix_upper) +
                std::string(getFuncName_upper<func>()) + ",";

            PRINT_nobreak(outFile, funcname);

            constexpr bool isOzaki2 = false;

            const double err = evaluate_error<isOzaki2, T, accu_t>(
                outFile,
                [&]() {
                    CHECK_CUBLAS(testTraits<T>::trsm(
                        handle,
                        side, uplo, trans, diag,
                        mi, ni,
                        &alpha,
                        A, static_cast<int64_t>(lda),
                        B, static_cast<int64_t>(ldb)));
                },
                m, n, B0, B, ldb, B_hi, ldb, stream);

            for (unsigned num_moduli = num_moduli_min + 1; num_moduli <= num_moduli_max; ++num_moduli) {
                PRINT_nobreak(outFile, err << ",");
            }
            PRINT(outFile, "");
            CHECK_CUDA(cudaStreamSynchronize(stream));
        }

        //-------------------------------
        // cuBLAS Ozaki-I
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
                    std::to_string(m) + "," +
                    std::to_string(n) + "," +
                    std::string("OS1-") + std::to_string(num_slice) + ",";

                PRINT_nobreak(outFile, funcname);

                constexpr bool isOzaki2 = false;

                const double err = evaluate_error<isOzaki2, T, accu_t>(
                    outFile,
                    [&]() {
                        CHECK_CUBLAS(testTraits<T>::trsm(
                            handle_emu,
                            side, uplo, trans, diag,
                            mi, ni,
                            &alpha,
                            A, static_cast<int64_t>(lda),
                            B, static_cast<int64_t>(ldb)));
                    },
                    m, n, B0, B, ldb, B_hi, ldb, stream);

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
        // cuBLAS BF16x9
        //-------------------------------
#if avail_BF16x9
        if (use_bf16x9) {
            const std::string funcname =
                std::to_string(phi) + "," +
                std::to_string(m) + "," +
                std::to_string(n) + "," +
                std::string("BF16x9") + ",";

            PRINT_nobreak(outFile, funcname);

            CHECK_CUBLAS(cublasSetMathMode(handle_emu, CUBLAS_FP32_EMULATED_BF16X9_MATH));

            constexpr bool isOzaki2 = false;

            const double err = evaluate_error<isOzaki2, T, accu_t>(
                outFile,
                [&]() {
                    CHECK_CUBLAS(testTraits<T>::trsm(
                        handle_emu,
                        side, uplo, trans, diag,
                        mi, ni,
                        &alpha,
                        A, static_cast<int64_t>(lda),
                        B, static_cast<int64_t>(ldb)));
                },
                m, n, B0, B, ldb, B_hi, ldb, stream);

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
                std::to_string(m) + "," +
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
                        gemmul8::trsm<T, backend, T>(
                            handle,
                            side, uplo, trans, diag,
                            m, n,
                            &alpha,
                            A, lda,
                            B, ldb,
                            num_moduli, fastmode,
                            work_emu);
                    },
                    m, n, B0, B, ldb, B_hi, ldb, stream);
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
                std::to_string(m) + "," +
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
                        gemmul8::trsm<T, backend, T>(
                            handle,
                            side, uplo, trans, diag,
                            m, n,
                            &alpha,
                            A, lda,
                            B, ldb,
                            num_moduli, fastmode,
                            work_emu);
                    },
                    m, n, B0, B, ldb, B_hi, ldb, stream);
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
                std::to_string(m) + "," +
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
                        gemmul8::trsmLt<T, backend, T>(
                            handleLt,
                            side, uplo, trans, diag,
                            m, n,
                            &alpha,
                            A, lda,
                            B, ldb,
                            num_moduli, fastmode,
                            work_emu,
                            stream);
                    },
                    m, n, B0, B, ldb, B_hi, ldb, stream);
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
                std::to_string(m) + "," +
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
                        gemmul8::trsmLt<T, backend, T>(
                            handleLt,
                            side, uplo, trans, diag,
                            m, n,
                            &alpha,
                            A, lda,
                            B, ldb,
                            num_moduli, fastmode,
                            work_emu,
                            stream);
                    },
                    m, n, B0, B, ldb, B_hi, ldb, stream);
            }

            PRINT(outFile, "");
            CHECK_CUDA(cudaStreamSynchronize(stream));
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
    CHECK_CUDA(cudaFreeAsync(B_hi, stream));
    CHECK_CUDA(cudaFreeAsync(B, stream));
    CHECK_CUDA(cudaFreeAsync(B0, stream));
    CHECK_CUDA(cudaFreeAsync(A, stream));

    CHECK_CUBLAS(cublasLtDestroy(handleLt));
    CHECK_CUBLAS(cublasDestroy(handle_emu));
    CHECK_CUBLAS(cublasDestroy(handle));
    cudaStreamDestroy(stream);

    outFile.close();
}

template void check_accuracy<float>(
    std::string &, std::string &, size_t, cublasSideMode_t, cublasFillMode_t,
    cublasOperation_t, cublasDiagType_t,
    const bool, const bool, const bool, const bool);

template void check_accuracy<double>(
    std::string &, std::string &, size_t, cublasSideMode_t, cublasFillMode_t,
    cublasOperation_t, cublasDiagType_t,
    const bool, const bool, const bool, const bool);

template void check_accuracy<cuFloatComplex>(
    std::string &, std::string &, size_t, cublasSideMode_t, cublasFillMode_t,
    cublasOperation_t, cublasDiagType_t,
    const bool, const bool, const bool, const bool);

template void check_accuracy<cuDoubleComplex>(
    std::string &, std::string &, size_t, cublasSideMode_t, cublasFillMode_t,
    cublasOperation_t, cublasDiagType_t,
    const bool, const bool, const bool, const bool);

} // namespace bench::accuracy::trsm
