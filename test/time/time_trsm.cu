#include "../common/common.hpp"
#include "../common/ozaki1.hpp"
#include "time_trsm.hpp"

namespace bench::time::trsm {

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

template <bool isOzaki2 = false, class Reset, class F>
inline void evaluate_time(
    std::ofstream &outFile,
    const double computational_cost,
    const std::string &func_name,
    Reset &&reset,
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

        double elapsed_time = 0.0;
        for (unsigned i = 0; i < warmup_max; ++i) {
            reset();
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
            reset();
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
            time0_rep[rep]   = time0_med;
            time1_rep[rep]   = time1_med;
        }
    }

    double times_rep_med = calc_median(times_rep);
    double TFLOPS        = computational_cost / times_rep_med * 1.0e-12;
    if constexpr (isOzaki2) {
        double time0_rep_med = calc_median(time0_rep);
        double time1_rep_med = calc_median(time1_rep);
        PRINT(outFile, func_name << "," << TFLOPS << "," << times_rep_med << "," << time0_rep_med << "," << time1_rep_med << ",");
    } else {
        PRINT(outFile, func_name << "," << TFLOPS << "," << times_rep_med << "," << "," << ",");
    }

    CHECK_CUDA(cudaEventDestroy(stop));
    CHECK_CUDA(cudaEventDestroy(start));
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
void check_time(
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
                  "Unsupported type for TRSM time test.");

    if (side != CUBLAS_SIDE_LEFT && side != CUBLAS_SIDE_RIGHT) {
        assert(false && "TRSM requires side = LEFT or RIGHT.");
        return;
    }
    if (uplo != CUBLAS_FILL_MODE_UPPER && uplo != CUBLAS_FILL_MODE_LOWER) {
        assert(false && "TRSM requires uplo = UPPER or LOWER.");
        return;
    }
    if (!testTraits<T>::is_complex && trans == CUBLAS_OP_C) {
        assert(false && "real TRSM time test does not use trans = C.");
        return;
    }

    std::string square_tag = is_square ? std::string("square_") : std::string("");
    std::string fileName   = std::string("oz2_results_") +
                             testTraits<T>::prefix +
                             std::string(getFuncName<func>()) +
                             std::string("_time_") +
                             square_tag +
                             std::string(side_tag(side)) +
                             std::string(uplo_tag(uplo)) +
                             std::string(trans_tag(trans)) +
                             std::string(diag_tag(diag)) +
                             memTag(memory_limit) +
                             deviceName + "_" + dateTime + ".csv";

    std::ofstream outFile(fileName);

    PRINT(outFile,
          "phi,m,n,function,TFLOPS,total_time[sec],standard_trsm,gemmul8_gemm,");

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

    for (auto &mn : N_list) {
        std::vector<size_t> K_list;
        if (is_square) {
            K_list = std::vector<size_t>{mn};
        } else {
            K_list = N_list;
        }

        for (auto &k : K_list) {
            const size_t m = (side == CUBLAS_SIDE_LEFT) ? k : mn;
            const size_t n = (side == CUBLAS_SIDE_LEFT) ? mn : k;

            const int64_t mi = static_cast<int64_t>(m);
            const int64_t ni = static_cast<int64_t>(n);

            const size_t lda = k;
            const size_t ldb = m;

            T *A            = nullptr;
            T *B0           = nullptr;
            T *B            = nullptr;
            void *work_emu  = nullptr;
            void *work_blas = nullptr;

            bool run_oz2_i8 = run_Ozaki2_I8;
            bool run_oz2_f8 = run_Ozaki2_F8;
            bool run_oz1_i8 = run_Ozaki1_I8 && use_ozaki1;

            const double comp_cost            = getFuncCost<func>(m, n, k, testTraits<T>::is_complex);
            const size_t size_A               = k * k * sizeof(T);
            const size_t size_B               = m * n * sizeof(T);
            constexpr size_t lwork_blas       = size_t(32) << 20;
            constexpr size_t safety_margin    = size_t(256) << 20;
            constexpr size_t alignment_margin = size_t(1) << 20;
            size_t lwork_ABC                  = lwork_blas + size_A + size_B + size_B;

            size_t free_bytes  = 0;
            size_t total_bytes = 0;
            CHECK_CUDA(cudaStreamSynchronize(stream));
            CHECK_CUDA(cudaMemGetInfo(&free_bytes, &total_bytes));

            const size_t lwork_total = lwork_ABC + safety_margin + alignment_margin;
            bool alloc_ok            = (lwork_total <= free_bytes);
            CHECK_CUDA(cudaStreamSynchronize(stream));
            if (alloc_ok) alloc_ok = alloc_ok && (cudaMallocAsync(reinterpret_cast<void **>(&A), size_A, stream) == cudaSuccess);
            CHECK_CUDA(cudaStreamSynchronize(stream));
            if (alloc_ok) alloc_ok = alloc_ok && (cudaMallocAsync(reinterpret_cast<void **>(&B0), size_B, stream) == cudaSuccess);
            CHECK_CUDA(cudaStreamSynchronize(stream));
            if (alloc_ok) alloc_ok = alloc_ok && (cudaMallocAsync(reinterpret_cast<void **>(&B), size_B, stream) == cudaSuccess);
            CHECK_CUDA(cudaStreamSynchronize(stream));
            if (alloc_ok) alloc_ok = alloc_ok && (cudaMallocAsync(&work_blas, lwork_blas, stream) == cudaSuccess);
            CHECK_CUDA(cudaStreamSynchronize(stream));

            if (!alloc_ok) {
                free_async_if_needed(A, stream);
                free_async_if_needed(B0, stream);
                free_async_if_needed(B, stream);
                free_async_if_needed_void(work_blas, stream);
                CHECK_CUDA(cudaStreamSynchronize(stream));
                continue;
            }

            CHECK_CUBLAS(cublasSetWorkspace(handle, work_blas, lwork_blas));
            CHECK_CUDA(cudaStreamSynchronize(stream));

            // test matrices
            makemat::randmat<T>(k, k, A, phi, seedA, stream, uplo, diag);
            stabilize_trsm_triangular<T>(k, A, lda, uplo, diag, stream);

            makemat::set_ones<T>(m, n, B, ldb, stream);
            CHECK_CUDA(cudaMemsetAsync(B0, 0, size_B, stream));

            const T trmm_alpha = testTraits<T>::one();
            CHECK_CUBLAS(testTraits<T>::trmm(
                handle,
                side, uplo, trans, diag,
                mi, ni,
                &trmm_alpha,
                A, static_cast<int64_t>(lda),
                B, static_cast<int64_t>(ldb),
                B0, static_cast<int64_t>(ldb)));

            auto reset_rhs = [&]() {
                CHECK_CUDA(cudaMemcpyAsync(B, B0, size_B, cudaMemcpyDeviceToDevice, stream));
            };

            bool run_native = false;

            //-------------------------------
            // fast mode int8
            //-------------------------------
            if (run_oz2_i8) {
                for (unsigned num_moduli = num_moduli_min; num_moduli <= num_moduli_max; ++num_moduli) {

                    const size_t lwork_gemmul8_i8 = (memory_limit > 0) ? memory_limit : gemmul8::workSizeTrsm<T, gemmul8::Backend::INT8>(handle, side, m, n, num_moduli);
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
                        std::to_string(m) + "," +
                        std::to_string(n) + "," +
                        std::string("OS2-i8-fast-") + std::to_string(num_moduli);

                    constexpr bool isOzaki2            = true;
                    constexpr bool fastmode            = true;
                    constexpr gemmul8::Backend backend = gemmul8::Backend::INT8;

                    evaluate_time<isOzaki2>(
                        outFile,
                        comp_cost,
                        funcname,
                        reset_rhs,
                        [&]() { return gemmul8::trsm<T, backend, T>(
                                    handle,
                                    side, uplo, trans, diag,
                                    m, n,
                                    &alpha,
                                    A, lda,
                                    B, ldb,
                                    num_moduli, fastmode,
                                    work_emu); },
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

                    const size_t lwork_gemmul8_i8 = (memory_limit > 0) ? memory_limit : gemmul8::workSizeTrsm<T, gemmul8::Backend::INT8>(handle, side, m, n, num_moduli);
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
                        std::to_string(m) + "," +
                        std::to_string(n) + "," +
                        std::string("OS2-i8-accu-") + std::to_string(num_moduli);

                    constexpr bool isOzaki2            = true;
                    constexpr bool fastmode            = false;
                    constexpr gemmul8::Backend backend = gemmul8::Backend::INT8;

                    evaluate_time<isOzaki2>(
                        outFile,
                        comp_cost,
                        funcname,
                        reset_rhs,
                        [&]() { return gemmul8::trsm<T, backend, T>(
                                    handle,
                                    side, uplo, trans, diag,
                                    m, n,
                                    &alpha,
                                    A, lda,
                                    B, ldb,
                                    num_moduli, fastmode,
                                    work_emu); },
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

                    const size_t lwork_gemmul8_f8 = (memory_limit > 0) ? memory_limit : gemmul8::workSizeTrsmLt<T, gemmul8::Backend::FP8>(handleLt, side, m, n, num_moduli);
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
                        std::to_string(m) + "," +
                        std::to_string(n) + "," +
                        std::string("OS2-f8-fast-") + std::to_string(num_moduli);

                    constexpr bool isOzaki2            = true;
                    constexpr bool fastmode            = true;
                    constexpr gemmul8::Backend backend = gemmul8::Backend::FP8;

                    evaluate_time<isOzaki2>(
                        outFile,
                        comp_cost,
                        funcname,
                        reset_rhs,
                        [&]() { return gemmul8::trsmLt<T, backend, T>(
                                    handleLt,
                                    side, uplo, trans, diag,
                                    m, n,
                                    &alpha,
                                    A, lda,
                                    B, ldb,
                                    num_moduli, fastmode,
                                    work_emu,
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

                    const size_t lwork_gemmul8_f8 = (memory_limit > 0) ? memory_limit : gemmul8::workSizeTrsmLt<T, gemmul8::Backend::FP8>(handleLt, side, m, n, num_moduli);
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
                        std::to_string(m) + "," +
                        std::to_string(n) + "," +
                        std::string("OS2-f8-accu-") + std::to_string(num_moduli);

                    constexpr bool isOzaki2            = true;
                    constexpr bool fastmode            = false;
                    constexpr gemmul8::Backend backend = gemmul8::Backend::FP8;

                    evaluate_time<isOzaki2>(
                        outFile,
                        comp_cost,
                        funcname,
                        reset_rhs,
                        [&]() { return gemmul8::trsmLt<T, backend, T>(
                                    handleLt,
                                    side, uplo, trans, diag,
                                    m, n,
                                    &alpha,
                                    A, lda,
                                    B, ldb,
                                    num_moduli, fastmode,
                                    work_emu,
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

                    int mantissaBitCount = num_slice * 8 - 1;
                    const size_t lwork_ozaki1 =
                        std::max<size_t>(
                            lwork_blas,
                            ozaki1::workSize(m, n, k, 1,
                                             testTraits<T>::is_complex,
                                             mantissaBitCount));
                    bool alloc_ok = (lwork_total + lwork_ozaki1 <= free_bytes);
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
                        std::to_string(m) + "," +
                        std::to_string(n) + "," +
                        std::string("OS1-") + std::to_string(num_slice);

                    constexpr bool isOzaki2 = false;

                    evaluate_time<isOzaki2>(
                        outFile,
                        comp_cost,
                        funcname,
                        reset_rhs,
                        [&]() { CHECK_CUBLAS(testTraits<T>::trsm(
                                    handle_emu,
                                    side, uplo, trans, diag,
                                    mi, ni,
                                    &alpha,
                                    A, static_cast<int64_t>(lda),
                                    B, static_cast<int64_t>(ldb))); },
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
                    std::to_string(m) + "," +
                    std::to_string(n) + "," +
                    std::string("BF16x9");

                constexpr bool isOzaki2 = false;

                evaluate_time<isOzaki2>(
                    outFile,
                    comp_cost,
                    funcname,
                    reset_rhs,
                    [&]() { CHECK_CUBLAS(testTraits<T>::trsm(
                                handle_emu,
                                side, uplo, trans, diag,
                                mi, ni,
                                &alpha,
                                A, static_cast<int64_t>(lda),
                                B, static_cast<int64_t>(ldb))); },
                    stream);

                CHECK_CUDA(cudaStreamSynchronize(stream));
                CHECK_CUBLAS(cublasSetWorkspace(handle_emu, nullptr, 0));
                CHECK_CUBLAS(cublasSetMathMode(handle_emu, CUBLAS_DEFAULT_MATH));
            }
#endif

            if (!run_native) {
                CHECK_CUDA(cudaStreamSynchronize(stream));
                free_async_if_needed(A, stream);
                free_async_if_needed(B0, stream);
                free_async_if_needed(B, stream);
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
                    std::to_string(m) + "," +
                    std::to_string(n) + "," +
                    std::string(1, testTraits<T>::prefix_upper) + std::string(getFuncName_upper<func>());

                constexpr bool isOzaki2 = false;

                evaluate_time<isOzaki2>(
                    outFile,
                    comp_cost,
                    funcname,
                    reset_rhs,
                    [&]() { CHECK_CUBLAS(testTraits<T>::trsm(
                                handle,
                                side, uplo, trans, diag,
                                mi, ni,
                                &alpha,
                                A, static_cast<int64_t>(lda),
                                B, static_cast<int64_t>(ldb))); },
                    stream);
                CHECK_CUDA(cudaStreamSynchronize(stream));
            }

            CHECK_CUDA(cudaStreamSynchronize(stream));
            free_async_if_needed(A, stream);
            free_async_if_needed(B0, stream);
            free_async_if_needed(B, stream);
            free_async_if_needed_void(work_blas, stream);
            CHECK_CUDA(cudaStreamSynchronize(stream));
        }
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
    std::string &, std::string &, size_t, cublasSideMode_t, cublasFillMode_t,
    cublasOperation_t, cublasDiagType_t,
    const bool, const bool, const bool, const bool);

template void check_time<double>(
    std::string &, std::string &, size_t, cublasSideMode_t, cublasFillMode_t,
    cublasOperation_t, cublasDiagType_t,
    const bool, const bool, const bool, const bool);

template void check_time<cuFloatComplex>(
    std::string &, std::string &, size_t, cublasSideMode_t, cublasFillMode_t,
    cublasOperation_t, cublasDiagType_t,
    const bool, const bool, const bool, const bool);

template void check_time<cuDoubleComplex>(
    std::string &, std::string &, size_t, cublasSideMode_t, cublasFillMode_t,
    cublasOperation_t, cublasDiagType_t,
    const bool, const bool, const bool, const bool);

} // namespace bench::time::trsm
