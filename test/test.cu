#include "time/time_gemm.hpp"
#include "time/time_symm.hpp"
#include "time/time_syrk.hpp"
#include "time/time_syr2k.hpp"
#include "time/time_syrkx.hpp"
#include "time/time_hemm.hpp"
#include "time/time_herk.hpp"
#include "time/time_her2k.hpp"
#include "time/time_herkx.hpp"
#include "time/time_trmm.hpp"
#include "time/time_trsm.hpp"
#include "time/time_trtrmm.hpp"

#include "accuracy/accuracy_gemm.hpp"
#include "accuracy/accuracy_symm.hpp"
#include "accuracy/accuracy_syrk.hpp"
#include "accuracy/accuracy_syr2k.hpp"
#include "accuracy/accuracy_syrkx.hpp"
#include "accuracy/accuracy_hemm.hpp"
#include "accuracy/accuracy_herk.hpp"
#include "accuracy/accuracy_her2k.hpp"
#include "accuracy/accuracy_herkx.hpp"
#include "accuracy/accuracy_trmm.hpp"
#include "accuracy/accuracy_trtrmm.hpp"
#include "accuracy/accuracy_trsm.hpp"

#include <optional>
#include <cctype>
#include <stdexcept>

inline std::string upper_string(std::string s) {
    for (char &c : s) {
        c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    }
    return s;
}

inline bool split_key_value(
    const std::string &arg,
    std::string &key,
    std::string &value //
) {
    const auto pos = arg.find('=');
    if (pos == std::string::npos) return false;

    key   = arg.substr(0, pos);
    value = arg.substr(pos + 1);
    key   = upper_string(key);
    value = upper_string(value);
    return true;
}

inline std::vector<cublasFillMode_t> parse_uplo_list(const std::string &value) {
    if (value == "ALL" || value == "A") {
        return {CUBLAS_FILL_MODE_UPPER, CUBLAS_FILL_MODE_LOWER};
    }
    if (value == "UPPER" || value == "U") {
        return {CUBLAS_FILL_MODE_UPPER};
    }
    if (value == "LOWER" || value == "L") {
        return {CUBLAS_FILL_MODE_LOWER};
    }

    throw std::runtime_error("invalid uplo option: " + value);
}

inline std::vector<cublasSideMode_t> parse_side_list(const std::string &value) {
    if (value == "ALL" || value == "A") {
        return {CUBLAS_SIDE_LEFT, CUBLAS_SIDE_RIGHT};
    }
    if (value == "LEFT" || value == "L") {
        return {CUBLAS_SIDE_LEFT};
    }
    if (value == "RIGHT" || value == "R") {
        return {CUBLAS_SIDE_RIGHT};
    }

    throw std::runtime_error("invalid side option: " + value);
}

inline std::vector<cublasOperation_t> parse_trans_list(const std::string &value) {
    if (value == "ALL" || value == "A") {
        return {CUBLAS_OP_N, CUBLAS_OP_T, CUBLAS_OP_C};
    }
    if (value == "N" || value == "NO" || value == "NONTRANS") {
        return {CUBLAS_OP_N};
    }
    if (value == "T" || value == "TRANS") {
        return {CUBLAS_OP_T};
    }
    if (value == "C" || value == "CONJ" || value == "CONJTRANS") {
        return {CUBLAS_OP_C};
    }

    throw std::runtime_error("invalid trans option: " + value);
}

inline std::vector<cublasDiagType_t> parse_diag_list(const std::string &value) {
    if (value == "ALL" || value == "A") {
        return {CUBLAS_DIAG_NON_UNIT, CUBLAS_DIAG_UNIT};
    }
    if (value == "NONUNIT" || value == "NON_UNIT" || value == "N") {
        return {CUBLAS_DIAG_NON_UNIT};
    }
    if (value == "UNIT" || value == "U") {
        return {CUBLAS_DIAG_UNIT};
    }

    throw std::runtime_error("invalid diag option: " + value);
}

template <typename F>
inline void for_each_gemm_param(
    const std::vector<cublasOperation_t> &trans_A_list,
    const std::vector<cublasOperation_t> &trans_B_list,
    F f //
) {
    for (auto trans_A : trans_A_list) {
        for (auto trans_B : trans_B_list) {
            f(trans_A, trans_B);
        }
    }
}

template <typename F>
inline void for_each_side_uplo_param(
    const std::vector<cublasSideMode_t> &side_list,
    const std::vector<cublasFillMode_t> &uplo_list,
    F f //
) {
    for (auto side : side_list) {
        for (auto uplo : uplo_list) {
            f(uplo, side);
        }
    }
}

template <typename F>
inline void for_each_syr_param(
    const std::vector<cublasFillMode_t> &uplo_list,
    const std::vector<cublasOperation_t> &trans_list,
    F f //
) {
    for (auto uplo : uplo_list) {
        for (auto trans : trans_list) {
            // SYRK/SYR2K/SYRKX use N or T in the current tests.
            if (trans == CUBLAS_OP_C) continue;
            f(uplo, trans);
        }
    }
}

template <typename F>
inline void for_each_her_param(
    const std::vector<cublasFillMode_t> &uplo_list,
    const std::vector<cublasOperation_t> &trans_list,
    F f //
) {
    for (auto uplo : uplo_list) {
        for (auto trans : trans_list) {
            // HERK/HER2K/HERKX use N or C in the current tests.
            if (trans == CUBLAS_OP_T) continue;
            f(uplo, trans);
        }
    }
}

template <typename F>
inline void for_each_tri_param(
    const std::vector<cublasSideMode_t> &side_list,
    const std::vector<cublasFillMode_t> &uplo_list,
    const std::vector<cublasOperation_t> &trans_list,
    const std::vector<cublasDiagType_t> &diag_list,
    F f //
) {
    for (auto side : side_list) {
        for (auto uplo : uplo_list) {
            for (auto trans : trans_list) {
                for (auto diag : diag_list) {
                    f(side, uplo, trans, diag);
                }
            }
        }
    }
}

template <typename F>
inline void for_each_trtrmm_param(
    const std::vector<cublasFillMode_t> &uplo_A_list,
    const std::vector<cublasFillMode_t> &uplo_B_list,
    const std::vector<cublasOperation_t> &trans_A_list,
    const std::vector<cublasOperation_t> &trans_B_list,
    const std::vector<cublasDiagType_t> &diag_A_list,
    const std::vector<cublasDiagType_t> &diag_B_list,
    F f //
) {
    for (auto uplo_A : uplo_A_list) {
        for (auto uplo_B : uplo_B_list) {
            for (auto trans_A : trans_A_list) {
                for (auto trans_B : trans_B_list) {
                    for (auto diag_A : diag_A_list) {
                        for (auto diag_B : diag_B_list) {
                            f(uplo_A, uplo_B, trans_A, trans_B, diag_A, diag_B);
                        }
                    }
                }
            }
        }
    }
}

void print_options(const char *prog) {
    std::cout
        << "\nUsage:\n"
        << "  " << prog << " <test-option>... <routine-option>... [disable-option]...\n"
        << "\n"
        << "Test options:\n"
        << "  accuracy_square           Run accuracy tests for square matrices (n = 8192)\n"
        << "  accuracy_rectangle        Run accuracy tests for rectangular matrices\n"
        << "  time_square               Run timing tests for square matrices\n"
        << "  time_rectangle            Run timing tests for rectangular matrices\n"
        << "\n"
        << "Routine options:\n"
        << "  GEMM                      Run GEMM\n"
        << "  SYMM                      Run SYMM\n"
        << "  SYRK                      Run SYRK\n"
        << "  SYR2K                     Run SYR2K\n"
        << "  SYRKX                     Run SYRKX\n"
        << "  HEMM                      Run HEMM\n"
        << "  HERK                      Run HERK\n"
        << "  HER2K                     Run HER2K\n"
        << "  HERKX                     Run HERKX\n"
        << "  TRMM                      Run TRMM\n"
        << "  TRTRMM                    Run TRTRMM\n"
        << "  TRSM                      Run TRSM\n"
        << "\n"
        << "Precision options:\n"
        << "  S                         Run single-precision operations\n"
        << "  D                         Run double-precision operations\n"
        << "  C                         Run single-precision complex operations\n"
        << "  Z                         Run double-precision complex operations\n"
        << "\n"
        << "BLAS parameter options:\n"
        << "  trans=all|N|T|C           for SYRK/SYR2K/SYRKX/HERK/HER2K/HERKX/TRSM/TRMM\n"
        << "  transA=all|N|T|C          for GEMM/TRTRMM\n"
        << "  transB=all|N|T|C          for GEMM/TRTRMM\n"
        << "  uplo=all|upper|lower      for SYMM/SYRK/SYR2K/SYRKX/HEMM/HERK/HER2K/HERKX/TRSM/TRMM\n"
        << "  uploA=all|upper|lower     for TRTRMM\n"
        << "  uploB=all|upper|lower     for TRTRMM\n"
        << "  diag=all|nonunit|unit     for TRSM/TRMM\n"
        << "  diag_A=all|nonunit|unit   for TRTRMM\n"
        << "  diag_B=all|nonunit|unit   for TRTRMM\n"
        << "  side=all|left|right       for SYMM/HEMM/TRSM/TRMM\n"
        << "\n"
        << "Disable options:\n"
        << "  no_Ozaki2_INT8            Disable Ozaki-II INT8\n"
        << "  no_Ozaki2_FP8             Disable Ozaki-II FP8\n"
        << "  no_Ozaki1_INT8            Disable Ozaki-I INT8\n"
        << "\n"
        << "Help options:\n"
        << "  -h, --help, help Show this message and exit\n"
        << "\n"
        << "Examples:\n"
        << "  " << prog << " accuracy_rectangle GEMM C Z transA=N transB=N\n"
        << "  " << prog << " time_square GEMM S D no_Ozaki1_INT8\n"
        << "\n";
}

int main(int argc, char **argv) {
    if (argc == 1) {
        print_options(argv[0]);
        return 0;
    }

    std::chrono::system_clock::time_point start, stop;
    std::string startTime  = getCurrentDateTime(start);
    std::string deviceName = printEnvironmentInfo(startTime);

    bool isHopper = false;
#if defined(__CUDACC__) || defined(__NVCC__)
    auto [cc_major, cc_minor] = getComputeCapability(-1);
    isHopper                  = cc_major == 9;
#endif

    bool run_accuracy_rec = false;
    bool run_accuracy_sqr = false;
    bool run_time_rec     = false;
    bool run_time_sqr     = false;

    bool run_S = false;
    bool run_D = false;
    bool run_C = false;
    bool run_Z = false;

    bool run_Ozaki2_I8 = true;
    bool run_Ozaki2_F8 = true;
    bool run_Ozaki1_I8 = true;

    bool run_GEMM   = false;
    bool run_SYMM   = false;
    bool run_SYRK   = false;
    bool run_SYR2K  = false;
    bool run_SYRKX  = false;
    bool run_HEMM   = false;
    bool run_HERK   = false;
    bool run_HER2K  = false;
    bool run_HERKX  = false;
    bool run_TRMM   = false;
    bool run_TRTRMM = false;
    bool run_TRSM   = false;

    std::vector<cublasOperation_t> trans_list   = {CUBLAS_OP_N, CUBLAS_OP_T, CUBLAS_OP_C};
    std::vector<cublasOperation_t> trans_A_list = {CUBLAS_OP_N, CUBLAS_OP_T, CUBLAS_OP_C};
    std::vector<cublasOperation_t> trans_B_list = {CUBLAS_OP_N, CUBLAS_OP_T, CUBLAS_OP_C};
    std::vector<cublasFillMode_t> uplo_list     = {CUBLAS_FILL_MODE_UPPER, CUBLAS_FILL_MODE_LOWER};
    std::vector<cublasFillMode_t> uplo_A_list   = {CUBLAS_FILL_MODE_UPPER, CUBLAS_FILL_MODE_LOWER};
    std::vector<cublasFillMode_t> uplo_B_list   = {CUBLAS_FILL_MODE_UPPER, CUBLAS_FILL_MODE_LOWER};
    std::vector<cublasDiagType_t> diag_list     = {CUBLAS_DIAG_NON_UNIT, CUBLAS_DIAG_UNIT};
    std::vector<cublasDiagType_t> diag_A_list   = {CUBLAS_DIAG_NON_UNIT, CUBLAS_DIAG_UNIT};
    std::vector<cublasDiagType_t> diag_B_list   = {CUBLAS_DIAG_NON_UNIT, CUBLAS_DIAG_UNIT};
    std::vector<cublasSideMode_t> side_list     = {CUBLAS_SIDE_LEFT, CUBLAS_SIDE_RIGHT};

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help" || arg == "help") {
            print_options(argv[0]);
            return 0;
        }

        std::string key;
        std::string value;

        if (split_key_value(arg, key, value)) {
            try {
                if (key == "UPLO") {
                    uplo_list   = parse_uplo_list(value);
                    uplo_A_list = uplo_list;
                    uplo_B_list = uplo_list;
                    continue;
                }
                if (key == "UPLO_A" || key == "UPLOA") {
                    uplo_A_list = parse_uplo_list(value);
                    continue;
                }
                if (key == "UPLO_B" || key == "UPLOB") {
                    uplo_B_list = parse_uplo_list(value);
                    continue;
                }

                if (key == "SIDE") {
                    side_list = parse_side_list(value);
                    continue;
                }

                if (key == "TRANS") {
                    trans_list   = parse_trans_list(value);
                    trans_A_list = trans_list;
                    trans_B_list = trans_list;
                    continue;
                }
                if (key == "TRANS_A" || key == "TRANSA") {
                    trans_A_list = parse_trans_list(value);
                    continue;
                }
                if (key == "TRANS_B" || key == "TRANSB") {
                    trans_B_list = parse_trans_list(value);
                    continue;
                }

                if (key == "DIAG") {
                    diag_list   = parse_diag_list(value);
                    diag_A_list = diag_list;
                    diag_B_list = diag_list;
                    continue;
                }
                if (key == "DIAG_A" || key == "DIAGA") {
                    diag_A_list = parse_diag_list(value);
                    continue;
                }
                if (key == "DIAG_B" || key == "DIAGB") {
                    diag_B_list = parse_diag_list(value);
                    continue;
                }
            } catch (const std::exception &e) {
                std::cerr << e.what() << "\n\n";
                print_options(argv[0]);
                return 1;
            }

            std::cerr << "Unknown option: " << arg << "\n\n";
            print_options(argv[0]);
            return 1;
        }

        if (arg == "accuracy_square") {
            run_accuracy_sqr = true;
            continue;
        }
        if (arg == "accuracy_rectangle") {
            run_accuracy_rec = true;
            continue;
        }
        if (arg == "time_square") {
            run_time_sqr = true;
            continue;
        }
        if (arg == "time_rectangle") {
            run_time_rec = true;
            continue;
        }

        if (arg == "S") {
            run_S = true;
            continue;
        }
        if (arg == "D") {
            run_D = true;
            continue;
        }
        if (arg == "C") {
            run_C = true;
            continue;
        }
        if (arg == "Z") {
            run_Z = true;
            continue;
        }

        if (arg == "no_Ozaki2_INT8") {
            run_Ozaki2_I8 = false;
            continue;
        }
        if (arg == "no_Ozaki2_FP8") {
            run_Ozaki2_F8 = false;
            continue;
        }
        if (arg == "no_Ozaki1_INT8") {
            run_Ozaki1_I8 = false;
            continue;
        }

        if (arg == "GEMM") {
            run_GEMM = true;
            continue;
        }
        if (arg == "SYMM") {
            run_SYMM = true;
            continue;
        }
        if (arg == "SYRK") {
            run_SYRK = true;
            continue;
        }
        if (arg == "SYR2K") {
            run_SYR2K = true;
            continue;
        }
        if (arg == "SYRKX") {
            run_SYRKX = true;
            continue;
        }
        if (arg == "HEMM") {
            run_HEMM = true;
            continue;
        }
        if (arg == "HERK") {
            run_HERK = true;
            continue;
        }
        if (arg == "HER2K") {
            run_HER2K = true;
            continue;
        }
        if (arg == "HERKX") {
            run_HERKX = true;
            continue;
        }
        if (arg == "TRMM") {
            run_TRMM = true;
            continue;
        }
        if (arg == "TRTRMM") {
            run_TRTRMM = true;
            continue;
        }
        if (arg == "TRSM") {
            run_TRSM = true;
            continue;
        }

        std::cerr << "Unknown option: " << arg << "\n\n";
        print_options(argv[0]);
        return 1;
    }

    run_Ozaki2_F8 = (isHopper) ? false : run_Ozaki2_F8;

    if (run_accuracy_sqr) {

        if (run_GEMM) {
            auto run_gemm_accuracy = [&](cublasOperation_t transa, cublasOperation_t transb) {
                if (transa != CUBLAS_OP_C && transb != CUBLAS_OP_C) {
                    if (run_S) bench::accuracy::gemm::check_accuracy<float>(deviceName, startTime, transa, transb, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, true);
                    if (run_D) bench::accuracy::gemm::check_accuracy<double>(deviceName, startTime, transa, transb, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, true);
                }
                if (run_C) bench::accuracy::gemm::check_accuracy<cuFloatComplex>(deviceName, startTime, transa, transb, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, true);
                if (run_Z) bench::accuracy::gemm::check_accuracy<cuDoubleComplex>(deviceName, startTime, transa, transb, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, true);
            };

            for_each_gemm_param(
                trans_A_list,
                trans_B_list,
                run_gemm_accuracy);
        }

        if (run_SYMM) {
            auto run_symm_accuracy = [&](cublasFillMode_t uplo, cublasSideMode_t side) {
                if (run_S) bench::accuracy::symm::check_accuracy<float>(deviceName, startTime, uplo, side, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, true);
                if (run_D) bench::accuracy::symm::check_accuracy<double>(deviceName, startTime, uplo, side, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, true);
                if (run_C) bench::accuracy::symm::check_accuracy<cuFloatComplex>(deviceName, startTime, uplo, side, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, true);
                if (run_Z) bench::accuracy::symm::check_accuracy<cuDoubleComplex>(deviceName, startTime, uplo, side, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, true);
            };

            for_each_side_uplo_param(
                side_list,
                uplo_list,
                run_symm_accuracy);
        }

        if (run_SYRK) {
            auto run_syrk_accuracy = [&](cublasFillMode_t uplo, cublasOperation_t trans) {
                if (run_S) bench::accuracy::syrk::check_accuracy<float>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, true);
                if (run_D) bench::accuracy::syrk::check_accuracy<double>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, true);
                if (run_C) bench::accuracy::syrk::check_accuracy<cuFloatComplex>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, true);
                if (run_Z) bench::accuracy::syrk::check_accuracy<cuDoubleComplex>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, true);
            };

            for_each_syr_param(
                uplo_list,
                trans_list,
                run_syrk_accuracy);
        }

        if (run_SYR2K) {
            auto run_syr2k_accuracy = [&](cublasFillMode_t uplo, cublasOperation_t trans) {
                if (run_S) bench::accuracy::syr2k::check_accuracy<float>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, true);
                if (run_D) bench::accuracy::syr2k::check_accuracy<double>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, true);
                if (run_C) bench::accuracy::syr2k::check_accuracy<cuFloatComplex>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, true);
                if (run_Z) bench::accuracy::syr2k::check_accuracy<cuDoubleComplex>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, true);
            };

            for_each_syr_param(
                uplo_list,
                trans_list,
                run_syr2k_accuracy);
        }

        if (run_SYRKX) {
            auto run_syrkx_accuracy = [&](cublasFillMode_t uplo, cublasOperation_t trans) {
                if (run_S) bench::accuracy::syrkx::check_accuracy<float>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, true);
                if (run_D) bench::accuracy::syrkx::check_accuracy<double>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, true);
                if (run_C) bench::accuracy::syrkx::check_accuracy<cuFloatComplex>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, true);
                if (run_Z) bench::accuracy::syrkx::check_accuracy<cuDoubleComplex>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, true);
            };

            for_each_syr_param(
                uplo_list,
                trans_list,
                run_syrkx_accuracy);
        }

        if (run_HEMM) {
            auto run_hemm_accuracy = [&](cublasFillMode_t uplo, cublasSideMode_t side) {
                if (run_C) bench::accuracy::hemm::check_accuracy<cuFloatComplex>(deviceName, startTime, uplo, side, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, true);
                if (run_Z) bench::accuracy::hemm::check_accuracy<cuDoubleComplex>(deviceName, startTime, uplo, side, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, true);
            };

            for_each_side_uplo_param(
                side_list,
                uplo_list,
                run_hemm_accuracy);
        }

        if (run_HERK) {
            auto run_herk_accuracy = [&](cublasFillMode_t uplo, cublasOperation_t trans) {
                if (run_C) bench::accuracy::herk::check_accuracy<cuFloatComplex>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, true);
                if (run_Z) bench::accuracy::herk::check_accuracy<cuDoubleComplex>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, true);
            };

            for_each_her_param(
                uplo_list,
                trans_list,
                run_herk_accuracy);
        }

        if (run_HER2K) {
            auto run_her2k_accuracy = [&](cublasFillMode_t uplo, cublasOperation_t trans) {
                if (run_C) bench::accuracy::her2k::check_accuracy<cuFloatComplex>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, true);
                if (run_Z) bench::accuracy::her2k::check_accuracy<cuDoubleComplex>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, true);
            };

            for_each_her_param(
                uplo_list,
                trans_list,
                run_her2k_accuracy);
        }

        if (run_HERKX) {
            auto run_herkx_accuracy = [&](cublasFillMode_t uplo, cublasOperation_t trans) {
                if (run_C) bench::accuracy::herkx::check_accuracy<cuFloatComplex>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, true);
                if (run_Z) bench::accuracy::herkx::check_accuracy<cuDoubleComplex>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, true);
            };

            for_each_her_param(
                uplo_list,
                trans_list,
                run_herkx_accuracy);
        }

        if (run_TRMM) {
            auto run_trmm_accuracy = [&](cublasSideMode_t side, cublasFillMode_t uplo, cublasOperation_t trans, cublasDiagType_t diag) {
                if (trans != CUBLAS_OP_C) {
                    if (run_S) bench::accuracy::trmm::check_accuracy<float>(deviceName, startTime, side, uplo, trans, diag, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, true);
                    if (run_D) bench::accuracy::trmm::check_accuracy<double>(deviceName, startTime, side, uplo, trans, diag, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, true);
                }
                if (run_C) bench::accuracy::trmm::check_accuracy<cuFloatComplex>(deviceName, startTime, side, uplo, trans, diag, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, true);
                if (run_Z) bench::accuracy::trmm::check_accuracy<cuDoubleComplex>(deviceName, startTime, side, uplo, trans, diag, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, true);
            };

            for_each_tri_param(
                side_list,
                uplo_list,
                trans_list,
                diag_list,
                run_trmm_accuracy);
        }

        if (run_TRSM) {
            auto run_trsm_accuracy = [&](cublasSideMode_t side, cublasFillMode_t uplo, cublasOperation_t trans, cublasDiagType_t diag, bool is_square) {
                if (trans != CUBLAS_OP_C) {
                    if (run_S) bench::accuracy::trsm::check_accuracy<float>(deviceName, startTime, side, uplo, trans, diag, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, is_square);
                    if (run_D) bench::accuracy::trsm::check_accuracy<double>(deviceName, startTime, side, uplo, trans, diag, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, is_square);
                }
                if (run_C) bench::accuracy::trsm::check_accuracy<cuFloatComplex>(deviceName, startTime, side, uplo, trans, diag, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, is_square);
                if (run_Z) bench::accuracy::trsm::check_accuracy<cuDoubleComplex>(deviceName, startTime, side, uplo, trans, diag, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, is_square);
            };

            for_each_tri_param(
                side_list,
                uplo_list,
                trans_list,
                diag_list,
                [&](cublasSideMode_t side,
                    cublasFillMode_t uplo,
                    cublasOperation_t trans,
                    cublasDiagType_t diag) {
                    run_trsm_accuracy(side, uplo, trans, diag, true);
                });
        }

        if (run_TRTRMM) {
            auto run_trtrmm_accuracy = [&](cublasFillMode_t uplo_A, cublasFillMode_t uplo_B, cublasOperation_t trans_A, cublasOperation_t trans_B, cublasDiagType_t diag_A, cublasDiagType_t diag_B) {
                if (trans_A != CUBLAS_OP_C && trans_B != CUBLAS_OP_C) {
                    if (run_S) bench::accuracy::trtrmm::check_accuracy<float>(deviceName, startTime, uplo_A, uplo_B, trans_A, trans_B, diag_A, diag_B, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
                    if (run_D) bench::accuracy::trtrmm::check_accuracy<double>(deviceName, startTime, uplo_A, uplo_B, trans_A, trans_B, diag_A, diag_B, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
                }
                if (run_C) bench::accuracy::trtrmm::check_accuracy<cuFloatComplex>(deviceName, startTime, uplo_A, uplo_B, trans_A, trans_B, diag_A, diag_B, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
                if (run_Z) bench::accuracy::trtrmm::check_accuracy<cuDoubleComplex>(deviceName, startTime, uplo_A, uplo_B, trans_A, trans_B, diag_A, diag_B, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
            };

            for_each_trtrmm_param(
                uplo_A_list,
                uplo_B_list,
                trans_A_list,
                trans_B_list,
                diag_A_list,
                diag_B_list,
                run_trtrmm_accuracy);
        }
    }

    if (run_accuracy_rec) {

        if (run_GEMM) {
            auto run_gemm_accuracy = [&](cublasOperation_t transa, cublasOperation_t transb) {
                if (transa != CUBLAS_OP_C && transb != CUBLAS_OP_C) {
                    if (run_S) bench::accuracy::gemm::check_accuracy<float>(deviceName, startTime, transa, transb, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
                    if (run_D) bench::accuracy::gemm::check_accuracy<double>(deviceName, startTime, transa, transb, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
                }
                if (run_C) bench::accuracy::gemm::check_accuracy<cuFloatComplex>(deviceName, startTime, transa, transb, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
                if (run_Z) bench::accuracy::gemm::check_accuracy<cuDoubleComplex>(deviceName, startTime, transa, transb, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
            };

            for_each_gemm_param(
                trans_A_list,
                trans_B_list,
                run_gemm_accuracy);
        }

        if (run_SYMM) {
            auto run_symm_accuracy = [&](cublasFillMode_t uplo, cublasSideMode_t side) {
                if (run_S) bench::accuracy::symm::check_accuracy<float>(deviceName, startTime, uplo, side, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
                if (run_D) bench::accuracy::symm::check_accuracy<double>(deviceName, startTime, uplo, side, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
                if (run_C) bench::accuracy::symm::check_accuracy<cuFloatComplex>(deviceName, startTime, uplo, side, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
                if (run_Z) bench::accuracy::symm::check_accuracy<cuDoubleComplex>(deviceName, startTime, uplo, side, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
            };

            for_each_side_uplo_param(
                side_list,
                uplo_list,
                run_symm_accuracy);
        }

        if (run_SYRK) {
            auto run_syrk_accuracy = [&](cublasFillMode_t uplo, cublasOperation_t trans) {
                if (run_S) bench::accuracy::syrk::check_accuracy<float>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
                if (run_D) bench::accuracy::syrk::check_accuracy<double>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
                if (run_C) bench::accuracy::syrk::check_accuracy<cuFloatComplex>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
                if (run_Z) bench::accuracy::syrk::check_accuracy<cuDoubleComplex>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
            };

            for_each_syr_param(
                uplo_list,
                trans_list,
                run_syrk_accuracy);
        }

        if (run_SYR2K) {
            auto run_syr2k_accuracy = [&](cublasFillMode_t uplo, cublasOperation_t trans) {
                if (run_S) bench::accuracy::syr2k::check_accuracy<float>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
                if (run_D) bench::accuracy::syr2k::check_accuracy<double>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
                if (run_C) bench::accuracy::syr2k::check_accuracy<cuFloatComplex>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
                if (run_Z) bench::accuracy::syr2k::check_accuracy<cuDoubleComplex>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
            };

            for_each_syr_param(
                uplo_list,
                trans_list,
                run_syr2k_accuracy);
        }

        if (run_SYRKX) {
            auto run_syrkx_accuracy = [&](cublasFillMode_t uplo, cublasOperation_t trans) {
                if (run_S) bench::accuracy::syrkx::check_accuracy<float>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
                if (run_D) bench::accuracy::syrkx::check_accuracy<double>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
                if (run_C) bench::accuracy::syrkx::check_accuracy<cuFloatComplex>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
                if (run_Z) bench::accuracy::syrkx::check_accuracy<cuDoubleComplex>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
            };

            for_each_syr_param(
                uplo_list,
                trans_list,
                run_syrkx_accuracy);
        }

        if (run_HEMM) {
            auto run_hemm_accuracy = [&](cublasFillMode_t uplo, cublasSideMode_t side) {
                if (run_C) bench::accuracy::hemm::check_accuracy<cuFloatComplex>(deviceName, startTime, uplo, side, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
                if (run_Z) bench::accuracy::hemm::check_accuracy<cuDoubleComplex>(deviceName, startTime, uplo, side, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
            };

            for_each_side_uplo_param(
                side_list,
                uplo_list,
                run_hemm_accuracy);
        }

        if (run_HERK) {
            auto run_herk_accuracy = [&](cublasFillMode_t uplo, cublasOperation_t trans) {
                if (run_C) bench::accuracy::herk::check_accuracy<cuFloatComplex>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
                if (run_Z) bench::accuracy::herk::check_accuracy<cuDoubleComplex>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
            };

            for_each_her_param(
                uplo_list,
                trans_list,
                run_herk_accuracy);
        }

        if (run_HER2K) {
            auto run_her2k_accuracy = [&](cublasFillMode_t uplo, cublasOperation_t trans) {
                if (run_C) bench::accuracy::her2k::check_accuracy<cuFloatComplex>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
                if (run_Z) bench::accuracy::her2k::check_accuracy<cuDoubleComplex>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
            };

            for_each_her_param(
                uplo_list,
                trans_list,
                run_her2k_accuracy);
        }

        if (run_HERKX) {
            auto run_herkx_accuracy = [&](cublasFillMode_t uplo, cublasOperation_t trans) {
                if (run_C) bench::accuracy::herkx::check_accuracy<cuFloatComplex>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
                if (run_Z) bench::accuracy::herkx::check_accuracy<cuDoubleComplex>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
            };

            for_each_her_param(
                uplo_list,
                trans_list,
                run_herkx_accuracy);
        }

        if (run_TRMM) {
            auto run_trmm_accuracy = [&](cublasSideMode_t side, cublasFillMode_t uplo, cublasOperation_t trans, cublasDiagType_t diag) {
                if (trans != CUBLAS_OP_C) {
                    if (run_S) bench::accuracy::trmm::check_accuracy<float>(deviceName, startTime, side, uplo, trans, diag, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
                    if (run_D) bench::accuracy::trmm::check_accuracy<double>(deviceName, startTime, side, uplo, trans, diag, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
                }
                if (run_C) bench::accuracy::trmm::check_accuracy<cuFloatComplex>(deviceName, startTime, side, uplo, trans, diag, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
                if (run_Z) bench::accuracy::trmm::check_accuracy<cuDoubleComplex>(deviceName, startTime, side, uplo, trans, diag, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
            };

            for_each_tri_param(
                side_list,
                uplo_list,
                trans_list,
                diag_list,
                run_trmm_accuracy);
        }

        if (run_TRSM) {
            auto run_trsm_accuracy = [&](cublasSideMode_t side, cublasFillMode_t uplo, cublasOperation_t trans, cublasDiagType_t diag, bool is_square) {
                if (trans != CUBLAS_OP_C) {
                    if (run_S) bench::accuracy::trsm::check_accuracy<float>(deviceName, startTime, side, uplo, trans, diag, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, is_square);
                    if (run_D) bench::accuracy::trsm::check_accuracy<double>(deviceName, startTime, side, uplo, trans, diag, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, is_square);
                }
                if (run_C) bench::accuracy::trsm::check_accuracy<cuFloatComplex>(deviceName, startTime, side, uplo, trans, diag, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, is_square);
                if (run_Z) bench::accuracy::trsm::check_accuracy<cuDoubleComplex>(deviceName, startTime, side, uplo, trans, diag, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, is_square);
            };

            for_each_tri_param(
                side_list,
                uplo_list,
                trans_list,
                diag_list,
                [&](cublasSideMode_t side,
                    cublasFillMode_t uplo,
                    cublasOperation_t trans,
                    cublasDiagType_t diag) {
                    run_trsm_accuracy(side, uplo, trans, diag, false);
                });
        }

        if (run_TRTRMM) {
            auto run_trtrmm_accuracy = [&](cublasFillMode_t uplo_A, cublasFillMode_t uplo_B, cublasOperation_t trans_A, cublasOperation_t trans_B, cublasDiagType_t diag_A, cublasDiagType_t diag_B) {
                if (trans_A != CUBLAS_OP_C && trans_B != CUBLAS_OP_C) {
                    if (run_S) bench::accuracy::trtrmm::check_accuracy<float>(deviceName, startTime, uplo_A, uplo_B, trans_A, trans_B, diag_A, diag_B, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
                    if (run_D) bench::accuracy::trtrmm::check_accuracy<double>(deviceName, startTime, uplo_A, uplo_B, trans_A, trans_B, diag_A, diag_B, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
                }
                if (run_C) bench::accuracy::trtrmm::check_accuracy<cuFloatComplex>(deviceName, startTime, uplo_A, uplo_B, trans_A, trans_B, diag_A, diag_B, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
                if (run_Z) bench::accuracy::trtrmm::check_accuracy<cuDoubleComplex>(deviceName, startTime, uplo_A, uplo_B, trans_A, trans_B, diag_A, diag_B, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
            };

            for_each_trtrmm_param(
                uplo_A_list,
                uplo_B_list,
                trans_A_list,
                trans_B_list,
                diag_A_list,
                diag_B_list,
                run_trtrmm_accuracy);
        }
    }

    if (run_time_sqr) {

        if (run_GEMM) {
            auto run_gemm_time = [&](cublasOperation_t transa, cublasOperation_t transb) {
                if (transa != CUBLAS_OP_C && transb != CUBLAS_OP_C) {
                    if (run_S) bench::time::gemm::check_time<float>(deviceName, startTime, transa, transb, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, true);
                    if (run_D) bench::time::gemm::check_time<double>(deviceName, startTime, transa, transb, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, true);
                }
                if (run_C) bench::time::gemm::check_time<cuFloatComplex>(deviceName, startTime, transa, transb, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, true);
                if (run_Z) bench::time::gemm::check_time<cuDoubleComplex>(deviceName, startTime, transa, transb, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, true);
            };

            for_each_gemm_param(
                trans_A_list,
                trans_B_list,
                run_gemm_time);
        }

        if (run_SYMM) {
            auto run_symm_time = [&](cublasFillMode_t uplo, cublasSideMode_t side) {
                if (run_S) bench::time::symm::check_time<float>(deviceName, startTime, uplo, side, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, true);
                if (run_D) bench::time::symm::check_time<double>(deviceName, startTime, uplo, side, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, true);
                if (run_C) bench::time::symm::check_time<cuFloatComplex>(deviceName, startTime, uplo, side, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, true);
                if (run_Z) bench::time::symm::check_time<cuDoubleComplex>(deviceName, startTime, uplo, side, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, true);
            };

            for (auto side : side_list) {
                for (auto uplo : uplo_list) {
                    run_symm_time(uplo, side);
                }
            }
        }

        if (run_SYRK) {
            auto run_syrk_time = [&](cublasFillMode_t uplo, cublasOperation_t trans) {
                if (run_S) bench::time::syrk::check_time<float>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, true);
                if (run_D) bench::time::syrk::check_time<double>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, true);
                if (run_C) bench::time::syrk::check_time<cuFloatComplex>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, true);
                if (run_Z) bench::time::syrk::check_time<cuDoubleComplex>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, true);
            };

            for_each_syr_param(
                uplo_list,
                trans_list,
                run_syrk_time);
        }

        if (run_SYR2K) {
            auto run_syr2k_time = [&](cublasFillMode_t uplo, cublasOperation_t trans) {
                if (run_S) bench::time::syr2k::check_time<float>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, true);
                if (run_D) bench::time::syr2k::check_time<double>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, true);
                if (run_C) bench::time::syr2k::check_time<cuFloatComplex>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, true);
                if (run_Z) bench::time::syr2k::check_time<cuDoubleComplex>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, true);
            };

            for_each_syr_param(
                uplo_list,
                trans_list,
                run_syr2k_time);
        }

        if (run_SYRKX) {
            auto run_syrkx_time = [&](cublasFillMode_t uplo, cublasOperation_t trans) {
                if (run_S) bench::time::syrkx::check_time<float>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, true);
                if (run_D) bench::time::syrkx::check_time<double>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, true);
                if (run_C) bench::time::syrkx::check_time<cuFloatComplex>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, true);
                if (run_Z) bench::time::syrkx::check_time<cuDoubleComplex>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, true);
            };

            for_each_syr_param(
                uplo_list,
                trans_list,
                run_syrkx_time);
        }

        if (run_HEMM) {
            auto run_hemm_time = [&](cublasFillMode_t uplo, cublasSideMode_t side) {
                if (run_C) bench::time::hemm::check_time<cuFloatComplex>(deviceName, startTime, uplo, side, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, true);
                if (run_Z) bench::time::hemm::check_time<cuDoubleComplex>(deviceName, startTime, uplo, side, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, true);
            };

            for_each_side_uplo_param(
                side_list,
                uplo_list,
                run_hemm_time);
        }

        if (run_HERK) {
            auto run_herk_time = [&](cublasFillMode_t uplo, cublasOperation_t trans) {
                if (run_C) bench::time::herk::check_time<cuFloatComplex>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, true);
                if (run_Z) bench::time::herk::check_time<cuDoubleComplex>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, true);
            };

            for_each_her_param(
                uplo_list,
                trans_list,
                run_herk_time);
        }

        if (run_HER2K) {
            auto run_her2k_time = [&](cublasFillMode_t uplo, cublasOperation_t trans) {
                if (run_C) bench::time::her2k::check_time<cuFloatComplex>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, true);
                if (run_Z) bench::time::her2k::check_time<cuDoubleComplex>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, true);
            };

            for_each_her_param(
                uplo_list,
                trans_list,
                run_her2k_time);
        }

        if (run_HERKX) {
            auto run_herkx_time = [&](cublasFillMode_t uplo, cublasOperation_t trans) {
                if (run_C) bench::time::herkx::check_time<cuFloatComplex>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, true);
                if (run_Z) bench::time::herkx::check_time<cuDoubleComplex>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, true);
            };

            for_each_her_param(
                uplo_list,
                trans_list,
                run_herkx_time);
        }

        if (run_TRMM) {
            auto run_trmm_time = [&](cublasSideMode_t side, cublasFillMode_t uplo, cublasOperation_t trans, cublasDiagType_t diag) {
                if (trans != CUBLAS_OP_C) {
                    if (run_S) bench::time::trmm::check_time<float>(deviceName, startTime, side, uplo, trans, diag, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, true);
                    if (run_D) bench::time::trmm::check_time<double>(deviceName, startTime, side, uplo, trans, diag, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, true);
                }
                if (run_C) bench::time::trmm::check_time<cuFloatComplex>(deviceName, startTime, side, uplo, trans, diag, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, true);
                if (run_Z) bench::time::trmm::check_time<cuDoubleComplex>(deviceName, startTime, side, uplo, trans, diag, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, true);
            };

            for_each_tri_param(
                side_list,
                uplo_list,
                trans_list,
                diag_list,
                run_trmm_time);
        }

        if (run_TRSM) {
            auto run_trsm_time = [&](cublasSideMode_t side, cublasFillMode_t uplo, cublasOperation_t trans, cublasDiagType_t diag) {
                if (trans != CUBLAS_OP_C) {
                    if (run_S) bench::time::trsm::check_time<float>(deviceName, startTime, side, uplo, trans, diag, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, true);
                    if (run_D) bench::time::trsm::check_time<double>(deviceName, startTime, side, uplo, trans, diag, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, true);
                }
                if (run_C) bench::time::trsm::check_time<cuFloatComplex>(deviceName, startTime, side, uplo, trans, diag, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, true);
                if (run_Z) bench::time::trsm::check_time<cuDoubleComplex>(deviceName, startTime, side, uplo, trans, diag, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8, true);
            };

            for_each_tri_param(
                side_list,
                uplo_list,
                trans_list,
                diag_list,
                run_trsm_time);
        }

        if (run_TRTRMM) {
            auto run_trtrmm_time = [&](cublasFillMode_t uplo_A, cublasFillMode_t uplo_B, cublasOperation_t trans_A, cublasOperation_t trans_B, cublasDiagType_t diag_A, cublasDiagType_t diag_B) {
                if (trans_A != CUBLAS_OP_C && trans_B != CUBLAS_OP_C) {
                    if (run_S) bench::time::trtrmm::check_time<float>(deviceName, startTime, uplo_A, uplo_B, trans_A, trans_B, diag_A, diag_B, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
                    if (run_D) bench::time::trtrmm::check_time<double>(deviceName, startTime, uplo_A, uplo_B, trans_A, trans_B, diag_A, diag_B, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
                }
                if (run_C) bench::time::trtrmm::check_time<cuFloatComplex>(deviceName, startTime, uplo_A, uplo_B, trans_A, trans_B, diag_A, diag_B, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
                if (run_Z) bench::time::trtrmm::check_time<cuDoubleComplex>(deviceName, startTime, uplo_A, uplo_B, trans_A, trans_B, diag_A, diag_B, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
            };

            for_each_trtrmm_param(
                uplo_A_list,
                uplo_B_list,
                trans_A_list,
                trans_B_list,
                diag_A_list,
                diag_B_list,
                run_trtrmm_time);
        }
    }

    if (run_time_rec) {

        if (run_GEMM) {
            auto run_gemm_time = [&](cublasOperation_t transa, cublasOperation_t transb) {
                if (transa != CUBLAS_OP_C && transb != CUBLAS_OP_C) {
                    if (run_S) bench::time::gemm::check_time<float>(deviceName, startTime, transa, transb, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
                    if (run_D) bench::time::gemm::check_time<double>(deviceName, startTime, transa, transb, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
                }
                if (run_C) bench::time::gemm::check_time<cuFloatComplex>(deviceName, startTime, transa, transb, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
                if (run_Z) bench::time::gemm::check_time<cuDoubleComplex>(deviceName, startTime, transa, transb, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
            };

            for_each_gemm_param(
                trans_A_list,
                trans_B_list,
                run_gemm_time);
        }

        if (run_SYMM) {
            auto run_symm_time = [&](cublasFillMode_t uplo, cublasSideMode_t side) {
                if (run_S) bench::time::symm::check_time<float>(deviceName, startTime, uplo, side, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
                if (run_D) bench::time::symm::check_time<double>(deviceName, startTime, uplo, side, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
                if (run_C) bench::time::symm::check_time<cuFloatComplex>(deviceName, startTime, uplo, side, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
                if (run_Z) bench::time::symm::check_time<cuDoubleComplex>(deviceName, startTime, uplo, side, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
            };

            for_each_side_uplo_param(
                side_list,
                uplo_list,
                run_symm_time);
        }

        if (run_SYRK) {
            auto run_syrk_time = [&](cublasFillMode_t uplo, cublasOperation_t trans) {
                if (run_S) bench::time::syrk::check_time<float>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
                if (run_D) bench::time::syrk::check_time<double>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
                if (run_C) bench::time::syrk::check_time<cuFloatComplex>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
                if (run_Z) bench::time::syrk::check_time<cuDoubleComplex>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
            };

            for_each_syr_param(
                uplo_list,
                trans_list,
                run_syrk_time);
        }

        if (run_SYR2K) {
            auto run_syr2k_time = [&](cublasFillMode_t uplo, cublasOperation_t trans) {
                if (run_S) bench::time::syr2k::check_time<float>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
                if (run_D) bench::time::syr2k::check_time<double>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
                if (run_C) bench::time::syr2k::check_time<cuFloatComplex>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
                if (run_Z) bench::time::syr2k::check_time<cuDoubleComplex>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
            };

            for_each_syr_param(
                uplo_list,
                trans_list,
                run_syr2k_time);
        }

        if (run_SYRKX) {
            auto run_syrkx_time = [&](cublasFillMode_t uplo, cublasOperation_t trans) {
                if (run_S) bench::time::syrkx::check_time<float>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
                if (run_D) bench::time::syrkx::check_time<double>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
                if (run_C) bench::time::syrkx::check_time<cuFloatComplex>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
                if (run_Z) bench::time::syrkx::check_time<cuDoubleComplex>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
            };

            for_each_syr_param(
                uplo_list,
                trans_list,
                run_syrkx_time);
        }

        if (run_HEMM) {
            auto run_hemm_time = [&](cublasFillMode_t uplo, cublasSideMode_t side) {
                if (run_C) bench::time::hemm::check_time<cuFloatComplex>(deviceName, startTime, uplo, side, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
                if (run_Z) bench::time::hemm::check_time<cuDoubleComplex>(deviceName, startTime, uplo, side, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
            };

            for_each_side_uplo_param(
                side_list,
                uplo_list,
                run_hemm_time);
        }

        if (run_HERK) {
            auto run_herk_time = [&](cublasFillMode_t uplo, cublasOperation_t trans) {
                if (run_C) bench::time::herk::check_time<cuFloatComplex>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
                if (run_Z) bench::time::herk::check_time<cuDoubleComplex>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
            };

            for_each_her_param(
                uplo_list,
                trans_list,
                run_herk_time);
        }

        if (run_HER2K) {
            auto run_her2k_time = [&](cublasFillMode_t uplo, cublasOperation_t trans) {
                if (run_C) bench::time::her2k::check_time<cuFloatComplex>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
                if (run_Z) bench::time::her2k::check_time<cuDoubleComplex>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
            };

            for_each_her_param(
                uplo_list,
                trans_list,
                run_her2k_time);
        }

        if (run_HERKX) {
            auto run_herkx_time = [&](cublasFillMode_t uplo, cublasOperation_t trans) {
                if (run_C) bench::time::herkx::check_time<cuFloatComplex>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
                if (run_Z) bench::time::herkx::check_time<cuDoubleComplex>(deviceName, startTime, uplo, trans, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
            };

            for_each_her_param(
                uplo_list,
                trans_list,
                run_herkx_time);
        }

        if (run_TRMM) {
            auto run_trmm_time = [&](cublasSideMode_t side, cublasFillMode_t uplo, cublasOperation_t trans, cublasDiagType_t diag) {
                if (trans != CUBLAS_OP_C) {
                    if (run_S) bench::time::trmm::check_time<float>(deviceName, startTime, side, uplo, trans, diag, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
                    if (run_D) bench::time::trmm::check_time<double>(deviceName, startTime, side, uplo, trans, diag, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
                }
                if (run_C) bench::time::trmm::check_time<cuFloatComplex>(deviceName, startTime, side, uplo, trans, diag, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
                if (run_Z) bench::time::trmm::check_time<cuDoubleComplex>(deviceName, startTime, side, uplo, trans, diag, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
            };

            for_each_tri_param(
                side_list,
                uplo_list,
                trans_list,
                diag_list,
                run_trmm_time);
        }

        if (run_TRSM) {
            auto run_trsm_time = [&](cublasSideMode_t side, cublasFillMode_t uplo, cublasOperation_t trans, cublasDiagType_t diag) {
                if (trans != CUBLAS_OP_C) {
                    if (run_S) bench::time::trsm::check_time<float>(deviceName, startTime, side, uplo, trans, diag, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
                    if (run_D) bench::time::trsm::check_time<double>(deviceName, startTime, side, uplo, trans, diag, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
                }
                if (run_C) bench::time::trsm::check_time<cuFloatComplex>(deviceName, startTime, side, uplo, trans, diag, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
                if (run_Z) bench::time::trsm::check_time<cuDoubleComplex>(deviceName, startTime, side, uplo, trans, diag, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
            };

            for_each_tri_param(
                side_list,
                uplo_list,
                trans_list,
                diag_list,
                run_trsm_time);
        }

        if (run_TRTRMM) {
            auto run_trtrmm_time = [&](cublasFillMode_t uplo_A, cublasFillMode_t uplo_B, cublasOperation_t trans_A, cublasOperation_t trans_B, cublasDiagType_t diag_A, cublasDiagType_t diag_B) {
                if (trans_A != CUBLAS_OP_C && trans_B != CUBLAS_OP_C) {
                    if (run_S) bench::time::trtrmm::check_time<float>(deviceName, startTime, uplo_A, uplo_B, trans_A, trans_B, diag_A, diag_B, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
                    if (run_D) bench::time::trtrmm::check_time<double>(deviceName, startTime, uplo_A, uplo_B, trans_A, trans_B, diag_A, diag_B, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
                }
                if (run_C) bench::time::trtrmm::check_time<cuFloatComplex>(deviceName, startTime, uplo_A, uplo_B, trans_A, trans_B, diag_A, diag_B, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
                if (run_Z) bench::time::trtrmm::check_time<cuDoubleComplex>(deviceName, startTime, uplo_A, uplo_B, trans_A, trans_B, diag_A, diag_B, run_Ozaki2_I8, run_Ozaki2_F8, run_Ozaki1_I8);
            };

            for_each_trtrmm_param(
                uplo_A_list,
                uplo_B_list,
                trans_A_list,
                trans_B_list,
                diag_A_list,
                diag_B_list,
                run_trtrmm_time);
        }
    }

    std::string endTime = getCurrentDateTime(stop);
    auto sec            = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count() * 1.e-9;
    std::cout << std::endl;
    std::cout << "=========================" << std::endl;
    std::cout << "start        : " << startTime << std::endl;
    std::cout << "end          : " << endTime << std::endl;
    std::cout << "elapsed time : " << sec << " [sec]" << " (" << sec / 60.0 << "[min])" << std::endl;
    std::cout << "=========================" << std::endl;
    std::cout << std::endl;

    CHECK_CUDA(cudaDeviceReset());
    return 0;
}
