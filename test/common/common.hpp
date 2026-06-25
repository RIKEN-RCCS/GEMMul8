#pragma once
#include "gemmul8.hpp"

#include <fstream>
#include <functional>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <string_view>
#include <chrono>
#include <type_traits>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cassert>
#include <utility>
#include <stdexcept>

#include "self_hipify.hpp"
#include "gpu_info.hpp"
#include "eval.hpp"
#include "make_matrix.hpp"

#ifdef PRINT
    #undef PRINT
#endif

#ifdef PRINT_nobreak
    #undef PRINT_nobreak
#endif

#ifdef avail_Ozaki1
    #undef avail_Ozaki1
#endif

#ifdef avail_BF16x9
    #undef avail_BF16x9
#endif

#if defined(__CUDACC__) && defined(CUBLAS_VER_MAJOR) && defined(CUBLAS_VER_MINOR)
    #if (CUBLAS_VER_MAJOR * 100 + CUBLAS_VER_MINOR) >= 1301
        #define avail_Ozaki1 1
    #else
        #define avail_Ozaki1 0
    #endif
    #if (CUBLAS_VER_MAJOR * 100 + CUBLAS_VER_MINOR) >= 1209
        #define avail_BF16x9 1
    #else
        #define avail_BF16x9 0
    #endif
#else
    #define avail_Ozaki1 0
    #define avail_BF16x9 0
#endif

#define PRINT(outFile, LINE)                               \
    do {                                                   \
        (outFile) << std::scientific << LINE << std::endl; \
        std::cout << std::scientific << LINE << std::endl; \
    } while (0)

#define PRINT_nobreak(outFile, LINE)          \
    do {                                      \
        (outFile) << std::scientific << LINE; \
        std::cout << std::scientific << LINE; \
    } while (0)

inline constexpr unsigned warmup_min  = 10;
inline constexpr unsigned warmup_max  = 100;
inline constexpr double warmup_ms_min = 5000.0;

inline constexpr unsigned mainloop_min  = 30;
inline constexpr unsigned mainloop_max  = 100;
inline constexpr double mainloop_ms_min = 15000.0;

inline constexpr unsigned repetitions = 1;

inline constexpr unsigned long long seedA = 12345;
inline constexpr unsigned long long seedB = 54321;

inline std::vector<int> oz1_slice_list = {7, 11};

inline std::vector<size_t> N_list   = {1024, 2048, 4096, 8192, 16384, 32768, 65536};
inline std::vector<double> phi_list = {-1.0, 0.0, 0.5, 1.0, 2.0, 4.0};

template <gemmul8::Backend backend> inline constexpr char backendType = 'f';
template <> inline constexpr char backendType<gemmul8::Backend::INT8> = 'i';

template <typename T> struct testTraits;

template <> struct testTraits<float> {
    static constexpr unsigned NUM_MODULI_MIN          = 6;
    static constexpr unsigned NUM_MODULI_MAX          = 9;
    static constexpr unsigned NUM_MODULI_MIN_accuracy = 3;
    static constexpr unsigned NUM_MODULI_MAX_accuracy = 12;
    static constexpr char prefix                      = 's';
    static constexpr char prefix_upper                = 'S';
    static constexpr bool is_complex                  = false;
    static constexpr bool is_double                   = false;
    using accu_t                                      = DD_real;

    static constexpr auto gemm   = &cublasSgemm;
    static constexpr auto symm   = &cublasSsymm;
    static constexpr auto syrk   = &cublasSsyrk;
    static constexpr auto syr2k  = &cublasSsyr2k;
    static constexpr auto syrkx  = &cublasSsyrkx;
    static constexpr auto trmm   = &cublasStrmm;
    static constexpr auto trsm   = &cublasStrsm;
    static constexpr auto trtrmm = &cublasSgemm;

    static constexpr float one() { return 1.0f; }
    static constexpr float mone() { return -1.0f; }
    static constexpr float zero() { return 0.0f; }
};

template <> struct testTraits<double> {
    static constexpr unsigned NUM_MODULI_MIN          = 12;
    static constexpr unsigned NUM_MODULI_MAX          = 16;
    static constexpr unsigned NUM_MODULI_MIN_accuracy = 9;
    static constexpr unsigned NUM_MODULI_MAX_accuracy = 20;
    static constexpr char prefix                      = 'd';
    static constexpr char prefix_upper                = 'D';
    static constexpr bool is_complex                  = false;
    static constexpr bool is_double                   = true;
    using accu_t                                      = DD_real;

    static constexpr auto gemm   = &cublasDgemm;
    static constexpr auto symm   = &cublasDsymm;
    static constexpr auto syrk   = &cublasDsyrk;
    static constexpr auto syr2k  = &cublasDsyr2k;
    static constexpr auto syrkx  = &cublasDsyrkx;
    static constexpr auto trmm   = &cublasDtrmm;
    static constexpr auto trsm   = &cublasDtrsm;
    static constexpr auto trtrmm = &cublasDgemm;

    static constexpr double one() { return 1.0; }
    static constexpr double mone() { return -1.0; }
    static constexpr double zero() { return 0.0; }
};

template <> struct testTraits<cuFloatComplex> {
    static constexpr unsigned NUM_MODULI_MIN          = 6;
    static constexpr unsigned NUM_MODULI_MAX          = 9;
    static constexpr unsigned NUM_MODULI_MIN_accuracy = 3;
    static constexpr unsigned NUM_MODULI_MAX_accuracy = 12;
    static constexpr char prefix                      = 'c';
    static constexpr char prefix_upper                = 'C';
    static constexpr bool is_complex                  = true;
    static constexpr bool is_double                   = false;
    using accu_t                                      = DD_complex;

    static constexpr auto gemm   = &cublasCgemm;
    static constexpr auto symm   = &cublasCsymm;
    static constexpr auto syrk   = &cublasCsyrk;
    static constexpr auto syr2k  = &cublasCsyr2k;
    static constexpr auto syrkx  = &cublasCsyrkx;
    static constexpr auto trmm   = &cublasCtrmm;
    static constexpr auto trsm   = &cublasCtrsm;
    static constexpr auto hemm   = &cublasChemm;
    static constexpr auto herk   = &cublasCherk;
    static constexpr auto her2k  = &cublasCher2k;
    static constexpr auto herkx  = &cublasCherkx;
    static constexpr auto trtrmm = &cublasCgemm;

    static constexpr cuFloatComplex one() { return cuFloatComplex{1.0f, 0.0f}; }
    static constexpr cuFloatComplex mone() { return cuFloatComplex{-1.0f, 0.0f}; }
    static constexpr cuFloatComplex zero() { return cuFloatComplex{0.0f, 0.0f}; }
};

template <> struct testTraits<cuDoubleComplex> {
    static constexpr unsigned NUM_MODULI_MIN          = 12;
    static constexpr unsigned NUM_MODULI_MAX          = 16;
    static constexpr unsigned NUM_MODULI_MIN_accuracy = 9;
    static constexpr unsigned NUM_MODULI_MAX_accuracy = 20;
    static constexpr char prefix                      = 'z';
    static constexpr char prefix_upper                = 'Z';
    static constexpr bool is_complex                  = true;
    static constexpr bool is_double                   = true;
    using accu_t                                      = DD_complex;

    static constexpr auto gemm   = &cublasZgemm;
    static constexpr auto symm   = &cublasZsymm;
    static constexpr auto syrk   = &cublasZsyrk;
    static constexpr auto syr2k  = &cublasZsyr2k;
    static constexpr auto syrkx  = &cublasZsyrkx;
    static constexpr auto trmm   = &cublasZtrmm;
    static constexpr auto trsm   = &cublasZtrsm;
    static constexpr auto hemm   = &cublasZhemm;
    static constexpr auto herk   = &cublasZherk;
    static constexpr auto her2k  = &cublasZher2k;
    static constexpr auto herkx  = &cublasZherkx;
    static constexpr auto trtrmm = &cublasZgemm;
#if defined(__CUDACC__)
    static constexpr auto gemm3m = &cublasZgemm3m;
#endif

    static constexpr cuDoubleComplex one() { return cuDoubleComplex{1.0, 0.0}; }
    static constexpr cuDoubleComplex mone() { return cuDoubleComplex{-1.0, 0.0}; }
    static constexpr cuDoubleComplex zero() { return cuDoubleComplex{0.0, 0.0}; }
};

template <gemmul8::Func>
inline constexpr bool dependent_false_func_v = false;

template <gemmul8::Func func, typename T>
constexpr auto getCublasFunc() {
    if constexpr (func == gemmul8::Func::gemm) {
        return testTraits<T>::gemm;
    } else if constexpr (func == gemmul8::Func::symm) {
        return testTraits<T>::symm;
    } else if constexpr (func == gemmul8::Func::syrk) {
        return testTraits<T>::syrk;
    } else if constexpr (func == gemmul8::Func::syr2k) {
        return testTraits<T>::syr2k;
    } else if constexpr (func == gemmul8::Func::syrkx) {
        return testTraits<T>::syrkx;
    } else if constexpr (func == gemmul8::Func::trmm) {
        return testTraits<T>::trmm;
    } else if constexpr (func == gemmul8::Func::trsm) {
        return testTraits<T>::trsm;
    } else if constexpr (func == gemmul8::Func::hemm) {
        return testTraits<T>::hemm;
    } else if constexpr (func == gemmul8::Func::herk) {
        return testTraits<T>::herk;
    } else if constexpr (func == gemmul8::Func::her2k) {
        return testTraits<T>::her2k;
    } else if constexpr (func == gemmul8::Func::herkx) {
        return testTraits<T>::herkx;
    } else if constexpr (func == gemmul8::Func::trtrmm) {
        return testTraits<T>::trtrmm;
    } else {
        static_assert(dependent_false_func_v<func>, "Unsupported gemmul8::Func.");
    }
}

template <gemmul8::Func func, typename T>
constexpr auto getDDFunc() {
    if constexpr (func == gemmul8::Func::gemm) {
        return &eval::DDgemm;
    } else if constexpr (func == gemmul8::Func::symm) {
        return &eval::DDsymm;
    } else if constexpr (func == gemmul8::Func::syrk) {
        return &eval::DDsyrk;
    } else if constexpr (func == gemmul8::Func::syr2k) {
        return &eval::DDsyr2k;
    } else if constexpr (func == gemmul8::Func::syrkx) {
        return &eval::DDsyrkx;
    } else if constexpr (func == gemmul8::Func::trmm) {
        return &eval::DDtrmm;
    } else if constexpr (func == gemmul8::Func::trsm) {
        return &eval::DDtrsm;
    } else if constexpr (func == gemmul8::Func::hemm) {
        return &eval::DDhemm;
    } else if constexpr (func == gemmul8::Func::herk) {
        return &eval::DDherk;
    } else if constexpr (func == gemmul8::Func::her2k) {
        return &eval::DDher2k;
    } else if constexpr (func == gemmul8::Func::herkx) {
        return &eval::DDherkx;
    } else if constexpr (func == gemmul8::Func::trtrmm) {
        return &eval::DDtrtrmm;
    } else {
        static_assert(dependent_false_func_v<func>, "Unsupported gemmul8::Func.");
    }
}

template <gemmul8::Func func>
constexpr std::string_view getFuncName() {
    if constexpr (func == gemmul8::Func::gemm) {
        return "gemm";
    } else if constexpr (func == gemmul8::Func::symm) {
        return "symm";
    } else if constexpr (func == gemmul8::Func::syrk) {
        return "syrk";
    } else if constexpr (func == gemmul8::Func::syr2k) {
        return "syr2k";
    } else if constexpr (func == gemmul8::Func::syrkx) {
        return "syrkx";
    } else if constexpr (func == gemmul8::Func::trmm) {
        return "trmm";
    } else if constexpr (func == gemmul8::Func::trsm) {
        return "trsm";
    } else if constexpr (func == gemmul8::Func::hemm) {
        return "hemm";
    } else if constexpr (func == gemmul8::Func::herk) {
        return "herk";
    } else if constexpr (func == gemmul8::Func::her2k) {
        return "her2k";
    } else if constexpr (func == gemmul8::Func::herkx) {
        return "herkx";
    } else if constexpr (func == gemmul8::Func::trtrmm) {
        return "trtrmm";
    } else {
        static_assert(dependent_false_func_v<func>, "Unsupported gemmul8::Func.");
    }
}

template <gemmul8::Func func>
constexpr std::string_view getFuncName_upper() {
    if constexpr (func == gemmul8::Func::gemm) {
        return "GEMM";
    } else if constexpr (func == gemmul8::Func::symm) {
        return "SYMM";
    } else if constexpr (func == gemmul8::Func::syrk) {
        return "SYRK";
    } else if constexpr (func == gemmul8::Func::syr2k) {
        return "SYR2K";
    } else if constexpr (func == gemmul8::Func::syrkx) {
        return "SYRKX";
    } else if constexpr (func == gemmul8::Func::trmm) {
        return "TRMM";
    } else if constexpr (func == gemmul8::Func::trsm) {
        return "TRSM";
    } else if constexpr (func == gemmul8::Func::hemm) {
        return "HEMM";
    } else if constexpr (func == gemmul8::Func::herk) {
        return "HERK";
    } else if constexpr (func == gemmul8::Func::her2k) {
        return "HER2K";
    } else if constexpr (func == gemmul8::Func::herkx) {
        return "HERKX";
    } else if constexpr (func == gemmul8::Func::trtrmm) {
        return "TRTRMM (GEMM)";
    } else {
        static_assert(dependent_false_func_v<func>, "Unsupported gemmul8::Func.");
    }
}

template <gemmul8::Func func>
constexpr double getFuncCost(
    const size_t m,
    const size_t n,
    const size_t k,
    const bool Complex //
) {
    const double dm   = static_cast<double>(m);
    const double dn   = static_cast<double>(n);
    const double dk   = static_cast<double>(k);
    const double cmul = (Complex) ? 4.0 : 1.0;

    if constexpr (func == gemmul8::Func::gemm) {
        // GEMM: C = A B
        return 2.0 * dm * dn * dk * cmul;

    } else if constexpr (func == gemmul8::Func::symm) {
        // SYMM: C = A B or C = B A
        // k should be m for left-side SYMM, and n for right-side SYMM.
        return 2.0 * dm * dn * dk * cmul;

    } else if constexpr (func == gemmul8::Func::syrk) {
        // SYRK: triangular part of C = A A^T
        // For square C, call with m == n.
        return dm * dn * dk * cmul;

    } else if constexpr (func == gemmul8::Func::syr2k) {
        // SYR2K: triangular part of C = A B^T + B A^T
        // For square C, call with m == n.
        return 2.0 * dm * dn * dk * cmul;

    } else if constexpr (func == gemmul8::Func::syrkx) {
        // SYRKX: triangular part of C = A B^T
        // For square C, call with m == n.
        return dm * dn * dk * cmul;

    } else if constexpr (func == gemmul8::Func::trmm) {
        // TRMM: triangular matrix times dense matrix.
        // k should be the triangular dimension.
        return dm * dn * dk * cmul;

    } else if constexpr (func == gemmul8::Func::trsm) {
        // TRSM: triangular solve with dense RHS.
        // Divisions and lower-order terms are ignored.
        return dm * dn * dk * cmul;

    } else if constexpr (func == gemmul8::Func::hemm) {
        // HEMM: C = A B or C = B A, with Hermitian A.
        // k should be m for left-side HEMM, and n for right-side HEMM.
        return 2.0 * dm * dn * dk * cmul;

    } else if constexpr (func == gemmul8::Func::herk) {
        // HERK: triangular part of C = A A^H
        // For square C, call with m == n.
        return dm * dn * dk * cmul;

    } else if constexpr (func == gemmul8::Func::her2k) {
        // HER2K: triangular part of C = A B^H + B A^H
        // For square C, call with m == n.
        return 2.0 * dm * dn * dk * cmul;

    } else if constexpr (func == gemmul8::Func::herkx) {
        // HERKX: triangular part of C = A B^H
        // For square C, call with m == n.
        return dm * dn * dk * cmul;

    } else if constexpr (func == gemmul8::Func::trtrmm) {
        // Triangular-by-triangular multiplication.
        // For an n-by-n triangular product, call with m == n == k.
        // Leading term is approximately n^3 / 3.
        return (dm * dn * dk / 3.0) * cmul;

    } else {
        static_assert(dependent_false_func_v<func>, "Unsupported gemmul8::Func.");
    }
}

inline std::string opTag(cublasOperation_t op) {
    if (op == CUBLAS_OP_N) return "n_";
    if (op == CUBLAS_OP_T) return "t_";
    if (op == CUBLAS_OP_C) return "c_";
    return "unknown_";
}

inline std::string diagTag(cublasDiagType_t diag) {
    if (diag == CUBLAS_DIAG_NON_UNIT) return "nonunit_";
    if (diag == CUBLAS_DIAG_UNIT) return "unit_";
    return "unknown_";
}

template <typename T>
inline double calc_median(std::vector<T> &times) {
    std::sort(times.begin(), times.end());
    const size_t num = times.size();
    if (num % 2 == 1) {
        return double(times[num / 2]);
    } else {
        return (double(times[num / 2]) + double(times[num / 2 - 1])) * 0.5;
    }
}
