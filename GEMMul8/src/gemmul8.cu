#include "../include/gemmul8.hpp"

#include <algorithm>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <random>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>
#ifndef _WIN32
    #include <dlfcn.h>
#endif

#if defined(__NVCC__)
    #include <cuda_runtime.h>
    #include <cublas_v2.h>
    #include <cublasLt.h>
    #include <cuComplex.h>
    #include <cuda_fp8.h>
#endif

#include "self_hipify.hpp"
#include "template_type.hpp"
#include "common.hpp"
#include "template_math.hpp"
#include "table.hpp"
#include "find_max.hpp"
#include "mod.hpp"
#include "matmult.hpp"
#include "scaling.hpp"
#include "scaling_fast_real.hpp"
#include "scaling_fast_complex.hpp"
#include "scaling_accu_real.hpp"
#include "scaling_accu_complex.hpp"
#include "conv_hi2mid_real.hpp"
#include "conv_hi2mid_complex.hpp"
#include "inverse_scaling_real.hpp"
#include "inverse_scaling_complex.hpp"
#include "gemmul8_real.hpp"
#include "gemmul8_complex.hpp"

#if !defined(GEMM_ARGS)
    #define GEMM_ARGS(T) cublasHandle_t handle,                                  \
                         cublasOperation_t op_A, cublasOperation_t op_B,         \
                         size_t m, size_t n, size_t k,                           \
                         const T *alpha, const T *const A, size_t lda,           \
                         const T *const B, size_t ldb,                           \
                         const T *beta, T *const C, size_t ldc,                  \
                         unsigned num_moduli, bool fastmode,                     \
                         void *const work, void *const workA, void *const workB, \
                         bool enable_skip_scalA, bool enable_skip_scalB,         \
                         bool skip_scalA, bool skip_scalB
#endif

#if !defined(GEMMLt_ARGS)
    #define GEMMLt_ARGS(T) cublasLtHandle_t handle,                                \
                           cublasOperation_t op_A, cublasOperation_t op_B,         \
                           size_t m, size_t n, size_t k,                           \
                           const T *alpha, const T *const A, size_t lda,           \
                           const T *const B, size_t ldb,                           \
                           const T *beta, T *const C, size_t ldc,                  \
                           unsigned num_moduli, bool fastmode,                     \
                           void *const work, void *const workA, void *const workB, \
                           bool enable_skip_scalA, bool enable_skip_scalB,         \
                           bool skip_scalA, bool skip_scalB,                       \
                           cudaStream_t stream
#endif

#if !defined(GEMM_CALL_ARGS)
    #define GEMM_CALL_ARGS Handle_t(handle), op_A, op_B, m, n, k,    \
                           alpha, A, lda, B, ldb, beta, C, ldc,      \
                           num_moduli, fastmode, work, workA, workB, \
                           enable_skip_scalA, enable_skip_scalB,     \
                           skip_scalA, skip_scalB, stream
#endif

namespace gemmul8 {

//------------------------------
// Calculate required work size (INT8)
//------------------------------
template <> size_t workSize<false, Backend::INT8>(size_t m, size_t n, size_t k, unsigned num_moduli, bool enable_skip_scalA, bool enable_skip_scalB, size_t *workSizeA, size_t *workSizeB) {
    return real::workSize<Backend::INT8>(m, n, k, num_moduli, enable_skip_scalA, enable_skip_scalB, workSizeA, workSizeB);
}
template <> size_t workSize<true, Backend::INT8>(size_t m, size_t n, size_t k, unsigned num_moduli, bool enable_skip_scalA, bool enable_skip_scalB, size_t *workSizeA, size_t *workSizeB) {
    return complex::workSize<Backend::INT8>(m, n, k, num_moduli, enable_skip_scalA, enable_skip_scalB, workSizeA, workSizeB);
}

//------------------------------
// Calculate required work size (FP8)
//------------------------------
template <> size_t workSize<true, Backend::FP8>(size_t m, size_t n, size_t k, unsigned num_moduli, bool enable_skip_scalA, bool enable_skip_scalB, size_t *workSizeA, size_t *workSizeB) {
    return complex::workSize<Backend::FP8>(m, n, k, num_moduli, enable_skip_scalA, enable_skip_scalB, workSizeA, workSizeB);
}
template <> size_t workSize<false, Backend::FP8>(size_t m, size_t n, size_t k, unsigned num_moduli, bool enable_skip_scalA, bool enable_skip_scalB, size_t *workSizeA, size_t *workSizeB) {
    return real::workSize<Backend::FP8>(m, n, k, num_moduli, enable_skip_scalA, enable_skip_scalB, workSizeA, workSizeB);
}

//------------------------------
// GEMM emulation using INT8 Tensor Cores (cuBLAS)
//------------------------------
template <> std::vector<double> gemm<double, Backend::INT8>(GEMM_ARGS(double)) {
    cudaStream_t stream;
    cublasGetStream(handle, &stream);
    return real::gemm<double, Backend::INT8>(GEMM_CALL_ARGS);
}
template <> std::vector<double> gemm<float, Backend::INT8>(GEMM_ARGS(float)) {
    cudaStream_t stream;
    cublasGetStream(handle, &stream);
    return real::gemm<float, Backend::INT8>(GEMM_CALL_ARGS);
}
template <> std::vector<double> gemm<cuFloatComplex, Backend::INT8>(GEMM_ARGS(cuFloatComplex)) {
    cudaStream_t stream;
    cublasGetStream(handle, &stream);
    return complex::gemm<cuFloatComplex, Backend::INT8>(GEMM_CALL_ARGS);
}
template <> std::vector<double> gemm<cuDoubleComplex, Backend::INT8>(GEMM_ARGS(cuDoubleComplex)) {
    cudaStream_t stream;
    cublasGetStream(handle, &stream);
    return complex::gemm<cuDoubleComplex, Backend::INT8>(GEMM_CALL_ARGS);
}

//------------------------------
// GEMM emulation using FP8 Tensor Cores (cuBLAS)
//------------------------------
// NOT SUPPORTED

//------------------------------
// GEMM emulation using INT8 Tensor Cores (cuBLASLt)
//------------------------------
template <> std::vector<double> gemm<double, Backend::INT8>(GEMMLt_ARGS(double)) { return real::gemm<double, Backend::INT8>(GEMM_CALL_ARGS); }
template <> std::vector<double> gemm<float, Backend::INT8>(GEMMLt_ARGS(float)) { return real::gemm<float, Backend::INT8>(GEMM_CALL_ARGS); }
template <> std::vector<double> gemm<cuFloatComplex, Backend::INT8>(GEMMLt_ARGS(cuFloatComplex)) { return complex::gemm<cuFloatComplex, Backend::INT8>(GEMM_CALL_ARGS); }
template <> std::vector<double> gemm<cuDoubleComplex, Backend::INT8>(GEMMLt_ARGS(cuDoubleComplex)) { return complex::gemm<cuDoubleComplex, Backend::INT8>(GEMM_CALL_ARGS); }

//------------------------------
// GEMM emulation using FP8 Tensor Cores (cuBLASLt)
//------------------------------
template <> std::vector<double> gemm<double, Backend::FP8>(GEMMLt_ARGS(double)) { return real::gemm<double, Backend::FP8>(GEMM_CALL_ARGS); }
template <> std::vector<double> gemm<float, Backend::FP8>(GEMMLt_ARGS(float)) { return real::gemm<float, Backend::FP8>(GEMM_CALL_ARGS); }
template <> std::vector<double> gemm<cuFloatComplex, Backend::FP8>(GEMMLt_ARGS(cuFloatComplex)) { return complex::gemm<cuFloatComplex, Backend::FP8>(GEMM_CALL_ARGS); }
template <> std::vector<double> gemm<cuDoubleComplex, Backend::FP8>(GEMMLt_ARGS(cuDoubleComplex)) { return complex::gemm<cuDoubleComplex, Backend::FP8>(GEMM_CALL_ARGS); }

} // namespace gemmul8
