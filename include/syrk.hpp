/**
 * SYRK
 * ----
 * Computes a symmetric rank-k update:
 *
 *   If trans == CUBLAS_OP_N / HIPBLAS_OP_N:
 *     C := alpha * A * A^T + beta * C.
 *
 *   Otherwise:
 *     C := alpha * A^T * A + beta * C.
 *
 * C is symmetric.  Only the triangle specified by uplo is updated.
 *
 * This interface follows the cublasXsyrk / hipblasXsyrk argument convention,
 * with additional GEMMul8-specific arguments.
 *
 * Matrix dimensions:
 *
 *   C is n-by-n.
 *
 *   If trans == CUBLAS_OP_N / HIPBLAS_OP_N:
 *     A is n-by-k.
 *
 *   Otherwise:
 *     A is k-by-n.
 *
 * Template parameters:
 *
 *   TA:
 *     Element type of A.
 *
 *   BACKEND:
 *     Low-precision backend used internally.
 *     Select one of the following:
 *       - INT8-based emulation: gemmul8::Backend::INT8
 *       - FP8-based emulation : gemmul8::Backend::FP8
 *     Defaults to gemmul8::Backend::INT8.
 *
 *   TC:
 *     Element type of C, alpha, and beta.  Defaults to TA.
 *
 * Template-argument order:
 *
 *   Both template-argument orders are supported:
 *
 *     gemmul8::syrk<TA, gemmul8::Backend::INT8, TC>(...);
 *     gemmul8::syrk<gemmul8::Backend::INT8, TA, TC>(...);
 *
 *   If TA and TC can be deduced from the function arguments, the
 *   alternative order can be shortened to:
 *
 *     gemmul8::syrk<gemmul8::Backend::INT8>(...);
 *
 * GEMMul8-specific arguments:
 *
 *   num_moduli:
 *     Number of moduli used in the Ozaki-scheme decomposition.
 *     Valid range:
 *       2 <= num_moduli <= 20 for FP64 output,
 *       2 <= num_moduli <= 13 for FP32 output.
 *
 *   fastmode:
 *     If true, use the fast scaling path.
 *     If false, use the more accurate scaling path.
 *
 *   work:
 *     Workspace pointer.  If work == nullptr, the SYRK operation is not
 *     executed and the required workspace sizes are returned.
 *
 *   workA:
 *     Optional workspace pointer for A.  If nullptr, the corresponding region
 *     is taken from work.
 *
 *   enable_skip_scalA:
 *     Reserve additional workspace that allows reuse of the previously scaled
 *     matrix A.
 *
 *   skip_scalA:
 *     If enabled, skip scaling A and reuse the corresponding precomputed
 *     scaled matrix.
 *
 * Return value:
 *
 *   If work != nullptr:
 *     Executes SYRK and returns elapsed times (in seconds) of the four internal phases:
 *       t[0]: time for scaling and quantization
 *       t[1]: time for low-precision matrix multiplication
 *       t[2]: time for re-quantization of matrix products
 *       t[3]: time for final reduction of CRT and undo scaling
 *
 *   If work == nullptr:
 *     Does not execute SYRK and returns required workspace sizes in bytes:
 *       t[0]: total workspace size (including size of workA)
 *       t[1]: workspace size associated with A (size of workA)
 *       t[2]: zero or unused
 *
 * Note:
 * 
 *   When BACKEND = Backend::INT8 and fastmode = false, currently supports 
 *   the inner dimension of the underlying matrix multiplication of at most 2^17.
 *
 * Lt variant:
 *
 *   syrkLt() uses a cublasLtHandle_t / hipblasLtHandle_t and appends a stream
 *   argument at the end.
 */
#pragma once
#include "types.hpp"

namespace gemmul8 {

//------------------------------
// CUDA
//------------------------------
#if defined(__CUDACC__)

template <typename TA, Backend BACKEND = Backend::INT8, typename TC = TA>
std::vector<double> syrk(
    cublasHandle_t handle,
    cublasFillMode_t uplo, cublasOperation_t trans,
    size_t n, size_t k,
    const TC *alpha, const TA *const A, size_t lda,
    const TC *beta, TC *const C, size_t ldc,
    int num_moduli, bool fastmode,
    void *const work,
    void *const workA      = nullptr,
    bool enable_skip_scalA = false,
    bool skip_scalA        = false);

template <typename TA, Backend BACKEND = Backend::INT8, typename TC = TA>
std::vector<double> syrkLt(
    cublasLtHandle_t handle,
    cublasFillMode_t uplo, cublasOperation_t trans,
    size_t n, size_t k,
    const TC *alpha, const TA *const A, size_t lda,
    const TC *beta, TC *const C, size_t ldc,
    int num_moduli, bool fastmode,
    void *const work,
    void *const workA      = nullptr,
    bool enable_skip_scalA = false,
    bool skip_scalA        = false,
    cudaStream_t stream    = 0);

// Alternative template-argument order:
template <Backend BACKEND, typename TA, typename TC = TA>
inline std::vector<double> syrk(
    cublasHandle_t handle,
    cublasFillMode_t uplo, cublasOperation_t trans,
    size_t n, size_t k,
    const TC *alpha, const TA *const A, size_t lda,
    const TC *beta, TC *const C, size_t ldc,
    int num_moduli, bool fastmode,
    void *const work,
    void *const workA      = nullptr,
    bool enable_skip_scalA = false,
    bool skip_scalA        = false //
) {
    return syrk<TA, BACKEND, TC>(
        handle, uplo, trans, n, k,
        alpha, A, lda, beta, C, ldc,
        num_moduli, fastmode,
        work, workA,
        enable_skip_scalA, skip_scalA);
}

// Alternative template-argument order:
template <Backend BACKEND, typename TA, typename TC = TA>
inline std::vector<double> syrkLt(
    cublasLtHandle_t handle,
    cublasFillMode_t uplo, cublasOperation_t trans,
    size_t n, size_t k,
    const TC *alpha, const TA *const A, size_t lda,
    const TC *beta, TC *const C, size_t ldc,
    int num_moduli, bool fastmode,
    void *const work,
    void *const workA      = nullptr,
    bool enable_skip_scalA = false,
    bool skip_scalA        = false,
    cudaStream_t stream    = 0 //
) {
    return syrkLt<TA, BACKEND, TC>(
        handle, uplo, trans, n, k,
        alpha, A, lda, beta, C, ldc,
        num_moduli, fastmode,
        work, workA,
        enable_skip_scalA, skip_scalA, stream);
}

#endif

//------------------------------
// HIP
//------------------------------
#if defined(__HIPCC__)

template <typename TA, Backend BACKEND = Backend::INT8, typename TC = TA>
std::vector<double> syrk(
    hipblasHandle_t handle,
    hipblasFillMode_t uplo, hipblasOperation_t trans,
    size_t n, size_t k,
    const TC *alpha, const TA *const A, size_t lda,
    const TC *beta, TC *const C, size_t ldc,
    int num_moduli, bool fastmode,
    void *const work,
    void *const workA      = nullptr,
    bool enable_skip_scalA = false,
    bool skip_scalA        = false);

template <typename TA, Backend BACKEND = Backend::INT8, typename TC = TA>
std::vector<double> syrkLt(
    hipblasLtHandle_t handle,
    hipblasFillMode_t uplo, hipblasOperation_t trans,
    size_t n, size_t k,
    const TC *alpha, const TA *const A, size_t lda,
    const TC *beta, TC *const C, size_t ldc,
    int num_moduli, bool fastmode,
    void *const work,
    void *const workA      = nullptr,
    bool enable_skip_scalA = false,
    bool skip_scalA        = false,
    hipStream_t stream     = 0);

// Alternative template-argument order:
template <Backend BACKEND, typename TA, typename TC = TA>
inline std::vector<double> syrk(
    hipblasHandle_t handle,
    hipblasFillMode_t uplo, hipblasOperation_t trans,
    size_t n, size_t k,
    const TC *alpha, const TA *const A, size_t lda,
    const TC *beta, TC *const C, size_t ldc,
    int num_moduli, bool fastmode,
    void *const work,
    void *const workA      = nullptr,
    bool enable_skip_scalA = false,
    bool skip_scalA        = false //
) {
    return syrk<TA, BACKEND, TC>(
        handle, uplo, trans, n, k,
        alpha, A, lda, beta, C, ldc,
        num_moduli, fastmode,
        work, workA,
        enable_skip_scalA, skip_scalA);
}

// Alternative template-argument order:
template <Backend BACKEND, typename TA, typename TC = TA>
inline std::vector<double> syrkLt(
    hipblasLtHandle_t handle,
    hipblasFillMode_t uplo, hipblasOperation_t trans,
    size_t n, size_t k,
    const TC *alpha, const TA *const A, size_t lda,
    const TC *beta, TC *const C, size_t ldc,
    int num_moduli, bool fastmode,
    void *const work,
    void *const workA      = nullptr,
    bool enable_skip_scalA = false,
    bool skip_scalA        = false,
    hipStream_t stream     = 0 //
) {
    return syrkLt<TA, BACKEND, TC>(
        handle, uplo, trans, n, k,
        alpha, A, lda, beta, C, ldc,
        num_moduli, fastmode,
        work, workA,
        enable_skip_scalA, skip_scalA, stream);
}

#endif

} // namespace gemmul8
