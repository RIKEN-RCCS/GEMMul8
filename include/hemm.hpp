/**
 * HEMM
 * ----
 * Computes a Hermitian matrix-matrix product:
 *
 *   If side == CUBLAS_SIDE_LEFT / HIPBLAS_SIDE_LEFT:
 *     C := alpha * A * B + beta * C.
 *
 *   If side == CUBLAS_SIDE_RIGHT / HIPBLAS_SIDE_RIGHT:
 *     C := alpha * B * A + beta * C.
 *
 * A is Hermitian.  The triangle specified by uplo is referenced.
 * The diagonal entries of A are assumed to be real.
 *
 * This interface follows the cublasXhemm / hipblasXhemm argument convention,
 * with additional GEMMul8-specific arguments.
 *
 * Matrix dimensions:
 *
 *   B is m-by-n.
 *   C is m-by-n.
 *
 *   If side == CUBLAS_SIDE_LEFT / HIPBLAS_SIDE_LEFT:
 *     A is m-by-m.
 *
 *   If side == CUBLAS_SIDE_RIGHT / HIPBLAS_SIDE_RIGHT:
 *     A is n-by-n.
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
 *   TB:
 *     Element type of B.  Defaults to TA.
 *
 *   TC:
 *     Element type of C, alpha, and beta.  Defaults to TA.
 *
 * Template-argument order:
 *
 *   Both template-argument orders are supported:
 *
 *     gemmul8::hemm<TA, gemmul8::Backend::INT8, TB, TC>(...);
 *     gemmul8::hemm<gemmul8::Backend::INT8, TA, TB, TC>(...);
 *
 *   If TA, TB, and TC can be deduced from the function arguments, the
 *   alternative order can be shortened to:
 *
 *     gemmul8::hemm<gemmul8::Backend::INT8>(...);
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
 *     Workspace pointer.  If work == nullptr, the HEMM operation is not
 *     executed and the required workspace sizes are returned.
 *
 *   workA, workB:
 *     Optional workspace pointers for A and B.  If nullptr, the corresponding
 *     regions are taken from work.
 *
 *   enable_skip_scalA, enable_skip_scalB:
 *     Reserve additional workspace that allows reuse of previously scaled
 *     matrices.
 *
 *   skip_scalA, skip_scalB:
 *     If enabled, skip scaling A and/or B and reuse the corresponding
 *     precomputed scaled matrices.
 *
 * Return value:
 *
 *   If work != nullptr:
 *     Executes HEMM and returns elapsed times (in seconds) of the four internal phases:
 *       t[0]: time for scaling and quantization
 *       t[1]: time for low-precision matrix multiplication
 *       t[2]: time for re-quantization of matrix products
 *       t[3]: time for final reduction of CRT and undo scaling
 *
 *   If work == nullptr:
 *     Does not execute HEMM and returns required workspace sizes in bytes:
 *       t[0]: total workspace size (including sizes of workA and workB)
 *       t[1]: workspace size associated with A (size of workA)
 *       t[2]: workspace size associated with B (size of workB)
 *
 * Note:
 * 
 *   When BACKEND = Backend::INT8 and fastmode = false, currently supports 
 *   the inner dimension of the underlying matrix multiplication of at most 2^17.
 *
 * Lt variant:
 *
 *   hemmLt() uses a cublasLtHandle_t / hipblasLtHandle_t and appends a stream
 *   argument at the end.
 */
#pragma once
#include "types.hpp"

namespace gemmul8 {

//------------------------------
// CUDA
//------------------------------
#if defined(__CUDACC__)

template <typename TA, Backend BACKEND = Backend::INT8, typename TB = TA, typename TC = TA>
std::vector<double> hemm(
    cublasHandle_t handle,
    cublasSideMode_t side, cublasFillMode_t uplo,
    size_t m, size_t n,
    const TC *alpha,
    const TA *const A, size_t lda,
    const TB *const B, size_t ldb,
    const TC *beta,
    TC *const C, size_t ldc,
    int num_moduli, bool fastmode,
    void *const work,
    void *const workA      = nullptr,
    void *const workB      = nullptr,
    bool enable_skip_scalA = false,
    bool enable_skip_scalB = false,
    bool skip_scalA        = false,
    bool skip_scalB        = false);

template <typename TA, Backend BACKEND = Backend::INT8, typename TB = TA, typename TC = TA>
std::vector<double> hemmLt(
    cublasLtHandle_t handle,
    cublasSideMode_t side, cublasFillMode_t uplo,
    size_t m, size_t n,
    const TC *alpha,
    const TA *const A, size_t lda,
    const TB *const B, size_t ldb,
    const TC *beta,
    TC *const C, size_t ldc,
    int num_moduli, bool fastmode,
    void *const work,
    void *const workA      = nullptr,
    void *const workB      = nullptr,
    bool enable_skip_scalA = false,
    bool enable_skip_scalB = false,
    bool skip_scalA        = false,
    bool skip_scalB        = false,
    cudaStream_t stream    = 0);

// Alternative template-argument order:
template <Backend BACKEND, typename TA, typename TB = TA, typename TC = TA>
inline std::vector<double> hemm(
    cublasHandle_t handle,
    cublasSideMode_t side, cublasFillMode_t uplo,
    size_t m, size_t n,
    const TC *alpha,
    const TA *const A, size_t lda,
    const TB *const B, size_t ldb,
    const TC *beta,
    TC *const C, size_t ldc,
    int num_moduli, bool fastmode,
    void *const work,
    void *const workA      = nullptr,
    void *const workB      = nullptr,
    bool enable_skip_scalA = false,
    bool enable_skip_scalB = false,
    bool skip_scalA        = false,
    bool skip_scalB        = false //
) {
    return hemm<TA, BACKEND, TB, TC>(
        handle, side, uplo, m, n,
        alpha, A, lda, B, ldb, beta, C, ldc,
        num_moduli, fastmode,
        work, workA, workB,
        enable_skip_scalA, enable_skip_scalB,
        skip_scalA, skip_scalB);
}

// Alternative template-argument order:
template <Backend BACKEND, typename TA, typename TB = TA, typename TC = TA>
inline std::vector<double> hemmLt(
    cublasLtHandle_t handle,
    cublasSideMode_t side, cublasFillMode_t uplo,
    size_t m, size_t n,
    const TC *alpha,
    const TA *const A, size_t lda,
    const TB *const B, size_t ldb,
    const TC *beta,
    TC *const C, size_t ldc,
    int num_moduli, bool fastmode,
    void *const work,
    void *const workA      = nullptr,
    void *const workB      = nullptr,
    bool enable_skip_scalA = false,
    bool enable_skip_scalB = false,
    bool skip_scalA        = false,
    bool skip_scalB        = false,
    cudaStream_t stream    = 0 //
) {
    return hemmLt<TA, BACKEND, TB, TC>(
        handle, side, uplo, m, n,
        alpha, A, lda, B, ldb, beta, C, ldc,
        num_moduli, fastmode,
        work, workA, workB,
        enable_skip_scalA, enable_skip_scalB,
        skip_scalA, skip_scalB, stream);
}

#endif

//------------------------------
// HIP
//------------------------------
#if defined(__HIPCC__)

template <typename TA, Backend BACKEND = Backend::INT8, typename TB = TA, typename TC = TA>
std::vector<double> hemm(
    hipblasHandle_t handle,
    hipblasSideMode_t side, hipblasFillMode_t uplo,
    size_t m, size_t n,
    const TC *alpha,
    const TA *const A, size_t lda,
    const TB *const B, size_t ldb,
    const TC *beta,
    TC *const C, size_t ldc,
    int num_moduli, bool fastmode,
    void *const work,
    void *const workA      = nullptr,
    void *const workB      = nullptr,
    bool enable_skip_scalA = false,
    bool enable_skip_scalB = false,
    bool skip_scalA        = false,
    bool skip_scalB        = false);

template <typename TA, Backend BACKEND = Backend::INT8, typename TB = TA, typename TC = TA>
std::vector<double> hemmLt(
    hipblasLtHandle_t handle,
    hipblasSideMode_t side, hipblasFillMode_t uplo,
    size_t m, size_t n,
    const TC *alpha,
    const TA *const A, size_t lda,
    const TB *const B, size_t ldb,
    const TC *beta,
    TC *const C, size_t ldc,
    int num_moduli, bool fastmode,
    void *const work,
    void *const workA      = nullptr,
    void *const workB      = nullptr,
    bool enable_skip_scalA = false,
    bool enable_skip_scalB = false,
    bool skip_scalA        = false,
    bool skip_scalB        = false,
    hipStream_t stream     = 0);

// Alternative template-argument order:
template <Backend BACKEND, typename TA, typename TB = TA, typename TC = TA>
inline std::vector<double> hemm(
    hipblasHandle_t handle,
    hipblasSideMode_t side, hipblasFillMode_t uplo,
    size_t m, size_t n,
    const TC *alpha,
    const TA *const A, size_t lda,
    const TB *const B, size_t ldb,
    const TC *beta,
    TC *const C, size_t ldc,
    int num_moduli, bool fastmode,
    void *const work,
    void *const workA      = nullptr,
    void *const workB      = nullptr,
    bool enable_skip_scalA = false,
    bool enable_skip_scalB = false,
    bool skip_scalA        = false,
    bool skip_scalB        = false //
) {
    return hemm<TA, BACKEND, TB, TC>(
        handle, side, uplo, m, n,
        alpha, A, lda, B, ldb, beta, C, ldc,
        num_moduli, fastmode,
        work, workA, workB,
        enable_skip_scalA, enable_skip_scalB,
        skip_scalA, skip_scalB);
}

// Alternative template-argument order:
template <Backend BACKEND, typename TA, typename TB = TA, typename TC = TA>
inline std::vector<double> hemmLt(
    hipblasLtHandle_t handle,
    hipblasSideMode_t side, hipblasFillMode_t uplo,
    size_t m, size_t n,
    const TC *alpha,
    const TA *const A, size_t lda,
    const TB *const B, size_t ldb,
    const TC *beta,
    TC *const C, size_t ldc,
    int num_moduli, bool fastmode,
    void *const work,
    void *const workA      = nullptr,
    void *const workB      = nullptr,
    bool enable_skip_scalA = false,
    bool enable_skip_scalB = false,
    bool skip_scalA        = false,
    bool skip_scalB        = false,
    hipStream_t stream     = 0 //
) {
    return hemmLt<TA, BACKEND, TB, TC>(
        handle, side, uplo, m, n,
        alpha, A, lda, B, ldb, beta, C, ldc,
        num_moduli, fastmode,
        work, workA, workB,
        enable_skip_scalA, enable_skip_scalB,
        skip_scalA, skip_scalB, stream);
}

#endif

} // namespace gemmul8
