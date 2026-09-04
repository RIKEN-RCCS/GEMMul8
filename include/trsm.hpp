/**
 * TRSM
 * ----
 * Solves a triangular system with multiple right-hand sides:
 *
 *   If side == CUBLAS_SIDE_LEFT / HIPBLAS_SIDE_LEFT:
 *     op(A) * X = alpha * B.
 *
 *   If side == CUBLAS_SIDE_RIGHT / HIPBLAS_SIDE_RIGHT:
 *     X * op(A) = alpha * B.
 *
 * A is triangular.  The triangle specified by uplo is referenced.
 * If diag == CUBLAS_DIAG_UNIT / HIPBLAS_DIAG_UNIT, the diagonal entries
 * of A are treated as one.
 *
 * This interface follows the cublasXtrsm / hipblasXtrsm argument convention,
 * with additional GEMMul8-specific arguments.
 *
 * Matrix dimensions:
 *
 *   B is m-by-n.
 *   X is m-by-n.
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
 *     Element type of B and alpha.  Defaults to TA.
 *
 * Template-argument order:
 *
 *   Both template-argument orders are supported:
 *
 *     gemmul8::trsm<TA, gemmul8::Backend::INT8, TB>(...);
 *     gemmul8::trsm<gemmul8::Backend::INT8, TA, TB>(...);
 *
 *   If TA and TB can be deduced from the function arguments, the
 *   alternative order can be shortened to:
 *
 *     gemmul8::trsm<gemmul8::Backend::INT8>(...);
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
 *     Workspace pointer.  If work == nullptr, the TRSM operation is not
 *     executed and the required workspace sizes are returned.
 *
 * Return value:
 *
 *   If work != nullptr:
 *     Executes TRSM and returns elapsed times (in seconds) of the four internal phases:
 *       t[0]: time for cuBLAS/hipBLAS trsm
 *       t[1]: time for gemm
 *
 *   If work == nullptr:
 *     Does not execute TRSM and returns required workspace sizes in bytes:
 *       t[0]: total workspace size
 *
 * Lt variant:
 *
 *   trsmLt() uses a cublasLtHandle_t / hipblasLtHandle_t and appends a stream
 *   argument at the end.
 */
#pragma once
#include "types.hpp"
#include "config.hpp"

namespace gemmul8 {

//------------------------------
// CUDA
//------------------------------
#if defined(__CUDACC__)

template <typename TA, Backend BACKEND = Backend::INT8, typename TB = TA>
std::vector<double> trsm(
    cublasHandle_t handle,
    cublasSideMode_t side, cublasFillMode_t uplo,
    cublasOperation_t trans, cublasDiagType_t diag,
    size_t m, size_t n,
    const TB *alpha,
    const TA *const A, size_t lda,
    TB *const B, size_t ldb,
    int num_moduli, bool fastmode,
    void *const work);

template <typename TA, Backend BACKEND = Backend::INT8, typename TB = TA>
std::vector<double> trsmLt(
    cublasLtHandle_t handle,
    cublasSideMode_t side, cublasFillMode_t uplo,
    cublasOperation_t trans, cublasDiagType_t diag,
    size_t m, size_t n,
    const TB *alpha,
    const TA *const A, size_t lda,
    TB *const B, size_t ldb,
    int num_moduli, bool fastmode,
    void *const work,
    cudaStream_t stream = 0);

// Alternative template-argument order:
template <Backend BACKEND, typename TA, typename TB = TA>
inline std::vector<double> trsm(
    cublasHandle_t handle,
    cublasSideMode_t side, cublasFillMode_t uplo,
    cublasOperation_t trans, cublasDiagType_t diag,
    size_t m, size_t n,
    const TB *alpha,
    const TA *const A, size_t lda,
    TB *const B, size_t ldb,
    int num_moduli, bool fastmode,
    void *const work //
) {
    return trsm<TA, BACKEND, TB>(
        handle, side, uplo, trans, diag, m, n,
        alpha, A, lda, B, ldb,
        num_moduli, fastmode, work);
}

// Alternative template-argument order:
template <Backend BACKEND, typename TA, typename TB = TA>
inline std::vector<double> trsmLt(
    cublasLtHandle_t handle,
    cublasSideMode_t side, cublasFillMode_t uplo,
    cublasOperation_t trans, cublasDiagType_t diag,
    size_t m, size_t n,
    const TB *alpha,
    const TA *const A, size_t lda,
    TB *const B, size_t ldb,
    int num_moduli, bool fastmode,
    void *const work,
    cudaStream_t stream = 0 //
) {
    return trsmLt<TA, BACKEND, TB>(
        handle, side, uplo, trans, diag, m, n,
        alpha, A, lda, B, ldb,
        num_moduli, fastmode, work, stream);
}

#endif

//------------------------------
// HIP
//------------------------------
#if defined(__HIPCC__)

template <typename TA, Backend BACKEND = Backend::INT8, typename TB = TA>
std::vector<double> trsm(
    hipblasHandle_t handle,
    hipblasSideMode_t side, hipblasFillMode_t uplo,
    hipblasOperation_t trans, hipblasDiagType_t diag,
    size_t m, size_t n,
    const TB *alpha,
    const TA *const A, size_t lda,
    TB *const B, size_t ldb,
    int num_moduli, bool fastmode,
    void *const work);

template <typename TA, Backend BACKEND = Backend::INT8, typename TB = TA>
std::vector<double> trsmLt(
    hipblasLtHandle_t handle,
    hipblasSideMode_t side, hipblasFillMode_t uplo,
    hipblasOperation_t trans, hipblasDiagType_t diag,
    size_t m, size_t n,
    const TB *alpha,
    const TA *const A, size_t lda,
    TB *const B, size_t ldb,
    int num_moduli, bool fastmode,
    void *const work,
    hipStream_t stream = 0);

// Alternative template-argument order:
template <Backend BACKEND, typename TA, typename TB = TA>
inline std::vector<double> trsm(
    hipblasHandle_t handle,
    hipblasSideMode_t side, hipblasFillMode_t uplo,
    hipblasOperation_t trans, hipblasDiagType_t diag,
    size_t m, size_t n,
    const TB *alpha,
    const TA *const A, size_t lda,
    TB *const B, size_t ldb,
    int num_moduli, bool fastmode,
    void *const work //
) {
    return trsm<TA, BACKEND, TB>(
        handle, side, uplo, trans, diag, m, n,
        alpha, A, lda, B, ldb,
        num_moduli, fastmode, work);
}

// Alternative template-argument order:
template <Backend BACKEND, typename TA, typename TB = TA>
inline std::vector<double> trsmLt(
    hipblasLtHandle_t handle,
    hipblasSideMode_t side, hipblasFillMode_t uplo,
    hipblasOperation_t trans, hipblasDiagType_t diag,
    size_t m, size_t n,
    const TB *alpha,
    const TA *const A, size_t lda,
    TB *const B, size_t ldb,
    int num_moduli, bool fastmode,
    void *const work,
    hipStream_t stream = 0 //
) {
    return trsmLt<TA, BACKEND, TB>(
        handle, side, uplo, trans, diag, m, n,
        alpha, A, lda, B, ldb,
        num_moduli, fastmode, work, stream);
}

#endif

} // namespace gemmul8
