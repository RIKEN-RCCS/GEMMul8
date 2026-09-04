/**
 * Workspace size query
 * --------------------
 * Returns the required workspace size in bytes for the specified GEMMul8
 * operation.
 *
 * This function provides a lightweight workspace query interface without
 * calling the corresponding BLAS-like routine.  The returned value is the total
 * required workspace size and corresponds to
 *
 *   t[0]
 *
 * in the return value of each routine when that routine is called with
 * work == nullptr.
 *
 * Template parameters:
 *
 *   is_Complex:
 *     If true, query the workspace size for complex-valued routines.
 *     If false, query the workspace size for real-valued routines.
 *
 *   BACKEND:
 *     Low-precision backend used internally.  Defaults to Backend::INT8.
 *
 *   FUNC:
 *     Operation for which the workspace size is queried.  Defaults to
 *     Func::gemm.
 *
 *   T:
 *     Element type of A & B for trsm.
 *
 * Template-argument order:
 *
 *   Both template-argument orders are supported:
 *
 *     gemmul8::workSize<is_Complex, gemmul8::Backend::INT8, FUNC>(...);
 *     gemmul8::workSize<gemmul8::Backend::INT8, is_Complex, FUNC>(...);
 *     gemmul8::workSizeTrsm<T, gemmul8::Backend::INT8>(handle, ...);
 *     gemmul8::workSizeTrsm<gemmul8::Backend::INT8, T>(handle, ...);
 *     gemmul8::workSizeTrsmLt<T, gemmul8::Backend::INT8>(ltHandle, ...);
 *     gemmul8::workSizeTrsmLt<gemmul8::Backend::INT8, T>(ltHandle, ...);
 *
 *   If is_Complex and FUNC use their defaults, the alternative order can be
 *   shortened to:
 *
 *     gemmul8::workSize<gemmul8::Backend::INT8>(...);
 *     gemmul8::workSizeTrsm<gemmul8::Backend::INT8>(handle, ...);
 *     gemmul8::workSizeTrsmLt<gemmul8::Backend::INT8>(ltHandle, ...);
 *
 * Arguments:
 *
 *   m, n, k:
 *     Compact matrix-size parameters used for workspace estimation.
 *
 *     m and n are always the dimensions of the output matrix:
 *
 *       output is m-by-n.
 *
 *     k is the inner dimension for matrix-matrix products, or the effective
 *     order of the structured input matrix when the operation does not have a
 *     separate inner dimension.
 *
 *   num_moduli:
 *     Number of moduli used in the Ozaki-scheme decomposition.
 *     Valid range:
 *       2 <= num_moduli <= 20 for FP64 output,
 *       2 <= num_moduli <= 13 for FP32 output.
 *
 *   enable_skip_scalA, enable_skip_scalB:
 *     If true, include the additional workspace required to allow reuse of
 *     previously scaled matrices A and/or B.
 *
 *   workSizeA:
 *     Optional output pointer.  If not nullptr, the workspace size associated
 *     with A is written to *workSizeA.  This value corresponds to t[1] in the
 *     return value of each routine when called with work == nullptr.
 *
 *   workSizeB:
 *     Optional output pointer.  If not nullptr, the workspace size associated
 *     with B is written to *workSizeB.  This value corresponds to t[2] in the
 *     return value of each routine when called with work == nullptr.
 *
 *   fastmode:
 *     If true, use the fast scaling path.
 *     If false, use the more accurate scaling path.
 *
 *   Examples:
 *
 *     Func::gemm:
 *       op(A) is m-by-k,
 *       op(B) is k-by-n,
 *       C     is m-by-n.
 *
 *     Func::symm, Func::hemm:
 *       C is m-by-n.
 *       If side == LEFT,  use workSize(m, n, m, ...).
 *       If side == RIGHT, use workSize(m, n, n, ...).
 *
 *     Func::syrk, Func::herk:
 *       C is n-by-n in the BLAS routine.
 *       Use workSize(n, n, k, ...).
 *
 *     Func::syr2k, Func::syrkx, Func::her2k, Func::herkx:
 *       C is n-by-n in the BLAS routine.
 *       Use workSize(n, n, k, ...).
 *
 *     Func::trmm:
 *       C is m-by-n.
 *       If side == LEFT,  use workSize(m, n, m, ...).
 *       If side == RIGHT, use workSize(m, n, n, ...).
 *
 *     Func::trtrmm:
 *       A, B, and C are n-by-n.
 *       Use workSize(n, n, n, ...).
 *
 *     Func::trsm:
 *       Use workSizeTrsm(handle, side, m, n, ...) for a BLAS handle, or
 *       workSizeTrsmLt(handle, side, m, n, ...) for an Lt handle.
 *       The required workspace depends on side and the handle-local TRSM block size.
 *
 * Return value:
 *
 *   Total required workspace size in bytes.
 *
 * Notes:
 *
 *   This function does not allocate memory.
 *
 *   For routines with only one input matrix requiring scaling, such as syrk()
 *   and herk(), the workspace size associated with B may be zero or unused.
 *
 *   The side, uplo, trans, and diag arguments of the corresponding BLAS-like
 *   routines are not passed to workSize().  When they affect the effective
 *   size of the structured operand, choose k according to the examples above.
 */
#pragma once
#include "types.hpp"
#include "config.hpp"

namespace gemmul8 {

// for gemm, symm, syr2k, syrkx, syrk, hemm, her2k, herkx, herk, trmm, trtrmm
template <bool is_Complex = false, Backend BACKEND = Backend::INT8, Func FUNC = Func::gemm>
size_t workSize(
    size_t m, size_t n, size_t k,
    int num_moduli,
    bool enable_skip_scalA = false,   // [optional] Input: Reserve extra space for A to allow skip_scalA
    bool enable_skip_scalB = false,   // [optional] Input: Reserve extra space for B to allow skip_scalB
    size_t *workSizeA      = nullptr, // [optional] Output: workspace size used for A8i and sftA
    size_t *workSizeB      = nullptr, // [optional] Output: workspace size used for B8i and sftB
    bool fastmode          = false    // [optioanl] Input: Additional workspace is required for accurate mode (fastmode=false) when enable_skip_scal* = true except syrk & herk
);

// Alternative template-argument order:
template <Backend BACKEND = Backend::INT8, bool is_Complex = false, Func FUNC = Func::gemm>
inline size_t workSize(
    size_t m, size_t n, size_t k,
    int num_moduli,
    bool enable_skip_scalA = false,
    bool enable_skip_scalB = false,
    size_t *workSizeA      = nullptr,
    size_t *workSizeB      = nullptr,
    bool fastmode          = false //
) {
    return workSize<is_Complex, BACKEND, FUNC>(
        m, n, k,
        num_moduli,
        enable_skip_scalA, enable_skip_scalB,
        workSizeA, workSizeB, fastmode);
}

// for trsm
//------------------------------
// CUDA
//------------------------------
#if defined(__CUDACC__)

template <typename T = double, Backend BACKEND = Backend::INT8>
size_t workSizeTrsm(
    cublasHandle_t handle,
    cublasSideMode_t side,
    size_t m, size_t n,
    int num_moduli);

template <typename T = double, Backend BACKEND = Backend::INT8>
size_t workSizeTrsmLt(
    cublasLtHandle_t handle,
    cublasSideMode_t side,
    size_t m, size_t n,
    int num_moduli);

// Alternative template-argument order:
template <Backend BACKEND = Backend::INT8, typename T = double>
inline size_t workSizeTrsm(
    cublasHandle_t handle,
    cublasSideMode_t side,
    size_t m, size_t n,
    int num_moduli //
) {
    return workSizeTrsm<T, BACKEND>(handle, side, m, n, num_moduli);
}

template <Backend BACKEND = Backend::INT8, typename T = double>
inline size_t workSizeTrsmLt(
    cublasLtHandle_t handle,
    cublasSideMode_t side,
    size_t m, size_t n,
    int num_moduli //
) {
    return workSizeTrsmLt<T, BACKEND>(handle, side, m, n, num_moduli);
}

#endif

//------------------------------
// HIP
//------------------------------
#if defined(__HIPCC__)

template <typename T = double, Backend BACKEND = Backend::INT8>
size_t workSizeTrsm(
    hipblasHandle_t handle,
    hipblasSideMode_t side,
    size_t m, size_t n,
    int num_moduli);

template <typename T = double, Backend BACKEND = Backend::INT8>
size_t workSizeTrsmLt(
    hipblasLtHandle_t handle,
    hipblasSideMode_t side,
    size_t m, size_t n,
    int num_moduli);

// Alternative template-argument order:
template <Backend BACKEND = Backend::INT8, typename T = double>
inline size_t workSizeTrsm(
    hipblasHandle_t handle,
    hipblasSideMode_t side,
    size_t m, size_t n,
    int num_moduli //
) {
    return workSizeTrsm<T, BACKEND>(handle, side, m, n, num_moduli);
}

template <Backend BACKEND = Backend::INT8, typename T = double>
inline size_t workSizeTrsmLt(
    hipblasLtHandle_t handle,
    hipblasSideMode_t side,
    size_t m, size_t n,
    int num_moduli //
) {
    return workSizeTrsmLt<T, BACKEND>(handle, side, m, n, num_moduli);
}

#endif

} // namespace gemmul8
