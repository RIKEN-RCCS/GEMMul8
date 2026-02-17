#pragma once
#if defined(__NVCC__)
    #include <cuda_runtime.h>
    #include <cublas_v2.h>
    #include <cublasLt.h>
    #include <cuComplex.h>
#endif
#if defined(__HIPCC__)
    #include <hip/hip_runtime.h>
    #include <hipblas/hipblas.h>
    #include <hipblaslt/hipblaslt.h>
    #include <hip/hip_complex.h>
#endif
#include <cstddef>
#include <vector>

namespace gemmul8 {

enum class Backend { INT8,
                     FP8 };

/***
 * workSize returns the required workspace size in bytes.
 */
template <bool is_Complex = false, Backend backend = Backend::INT8>
size_t workSize(
    size_t m,                         // Number of rows of C
    size_t n,                         // Number of columns of C
    size_t k,                         // Inner dimension <= 2^17
    unsigned num_moduli,              // #moduli, 2 <= num_moduli <= 20
    bool enable_skip_scalA = false,   // [optional] Reserve extra space for A to allow skip_scalA
    bool enable_skip_scalB = false,   // [optional] Reserve extra space for B to allow skip_scalB
    size_t *workSizeA      = nullptr, // [optional] Output: workspace size used for A8i and sftA
    size_t *workSizeB      = nullptr  // [optional] Output: workspace size used for B8i and sftB
);

/***
 * GEMM emulation using INT8 Tensor Cores
 */
#if defined(__NVCC__)
template <typename T, Backend backend = Backend::INT8>
std::vector<double> gemm(
    cublasHandle_t handle,            // Handle to the cuBLAS library context
    cublasOperation_t op_A,           // CUBLAS_OP_N, CUBLAS_OP_T, or CUBLAS_OP_C
    cublasOperation_t op_B,           // CUBLAS_OP_N, CUBLAS_OP_T, or CUBLAS_OP_C
    size_t m,                         // Number of rows of C
    size_t n,                         // Number of columns of C
    size_t k,                         // Inner dimension <= 2^17
    const T *alpha,                   // Scaling factor for op(A)*op(B)
    const T *const A,                 // 1-D device array of dimensions lda*k (CUBLAS_OP_N) or lda*m (CUBLAS_OP_T)
    size_t lda,                       // Leading dimension of A
    const T *const B,                 // 1-D device array of dimensions ldb*n (CUBLAS_OP_N) or ldb*k (CUBLAS_OP_T)
    size_t ldb,                       // Leading dimension of B
    const T *beta,                    // Scaling factor for C
    T *const C,                       // 1-D device array of dimensions ldc*n
    size_t ldc,                       // Leading dimension of C
    unsigned num_moduli,              // #moduli, 2 <= num_moduli <= 20
    bool fastmode,                    // false (accurate mode) or true (fast mode)
    void *const work,                 // Preallocated workspace
    void *const workA      = nullptr, // [optional] Separate workspace for A (if nullptr, uses work)
    void *const workB      = nullptr, // [optional] Separate workspace for B (if nullptr, uses work)
    bool enable_skip_scalA = false,   // [optional] Enables scaling-skip mechanism for A
    bool enable_skip_scalB = false,   // [optional] Enables scaling-skip mechanism for B
    bool skip_scalA        = false,   // [optional] If true, skip preprocessing for A
    bool skip_scalB        = false    // [optional] If true, skip preprocessing for B
);

template <typename T, Backend backend = Backend::INT8>
std::vector<double> gemm(
    cublasLtHandle_t handle,          // Handle to the cuBLASLt library context
    cublasOperation_t op_A,           // CUBLAS_OP_N, CUBLAS_OP_T, or CUBLAS_OP_C
    cublasOperation_t op_B,           // CUBLAS_OP_N, CUBLAS_OP_T, or CUBLAS_OP_C
    size_t m,                         // Number of rows of C
    size_t n,                         // Number of columns of C
    size_t k,                         // Inner dimension <= 2^17
    const T *alpha,                   // Scaling factor for op(A)*op(B)
    const T *const A,                 // 1-D device array of dimensions lda*k (CUBLAS_OP_N) or lda*m (CUBLAS_OP_T)
    size_t lda,                       // Leading dimension of A
    const T *const B,                 // 1-D device array of dimensions ldb*n (CUBLAS_OP_N) or ldb*k (CUBLAS_OP_T)
    size_t ldb,                       // Leading dimension of B
    const T *beta,                    // Scaling factor for C
    T *const C,                       // 1-D device array of dimensions ldc*n
    size_t ldc,                       // Leading dimension of C
    unsigned num_moduli,              // #moduli, 2 <= num_moduli <= 20
    bool fastmode,                    // false (accurate mode) or true (fast mode)
    void *const work,                 // Preallocated workspace
    void *const workA      = nullptr, // [optional] Separate workspace for A (if nullptr, uses work)
    void *const workB      = nullptr, // [optional] Separate workspace for B (if nullptr, uses work)
    bool enable_skip_scalA = false,   // [optional] Enables scaling-skip mechanism for A
    bool enable_skip_scalB = false,   // [optional] Enables scaling-skip mechanism for B
    bool skip_scalA        = false,   // [optional] If true, skip preprocessing for A
    bool skip_scalB        = false,   // [optional] If true, skip preprocessing for B
    cudaStream_t stream    = 0        // [optional] stream identifier
);
#endif

#if defined(__HIPCC__)
template <typename T, Backend backend = Backend::INT8>
std::vector<double> gemm(
    hipblasHandle_t handle,           // Handle to the hipBLAS library context
    hipblasOperation_t op_A,          // HIPBLAS_OP_N, HIPBLAS_OP_T, or HIPBLAS_OP_C
    hipblasOperation_t op_B,          // HIPBLAS_OP_N, HIPBLAS_OP_T, or HIPBLAS_OP_C
    size_t m,                         // Number of rows of C
    size_t n,                         // Number of columns of C
    size_t k,                         // Inner dimension <= 2^17
    const T *alpha,                   // Scaling factor for op(A)*op(B)
    const T *const A,                 // 1-D device array of dimensions lda*k (HIPBLAS_OP_N) or lda*m (HIPBLAS_OP_T)
    size_t lda,                       // Leading dimension of A
    const T *const B,                 // 1-D device array of dimensions ldb*n (HIPBLAS_OP_N) or ldb*k (HIPBLAS_OP_T)
    size_t ldb,                       // Leading dimension of B
    const T *beta,                    // Scaling factor for C
    T *const C,                       // 1-D device array of dimensions ldc*n
    size_t ldc,                       // Leading dimension of C
    unsigned num_moduli,              // #moduli, 2 <= num_moduli <= 20
    bool fastmode,                    // false (accurate mode) or true (fast mode)
    void *const work,                 // Preallocated workspace
    void *const workA      = nullptr, // [optional] Separate workspace for A (if nullptr, uses work)
    void *const workB      = nullptr, // [optional] Separate workspace for B (if nullptr, uses work)
    bool enable_skip_scalA = false,   // [optional] Enables scaling-skip mechanism for A
    bool enable_skip_scalB = false,   // [optional] Enables scaling-skip mechanism for B
    bool skip_scalA        = false,   // [optional] If true, skip preprocessing for A
    bool skip_scalB        = false    // [optional] If true, skip preprocessing for B
);

template <typename T, Backend backend = Backend::INT8>
std::vector<double> gemm(
    hipblasLtHandle_t handle,         // Handle to the hipBLAS library context
    hipblasOperation_t op_A,          // HIPBLAS_OP_N, HIPBLAS_OP_T, or HIPBLAS_OP_C
    hipblasOperation_t op_B,          // HIPBLAS_OP_N, HIPBLAS_OP_T, or HIPBLAS_OP_C
    size_t m,                         // Number of rows of C
    size_t n,                         // Number of columns of C
    size_t k,                         // Inner dimension <= 2^17
    const T *alpha,                   // Scaling factor for op(A)*op(B)
    const T *const A,                 // 1-D device array of dimensions lda*k (HIPBLAS_OP_N) or lda*m (HIPBLAS_OP_T)
    size_t lda,                       // Leading dimension of A
    const T *const B,                 // 1-D device array of dimensions ldb*n (HIPBLAS_OP_N) or ldb*k (HIPBLAS_OP_T)
    size_t ldb,                       // Leading dimension of B
    const T *beta,                    // Scaling factor for C
    T *const C,                       // 1-D device array of dimensions ldc*n
    size_t ldc,                       // Leading dimension of C
    unsigned num_moduli,              // #moduli, 2 <= num_moduli <= 20
    bool fastmode,                    // false (accurate mode) or true (fast mode)
    void *const work,                 // Preallocated workspace
    void *const workA      = nullptr, // [optional] Separate workspace for A (if nullptr, uses work)
    void *const workB      = nullptr, // [optional] Separate workspace for B (if nullptr, uses work)
    bool enable_skip_scalA = false,   // [optional] Enables scaling-skip mechanism for A
    bool enable_skip_scalB = false,   // [optional] Enables scaling-skip mechanism for B
    bool skip_scalA        = false,   // [optional] If true, skip preprocessing for A
    bool skip_scalB        = false,   // [optional] If true, skip preprocessing for B
    cudaStream_t stream    = 0        // [optional] stream identifier
);
#endif

} // namespace gemmul8
