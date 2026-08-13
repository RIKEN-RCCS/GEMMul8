#pragma once
#include "common.hpp"
#include "../mod/mod_f32.hpp"

namespace gemmul8::common {

template <Backend BACKEND>
inline void call_gemm_tn_raw(
    const cudaStream_t stream, Handle_t &h,
    int m, int n, int k,
    const void *alpha,
    const void *A, size_t lda,
    const void *B, size_t ldb,
    const void *beta,
    void *C, size_t ldc //
) {
    constexpr auto CUDA_R_LOW  = (BACKEND == Backend::INT8) ? CUDA_R_8I : CUDA_R_8F_E4M3;
    constexpr auto CUDA_R_HIGH = (BACKEND == Backend::INT8) ? CUDA_R_32I : CUDA_R_32F;
    constexpr auto COMP_TYPE   = (BACKEND == Backend::INT8) ? CUBLAS_COMPUTE_32I : CUBLAS_COMPUTE_32F;

    if (h.kind == HandleKind::cuBLAS) {
        cublasGemmEx(h.cublas, CUBLAS_OP_T, CUBLAS_OP_N, m, n, k,
                     alpha, A, CUDA_R_LOW, lda, B, CUDA_R_LOW, ldb,
                     beta, C, CUDA_R_HIGH, ldc,
                     COMP_TYPE, CUBLAS_GEMM_DEFAULT);
        return;
    }

    update_lt_layouts(h, m, n, k, lda, ldb, ldc, 1, 0, 0, 0);
    auto &plan = get_or_create_plan<BACKEND>(
        h, m, n, k, 1, lda, ldb, ldc, 0, 0, 0, LtBatchKind::None);

    const size_t ws = std::min(h.workspaceSizeInBytes, plan.ws);

    cublasLtMatmul(h.cublasLt, h.opDesc,
                   alpha, A, h.Adesc, B, h.Bdesc,
                   beta, C, h.Cdesc, C, h.Cdesc,
                   &plan.heur.algo, h.workspace, ws, stream);
}

template <Backend BACKEND>
inline void call_gemm_tn(
    const cudaStream_t stream, Handle_t &h,
    int m, int n, int k,
    const void *alpha,
    const void *A, size_t lda,
    const void *B, size_t ldb,
    const void *beta,
    void *C, size_t ldc //
) {
    if constexpr (BACKEND != Backend::FP8) {
        call_gemm_tn_raw<BACKEND>(
            stream, h,
            m, n, k,
            alpha, A, lda,
            B, ldb,
            beta, C, ldc);
        return;
    }

    if (!h.fp8_k_blocking || k <= h.fp8_k_block_first) {
        call_gemm_tn_raw<BACKEND>(
            stream, h,
            m, n, k,
            alpha, A, lda,
            B, ldb,
            beta, C, ldc);
        return;
    }

    using LowT = low_t<Backend::FP8>;
    using HiT  = hi_t<Backend::FP8>;

    constexpr HiT one = 1.0f;

    const LowT *A0 = reinterpret_cast<const LowT *>(A);
    const LowT *B0 = reinterpret_cast<const LowT *>(B);
    HiT *C0        = reinterpret_cast<HiT *>(C);

    const int KB0 = h.fp8_k_block_first;
    const int KB  = h.fp8_k_block_next;

    int k0 = 0;
    {
        const int kk = std::min(KB0, k);

        call_gemm_tn_raw<Backend::FP8>(
            stream, h,
            m, n, kk,
            alpha,
            A0, lda,
            B0, ldb,
            beta,
            C0, ldc);

        k0 = kk;
    }

    while (k0 < k) {
        mod_f32(stream, C0, m, n, ldc, h);

        const int kk = std::min(KB, k - k0);

        call_gemm_tn_raw<Backend::FP8>(
            stream, h,
            m, n, kk,
            alpha,
            A0 + k0, lda,
            B0 + k0, ldb,
            &one,
            C0, ldc);

        k0 += kk;
    }
}

template <Backend BACKEND>
inline void call_gemm_tn_strided_batched_raw(
    const cudaStream_t stream, Handle_t &h,
    int m, int n, int k, int batchCount,
    const void *alpha,
    const void *A, size_t lda, int64_t strideA,
    const void *B, size_t ldb, int64_t strideB,
    const void *beta,
    void *C, size_t ldc, int64_t strideC //
) {
    constexpr auto CUDA_R_LOW  = (BACKEND == Backend::INT8) ? CUDA_R_8I : CUDA_R_8F_E4M3;
    constexpr auto CUDA_R_HIGH = (BACKEND == Backend::INT8) ? CUDA_R_32I : CUDA_R_32F;
    constexpr auto COMP_TYPE   = (BACKEND == Backend::INT8) ? CUBLAS_COMPUTE_32I : CUBLAS_COMPUTE_32F;

    if (h.kind == HandleKind::cuBLAS) {
        cublasGemmStridedBatchedEx(h.cublas, CUBLAS_OP_T, CUBLAS_OP_N, m, n, k,
                                   alpha, A, CUDA_R_LOW, lda, strideA, B, CUDA_R_LOW, ldb, strideB,
                                   beta, C, CUDA_R_HIGH, ldc, strideC,
                                   batchCount, COMP_TYPE, CUBLAS_GEMM_DEFAULT);
        return;
    }

    update_lt_layouts(h, m, n, k, lda, ldb, ldc, batchCount, strideA, strideB, strideC);
    auto &plan = get_or_create_plan<BACKEND>(
        h, m, n, k, batchCount, lda, ldb, ldc, strideA, strideB, strideC, LtBatchKind::Strided);

    const size_t ws = std::min(h.workspaceSizeInBytes, plan.ws);

    cublasLtMatmul(h.cublasLt, h.opDesc,
                   alpha, A, h.Adesc, B, h.Bdesc,
                   beta, C, h.Cdesc, C, h.Cdesc,
                   &plan.heur.algo, h.workspace, ws, stream);
}

template <Backend BACKEND>
inline void call_gemm_tn_strided_batched(
    const cudaStream_t stream, Handle_t &h,
    int m, int n, int k, int batchCount,
    const void *alpha,
    const void *A, size_t lda, int64_t strideA,
    const void *B, size_t ldb, int64_t strideB,
    const void *beta,
    void *C, size_t ldc, int64_t strideC //
) {
    if constexpr (BACKEND != Backend::FP8) {
        call_gemm_tn_strided_batched_raw<BACKEND>(
            stream, h,
            m, n, k, batchCount,
            alpha,
            A, lda, strideA,
            B, ldb, strideB,
            beta,
            C, ldc, strideC);
        return;
    }

    if (!h.fp8_k_blocking || k <= h.fp8_k_block_first) {
        call_gemm_tn_strided_batched_raw<BACKEND>(
            stream, h,
            m, n, k, batchCount,
            alpha,
            A, lda, strideA,
            B, ldb, strideB,
            beta,
            C, ldc, strideC);
        return;
    }

    using LowT = low_t<Backend::FP8>;
    using HiT  = hi_t<Backend::FP8>;

    constexpr HiT one = 1.0f;

    const LowT *A0 = reinterpret_cast<const LowT *>(A);
    const LowT *B0 = reinterpret_cast<const LowT *>(B);
    HiT *C0        = reinterpret_cast<HiT *>(C);

    const int KB0 = h.fp8_k_block_first;
    const int KB  = h.fp8_k_block_next;

    int k0 = 0;
    {
        const int kk = std::min(KB0, k);

        call_gemm_tn_strided_batched_raw<Backend::FP8>(
            stream, h,
            m, n, kk, batchCount,
            alpha,
            A0 + k0, lda, strideA,
            B0 + k0, ldb, strideB,
            beta,
            C0, ldc, strideC);

        k0 = kk;
    }

    while (k0 < k) {
        mod_f32_strided(stream, C0, m, n, ldc, strideC, batchCount, h);

        const int kk = std::min(KB, k - k0);

        call_gemm_tn_strided_batched_raw<Backend::FP8>(
            stream, h,
            m, n, kk, batchCount,
            alpha,
            A0 + k0, lda, strideA,
            B0 + k0, ldb, strideB,
            &one,
            C0, ldc, strideC);

        k0 += kk;
    }
}

template <Backend BACKEND>
inline void call_gemm_tn_pointer_batched_raw(
    const cudaStream_t stream, Handle_t &h,
    int m, int n, int k, int batchCount,
    const void *alpha,
    void **Aarray, size_t lda,
    void **Barray, size_t ldb,
    const void *beta,
    void **Carray, size_t ldc //
) {
#if defined(__CUDACC__)
    constexpr auto CUDA_R_LOW  = (BACKEND == Backend::INT8) ? CUDA_R_8I : CUDA_R_8F_E4M3;
    constexpr auto CUDA_R_HIGH = (BACKEND == Backend::INT8) ? CUDA_R_32I : CUDA_R_32F;
    constexpr auto COMP_TYPE   = (BACKEND == Backend::INT8) ? CUBLAS_COMPUTE_32I : CUBLAS_COMPUTE_32F;

    if (h.kind == HandleKind::cuBLAS) {
        cublasGemmBatchedEx(h.cublas, CUBLAS_OP_T, CUBLAS_OP_N,
                            m, n, k,
                            alpha,
                            reinterpret_cast<const void *const *>(Aarray), CUDA_R_LOW, lda,
                            reinterpret_cast<const void *const *>(Barray), CUDA_R_LOW, ldb,
                            beta,
                            reinterpret_cast<void *const *>(Carray), CUDA_R_HIGH, ldc,
                            batchCount,
                            COMP_TYPE,
                            CUBLAS_GEMM_DEFAULT);
        return;
    }

    update_lt_layouts_pointer_array(h, m, n, k, lda, ldb, ldc, batchCount);
    auto &plan = get_or_create_plan<BACKEND>(
        h, m, n, k, batchCount, lda, ldb, ldc, 0, 0, 0, LtBatchKind::PointerArray);

    const size_t ws = std::min(h.workspaceSizeInBytes, plan.ws);

    cublasLtMatmul(h.cublasLt, h.opDesc,
                   alpha, Aarray, h.Adesc, Barray, h.Bdesc,
                   beta, Carray, h.Cdesc, Carray, h.Cdesc,
                   &plan.heur.algo, h.workspace, ws, stream);
#endif
}

template <Backend BACKEND>
inline void call_gemm_tn_pointer_batched(
    const cudaStream_t stream, Handle_t &h,
    int m, int n, int k, int batchCount,
    const void *alpha,
    void **Aarray, size_t lda,
    void **Barray, size_t ldb,
    const void *beta,
    void **Carray, size_t ldc //
) {
    if constexpr (BACKEND != Backend::FP8) {
        call_gemm_tn_pointer_batched_raw<BACKEND>(
            stream, h,
            m, n, k, batchCount,
            alpha,
            Aarray, lda,
            Barray, ldb,
            beta,
            Carray, ldc);
        return;
    }

    if (!h.fp8_k_blocking || k <= h.fp8_k_block_first) {
        call_gemm_tn_pointer_batched_raw<BACKEND>(
            stream, h,
            m, n, k, batchCount,
            alpha,
            Aarray, lda,
            Barray, ldb,
            beta,
            Carray, ldc);
        return;
    }

    using HiT = hi_t<Backend::FP8>;

    constexpr HiT one = 1.0f;

    const int KB0 = h.fp8_k_block_first;
    const int KB  = h.fp8_k_block_next;

    int k0 = 0;
    int kk = std::min(KB0, k);
    {
        call_gemm_tn_pointer_batched_raw<Backend::FP8>(
            stream, h,
            m, n, kk, batchCount,
            alpha,
            Aarray, lda,
            Barray, ldb,
            beta,
            Carray, ldc);

        k0 = kk;
    }

    while (k0 < k) {
        mod_f32_pointer_and_advance(stream, Aarray, Barray, Carray, m, n, ldc, batchCount, kk, h);

        kk = std::min(KB, k - k0);

        call_gemm_tn_pointer_batched_raw<Backend::FP8>(
            stream, h,
            m, n, kk, batchCount,
            alpha,
            Aarray, lda,
            Barray, ldb,
            &one,
            Carray, ldc);

        k0 += kk;
    }
}

//------------------------------
// ATxB
//------------------------------
template <Backend BACKEND>
inline void block_gemm_1(
    const cudaStream_t stream, Handle_t &handle, int m, int n, int k,
    const hi_t<BACKEND> *alpha, const low_t<BACKEND> *A, size_t lda, const low_t<BACKEND> *B, size_t ldb,
    const hi_t<BACKEND> *beta, hi_t<BACKEND> *C, size_t ldc //
) {
    const int nB      = handle.nB;
    const int nBlocks = (n + nB - 1) / nB;
    for (int ib = 0; ib < nBlocks; ++ib) {
        const int s   = ib * nB;
        const int rem = n - s;
        const int nn  = std::min(nB, rem);
        call_gemm_tn<BACKEND>(stream, handle, m, nn, k, alpha, A, lda, B + s * ldb, ldb, beta, C + s * ldc, ldc);
    }
}

template <Backend BACKEND>
inline void block_gemm_3(
    const cudaStream_t stream, Handle_t &handle, int m, int n, int k,
    const hi_t<BACKEND> *alpha1, const hi_t<BACKEND> *alpha2, const hi_t<BACKEND> *alpha3,
    const low_t<BACKEND> *A1, const low_t<BACKEND> *A2, const low_t<BACKEND> *A3, size_t lda,
    const low_t<BACKEND> *B1, const low_t<BACKEND> *B2, const low_t<BACKEND> *B3, size_t ldb,
    const hi_t<BACKEND> *beta1, const hi_t<BACKEND> *beta2, const hi_t<BACKEND> *beta3,
    hi_t<BACKEND> *C1, hi_t<BACKEND> *C2, hi_t<BACKEND> *C3, size_t ldc //
) {
    const int nB      = handle.nB;
    const int nBlocks = (n + nB - 1) / nB;
    for (int ib = 0; ib < nBlocks; ++ib) {
        const int s   = ib * nB;
        const int rem = n - s;
        const int nn  = std::min(nB, rem);
        call_gemm_tn<BACKEND>(stream, handle, m, nn, k, alpha1, A1, lda, B1 + s * ldb, ldb, beta1, C1 + s * ldc, ldc);
        call_gemm_tn<BACKEND>(stream, handle, m, nn, k, alpha2, A2, lda, B2 + s * ldb, ldb, beta2, C2 + s * ldc, ldc);
        call_gemm_tn<BACKEND>(stream, handle, m, nn, k, alpha3, A3, lda, B3 + s * ldb, ldb, beta3, C3 + s * ldc, ldc);
    }
}

template <Backend BACKEND>
inline void block_gemm_1_strided_batched(
    const cudaStream_t stream, Handle_t &handle,
    int m, int n, int k, int batchCount,
    const hi_t<BACKEND> *alpha,
    const low_t<BACKEND> *A, size_t lda, int64_t strideA,
    const low_t<BACKEND> *B, size_t ldb, int64_t strideB,
    const hi_t<BACKEND> *beta,
    hi_t<BACKEND> *C, size_t ldc, int64_t strideC //
) {
    const int nB      = handle.nB;
    const int nBlocks = (n + nB - 1) / nB;

    for (int ib = 0; ib < nBlocks; ++ib) {
        const int s   = ib * nB;
        const int rem = n - s;
        const int nn  = std::min(nB, rem);

        call_gemm_tn_strided_batched<BACKEND>(
            stream, handle,
            m, nn, k, batchCount,
            alpha,
            A, lda, strideA,
            B + s * ldb, ldb, strideB,
            beta,
            C + s * ldc, ldc, strideC);
    }
}

template <Backend BACKEND>
inline void block_gemm_3_strided_batched(
    const cudaStream_t stream, Handle_t &handle,
    int m, int n, int k, int batchCount,
    const hi_t<BACKEND> *alpha1, const hi_t<BACKEND> *alpha2, const hi_t<BACKEND> *alpha3,
    const low_t<BACKEND> *A1, const low_t<BACKEND> *A2, const low_t<BACKEND> *A3,
    size_t lda, int64_t strideA,
    const low_t<BACKEND> *B1, const low_t<BACKEND> *B2, const low_t<BACKEND> *B3,
    size_t ldb, int64_t strideB,
    const hi_t<BACKEND> *beta1, const hi_t<BACKEND> *beta2, const hi_t<BACKEND> *beta3,
    hi_t<BACKEND> *C1, hi_t<BACKEND> *C2, hi_t<BACKEND> *C3,
    size_t ldc, int64_t strideC //
) {
    const int nB      = handle.nB;
    const int nBlocks = (n + nB - 1) / nB;

    for (int ib = 0; ib < nBlocks; ++ib) {
        const int s   = ib * nB;
        const int rem = n - s;
        const int nn  = std::min(nB, rem);

        call_gemm_tn_strided_batched<BACKEND>(
            stream, handle,
            m, nn, k, batchCount,
            alpha1,
            A1, lda, strideA,
            B1 + s * ldb, ldb, strideB,
            beta1,
            C1 + s * ldc, ldc, strideC);

        call_gemm_tn_strided_batched<BACKEND>(
            stream, handle,
            m, nn, k, batchCount,
            alpha2,
            A2, lda, strideA,
            B2 + s * ldb, ldb, strideB,
            beta2,
            C2 + s * ldc, ldc, strideC);

        call_gemm_tn_strided_batched<BACKEND>(
            stream, handle,
            m, nn, k, batchCount,
            alpha3,
            A3, lda, strideA,
            B3 + s * ldb, ldb, strideB,
            beta3,
            C3 + s * ldc, ldc, strideC);
    }
}

//------------------------------
// ATxB (only uplo part)
//------------------------------
template <Backend BACKEND, cublasFillMode_t UPLO>
inline void block_ATxB_1(
    const cudaStream_t stream, Handle_t &handle, int n, int k,
    const hi_t<BACKEND> *alpha,
    const low_t<BACKEND> *A, size_t lda,
    const low_t<BACKEND> *B, size_t ldb,
    const hi_t<BACKEND> *beta,
    hi_t<BACKEND> *C, size_t ldc //
) {
    static_assert(UPLO == CUBLAS_FILL_MODE_UPPER || UPLO == CUBLAS_FILL_MODE_LOWER,
                  "block_ATxB requires UPPER or LOWER output UPLO.");

    const int nB    = handle.nB;
    const int n_pad = int(common::padding(size_t(n)));

    if (n <= nB) {
        call_gemm_tn<BACKEND>(stream, handle, n_pad, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
        return;
    }

    const int nBlocks = (n + nB - 1) / nB;
    for (int ib = 0; ib < nBlocks; ++ib) {
        const int s   = ib * nB;
        const int rem = n - s;
        const int mm  = std::min(nB, rem);

        const bool last_block = (ib + 1 == nBlocks);
        const int m           = (UPLO == CUBLAS_FILL_MODE_UPPER) ? (last_block ? (n_pad - s) : mm) : (n_pad - s);
        const int ncol        = (UPLO == CUBLAS_FILL_MODE_UPPER) ? rem : mm;

        call_gemm_tn<BACKEND>(
            stream, handle, m, ncol, k,
            alpha,
            A + s * lda, lda,
            B + s * ldb, ldb,
            beta,
            C + s * ldc + s, ldc);
    }
}

template <Backend BACKEND, cublasFillMode_t UPLO>
inline void block_ATxB_3(
    const cudaStream_t stream, Handle_t &handle, int n, int k,
    const hi_t<BACKEND> *alpha1, const hi_t<BACKEND> *alpha2, const hi_t<BACKEND> *alpha3,
    const low_t<BACKEND> *A1, const low_t<BACKEND> *A2, const low_t<BACKEND> *A3, size_t lda,
    const low_t<BACKEND> *B1, const low_t<BACKEND> *B2, const low_t<BACKEND> *B3, size_t ldb,
    const hi_t<BACKEND> *beta1, const hi_t<BACKEND> *beta2, const hi_t<BACKEND> *beta3,
    hi_t<BACKEND> *C1, hi_t<BACKEND> *C2, hi_t<BACKEND> *C3, size_t ldc //
) {
    static_assert(UPLO == CUBLAS_FILL_MODE_UPPER || UPLO == CUBLAS_FILL_MODE_LOWER,
                  "block_ATxB requires UPPER or LOWER output UPLO.");

    const int nB    = handle.nB;
    const int n_pad = int(common::padding(size_t(n)));

    if (n <= nB) {
        call_gemm_tn<BACKEND>(stream, handle, n_pad, n, k, alpha1, A1, lda, B1, ldb, beta1, C1, ldc);
        call_gemm_tn<BACKEND>(stream, handle, n_pad, n, k, alpha2, A2, lda, B2, ldb, beta2, C2, ldc);
        call_gemm_tn<BACKEND>(stream, handle, n_pad, n, k, alpha3, A3, lda, B3, ldb, beta3, C3, ldc);
        return;
    }

    const int nBlocks = (n + nB - 1) / nB;
    for (int ib = 0; ib < nBlocks; ++ib) {
        const int s   = ib * nB;
        const int rem = n - s;
        const int mm  = std::min(nB, rem);

        const bool last_block = (ib + 1 == nBlocks);
        const int m           = (UPLO == CUBLAS_FILL_MODE_UPPER) ? (last_block ? (n_pad - s) : mm) : (n_pad - s);
        const int ncol        = (UPLO == CUBLAS_FILL_MODE_UPPER) ? rem : mm;

        call_gemm_tn<BACKEND>(stream, handle, m, ncol, k, alpha1, A1 + s * lda, lda, B1 + s * ldb, ldb, beta1, C1 + s * ldc + s, ldc);
        call_gemm_tn<BACKEND>(stream, handle, m, ncol, k, alpha2, A2 + s * lda, lda, B2 + s * ldb, ldb, beta2, C2 + s * ldc + s, ldc);
        call_gemm_tn<BACKEND>(stream, handle, m, ncol, k, alpha3, A3 + s * lda, lda, B3 + s * ldb, ldb, beta3, C3 + s * ldc + s, ldc);
    }
}

template <Backend BACKEND, cublasFillMode_t UPLO>
inline void block_ATxB_1_strided_batched(
    const cudaStream_t stream, Handle_t &handle, int n, int k, int batchCount,
    const hi_t<BACKEND> *alpha,
    const low_t<BACKEND> *A, size_t lda, int64_t strideA,
    const low_t<BACKEND> *B, size_t ldb, int64_t strideB,
    const hi_t<BACKEND> *beta,
    hi_t<BACKEND> *C, size_t ldc, int64_t strideC //
) {
    static_assert(UPLO == CUBLAS_FILL_MODE_UPPER || UPLO == CUBLAS_FILL_MODE_LOWER,
                  "block_ATxB requires UPPER or LOWER output UPLO.");
    if (batchCount <= 0) return;

    const int nB    = handle.nB;
    const int n_pad = int(common::padding(size_t(n)));

    if (n <= nB) {
        call_gemm_tn_strided_batched<BACKEND>(
            stream, handle, n_pad, n, k, batchCount,
            alpha,
            A, lda, strideA,
            B, ldb, strideB,
            beta,
            C, ldc, strideC);
        return;
    }

    const int nBlocks = (n + nB - 1) / nB;
    for (int ib = 0; ib < nBlocks; ++ib) {
        const int s   = ib * nB;
        const int rem = n - s;
        const int mm  = std::min(nB, rem);

        const bool last_block = (ib + 1 == nBlocks);
        const int m           = (UPLO == CUBLAS_FILL_MODE_UPPER) ? (last_block ? (n_pad - s) : mm) : (n_pad - s);
        const int ncol        = (UPLO == CUBLAS_FILL_MODE_UPPER) ? rem : mm;

        call_gemm_tn_strided_batched<BACKEND>(
            stream, handle, m, ncol, k, batchCount,
            alpha,
            A + s * lda, lda, strideA,
            B + s * ldb, ldb, strideB,
            beta,
            C + s * ldc + s, ldc, strideC);
    }
}

template <Backend BACKEND, cublasFillMode_t UPLO>
inline void block_ATxB_3_strided_batched(
    const cudaStream_t stream, Handle_t &handle, int n, int k, int batchCount,
    const hi_t<BACKEND> *alpha1, const hi_t<BACKEND> *alpha2, const hi_t<BACKEND> *alpha3,
    const low_t<BACKEND> *A1, const low_t<BACKEND> *A2, const low_t<BACKEND> *A3, size_t lda, int64_t strideA,
    const low_t<BACKEND> *B1, const low_t<BACKEND> *B2, const low_t<BACKEND> *B3, size_t ldb, int64_t strideB,
    const hi_t<BACKEND> *beta1, const hi_t<BACKEND> *beta2, const hi_t<BACKEND> *beta3,
    hi_t<BACKEND> *C1, hi_t<BACKEND> *C2, hi_t<BACKEND> *C3, size_t ldc, int64_t strideC //
) {
    static_assert(UPLO == CUBLAS_FILL_MODE_UPPER || UPLO == CUBLAS_FILL_MODE_LOWER,
                  "block_ATxB requires UPPER or LOWER output UPLO.");
    if (batchCount <= 0) return;

    const int nB    = handle.nB;
    const int n_pad = int(common::padding(size_t(n)));

    if (n <= nB) {
        call_gemm_tn_strided_batched<BACKEND>(stream, handle, n_pad, n, k, batchCount, alpha1, A1, lda, strideA, B1, ldb, strideB, beta1, C1, ldc, strideC);
        call_gemm_tn_strided_batched<BACKEND>(stream, handle, n_pad, n, k, batchCount, alpha2, A2, lda, strideA, B2, ldb, strideB, beta2, C2, ldc, strideC);
        call_gemm_tn_strided_batched<BACKEND>(stream, handle, n_pad, n, k, batchCount, alpha3, A3, lda, strideA, B3, ldb, strideB, beta3, C3, ldc, strideC);
        return;
    }

    const int nBlocks = (n + nB - 1) / nB;
    for (int ib = 0; ib < nBlocks; ++ib) {
        const int s   = ib * nB;
        const int rem = n - s;
        const int mm  = std::min(nB, rem);

        const bool last_block = (ib + 1 == nBlocks);
        const int m           = (UPLO == CUBLAS_FILL_MODE_UPPER) ? (last_block ? (n_pad - s) : mm) : (n_pad - s);
        const int ncol        = (UPLO == CUBLAS_FILL_MODE_UPPER) ? rem : mm;

        call_gemm_tn_strided_batched<BACKEND>(stream, handle, m, ncol, k, batchCount, alpha1, A1 + s * lda, lda, strideA, B1 + s * ldb, ldb, strideB, beta1, C1 + s * ldc + s, ldc, strideC);
        call_gemm_tn_strided_batched<BACKEND>(stream, handle, m, ncol, k, batchCount, alpha2, A2 + s * lda, lda, strideA, B2 + s * ldb, ldb, strideB, beta2, C2 + s * ldc + s, ldc, strideC);
        call_gemm_tn_strided_batched<BACKEND>(stream, handle, m, ncol, k, batchCount, alpha3, A3 + s * lda, lda, strideA, B3 + s * ldb, ldb, strideB, beta3, C3 + s * ldc + s, ldc, strideC);
    }
}

//------------------------------
// trmm
//------------------------------
template <Backend BACKEND, cublasFillMode_t EFF_UPLO>
inline void block_trmm_left_1(
    const cudaStream_t stream, Handle_t &handle, int m, int n,
    const hi_t<BACKEND> *alpha,
    const low_t<BACKEND> *A, size_t lda,
    const low_t<BACKEND> *B, size_t ldb,
    const hi_t<BACKEND> *beta,
    hi_t<BACKEND> *C, size_t ldc //
) {
    static_assert(EFF_UPLO == CUBLAS_FILL_MODE_UPPER || EFF_UPLO == CUBLAS_FILL_MODE_LOWER,
                  "block_trmm_left requires effective UPPER or LOWER UPLO.");
    const int nB = handle.nB;

    if (m <= nB) {
        call_gemm_tn<BACKEND>(stream, handle, m, n, m, alpha, A, lda, B, ldb, beta, C, ldc);
        return;
    }

    const int nBlocks = (m + nB - 1) / nB;
    for (int ib = 0; ib < nBlocks; ++ib) {
        const int s  = ib * nB;
        const int mm = std::min(nB, m - s);

        if constexpr (EFF_UPLO == CUBLAS_FILL_MODE_UPPER) {
            const int rem = m - s;
            call_gemm_tn<BACKEND>(stream, handle, mm, n, rem, alpha, A + s * lda + s, lda, B + s, ldb, beta, C + s, ldc);
        } else {
            const int lead = s + mm;
            call_gemm_tn<BACKEND>(stream, handle, mm, n, lead, alpha, A + s * lda, lda, B, ldb, beta, C + s, ldc);
        }
    }
}

template <Backend BACKEND, cublasFillMode_t EFF_UPLO>
inline void block_trmm_left_3(
    const cudaStream_t stream, Handle_t &handle, int m, int n,
    const hi_t<BACKEND> *alpha1, const hi_t<BACKEND> *alpha2, const hi_t<BACKEND> *alpha3,
    const low_t<BACKEND> *A1, const low_t<BACKEND> *A2, const low_t<BACKEND> *A3, size_t lda,
    const low_t<BACKEND> *B1, const low_t<BACKEND> *B2, const low_t<BACKEND> *B3, size_t ldb,
    const hi_t<BACKEND> *beta1, const hi_t<BACKEND> *beta2, const hi_t<BACKEND> *beta3,
    hi_t<BACKEND> *C1, hi_t<BACKEND> *C2, hi_t<BACKEND> *C3, size_t ldc //
) {
    static_assert(EFF_UPLO == CUBLAS_FILL_MODE_UPPER || EFF_UPLO == CUBLAS_FILL_MODE_LOWER,
                  "block_trmm_left requires effective UPPER or LOWER UPLO.");
    const int nB = handle.nB;

    if (m <= nB) {
        call_gemm_tn<BACKEND>(stream, handle, m, n, m, alpha1, A1, lda, B1, ldb, beta1, C1, ldc);
        call_gemm_tn<BACKEND>(stream, handle, m, n, m, alpha2, A2, lda, B2, ldb, beta2, C2, ldc);
        call_gemm_tn<BACKEND>(stream, handle, m, n, m, alpha3, A3, lda, B3, ldb, beta3, C3, ldc);
        return;
    }

    const int nBlocks = (m + nB - 1) / nB;
    for (int ib = 0; ib < nBlocks; ++ib) {
        const int s  = ib * nB;
        const int mm = std::min(nB, m - s);

        if constexpr (EFF_UPLO == CUBLAS_FILL_MODE_UPPER) {
            const int rem = m - s;
            call_gemm_tn<BACKEND>(stream, handle, mm, n, rem, alpha1, A1 + s * lda + s, lda, B1 + s, ldb, beta1, C1 + s, ldc);
            call_gemm_tn<BACKEND>(stream, handle, mm, n, rem, alpha2, A2 + s * lda + s, lda, B2 + s, ldb, beta2, C2 + s, ldc);
            call_gemm_tn<BACKEND>(stream, handle, mm, n, rem, alpha3, A3 + s * lda + s, lda, B3 + s, ldb, beta3, C3 + s, ldc);
        } else {
            const int lead = s + mm;
            call_gemm_tn<BACKEND>(stream, handle, mm, n, lead, alpha1, A1 + s * lda, lda, B1, ldb, beta1, C1 + s, ldc);
            call_gemm_tn<BACKEND>(stream, handle, mm, n, lead, alpha2, A2 + s * lda, lda, B2, ldb, beta2, C2 + s, ldc);
            call_gemm_tn<BACKEND>(stream, handle, mm, n, lead, alpha3, A3 + s * lda, lda, B3, ldb, beta3, C3 + s, ldc);
        }
    }
}

template <Backend BACKEND, cublasFillMode_t EFF_UPLO>
inline void block_trmm_right_1(
    const cudaStream_t stream, Handle_t &handle, int m, int n,
    const hi_t<BACKEND> *alpha,
    const low_t<BACKEND> *A, size_t lda, // triangular n x n
    const low_t<BACKEND> *B, size_t ldb, // full operand, stored as B^T for call_gemm_tn
    const hi_t<BACKEND> *beta,
    hi_t<BACKEND> *C, size_t ldc //
) {
    static_assert(EFF_UPLO == CUBLAS_FILL_MODE_UPPER || EFF_UPLO == CUBLAS_FILL_MODE_LOWER,
                  "block_trmm_right requires effective UPPER or LOWER UPLO.");
    const int nB    = handle.nB;
    const int n_pad = int(common::padding(size_t(n)));

    if (n <= nB) {
        call_gemm_tn<BACKEND>(stream, handle, m, n, n_pad, alpha, B, ldb, A, lda, beta, C, ldc);
        return;
    }

    const int nBlocks = (n + nB - 1) / nB;
    for (int ib = 0; ib < nBlocks; ++ib) {
        const int s  = ib * nB;
        const int nn = std::min(nB, n - s);

        if constexpr (EFF_UPLO == CUBLAS_FILL_MODE_UPPER) {
            const bool last_block = (ib + 1 == nBlocks);
            const int lead        = last_block ? n_pad : (s + nn);
            call_gemm_tn<BACKEND>(stream, handle, m, nn, lead, alpha, B, ldb, A + s * lda, lda, beta, C + s * ldc, ldc);
        } else {
            const int rem = n_pad - s;
            // call_gemm_tn<BACKEND>(stream, handle, m, nn, rem, alpha, B, ldb, A + s * lda + s, lda, beta, C + s * ldc, ldc);

            call_gemm_tn<BACKEND>(stream, handle, m, nn, rem, alpha, B + s, ldb, A + s * lda + s, lda, beta, C + s * ldc, ldc);
        }
    }
}

template <Backend BACKEND, cublasFillMode_t EFF_UPLO>
inline void block_trmm_right_3(
    const cudaStream_t stream, Handle_t &handle, int m, int n,
    const hi_t<BACKEND> *alpha1, const hi_t<BACKEND> *alpha2, const hi_t<BACKEND> *alpha3,
    const low_t<BACKEND> *A1, const low_t<BACKEND> *A2, const low_t<BACKEND> *A3, size_t lda,
    const low_t<BACKEND> *B1, const low_t<BACKEND> *B2, const low_t<BACKEND> *B3, size_t ldb,
    const hi_t<BACKEND> *beta1, const hi_t<BACKEND> *beta2, const hi_t<BACKEND> *beta3,
    hi_t<BACKEND> *C1, hi_t<BACKEND> *C2, hi_t<BACKEND> *C3, size_t ldc //
) {
    static_assert(EFF_UPLO == CUBLAS_FILL_MODE_UPPER || EFF_UPLO == CUBLAS_FILL_MODE_LOWER,
                  "block_trmm_right requires effective UPPER or LOWER UPLO.");
    const int nB    = handle.nB;
    const int n_pad = int(common::padding(size_t(n)));

    if (n <= nB) {
        call_gemm_tn<BACKEND>(stream, handle, m, n, n_pad, alpha1, B1, ldb, A1, lda, beta1, C1, ldc);
        call_gemm_tn<BACKEND>(stream, handle, m, n, n_pad, alpha2, B2, ldb, A2, lda, beta2, C2, ldc);
        call_gemm_tn<BACKEND>(stream, handle, m, n, n_pad, alpha3, B3, ldb, A3, lda, beta3, C3, ldc);
        return;
    }

    const int nBlocks = (n + nB - 1) / nB;
    for (int ib = 0; ib < nBlocks; ++ib) {
        const int s  = ib * nB;
        const int nn = std::min(nB, n - s);

        if constexpr (EFF_UPLO == CUBLAS_FILL_MODE_UPPER) {
            const bool last_block = (ib + 1 == nBlocks);
            const int lead        = last_block ? n_pad : (s + nn);
            call_gemm_tn<BACKEND>(stream, handle, m, nn, lead, alpha1, B1, ldb, A1 + s * lda, lda, beta1, C1 + s * ldc, ldc);
            call_gemm_tn<BACKEND>(stream, handle, m, nn, lead, alpha2, B2, ldb, A2 + s * lda, lda, beta2, C2 + s * ldc, ldc);
            call_gemm_tn<BACKEND>(stream, handle, m, nn, lead, alpha3, B3, ldb, A3 + s * lda, lda, beta3, C3 + s * ldc, ldc);
        } else {
            const int rem = n_pad - s;
            // call_gemm_tn<BACKEND>(stream, handle, m, nn, rem, alpha1, B1, ldb, A1 + s * lda + s, lda, beta1, C1 + s * ldc, ldc);
            // call_gemm_tn<BACKEND>(stream, handle, m, nn, rem, alpha2, B2, ldb, A2 + s * lda + s, lda, beta2, C2 + s * ldc, ldc);
            // call_gemm_tn<BACKEND>(stream, handle, m, nn, rem, alpha3, B3, ldb, A3 + s * lda + s, lda, beta3, C3 + s * ldc, ldc);

            call_gemm_tn<BACKEND>(stream, handle, m, nn, rem, alpha1, B1 + s, ldb, A1 + s * lda + s, lda, beta1, C1 + s * ldc, ldc);
            call_gemm_tn<BACKEND>(stream, handle, m, nn, rem, alpha2, B2 + s, ldb, A2 + s * lda + s, lda, beta2, C2 + s * ldc, ldc);
            call_gemm_tn<BACKEND>(stream, handle, m, nn, rem, alpha3, B3 + s, ldb, A3 + s * lda + s, lda, beta3, C3 + s * ldc, ldc);
        }
    }
}

template <Backend BACKEND, cublasFillMode_t EFF_UPLO>
inline void block_trmm_left_1_strided_batched(
    const cudaStream_t stream, Handle_t &handle, int m, int n, int batchCount,
    const hi_t<BACKEND> *alpha,
    const low_t<BACKEND> *A, size_t lda, int64_t strideA,
    const low_t<BACKEND> *B, size_t ldb, int64_t strideB,
    const hi_t<BACKEND> *beta,
    hi_t<BACKEND> *C, size_t ldc, int64_t strideC //
) {
    static_assert(EFF_UPLO == CUBLAS_FILL_MODE_UPPER || EFF_UPLO == CUBLAS_FILL_MODE_LOWER,
                  "block_trmm_left requires effective UPPER or LOWER UPLO.");
    if (batchCount <= 0) return;

    const int nB = handle.nB;

    if (m <= nB) {
        call_gemm_tn_strided_batched<BACKEND>(stream, handle, m, n, m, batchCount, alpha, A, lda, strideA, B, ldb, strideB, beta, C, ldc, strideC);
        return;
    }

    const int nBlocks = (m + nB - 1) / nB;
    for (int ib = 0; ib < nBlocks; ++ib) {
        const int s  = ib * nB;
        const int mm = std::min(nB, m - s);

        if constexpr (EFF_UPLO == CUBLAS_FILL_MODE_UPPER) {
            const int rem = m - s;
            call_gemm_tn_strided_batched<BACKEND>(stream, handle, mm, n, rem, batchCount, alpha, A + s * lda + s, lda, strideA, B + s, ldb, strideB, beta, C + s, ldc, strideC);
        } else {
            const int lead = s + mm;
            call_gemm_tn_strided_batched<BACKEND>(stream, handle, mm, n, lead, batchCount, alpha, A + s * lda, lda, strideA, B, ldb, strideB, beta, C + s, ldc, strideC);
        }
    }
}

template <Backend BACKEND, cublasFillMode_t EFF_UPLO>
inline void block_trmm_left_3_strided_batched(
    const cudaStream_t stream, Handle_t &handle, int m, int n, int batchCount,
    const hi_t<BACKEND> *alpha1, const hi_t<BACKEND> *alpha2, const hi_t<BACKEND> *alpha3,
    const low_t<BACKEND> *A1, const low_t<BACKEND> *A2, const low_t<BACKEND> *A3, size_t lda, int64_t strideA,
    const low_t<BACKEND> *B1, const low_t<BACKEND> *B2, const low_t<BACKEND> *B3, size_t ldb, int64_t strideB,
    const hi_t<BACKEND> *beta1, const hi_t<BACKEND> *beta2, const hi_t<BACKEND> *beta3,
    hi_t<BACKEND> *C1, hi_t<BACKEND> *C2, hi_t<BACKEND> *C3, size_t ldc, int64_t strideC //
) {
    static_assert(EFF_UPLO == CUBLAS_FILL_MODE_UPPER || EFF_UPLO == CUBLAS_FILL_MODE_LOWER,
                  "block_trmm_left requires effective UPPER or LOWER UPLO.");
    if (batchCount <= 0) return;

    const int nB = handle.nB;

    auto gemm3 = [&](int mm, int nn, int kk,
                     const low_t<BACKEND> *a1, const low_t<BACKEND> *b1, hi_t<BACKEND> *c1,
                     const low_t<BACKEND> *a2, const low_t<BACKEND> *b2, hi_t<BACKEND> *c2,
                     const low_t<BACKEND> *a3, const low_t<BACKEND> *b3, hi_t<BACKEND> *c3) {
        call_gemm_tn_strided_batched<BACKEND>(stream, handle, mm, nn, kk, batchCount, alpha1, a1, lda, strideA, b1, ldb, strideB, beta1, c1, ldc, strideC);
        call_gemm_tn_strided_batched<BACKEND>(stream, handle, mm, nn, kk, batchCount, alpha2, a2, lda, strideA, b2, ldb, strideB, beta2, c2, ldc, strideC);
        call_gemm_tn_strided_batched<BACKEND>(stream, handle, mm, nn, kk, batchCount, alpha3, a3, lda, strideA, b3, ldb, strideB, beta3, c3, ldc, strideC);
    };

    if (m <= nB) {
        gemm3(m, n, m, A1, B1, C1, A2, B2, C2, A3, B3, C3);
        return;
    }

    const int nBlocks = (m + nB - 1) / nB;
    for (int ib = 0; ib < nBlocks; ++ib) {
        const int s  = ib * nB;
        const int mm = std::min(nB, m - s);

        if constexpr (EFF_UPLO == CUBLAS_FILL_MODE_UPPER) {
            const int rem = m - s;
            gemm3(mm, n, rem,
                  A1 + s * lda + s, B1 + s, C1 + s,
                  A2 + s * lda + s, B2 + s, C2 + s,
                  A3 + s * lda + s, B3 + s, C3 + s);
        } else {
            const int lead = s + mm;
            gemm3(mm, n, lead,
                  A1 + s * lda, B1, C1 + s,
                  A2 + s * lda, B2, C2 + s,
                  A3 + s * lda, B3, C3 + s);
        }
    }
}

template <Backend BACKEND, cublasFillMode_t EFF_UPLO>
inline void block_trmm_right_1_strided_batched(
    const cudaStream_t stream, Handle_t &handle, int m, int n, int batchCount,
    const hi_t<BACKEND> *alpha,
    const low_t<BACKEND> *A, size_t lda, int64_t strideA,
    const low_t<BACKEND> *B, size_t ldb, int64_t strideB,
    const hi_t<BACKEND> *beta,
    hi_t<BACKEND> *C, size_t ldc, int64_t strideC //
) {
    static_assert(EFF_UPLO == CUBLAS_FILL_MODE_UPPER || EFF_UPLO == CUBLAS_FILL_MODE_LOWER,
                  "block_trmm_right requires effective UPPER or LOWER UPLO.");
    if (batchCount <= 0) return;

    const int nB    = handle.nB;
    const int n_pad = int(common::padding(size_t(n)));

    if (n <= nB) {
        call_gemm_tn_strided_batched<BACKEND>(stream, handle, m, n, n_pad, batchCount, alpha, B, ldb, strideB, A, lda, strideA, beta, C, ldc, strideC);
        return;
    }

    const int nBlocks = (n + nB - 1) / nB;
    for (int ib = 0; ib < nBlocks; ++ib) {
        const int s  = ib * nB;
        const int nn = std::min(nB, n - s);

        if constexpr (EFF_UPLO == CUBLAS_FILL_MODE_UPPER) {
            const bool last_block = (ib + 1 == nBlocks);
            const int lead        = last_block ? n_pad : (s + nn);
            call_gemm_tn_strided_batched<BACKEND>(stream, handle, m, nn, lead, batchCount, alpha, B, ldb, strideB, A + s * lda, lda, strideA, beta, C + s * ldc, ldc, strideC);
        } else {
            const int rem = n_pad - s;
            // call_gemm_tn_strided_batched<BACKEND>(stream, handle, m, nn, rem, batchCount, alpha, B, ldb, strideB, A + s * lda + s, lda, strideA, beta, C + s * ldc, ldc, strideC);

            call_gemm_tn_strided_batched<BACKEND>(stream, handle, m, nn, rem, batchCount, alpha, B + s, ldb, strideB, A + s * lda + s, lda, strideA, beta, C + s * ldc, ldc, strideC);
        }
    }
}

template <Backend BACKEND, cublasFillMode_t EFF_UPLO>
inline void block_trmm_right_3_strided_batched(
    const cudaStream_t stream, Handle_t &handle, int m, int n, int batchCount,
    const hi_t<BACKEND> *alpha1, const hi_t<BACKEND> *alpha2, const hi_t<BACKEND> *alpha3,
    const low_t<BACKEND> *A1, const low_t<BACKEND> *A2, const low_t<BACKEND> *A3, size_t lda, int64_t strideA,
    const low_t<BACKEND> *B1, const low_t<BACKEND> *B2, const low_t<BACKEND> *B3, size_t ldb, int64_t strideB,
    const hi_t<BACKEND> *beta1, const hi_t<BACKEND> *beta2, const hi_t<BACKEND> *beta3,
    hi_t<BACKEND> *C1, hi_t<BACKEND> *C2, hi_t<BACKEND> *C3, size_t ldc, int64_t strideC //
) {
    static_assert(EFF_UPLO == CUBLAS_FILL_MODE_UPPER || EFF_UPLO == CUBLAS_FILL_MODE_LOWER,
                  "block_trmm_right requires effective UPPER or LOWER UPLO.");
    if (batchCount <= 0) return;

    const int nB    = handle.nB;
    const int n_pad = int(common::padding(size_t(n)));

    auto gemm3 = [&](int mm, int nn, int kk,
                     const low_t<BACKEND> *a1, const low_t<BACKEND> *b1, hi_t<BACKEND> *c1,
                     const low_t<BACKEND> *a2, const low_t<BACKEND> *b2, hi_t<BACKEND> *c2,
                     const low_t<BACKEND> *a3, const low_t<BACKEND> *b3, hi_t<BACKEND> *c3) {
        call_gemm_tn_strided_batched<BACKEND>(stream, handle, mm, nn, kk, batchCount, alpha1, b1, ldb, strideB, a1, lda, strideA, beta1, c1, ldc, strideC);
        call_gemm_tn_strided_batched<BACKEND>(stream, handle, mm, nn, kk, batchCount, alpha2, b2, ldb, strideB, a2, lda, strideA, beta2, c2, ldc, strideC);
        call_gemm_tn_strided_batched<BACKEND>(stream, handle, mm, nn, kk, batchCount, alpha3, b3, ldb, strideB, a3, lda, strideA, beta3, c3, ldc, strideC);
    };

    if (n <= nB) {
        gemm3(m, n, n_pad, A1, B1, C1, A2, B2, C2, A3, B3, C3);
        return;
    }

    const int nBlocks = (n + nB - 1) / nB;
    for (int ib = 0; ib < nBlocks; ++ib) {
        const int s  = ib * nB;
        const int nn = std::min(nB, n - s);

        if constexpr (EFF_UPLO == CUBLAS_FILL_MODE_UPPER) {
            const bool last_block = (ib + 1 == nBlocks);
            const int lead        = last_block ? n_pad : (s + nn);
            gemm3(m, nn, lead,
                  A1 + s * lda, B1, C1 + s * ldc,
                  A2 + s * lda, B2, C2 + s * ldc,
                  A3 + s * lda, B3, C3 + s * ldc);
        } else {
            const int rem = n_pad - s;
            // gemm3(m, nn, rem,
            //       A1 + s * lda + s, B1, C1 + s * ldc,
            //       A2 + s * lda + s, B2, C2 + s * ldc,
            //       A3 + s * lda + s, B3, C3 + s * ldc);

            gemm3(m, nn, rem,
                  A1 + s * lda + s, B1 + s, C1 + s * ldc,
                  A2 + s * lda + s, B2 + s, C2 + s * ldc,
                  A3 + s * lda + s, B3 + s, C3 + s * ldc);
        }
    }
}

template <Backend BACKEND, cublasSideMode_t SIDE, cublasFillMode_t EFF_UPLO>
inline void block_trmm_1(
    const cudaStream_t stream, Handle_t &handle,
    int m, int n,
    const hi_t<BACKEND> *alpha,
    const low_t<BACKEND> *A, size_t lda,
    const low_t<BACKEND> *B, size_t ldb,
    const hi_t<BACKEND> *beta,
    hi_t<BACKEND> *C, size_t ldc //
) {
    if constexpr (SIDE == CUBLAS_SIDE_LEFT) {
        block_trmm_left_1<BACKEND, EFF_UPLO>(stream, handle, m, n, alpha, A, lda, B, ldb, beta, C, ldc);
    } else {
        block_trmm_right_1<BACKEND, EFF_UPLO>(stream, handle, m, n, alpha, A, lda, B, ldb, beta, C, ldc);
    }
}

template <Backend BACKEND, cublasSideMode_t SIDE, cublasFillMode_t EFF_UPLO>
inline void block_trmm_3(
    const cudaStream_t stream, Handle_t &handle,
    int m, int n,
    const hi_t<BACKEND> *alpha1, const hi_t<BACKEND> *alpha2, const hi_t<BACKEND> *alpha3,
    const low_t<BACKEND> *A1, const low_t<BACKEND> *A2, const low_t<BACKEND> *A3, size_t lda,
    const low_t<BACKEND> *B1, const low_t<BACKEND> *B2, const low_t<BACKEND> *B3, size_t ldb,
    const hi_t<BACKEND> *beta1, const hi_t<BACKEND> *beta2, const hi_t<BACKEND> *beta3,
    hi_t<BACKEND> *C1, hi_t<BACKEND> *C2, hi_t<BACKEND> *C3, size_t ldc //
) {
    if constexpr (SIDE == CUBLAS_SIDE_LEFT) {
        block_trmm_left_3<BACKEND, EFF_UPLO>(stream, handle, m, n, alpha1, alpha2, alpha3, A1, A2, A3, lda, B1, B2, B3, ldb, beta1, beta2, beta3, C1, C2, C3, ldc);
    } else {
        block_trmm_right_3<BACKEND, EFF_UPLO>(stream, handle, m, n, alpha1, alpha2, alpha3, A1, A2, A3, lda, B1, B2, B3, ldb, beta1, beta2, beta3, C1, C2, C3, ldc);
    }
}

template <Backend BACKEND, cublasSideMode_t SIDE, cublasFillMode_t EFF_UPLO>
inline void block_trmm_1_strided_batched(
    const cudaStream_t stream, Handle_t &handle,
    int m, int n, int batchCount,
    const hi_t<BACKEND> *alpha,
    const low_t<BACKEND> *A, size_t lda, int64_t strideA,
    const low_t<BACKEND> *B, size_t ldb, int64_t strideB,
    const hi_t<BACKEND> *beta,
    hi_t<BACKEND> *C, size_t ldc, int64_t strideC //
) {
    if constexpr (SIDE == CUBLAS_SIDE_LEFT) {
        block_trmm_left_1_strided_batched<BACKEND, EFF_UPLO>(stream, handle, m, n, batchCount, alpha, A, lda, strideA, B, ldb, strideB, beta, C, ldc, strideC);
    } else {
        block_trmm_right_1_strided_batched<BACKEND, EFF_UPLO>(stream, handle, m, n, batchCount, alpha, A, lda, strideA, B, ldb, strideB, beta, C, ldc, strideC);
    }
}

template <Backend BACKEND, cublasSideMode_t SIDE, cublasFillMode_t EFF_UPLO>
inline void block_trmm_3_strided_batched(
    const cudaStream_t stream, Handle_t &handle,
    int m, int n, int batchCount,
    const hi_t<BACKEND> *alpha1, const hi_t<BACKEND> *alpha2, const hi_t<BACKEND> *alpha3,
    const low_t<BACKEND> *A1, const low_t<BACKEND> *A2, const low_t<BACKEND> *A3, size_t lda, int64_t strideA,
    const low_t<BACKEND> *B1, const low_t<BACKEND> *B2, const low_t<BACKEND> *B3, size_t ldb, int64_t strideB,
    const hi_t<BACKEND> *beta1, const hi_t<BACKEND> *beta2, const hi_t<BACKEND> *beta3,
    hi_t<BACKEND> *C1, hi_t<BACKEND> *C2, hi_t<BACKEND> *C3, size_t ldc, int64_t strideC //
) {
    if constexpr (SIDE == CUBLAS_SIDE_LEFT) {
        block_trmm_left_3_strided_batched<BACKEND, EFF_UPLO>(stream, handle, m, n, batchCount, alpha1, alpha2, alpha3, A1, A2, A3, lda, strideA, B1, B2, B3, ldb, strideB, beta1, beta2, beta3, C1, C2, C3, ldc, strideC);
    } else {
        block_trmm_right_3_strided_batched<BACKEND, EFF_UPLO>(stream, handle, m, n, batchCount, alpha1, alpha2, alpha3, A1, A2, A3, lda, strideA, B1, B2, B3, ldb, strideB, beta1, beta2, beta3, C1, C2, C3, ldc, strideC);
    }
}

//------------------------------
// trtrmm
//------------------------------
template <Backend BACKEND, cublasFillMode_t EFF_UPLO_A, cublasFillMode_t EFF_UPLO_B>
inline void block_trtrmm_1(
    const cudaStream_t stream, Handle_t &handle,
    int n,
    const hi_t<BACKEND> *alpha,
    const low_t<BACKEND> *A, size_t lda,
    const low_t<BACKEND> *B, size_t ldb,
    const hi_t<BACKEND> *beta,
    hi_t<BACKEND> *C, size_t ldc //
) {
    static_assert((EFF_UPLO_A == CUBLAS_FILL_MODE_UPPER || EFF_UPLO_A == CUBLAS_FILL_MODE_LOWER) &&
                      (EFF_UPLO_B == CUBLAS_FILL_MODE_UPPER || EFF_UPLO_B == CUBLAS_FILL_MODE_LOWER),
                  "block_trtrmm requires effective UPPER or LOWER UPLOs.");

    const int nB = handle.nB;

    if constexpr (EFF_UPLO_A == CUBLAS_FILL_MODE_UPPER && EFF_UPLO_B == CUBLAS_FILL_MODE_UPPER) {
        const int fullBlocks = n / nB;
        const int tail       = n - fullBlocks * nB;

        const int64_t diag_stride_A = int64_t(nB) * lda + nB;
        const int64_t diag_stride_B = int64_t(nB) * ldb + nB;
        const int64_t diag_stride_C = int64_t(nB) * ldc + nB;

        for (int d = 0; d < fullBlocks; ++d) {
            const int kk         = (d + 1) * nB;
            const int innerBatch = fullBlocks - d;
            if (innerBatch <= 0) continue;

            const low_t<BACKEND> *A0 = A;
            const low_t<BACKEND> *B0 = B + d * nB * ldb;
            hi_t<BACKEND> *C0        = C + d * nB * ldc;

            if (innerBatch == 1) {
                call_gemm_tn<BACKEND>(stream, handle, nB, nB, kk, alpha, A0, lda, B0, ldb, beta, C0, ldc);
            } else {
                call_gemm_tn_strided_batched<BACKEND>(stream, handle, nB, nB, kk, innerBatch, alpha, A0, lda, diag_stride_A, B0, ldb, diag_stride_B, beta, C0, ldc, diag_stride_C);
            }
        }

        if (tail > 0) {
            const int cj = fullBlocks * nB;
            const int nj = tail;

            for (int i = 0; i < fullBlocks; ++i) {
                const int ri = i * nB;
                const int kk = n - ri;
                call_gemm_tn<BACKEND>(stream, handle, nB, nj, kk, alpha, A + ri * lda + ri, lda, B + cj * ldb + ri, ldb, beta, C + cj * ldc + ri, ldc);
            }

            const int s = fullBlocks * nB;
            call_gemm_tn<BACKEND>(stream, handle, tail, tail, tail, alpha, A + s * lda + s, lda, B + s * ldb + s, ldb, beta, C + s * ldc + s, ldc);
        }
        return;
    }

    if constexpr (EFF_UPLO_A == CUBLAS_FILL_MODE_LOWER && EFF_UPLO_B == CUBLAS_FILL_MODE_LOWER) {
        const int fullBlocks = n / nB;
        const int tail       = n - fullBlocks * nB;

        const int64_t diag_stride_A = int64_t(nB) * lda + nB;
        const int64_t diag_stride_B = int64_t(nB) * ldb + nB;
        const int64_t diag_stride_C = int64_t(nB) * ldc + nB;

        for (int d = 0; d < fullBlocks; ++d) {
            const int kk         = (d + 1) * nB;
            const int innerBatch = fullBlocks - d;
            if (innerBatch <= 0) continue;

            const low_t<BACKEND> *A0 = A + d * nB * lda;
            const low_t<BACKEND> *B0 = B;
            hi_t<BACKEND> *C0        = C + d * nB;

            if (innerBatch == 1) {
                call_gemm_tn<BACKEND>(stream, handle, nB, nB, kk, alpha, A0, lda, B0, ldb, beta, C0, ldc);
            } else {
                call_gemm_tn_strided_batched<BACKEND>(stream, handle, nB, nB, kk, innerBatch, alpha, A0, lda, diag_stride_A, B0, ldb, diag_stride_B, beta, C0, ldc, diag_stride_C);
            }
        }

        if (tail > 0) {
            const int ri = fullBlocks * nB;
            const int mi = tail;

            for (int j = 0; j < fullBlocks; ++j) {
                const int cj = j * nB;
                const int kk = n - cj;
                call_gemm_tn<BACKEND>(stream, handle, mi, nB, kk, alpha, A + ri * lda + cj, lda, B + cj * ldb + cj, ldb, beta, C + cj * ldc + ri, ldc);
            }

            const int s = fullBlocks * nB;
            call_gemm_tn<BACKEND>(stream, handle, tail, tail, tail, alpha, A + s * lda + s, lda, B + s * ldb + s, ldb, beta, C + s * ldc + s, ldc);
        }
        return;
    }

    if constexpr (EFF_UPLO_A == CUBLAS_FILL_MODE_LOWER && EFF_UPLO_B == CUBLAS_FILL_MODE_UPPER) {
        const int blocks = (n + nB - 1) / nB;
        for (int p = 0; p < blocks; ++p) {
            const int s     = p * nB;
            const int kk    = std::min(nB, n - s);
            const int e     = s + kk;
            const int rem   = n - s;
            const int rem_l = n - e;

            call_gemm_tn<BACKEND>(stream, handle, kk, rem, e, alpha, A + s * lda, lda, B + s * ldb, ldb, beta, C + s * ldc + s, ldc);

            if (rem_l > 0) {
                call_gemm_tn<BACKEND>(stream, handle, rem_l, kk, e, alpha, A + e * lda, lda, B + s * ldb, ldb, beta, C + s * ldc + e, ldc);
            }
        }
        return;
    }

    // EFF_UPLO_A == UPPER && EFF_UPLO_B == LOWER
    const int blocks = (n + nB - 1) / nB;
    for (int p = 0; p < blocks; ++p) {
        const int s   = p * nB;
        const int kk  = std::min(nB, n - s);
        const int e   = s + kk;
        const int rem = n - s;

        call_gemm_tn<BACKEND>(stream, handle, kk, e, rem, alpha, A + s * lda + s, lda, B + s, ldb, beta, C + s, ldc);

        for (int j = p + 1; j < blocks; ++j) {
            const int cj   = j * nB;
            const int nj   = std::min(nB, n - cj);
            const int remj = n - cj;
            call_gemm_tn<BACKEND>(stream, handle, kk, nj, remj, alpha, A + s * lda + cj, lda, B + cj * ldb + cj, ldb, beta, C + cj * ldc + s, ldc);
        }
    }
}

template <Backend BACKEND, cublasFillMode_t EFF_UPLO_A, cublasFillMode_t EFF_UPLO_B>
inline void block_trtrmm_3(
    const cudaStream_t stream, Handle_t &handle,
    int n,
    const hi_t<BACKEND> *alpha1, const hi_t<BACKEND> *alpha2, const hi_t<BACKEND> *alpha3,
    const low_t<BACKEND> *A1, const low_t<BACKEND> *A2, const low_t<BACKEND> *A3, size_t lda,
    const low_t<BACKEND> *B1, const low_t<BACKEND> *B2, const low_t<BACKEND> *B3, size_t ldb,
    const hi_t<BACKEND> *beta1, const hi_t<BACKEND> *beta2, const hi_t<BACKEND> *beta3,
    hi_t<BACKEND> *C1, hi_t<BACKEND> *C2, hi_t<BACKEND> *C3, size_t ldc //
) {
    static_assert((EFF_UPLO_A == CUBLAS_FILL_MODE_UPPER || EFF_UPLO_A == CUBLAS_FILL_MODE_LOWER) &&
                      (EFF_UPLO_B == CUBLAS_FILL_MODE_UPPER || EFF_UPLO_B == CUBLAS_FILL_MODE_LOWER),
                  "block_trtrmm requires effective UPPER or LOWER UPLOs.");

    const int nB = handle.nB;

    auto gemm3 = [&](int m, int ncol, int kk,
                     const low_t<BACKEND> *a1, const low_t<BACKEND> *b1, hi_t<BACKEND> *c1,
                     const low_t<BACKEND> *a2, const low_t<BACKEND> *b2, hi_t<BACKEND> *c2,
                     const low_t<BACKEND> *a3, const low_t<BACKEND> *b3, hi_t<BACKEND> *c3) {
        call_gemm_tn<BACKEND>(stream, handle, m, ncol, kk, alpha1, a1, lda, b1, ldb, beta1, c1, ldc);
        call_gemm_tn<BACKEND>(stream, handle, m, ncol, kk, alpha2, a2, lda, b2, ldb, beta2, c2, ldc);
        call_gemm_tn<BACKEND>(stream, handle, m, ncol, kk, alpha3, a3, lda, b3, ldb, beta3, c3, ldc);
    };

    auto gemm3_batched = [&](int m, int ncol, int kk, int bcnt,
                             const low_t<BACKEND> *a1, int64_t stride_a, const low_t<BACKEND> *b1, int64_t stride_b, hi_t<BACKEND> *c1, int64_t stride_c,
                             const low_t<BACKEND> *a2, const low_t<BACKEND> *b2, hi_t<BACKEND> *c2,
                             const low_t<BACKEND> *a3, const low_t<BACKEND> *b3, hi_t<BACKEND> *c3) {
        call_gemm_tn_strided_batched<BACKEND>(stream, handle, m, ncol, kk, bcnt, alpha1, a1, lda, stride_a, b1, ldb, stride_b, beta1, c1, ldc, stride_c);
        call_gemm_tn_strided_batched<BACKEND>(stream, handle, m, ncol, kk, bcnt, alpha2, a2, lda, stride_a, b2, ldb, stride_b, beta2, c2, ldc, stride_c);
        call_gemm_tn_strided_batched<BACKEND>(stream, handle, m, ncol, kk, bcnt, alpha3, a3, lda, stride_a, b3, ldb, stride_b, beta3, c3, ldc, stride_c);
    };

    if constexpr (EFF_UPLO_A == CUBLAS_FILL_MODE_UPPER && EFF_UPLO_B == CUBLAS_FILL_MODE_UPPER) {
        const int fullBlocks = n / nB;
        const int tail       = n - fullBlocks * nB;

        const int64_t diag_stride_A = int64_t(nB) * lda + nB;
        const int64_t diag_stride_B = int64_t(nB) * ldb + nB;
        const int64_t diag_stride_C = int64_t(nB) * ldc + nB;

        for (int d = 0; d < fullBlocks; ++d) {
            const int kk         = (d + 1) * nB;
            const int innerBatch = fullBlocks - d;
            if (innerBatch <= 0) continue;

            const low_t<BACKEND> *A10 = A1;
            const low_t<BACKEND> *B10 = B1 + d * nB * ldb;
            hi_t<BACKEND> *C10        = C1 + d * nB * ldc;
            const low_t<BACKEND> *A20 = A2;
            const low_t<BACKEND> *B20 = B2 + d * nB * ldb;
            hi_t<BACKEND> *C20        = C2 + d * nB * ldc;
            const low_t<BACKEND> *A30 = A3;
            const low_t<BACKEND> *B30 = B3 + d * nB * ldb;
            hi_t<BACKEND> *C30        = C3 + d * nB * ldc;

            if (innerBatch == 1) {
                gemm3(nB, nB, kk, A10, B10, C10, A20, B20, C20, A30, B30, C30);
            } else {
                gemm3_batched(nB, nB, kk, innerBatch,
                              A10, diag_stride_A, B10, diag_stride_B, C10, diag_stride_C,
                              A20, B20, C20,
                              A30, B30, C30);
            }
        }

        if (tail > 0) {
            const int cj = fullBlocks * nB;
            const int nj = tail;

            for (int i = 0; i < fullBlocks; ++i) {
                const int ri = i * nB;
                const int kk = n - ri;
                gemm3(nB, nj, kk,
                      A1 + ri * lda + ri, B1 + cj * ldb + ri, C1 + cj * ldc + ri,
                      A2 + ri * lda + ri, B2 + cj * ldb + ri, C2 + cj * ldc + ri,
                      A3 + ri * lda + ri, B3 + cj * ldb + ri, C3 + cj * ldc + ri);
            }

            const int s = fullBlocks * nB;
            gemm3(tail, tail, tail,
                  A1 + s * lda + s, B1 + s * ldb + s, C1 + s * ldc + s,
                  A2 + s * lda + s, B2 + s * ldb + s, C2 + s * ldc + s,
                  A3 + s * lda + s, B3 + s * ldb + s, C3 + s * ldc + s);
        }
        return;
    }

    if constexpr (EFF_UPLO_A == CUBLAS_FILL_MODE_LOWER && EFF_UPLO_B == CUBLAS_FILL_MODE_LOWER) {
        const int fullBlocks = n / nB;
        const int tail       = n - fullBlocks * nB;

        const int64_t diag_stride_A = int64_t(nB) * lda + nB;
        const int64_t diag_stride_B = int64_t(nB) * ldb + nB;
        const int64_t diag_stride_C = int64_t(nB) * ldc + nB;

        for (int d = 0; d < fullBlocks; ++d) {
            const int kk         = (d + 1) * nB;
            const int innerBatch = fullBlocks - d;
            if (innerBatch <= 0) continue;

            const low_t<BACKEND> *A10 = A1 + d * nB * lda;
            const low_t<BACKEND> *B10 = B1;
            hi_t<BACKEND> *C10        = C1 + d * nB;
            const low_t<BACKEND> *A20 = A2 + d * nB * lda;
            const low_t<BACKEND> *B20 = B2;
            hi_t<BACKEND> *C20        = C2 + d * nB;
            const low_t<BACKEND> *A30 = A3 + d * nB * lda;
            const low_t<BACKEND> *B30 = B3;
            hi_t<BACKEND> *C30        = C3 + d * nB;

            if (innerBatch == 1) {
                gemm3(nB, nB, kk, A10, B10, C10, A20, B20, C20, A30, B30, C30);
            } else {
                gemm3_batched(nB, nB, kk, innerBatch,
                              A10, diag_stride_A, B10, diag_stride_B, C10, diag_stride_C,
                              A20, B20, C20,
                              A30, B30, C30);
            }
        }

        if (tail > 0) {
            const int ri = fullBlocks * nB;
            const int mi = tail;

            for (int j = 0; j < fullBlocks; ++j) {
                const int cj = j * nB;
                const int kk = n - cj;
                gemm3(mi, nB, kk,
                      A1 + ri * lda + cj, B1 + cj * ldb + cj, C1 + cj * ldc + ri,
                      A2 + ri * lda + cj, B2 + cj * ldb + cj, C2 + cj * ldc + ri,
                      A3 + ri * lda + cj, B3 + cj * ldb + cj, C3 + cj * ldc + ri);
            }

            const int s = fullBlocks * nB;
            gemm3(tail, tail, tail,
                  A1 + s * lda + s, B1 + s * ldb + s, C1 + s * ldc + s,
                  A2 + s * lda + s, B2 + s * ldb + s, C2 + s * ldc + s,
                  A3 + s * lda + s, B3 + s * ldb + s, C3 + s * ldc + s);
        }
        return;
    }

    if constexpr (EFF_UPLO_A == CUBLAS_FILL_MODE_LOWER && EFF_UPLO_B == CUBLAS_FILL_MODE_UPPER) {
        const int blocks = (n + nB - 1) / nB;
        for (int p = 0; p < blocks; ++p) {
            const int s     = p * nB;
            const int kk    = std::min(nB, n - s);
            const int e     = s + kk;
            const int rem   = n - s;
            const int rem_l = n - e;

            gemm3(kk, rem, e,
                  A1 + s * lda, B1 + s * ldb, C1 + s * ldc + s,
                  A2 + s * lda, B2 + s * ldb, C2 + s * ldc + s,
                  A3 + s * lda, B3 + s * ldb, C3 + s * ldc + s);

            if (rem_l > 0) {
                gemm3(rem_l, kk, e,
                      A1 + e * lda, B1 + s * ldb, C1 + s * ldc + e,
                      A2 + e * lda, B2 + s * ldb, C2 + s * ldc + e,
                      A3 + e * lda, B3 + s * ldb, C3 + s * ldc + e);
            }
        }
        return;
    }

    const int blocks = (n + nB - 1) / nB;
    for (int p = 0; p < blocks; ++p) {
        const int s   = p * nB;
        const int kk  = std::min(nB, n - s);
        const int e   = s + kk;
        const int rem = n - s;

        gemm3(kk, e, rem,
              A1 + s * lda + s, B1 + s, C1 + s,
              A2 + s * lda + s, B2 + s, C2 + s,
              A3 + s * lda + s, B3 + s, C3 + s);

        for (int j = p + 1; j < blocks; ++j) {
            const int cj   = j * nB;
            const int nj   = std::min(nB, n - cj);
            const int remj = n - cj;

            gemm3(kk, nj, remj,
                  A1 + s * lda + cj, B1 + cj * ldb + cj, C1 + cj * ldc + s,
                  A2 + s * lda + cj, B2 + cj * ldb + cj, C2 + cj * ldc + s,
                  A3 + s * lda + cj, B3 + cj * ldb + cj, C3 + cj * ldc + s);
        }
    }
}

template <Backend BACKEND, cublasFillMode_t EFF_UPLO_A, cublasFillMode_t EFF_UPLO_B>
inline void block_trtrmm_1_strided_batched(
    const cudaStream_t stream, Handle_t &handle,
    int n, int batchCount,
    const hi_t<BACKEND> *alpha,
    const low_t<BACKEND> *A, size_t lda, int64_t strideA,
    const low_t<BACKEND> *B, size_t ldb, int64_t strideB,
    const hi_t<BACKEND> *beta,
    hi_t<BACKEND> *C, size_t ldc, int64_t strideC //
) {
    static_assert((EFF_UPLO_A == CUBLAS_FILL_MODE_UPPER || EFF_UPLO_A == CUBLAS_FILL_MODE_LOWER) &&
                      (EFF_UPLO_B == CUBLAS_FILL_MODE_UPPER || EFF_UPLO_B == CUBLAS_FILL_MODE_LOWER),
                  "block_trtrmm requires effective UPPER or LOWER UPLOs.");
    if (batchCount <= 0) return;

    const int nB = handle.nB;

    auto gemm_mod_batched = [&](int mm, int nn, int kk,
                                const low_t<BACKEND> *a, const low_t<BACKEND> *b, hi_t<BACKEND> *c) {
        call_gemm_tn_strided_batched<BACKEND>(stream, handle, mm, nn, kk, batchCount, alpha, a, lda, strideA, b, ldb, strideB, beta, c, ldc, strideC);
    };

    auto gemm_diag_batched = [&](int mm, int nn, int kk, int innerBatch,
                                 const low_t<BACKEND> *a, int64_t diag_stride_A,
                                 const low_t<BACKEND> *b, int64_t diag_stride_B,
                                 hi_t<BACKEND> *c, int64_t diag_stride_C) {
        call_gemm_tn_strided_batched<BACKEND>(stream, handle, mm, nn, kk, innerBatch, alpha, a, lda, diag_stride_A, b, ldb, diag_stride_B, beta, c, ldc, diag_stride_C);
    };

    if constexpr (EFF_UPLO_A == CUBLAS_FILL_MODE_UPPER && EFF_UPLO_B == CUBLAS_FILL_MODE_UPPER) {
        const int fullBlocks = n / nB;
        const int tail       = n - fullBlocks * nB;

        const int64_t diag_stride_A = int64_t(nB) * lda + nB;
        const int64_t diag_stride_B = int64_t(nB) * ldb + nB;
        const int64_t diag_stride_C = int64_t(nB) * ldc + nB;

        for (int d = 0; d < fullBlocks; ++d) {
            const int kk         = (d + 1) * nB;
            const int innerBatch = fullBlocks - d;
            if (innerBatch <= 0) continue;

            if (innerBatch == 1) {
                gemm_mod_batched(nB, nB, kk,
                                 A,
                                 B + d * nB * ldb,
                                 C + d * nB * ldc);
            } else {
                for (int bidx = 0; bidx < batchCount; ++bidx) {
                    const low_t<BACKEND> *Ab = A + bidx * strideA;
                    const low_t<BACKEND> *Bb = B + bidx * strideB;
                    hi_t<BACKEND> *Cb        = C + bidx * strideC;
                    gemm_diag_batched(nB, nB, kk, innerBatch,
                                      Ab, diag_stride_A,
                                      Bb + d * nB * ldb, diag_stride_B,
                                      Cb + d * nB * ldc, diag_stride_C);
                }
            }
        }

        if (tail > 0) {
            const int cj = fullBlocks * nB;
            const int nj = tail;

            for (int i = 0; i < fullBlocks; ++i) {
                const int ri = i * nB;
                const int kk = n - ri;
                gemm_mod_batched(nB, nj, kk,
                                 A + ri * lda + ri,
                                 B + cj * ldb + ri,
                                 C + cj * ldc + ri);
            }

            const int s = fullBlocks * nB;
            gemm_mod_batched(tail, tail, tail,
                             A + s * lda + s,
                             B + s * ldb + s,
                             C + s * ldc + s);
        }
        return;
    }

    if constexpr (EFF_UPLO_A == CUBLAS_FILL_MODE_LOWER && EFF_UPLO_B == CUBLAS_FILL_MODE_LOWER) {
        const int fullBlocks = n / nB;
        const int tail       = n - fullBlocks * nB;

        const int64_t diag_stride_A = int64_t(nB) * lda + nB;
        const int64_t diag_stride_B = int64_t(nB) * ldb + nB;
        const int64_t diag_stride_C = int64_t(nB) * ldc + nB;

        for (int d = 0; d < fullBlocks; ++d) {
            const int kk         = (d + 1) * nB;
            const int innerBatch = fullBlocks - d;
            if (innerBatch <= 0) continue;

            if (innerBatch == 1) {
                gemm_mod_batched(nB, nB, kk,
                                 A + d * nB * lda,
                                 B,
                                 C + d * nB);
            } else {
                for (int bidx = 0; bidx < batchCount; ++bidx) {
                    const low_t<BACKEND> *Ab = A + bidx * strideA;
                    const low_t<BACKEND> *Bb = B + bidx * strideB;
                    hi_t<BACKEND> *Cb        = C + bidx * strideC;
                    gemm_diag_batched(nB, nB, kk, innerBatch,
                                      Ab + d * nB * lda, diag_stride_A,
                                      Bb, diag_stride_B,
                                      Cb + d * nB, diag_stride_C);
                }
            }
        }

        if (tail > 0) {
            const int ri = fullBlocks * nB;
            const int mi = tail;

            for (int j = 0; j < fullBlocks; ++j) {
                const int cj = j * nB;
                const int kk = n - cj;
                gemm_mod_batched(mi, nB, kk,
                                 A + ri * lda + cj,
                                 B + cj * ldb + cj,
                                 C + cj * ldc + ri);
            }

            const int s = fullBlocks * nB;
            gemm_mod_batched(tail, tail, tail,
                             A + s * lda + s,
                             B + s * ldb + s,
                             C + s * ldc + s);
        }
        return;
    }

    if constexpr (EFF_UPLO_A == CUBLAS_FILL_MODE_LOWER && EFF_UPLO_B == CUBLAS_FILL_MODE_UPPER) {
        const int blocks = (n + nB - 1) / nB;
        for (int p = 0; p < blocks; ++p) {
            const int s     = p * nB;
            const int kk    = std::min(nB, n - s);
            const int e     = s + kk;
            const int rem   = n - s;
            const int rem_l = n - e;

            gemm_mod_batched(kk, rem, e,
                             A + s * lda,
                             B + s * ldb,
                             C + s * ldc + s);

            if (rem_l > 0) {
                gemm_mod_batched(rem_l, kk, e,
                                 A + e * lda,
                                 B + s * ldb,
                                 C + s * ldc + e);
            }
        }
        return;
    }

    // EFF_UPLO_A == UPPER && EFF_UPLO_B == LOWER
    const int blocks = (n + nB - 1) / nB;
    for (int p = 0; p < blocks; ++p) {
        const int s   = p * nB;
        const int kk  = std::min(nB, n - s);
        const int e   = s + kk;
        const int rem = n - s;

        gemm_mod_batched(kk, e, rem,
                         A + s * lda + s,
                         B + s,
                         C + s);

        for (int j = p + 1; j < blocks; ++j) {
            const int cj   = j * nB;
            const int nj   = std::min(nB, n - cj);
            const int remj = n - cj;
            gemm_mod_batched(kk, nj, remj,
                             A + s * lda + cj,
                             B + cj * ldb + cj,
                             C + cj * ldc + s);
        }
    }
}

template <Backend BACKEND, cublasFillMode_t EFF_UPLO_A, cublasFillMode_t EFF_UPLO_B>
inline void block_trtrmm_3_strided_batched(
    const cudaStream_t stream, Handle_t &handle,
    int n, int batchCount,
    const hi_t<BACKEND> *alpha1, const hi_t<BACKEND> *alpha2, const hi_t<BACKEND> *alpha3,
    const low_t<BACKEND> *A1, const low_t<BACKEND> *A2, const low_t<BACKEND> *A3, size_t lda, int64_t strideA,
    const low_t<BACKEND> *B1, const low_t<BACKEND> *B2, const low_t<BACKEND> *B3, size_t ldb, int64_t strideB,
    const hi_t<BACKEND> *beta1, const hi_t<BACKEND> *beta2, const hi_t<BACKEND> *beta3,
    hi_t<BACKEND> *C1, hi_t<BACKEND> *C2, hi_t<BACKEND> *C3, size_t ldc, int64_t strideC //
) {
    static_assert((EFF_UPLO_A == CUBLAS_FILL_MODE_UPPER || EFF_UPLO_A == CUBLAS_FILL_MODE_LOWER) &&
                      (EFF_UPLO_B == CUBLAS_FILL_MODE_UPPER || EFF_UPLO_B == CUBLAS_FILL_MODE_LOWER),
                  "block_trtrmm requires effective UPPER or LOWER UPLOs.");
    if (batchCount <= 0) return;

    const int nB = handle.nB;

    auto gemm3_mod_batched = [&](int mm, int nn, int kk,
                                 const low_t<BACKEND> *a1, const low_t<BACKEND> *b1, hi_t<BACKEND> *c1,
                                 const low_t<BACKEND> *a2, const low_t<BACKEND> *b2, hi_t<BACKEND> *c2,
                                 const low_t<BACKEND> *a3, const low_t<BACKEND> *b3, hi_t<BACKEND> *c3) {
        call_gemm_tn_strided_batched<BACKEND>(stream, handle, mm, nn, kk, batchCount, alpha1, a1, lda, strideA, b1, ldb, strideB, beta1, c1, ldc, strideC);
        call_gemm_tn_strided_batched<BACKEND>(stream, handle, mm, nn, kk, batchCount, alpha2, a2, lda, strideA, b2, ldb, strideB, beta2, c2, ldc, strideC);
        call_gemm_tn_strided_batched<BACKEND>(stream, handle, mm, nn, kk, batchCount, alpha3, a3, lda, strideA, b3, ldb, strideB, beta3, c3, ldc, strideC);
    };

    auto gemm3_diag_batched = [&](int mm, int nn, int kk, int innerBatch,
                                  const low_t<BACKEND> *a1, int64_t stride_a, const low_t<BACKEND> *b1, int64_t stride_b, hi_t<BACKEND> *c1, int64_t stride_c,
                                  const low_t<BACKEND> *a2, const low_t<BACKEND> *b2, hi_t<BACKEND> *c2,
                                  const low_t<BACKEND> *a3, const low_t<BACKEND> *b3, hi_t<BACKEND> *c3) {
        call_gemm_tn_strided_batched<BACKEND>(stream, handle, mm, nn, kk, innerBatch, alpha1, a1, lda, stride_a, b1, ldb, stride_b, beta1, c1, ldc, stride_c);
        call_gemm_tn_strided_batched<BACKEND>(stream, handle, mm, nn, kk, innerBatch, alpha2, a2, lda, stride_a, b2, ldb, stride_b, beta2, c2, ldc, stride_c);
        call_gemm_tn_strided_batched<BACKEND>(stream, handle, mm, nn, kk, innerBatch, alpha3, a3, lda, stride_a, b3, ldb, stride_b, beta3, c3, ldc, stride_c);
    };

    if constexpr (EFF_UPLO_A == CUBLAS_FILL_MODE_UPPER && EFF_UPLO_B == CUBLAS_FILL_MODE_UPPER) {
        const int fullBlocks = n / nB;
        const int tail       = n - fullBlocks * nB;

        const int64_t diag_stride_A = int64_t(nB) * lda + nB;
        const int64_t diag_stride_B = int64_t(nB) * ldb + nB;
        const int64_t diag_stride_C = int64_t(nB) * ldc + nB;

        for (int d = 0; d < fullBlocks; ++d) {
            const int kk         = (d + 1) * nB;
            const int innerBatch = fullBlocks - d;
            if (innerBatch <= 0) continue;

            if (innerBatch == 1) {
                gemm3_mod_batched(nB, nB, kk,
                                  A1, B1 + d * nB * ldb, C1 + d * nB * ldc,
                                  A2, B2 + d * nB * ldb, C2 + d * nB * ldc,
                                  A3, B3 + d * nB * ldb, C3 + d * nB * ldc);
            } else {
                for (int bidx = 0; bidx < batchCount; ++bidx) {
                    const low_t<BACKEND> *A1b = A1 + bidx * strideA;
                    const low_t<BACKEND> *B1b = B1 + bidx * strideB;
                    hi_t<BACKEND> *C1b        = C1 + bidx * strideC;
                    const low_t<BACKEND> *A2b = A2 + bidx * strideA;
                    const low_t<BACKEND> *B2b = B2 + bidx * strideB;
                    hi_t<BACKEND> *C2b        = C2 + bidx * strideC;
                    const low_t<BACKEND> *A3b = A3 + bidx * strideA;
                    const low_t<BACKEND> *B3b = B3 + bidx * strideB;
                    hi_t<BACKEND> *C3b        = C3 + bidx * strideC;
                    gemm3_diag_batched(nB, nB, kk, innerBatch,
                                       A1b, diag_stride_A, B1b + d * nB * ldb, diag_stride_B, C1b + d * nB * ldc, diag_stride_C,
                                       A2b, B2b + d * nB * ldb, C2b + d * nB * ldc,
                                       A3b, B3b + d * nB * ldb, C3b + d * nB * ldc);
                }
            }
        }

        if (tail > 0) {
            const int cj = fullBlocks * nB;
            const int nj = tail;

            for (int i = 0; i < fullBlocks; ++i) {
                const int ri = i * nB;
                const int kk = n - ri;
                gemm3_mod_batched(nB, nj, kk,
                                  A1 + ri * lda + ri, B1 + cj * ldb + ri, C1 + cj * ldc + ri,
                                  A2 + ri * lda + ri, B2 + cj * ldb + ri, C2 + cj * ldc + ri,
                                  A3 + ri * lda + ri, B3 + cj * ldb + ri, C3 + cj * ldc + ri);
            }

            const int s = fullBlocks * nB;
            gemm3_mod_batched(tail, tail, tail,
                              A1 + s * lda + s, B1 + s * ldb + s, C1 + s * ldc + s,
                              A2 + s * lda + s, B2 + s * ldb + s, C2 + s * ldc + s,
                              A3 + s * lda + s, B3 + s * ldb + s, C3 + s * ldc + s);
        }
        return;
    }

    if constexpr (EFF_UPLO_A == CUBLAS_FILL_MODE_LOWER && EFF_UPLO_B == CUBLAS_FILL_MODE_LOWER) {
        const int fullBlocks = n / nB;
        const int tail       = n - fullBlocks * nB;

        const int64_t diag_stride_A = int64_t(nB) * lda + nB;
        const int64_t diag_stride_B = int64_t(nB) * ldb + nB;
        const int64_t diag_stride_C = int64_t(nB) * ldc + nB;

        for (int d = 0; d < fullBlocks; ++d) {
            const int kk         = (d + 1) * nB;
            const int innerBatch = fullBlocks - d;
            if (innerBatch <= 0) continue;

            if (innerBatch == 1) {
                gemm3_mod_batched(nB, nB, kk,
                                  A1 + d * nB * lda, B1, C1 + d * nB,
                                  A2 + d * nB * lda, B2, C2 + d * nB,
                                  A3 + d * nB * lda, B3, C3 + d * nB);
            } else {
                for (int bidx = 0; bidx < batchCount; ++bidx) {
                    const low_t<BACKEND> *A1b = A1 + bidx * strideA;
                    const low_t<BACKEND> *B1b = B1 + bidx * strideB;
                    hi_t<BACKEND> *C1b        = C1 + bidx * strideC;
                    const low_t<BACKEND> *A2b = A2 + bidx * strideA;
                    const low_t<BACKEND> *B2b = B2 + bidx * strideB;
                    hi_t<BACKEND> *C2b        = C2 + bidx * strideC;
                    const low_t<BACKEND> *A3b = A3 + bidx * strideA;
                    const low_t<BACKEND> *B3b = B3 + bidx * strideB;
                    hi_t<BACKEND> *C3b        = C3 + bidx * strideC;
                    gemm3_diag_batched(nB, nB, kk, innerBatch,
                                       A1b + d * nB * lda, diag_stride_A, B1b, diag_stride_B, C1b + d * nB, diag_stride_C,
                                       A2b + d * nB * lda, B2b, C2b + d * nB,
                                       A3b + d * nB * lda, B3b, C3b + d * nB);
                }
            }
        }

        if (tail > 0) {
            const int ri = fullBlocks * nB;
            const int mi = tail;

            for (int j = 0; j < fullBlocks; ++j) {
                const int cj = j * nB;
                const int kk = n - cj;
                gemm3_mod_batched(mi, nB, kk,
                                  A1 + ri * lda + cj, B1 + cj * ldb + cj, C1 + cj * ldc + ri,
                                  A2 + ri * lda + cj, B2 + cj * ldb + cj, C2 + cj * ldc + ri,
                                  A3 + ri * lda + cj, B3 + cj * ldb + cj, C3 + cj * ldc + ri);
            }

            const int s = fullBlocks * nB;
            gemm3_mod_batched(tail, tail, tail,
                              A1 + s * lda + s, B1 + s * ldb + s, C1 + s * ldc + s,
                              A2 + s * lda + s, B2 + s * ldb + s, C2 + s * ldc + s,
                              A3 + s * lda + s, B3 + s * ldb + s, C3 + s * ldc + s);
        }
        return;
    }

    if constexpr (EFF_UPLO_A == CUBLAS_FILL_MODE_LOWER && EFF_UPLO_B == CUBLAS_FILL_MODE_UPPER) {
        const int blocks = (n + nB - 1) / nB;
        for (int p = 0; p < blocks; ++p) {
            const int s     = p * nB;
            const int kk    = std::min(nB, n - s);
            const int e     = s + kk;
            const int rem   = n - s;
            const int rem_l = n - e;

            gemm3_mod_batched(kk, rem, e,
                              A1 + s * lda, B1 + s * ldb, C1 + s * ldc + s,
                              A2 + s * lda, B2 + s * ldb, C2 + s * ldc + s,
                              A3 + s * lda, B3 + s * ldb, C3 + s * ldc + s);

            if (rem_l > 0) {
                gemm3_mod_batched(rem_l, kk, e,
                                  A1 + e * lda, B1 + s * ldb, C1 + s * ldc + e,
                                  A2 + e * lda, B2 + s * ldb, C2 + s * ldc + e,
                                  A3 + e * lda, B3 + s * ldb, C3 + s * ldc + e);
            }
        }
        return;
    }

    const int blocks = (n + nB - 1) / nB;
    for (int p = 0; p < blocks; ++p) {
        const int s   = p * nB;
        const int kk  = std::min(nB, n - s);
        const int e   = s + kk;
        const int rem = n - s;

        gemm3_mod_batched(kk, e, rem,
                          A1 + s * lda + s, B1 + s, C1 + s,
                          A2 + s * lda + s, B2 + s, C2 + s,
                          A3 + s * lda + s, B3 + s, C3 + s);

        for (int j = p + 1; j < blocks; ++j) {
            const int cj   = j * nB;
            const int nj   = std::min(nB, n - cj);
            const int remj = n - cj;
            gemm3_mod_batched(kk, nj, remj,
                              A1 + s * lda + cj, B1 + cj * ldb + cj, C1 + cj * ldc + s,
                              A2 + s * lda + cj, B2 + cj * ldb + cj, C2 + cj * ldc + s,
                              A3 + s * lda + cj, B3 + cj * ldb + cj, C3 + cj * ldc + s);
        }
    }
}

//------------------------------
// MatMulKind dispatch interface
//------------------------------
template <MatMulKind KIND, Backend BACKEND,
          cublasFillMode_t UPLO_A = CUBLAS_FILL_MODE_FULL,
          cublasFillMode_t UPLO_B = CUBLAS_FILL_MODE_FULL,
          cublasFillMode_t UPLO_C = CUBLAS_FILL_MODE_FULL>
inline void block_matmul_1(
    const cudaStream_t stream, Handle_t &handle,
    int m, int n, int k,
    const hi_t<BACKEND> *alpha,
    const low_t<BACKEND> *A, size_t lda,
    const low_t<BACKEND> *B, size_t ldb,
    const hi_t<BACKEND> *beta,
    hi_t<BACKEND> *C, size_t ldc //
) {
    if constexpr (KIND == MatMulKind::Gemm) {
        block_gemm_1<BACKEND>(stream, handle, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);

    } else if constexpr (KIND == MatMulKind::ATxB) {
        static_assert(UPLO_C != CUBLAS_FILL_MODE_FULL, "ATxB output requires UPLO_C.");
        block_ATxB_1<BACKEND, UPLO_C>(stream, handle, n, k, alpha, A, lda, B, ldb, beta, C, ldc);

    } else if constexpr (KIND == MatMulKind::ATxA) {
        static_assert(UPLO_C != CUBLAS_FILL_MODE_FULL, "ATxA output requires UPLO_C.");
        block_ATxB_1<BACKEND, UPLO_C>(stream, handle, n, k, alpha, A, lda, B, ldb, beta, C, ldc);

    } else if constexpr (KIND == MatMulKind::AHxA) {
        static_assert(UPLO_C != CUBLAS_FILL_MODE_FULL, "AHxA output requires UPLO_C.");
        block_ATxB_1<BACKEND, UPLO_C>(stream, handle, n, k, alpha, A, lda, B, ldb, beta, C, ldc);

    } else if constexpr (KIND == MatMulKind::TrmmLeft) {
        static_assert(UPLO_A != CUBLAS_FILL_MODE_FULL, "Left TRMM requires effective UPLO_A.");
        block_trmm_left_1<BACKEND, UPLO_A>(stream, handle, m, n, alpha, A, lda, B, ldb, beta, C, ldc);

    } else if constexpr (KIND == MatMulKind::TrmmRight) {
        static_assert(UPLO_B != CUBLAS_FILL_MODE_FULL, "Right TRMM requires effective UPLO_B.");
        block_trmm_right_1<BACKEND, UPLO_B>(stream, handle, m, n, alpha, A, lda, B, ldb, beta, C, ldc);

    } else if constexpr (KIND == MatMulKind::Trtrmm) {
        static_assert(UPLO_A != CUBLAS_FILL_MODE_FULL && UPLO_B != CUBLAS_FILL_MODE_FULL, "TRTRMM requires effective UPLO_A and UPLO_B.");
        block_trtrmm_1<BACKEND, UPLO_A, UPLO_B>(stream, handle, n, alpha, A, lda, B, ldb, beta, C, ldc);
    }
}

template <MatMulKind KIND, Backend BACKEND,
          cublasFillMode_t UPLO_A = CUBLAS_FILL_MODE_FULL,
          cublasFillMode_t UPLO_B = CUBLAS_FILL_MODE_FULL,
          cublasFillMode_t UPLO_C = CUBLAS_FILL_MODE_FULL>
inline void block_matmul_3(
    const cudaStream_t stream, Handle_t &handle,
    int m, int n, int k,
    const hi_t<BACKEND> *alpha1, const hi_t<BACKEND> *alpha2, const hi_t<BACKEND> *alpha3,
    const low_t<BACKEND> *A1, const low_t<BACKEND> *A2, const low_t<BACKEND> *A3, size_t lda,
    const low_t<BACKEND> *B1, const low_t<BACKEND> *B2, const low_t<BACKEND> *B3, size_t ldb,
    const hi_t<BACKEND> *beta1, const hi_t<BACKEND> *beta2, const hi_t<BACKEND> *beta3,
    hi_t<BACKEND> *C1, hi_t<BACKEND> *C2, hi_t<BACKEND> *C3, size_t ldc //
) {
    if constexpr (KIND == MatMulKind::Gemm) {
        block_gemm_3<BACKEND>(stream, handle, m, n, k, alpha1, alpha2, alpha3, A1, A2, A3, lda, B1, B2, B3, ldb, beta1, beta2, beta3, C1, C2, C3, ldc);

    } else if constexpr (KIND == MatMulKind::ATxB) {
        static_assert(UPLO_C != CUBLAS_FILL_MODE_FULL, "ATxB output requires UPLO_C.");
        block_ATxB_3<BACKEND, UPLO_C>(stream, handle, n, k, alpha1, alpha2, alpha3, A1, A2, A3, lda, B1, B2, B3, ldb, beta1, beta2, beta3, C1, C2, C3, ldc);

    } else if constexpr (KIND == MatMulKind::ATxA) {
        static_assert(UPLO_C != CUBLAS_FILL_MODE_FULL, "ATxA output requires UPLO_C.");
        block_ATxB_3<BACKEND, UPLO_C>(stream, handle, n, k, alpha1, alpha2, alpha3, A1, A2, A3, lda, B1, B2, B3, ldb, beta1, beta2, beta3, C1, C2, C3, ldc);

    } else if constexpr (KIND == MatMulKind::AHxA) {
        static_assert(UPLO_C != CUBLAS_FILL_MODE_FULL, "AHxA output requires UPLO_C.");
        block_ATxB_3<BACKEND, UPLO_C>(stream, handle, n, k, alpha1, alpha2, alpha3, A1, A2, A3, lda, B1, B2, B3, ldb, beta1, beta2, beta3, C1, C2, C3, ldc);

    } else if constexpr (KIND == MatMulKind::TrmmLeft) {
        static_assert(UPLO_A != CUBLAS_FILL_MODE_FULL, "Left TRMM requires effective UPLO_A.");
        block_trmm_left_3<BACKEND, UPLO_A>(stream, handle, m, n, alpha1, alpha2, alpha3, A1, A2, A3, lda, B1, B2, B3, ldb, beta1, beta2, beta3, C1, C2, C3, ldc);

    } else if constexpr (KIND == MatMulKind::TrmmRight) {
        static_assert(UPLO_B != CUBLAS_FILL_MODE_FULL, "Right TRMM requires effective UPLO_B.");
        block_trmm_right_3<BACKEND, UPLO_B>(stream, handle, m, n, alpha1, alpha2, alpha3, A1, A2, A3, lda, B1, B2, B3, ldb, beta1, beta2, beta3, C1, C2, C3, ldc);

    } else if constexpr (KIND == MatMulKind::Trtrmm) {
        static_assert(UPLO_A != CUBLAS_FILL_MODE_FULL && UPLO_B != CUBLAS_FILL_MODE_FULL, "TRTRMM requires effective UPLO_A and UPLO_B.");
        block_trtrmm_3<BACKEND, UPLO_A, UPLO_B>(stream, handle, n, alpha1, alpha2, alpha3, A1, A2, A3, lda, B1, B2, B3, ldb, beta1, beta2, beta3, C1, C2, C3, ldc);
    }
}

template <MatMulKind KIND, Backend BACKEND,
          cublasFillMode_t UPLO_A = CUBLAS_FILL_MODE_FULL,
          cublasFillMode_t UPLO_B = CUBLAS_FILL_MODE_FULL,
          cublasFillMode_t UPLO_C = CUBLAS_FILL_MODE_FULL>
inline void block_matmul_1_strided_batched(
    const cudaStream_t stream, Handle_t &handle,
    int m, int n, int k, int batchCount,
    const hi_t<BACKEND> *alpha,
    const low_t<BACKEND> *A, size_t lda, int64_t strideA,
    const low_t<BACKEND> *B, size_t ldb, int64_t strideB,
    const hi_t<BACKEND> *beta,
    hi_t<BACKEND> *C, size_t ldc, int64_t strideC //
) {
    if constexpr (KIND == MatMulKind::Gemm) {
        block_gemm_1_strided_batched<BACKEND>(stream, handle, m, n, k, batchCount, alpha, A, lda, strideA, B, ldb, strideB, beta, C, ldc, strideC);

    } else if constexpr (KIND == MatMulKind::ATxB) {
        static_assert(UPLO_C != CUBLAS_FILL_MODE_FULL, "ATxB output requires UPLO_C.");
        block_ATxB_1_strided_batched<BACKEND, UPLO_C>(stream, handle, n, k, batchCount, alpha, A, lda, strideA, B, ldb, strideB, beta, C, ldc, strideC);

    } else if constexpr (KIND == MatMulKind::ATxA) {
        static_assert(UPLO_C != CUBLAS_FILL_MODE_FULL, "ATxA output requires UPLO_C.");
        block_ATxB_1_strided_batched<BACKEND, UPLO_C>(stream, handle, n, k, batchCount, alpha, A, lda, strideA, B, ldb, strideB, beta, C, ldc, strideC);

    } else if constexpr (KIND == MatMulKind::AHxA) {
        static_assert(UPLO_C != CUBLAS_FILL_MODE_FULL, "AHxA output requires UPLO_C.");
        block_ATxB_1_strided_batched<BACKEND, UPLO_C>(stream, handle, n, k, batchCount, alpha, A, lda, strideA, B, ldb, strideB, beta, C, ldc, strideC);

    } else if constexpr (KIND == MatMulKind::TrmmLeft) {
        static_assert(UPLO_A != CUBLAS_FILL_MODE_FULL, "Left TRMM requires effective UPLO_A.");
        block_trmm_left_1_strided_batched<BACKEND, UPLO_A>(stream, handle, m, n, batchCount, alpha, A, lda, strideA, B, ldb, strideB, beta, C, ldc, strideC);

    } else if constexpr (KIND == MatMulKind::TrmmRight) {
        static_assert(UPLO_B != CUBLAS_FILL_MODE_FULL, "Right TRMM requires effective UPLO_B.");
        block_trmm_right_1_strided_batched<BACKEND, UPLO_B>(stream, handle, m, n, batchCount, alpha, A, lda, strideA, B, ldb, strideB, beta, C, ldc, strideC);

    } else if constexpr (KIND == MatMulKind::Trtrmm) {
        static_assert(UPLO_A != CUBLAS_FILL_MODE_FULL && UPLO_B != CUBLAS_FILL_MODE_FULL, "TRTRMM requires effective UPLO_A and UPLO_B.");
        block_trtrmm_1_strided_batched<BACKEND, UPLO_A, UPLO_B>(stream, handle, n, batchCount, alpha, A, lda, strideA, B, ldb, strideB, beta, C, ldc, strideC);
    }
}

template <MatMulKind KIND, Backend BACKEND,
          cublasFillMode_t UPLO_A = CUBLAS_FILL_MODE_FULL,
          cublasFillMode_t UPLO_B = CUBLAS_FILL_MODE_FULL,
          cublasFillMode_t UPLO_C = CUBLAS_FILL_MODE_FULL>
inline void block_matmul_3_strided_batched(
    const cudaStream_t stream, Handle_t &handle,
    int m, int n, int k, int batchCount,
    const hi_t<BACKEND> *alpha1, const hi_t<BACKEND> *alpha2, const hi_t<BACKEND> *alpha3,
    const low_t<BACKEND> *A1, const low_t<BACKEND> *A2, const low_t<BACKEND> *A3, size_t lda, int64_t strideA,
    const low_t<BACKEND> *B1, const low_t<BACKEND> *B2, const low_t<BACKEND> *B3, size_t ldb, int64_t strideB,
    const hi_t<BACKEND> *beta1, const hi_t<BACKEND> *beta2, const hi_t<BACKEND> *beta3,
    hi_t<BACKEND> *C1, hi_t<BACKEND> *C2, hi_t<BACKEND> *C3, size_t ldc, int64_t strideC //
) {
    if constexpr (KIND == MatMulKind::Gemm) {
        block_gemm_3_strided_batched<BACKEND>(stream, handle, m, n, k, batchCount, alpha1, alpha2, alpha3, A1, A2, A3, lda, strideA, B1, B2, B3, ldb, strideB, beta1, beta2, beta3, C1, C2, C3, ldc, strideC);

    } else if constexpr (KIND == MatMulKind::ATxB) {
        static_assert(UPLO_C != CUBLAS_FILL_MODE_FULL, "ATxB output requires UPLO_C.");
        block_ATxB_3_strided_batched<BACKEND, UPLO_C>(stream, handle, n, k, batchCount, alpha1, alpha2, alpha3, A1, A2, A3, lda, strideA, B1, B2, B3, ldb, strideB, beta1, beta2, beta3, C1, C2, C3, ldc, strideC);

    } else if constexpr (KIND == MatMulKind::ATxA) {
        static_assert(UPLO_C != CUBLAS_FILL_MODE_FULL, "ATxA output requires UPLO_C.");
        block_ATxB_3_strided_batched<BACKEND, UPLO_C>(stream, handle, n, k, batchCount, alpha1, alpha2, alpha3, A1, A2, A3, lda, strideA, B1, B2, B3, ldb, strideB, beta1, beta2, beta3, C1, C2, C3, ldc, strideC);

    } else if constexpr (KIND == MatMulKind::AHxA) {
        static_assert(UPLO_C != CUBLAS_FILL_MODE_FULL, "AHxA output requires UPLO_C.");
        block_ATxB_3_strided_batched<BACKEND, UPLO_C>(stream, handle, n, k, batchCount, alpha1, alpha2, alpha3, A1, A2, A3, lda, strideA, B1, B2, B3, ldb, strideB, beta1, beta2, beta3, C1, C2, C3, ldc, strideC);

    } else if constexpr (KIND == MatMulKind::TrmmLeft) {
        static_assert(UPLO_A != CUBLAS_FILL_MODE_FULL, "Left TRMM requires effective UPLO_A.");
        block_trmm_left_3_strided_batched<BACKEND, UPLO_A>(stream, handle, m, n, batchCount, alpha1, alpha2, alpha3, A1, A2, A3, lda, strideA, B1, B2, B3, ldb, strideB, beta1, beta2, beta3, C1, C2, C3, ldc, strideC);

    } else if constexpr (KIND == MatMulKind::TrmmRight) {
        static_assert(UPLO_B != CUBLAS_FILL_MODE_FULL, "Right TRMM requires effective UPLO_B.");
        block_trmm_right_3_strided_batched<BACKEND, UPLO_B>(stream, handle, m, n, batchCount, alpha1, alpha2, alpha3, A1, A2, A3, lda, strideA, B1, B2, B3, ldb, strideB, beta1, beta2, beta3, C1, C2, C3, ldc, strideC);

    } else if constexpr (KIND == MatMulKind::Trtrmm) {
        static_assert(UPLO_A != CUBLAS_FILL_MODE_FULL && UPLO_B != CUBLAS_FILL_MODE_FULL, "TRTRMM requires effective UPLO_A and UPLO_B.");
        block_trtrmm_3_strided_batched<BACKEND, UPLO_A, UPLO_B>(stream, handle, n, batchCount, alpha1, alpha2, alpha3, A1, A2, A3, lda, strideA, B1, B2, B3, ldb, strideB, beta1, beta2, beta3, C1, C2, C3, ldc, strideC);
    }
}

} // namespace gemmul8::common