#pragma once

enum class HandleKind : unsigned char {
    cuBLAS,
    cuBLASLt
};

struct Handle_t {
    HandleKind kind;
    cublasHandle_t cublas{};
    cublasLtHandle_t cublasLt{};

    void *workspace             = nullptr;
    size_t workspaceSizeInBytes = 0;

    cublasLtMatmulDesc_t opDesc       = nullptr;
    cublasLtMatmulPreference_t pref   = nullptr;
    cublasLtMatrixLayout_t Adesc      = nullptr;
    cublasLtMatrixLayout_t Bdesc_blk  = nullptr;
    cublasLtMatrixLayout_t Cdesc_blk  = nullptr;
    cublasLtMatrixLayout_t Bdesc_tail = nullptr;
    cublasLtMatrixLayout_t Cdesc_tail = nullptr;

    cublasLtMatmulHeuristicResult_t heur_blk{};
    cublasLtMatmulHeuristicResult_t heur_tail{};
    int got_blk  = 0;
    int got_tail = 0;
    int nn_blk   = 0;
    int nn_tail  = 0;

    Handle_t(cublasHandle_t h) : kind(HandleKind::cuBLAS), cublas(h) {}
    Handle_t(cublasLtHandle_t h) : kind(HandleKind::cuBLASLt), cublasLt(h) {}
};

template <gemmul8::Backend backend>
inline void set_handle(
    const cudaStream_t stream, Handle_t &h,
    int m, int n, int k,
    size_t lda, size_t ldb, size_t ldc,
    void *workspace, size_t workspaceSizeInBytes //
) {
    const cudaDataType_t CUDA_R_LOW               = (backend == gemmul8::Backend::INT8) ? CUDA_R_8I : CUDA_R_8F_E4M3;
    const cudaDataType_t CUDA_R_HIGH              = (backend == gemmul8::Backend::INT8) ? CUDA_R_32I : CUDA_R_32F;
    const cublasComputeType_t CUBLAS_COMPUTE_TYPE = (backend == gemmul8::Backend::INT8) ? CUBLAS_COMPUTE_32I : CUBLAS_COMPUTE_32F;
    const cublasOperation_t transa                = CUBLAS_OP_T;
    const cublasOperation_t transb                = CUBLAS_OP_N;

    h.workspace            = workspace;
    h.workspaceSizeInBytes = workspaceSizeInBytes;

    if (h.kind == HandleKind::cuBLAS && h.workspace && h.workspaceSizeInBytes) {
        cublasSetWorkspace(h.cublas, h.workspace, h.workspaceSizeInBytes);
    }

    if (h.kind == HandleKind::cuBLASLt) {
        cublasLtMatmulDescCreate(&h.opDesc, CUBLAS_COMPUTE_TYPE, CUDA_R_HIGH);
        cublasLtMatmulDescSetAttribute(h.opDesc, CUBLASLT_MATMUL_DESC_TRANSA, &transa, sizeof(transa));
        cublasLtMatmulDescSetAttribute(h.opDesc, CUBLASLT_MATMUL_DESC_TRANSB, &transb, sizeof(transb));

        cublasLtMatmulPreferenceCreate(&h.pref);
        cublasLtMatmulPreferenceSetAttribute(h.pref, CUBLASLT_MATMUL_PREF_MAX_WORKSPACE_BYTES, &workspaceSizeInBytes, sizeof(workspaceSizeInBytes));

        if (n <= 12288) {
            h.nn_blk  = n;
            h.nn_tail = n;
        } else {
            h.nn_blk  = 8192;
            int r     = n % 8192;
            h.nn_tail = (r <= 4096) ? (r + 8192) : r;
        }

        cublasLtMatrixLayoutCreate(&h.Adesc, CUDA_R_LOW, k, m, (int64_t)lda);
        cublasLtMatrixLayoutCreate(&h.Bdesc_blk, CUDA_R_LOW, k, h.nn_blk, (int64_t)ldb);
        cublasLtMatrixLayoutCreate(&h.Cdesc_blk, CUDA_R_HIGH, m, h.nn_blk, (int64_t)ldc);
        cublasLtMatrixLayoutCreate(&h.Bdesc_tail, CUDA_R_LOW, k, h.nn_tail, (int64_t)ldb);
        cublasLtMatrixLayoutCreate(&h.Cdesc_tail, CUDA_R_HIGH, m, h.nn_tail, (int64_t)ldc);

        {
            int returned = 0;
            cublasLtMatmulAlgoGetHeuristic(
                h.cublasLt, h.opDesc,
                h.Adesc, h.Bdesc_blk, h.Cdesc_blk, h.Cdesc_blk,
                h.pref, 1, &h.heur_blk, &returned);
            h.got_blk = returned;
        }
        {
            int returned = 0;
            cublasLtMatmulAlgoGetHeuristic(
                h.cublasLt, h.opDesc,
                h.Adesc, h.Bdesc_tail, h.Cdesc_tail, h.Cdesc_tail,
                h.pref, 1, &h.heur_tail, &returned);
            h.got_tail = returned;
        }
    }
}

inline void cleanup_handle(Handle_t &h) {
    if (h.kind == HandleKind::cuBLASLt) {
        if (h.Cdesc_tail) cublasLtMatrixLayoutDestroy(h.Cdesc_tail), h.Cdesc_tail = nullptr;
        if (h.Bdesc_tail) cublasLtMatrixLayoutDestroy(h.Bdesc_tail), h.Bdesc_tail = nullptr;
        if (h.Cdesc_blk) cublasLtMatrixLayoutDestroy(h.Cdesc_blk), h.Cdesc_blk = nullptr;
        if (h.Bdesc_blk) cublasLtMatrixLayoutDestroy(h.Bdesc_blk), h.Bdesc_blk = nullptr;
        if (h.Adesc) cublasLtMatrixLayoutDestroy(h.Adesc), h.Adesc = nullptr;
        if (h.pref) cublasLtMatmulPreferenceDestroy(h.pref), h.pref = nullptr;
        if (h.opDesc) cublasLtMatmulDescDestroy(h.opDesc), h.opDesc = nullptr;
        h.got_blk = h.got_tail = 0;
    }
}

//------------------------------
// 1x INT8-GEMM
//------------------------------
inline void gemm_low_prec_i8x1(
    const cudaStream_t stream, Handle_t &handle,
    int m, int n, int k,
    const int32_t *alpha, const int8_t *A, size_t lda, const int8_t *B, size_t ldb,
    const int32_t *beta, int32_t *C, size_t ldc //
) {
    if (handle.kind == HandleKind::cuBLAS) {

        // cuBLAS
        const int blk = 8192;
        int rem       = n;
        int offset    = 0;

        while (rem > 0) {
            int nn = (rem <= 12288) ? rem : blk;

            cublasGemmEx(handle.cublas, CUBLAS_OP_T, CUBLAS_OP_N,
                         m, nn, k,
                         alpha,
                         A, CUDA_R_8I, lda,
                         B + offset * ldb, CUDA_R_8I, ldb,
                         beta,
                         C + offset * ldc, CUDA_R_32I, ldc,
                         CUBLAS_COMPUTE_32I, CUBLAS_GEMM_DEFAULT);

            offset += nn;
            rem -= nn;
        }

    } else {

        // cuBLASLt
        const int blk = 8192;
        int rem       = n;
        int offset    = 0;

        while (rem > 0) {
            int nn = (rem <= 12288) ? rem : blk;

            const bool is_blk = (nn == handle.nn_blk);
            auto Bdesc        = is_blk ? handle.Bdesc_blk : handle.Bdesc_tail;
            auto Cdesc        = is_blk ? handle.Cdesc_blk : handle.Cdesc_tail;
            auto heur         = is_blk ? &handle.heur_blk : &handle.heur_tail;

            cublasLtMatmul(
                handle.cublasLt, handle.opDesc,
                alpha, A, handle.Adesc, B + offset * ldb, Bdesc,
                beta, C + offset * ldc, Cdesc, C + offset * ldc, Cdesc,
                &heur->algo, handle.workspace, handle.workspaceSizeInBytes, stream);

            offset += nn;
            rem -= nn;
        }
    }
}

//------------------------------
// 1x FP8-GEMM (cuBLASLt only)
//------------------------------
inline void gemm_low_prec_f8x1(
    const cudaStream_t stream, Handle_t &handle,
    int m, int n, int k,
    const float *alpha, const __nv_fp8_e4m3 *A, size_t lda, const __nv_fp8_e4m3 *B, size_t ldb,
    const float *beta, float *C, size_t ldc //
) {
    const int blk = 8192;
    int rem       = n;
    int offset    = 0;

    while (rem > 0) {
        int nn = (rem <= 12288) ? rem : blk;

        const bool is_blk = (nn == handle.nn_blk);
        auto Bdesc        = is_blk ? handle.Bdesc_blk : handle.Bdesc_tail;
        auto Cdesc        = is_blk ? handle.Cdesc_blk : handle.Cdesc_tail;
        auto heur         = is_blk ? &handle.heur_blk : &handle.heur_tail;

        cublasLtMatmul(
            handle.cublasLt, handle.opDesc,
            alpha, A, handle.Adesc, B + offset * ldb, Bdesc,
            beta, C + offset * ldc, Cdesc, C + offset * ldc, Cdesc,
            &heur->algo, handle.workspace, handle.workspaceSizeInBytes, stream);

        offset += nn;
        rem -= nn;
    }
}

//------------------------------
// 3x INT8-GEMM
//------------------------------
inline void gemm_low_prec_i8x3(
    const cudaStream_t stream, Handle_t &handle,
    int m, int n, int k,
    const int32_t *alpha1, const int32_t *alpha2, const int32_t *alpha3,
    const int8_t *A1, const int8_t *A2, const int8_t *A3, size_t lda,
    const int8_t *B1, const int8_t *B2, const int8_t *B3, size_t ldb,
    const int32_t *beta1, const int32_t *beta2, const int32_t *beta3,
    int32_t *C1, int32_t *C2, int32_t *C3, size_t ldc //
) {
    if (handle.kind == HandleKind::cuBLAS) {

        // cuBLAS

        const int blk = 8192;
        int rem       = n;
        int offset    = 0;

        while (rem > 0) {
            int nn = (rem <= 12288) ? rem : blk;

            cublasGemmEx(handle.cublas, CUBLAS_OP_T, CUBLAS_OP_N,
                         m, nn, k,
                         alpha1,
                         A1, CUDA_R_8I, lda,
                         B1 + offset * ldb, CUDA_R_8I, ldb,
                         beta1,
                         C1 + offset * ldc, CUDA_R_32I, ldc,
                         CUBLAS_COMPUTE_32I, CUBLAS_GEMM_DEFAULT);

            cublasGemmEx(handle.cublas, CUBLAS_OP_T, CUBLAS_OP_N,
                         m, nn, k,
                         alpha2,
                         A2, CUDA_R_8I, lda,
                         B2 + offset * ldb, CUDA_R_8I, ldb,
                         beta2,
                         C2 + offset * ldc, CUDA_R_32I, ldc,
                         CUBLAS_COMPUTE_32I, CUBLAS_GEMM_DEFAULT);

            cublasGemmEx(handle.cublas, CUBLAS_OP_T, CUBLAS_OP_N,
                         m, nn, k,
                         alpha3,
                         A3, CUDA_R_8I, lda,
                         B3 + offset * ldb, CUDA_R_8I, ldb,
                         beta3,
                         C3 + offset * ldc, CUDA_R_32I, ldc,
                         CUBLAS_COMPUTE_32I, CUBLAS_GEMM_DEFAULT);

            offset += nn;
            rem -= nn;
        }

    } else {

        // cuBLASLt
        const int blk = 8192;
        int rem       = n;
        int offset    = 0;

        while (rem > 0) {
            int nn = (rem <= 12288) ? rem : blk;

            const bool is_blk = (nn == handle.nn_blk);
            auto Bdesc        = is_blk ? handle.Bdesc_blk : handle.Bdesc_tail;
            auto Cdesc        = is_blk ? handle.Cdesc_blk : handle.Cdesc_tail;
            auto heur         = is_blk ? &handle.heur_blk : &handle.heur_tail;

            cublasLtMatmul(
                handle.cublasLt, handle.opDesc,
                alpha1, A1, handle.Adesc, B1 + offset * ldb, Bdesc,
                beta1, C1 + offset * ldc, Cdesc, C1 + offset * ldc, Cdesc,
                &heur->algo, handle.workspace, handle.workspaceSizeInBytes, stream);

            cublasLtMatmul(
                handle.cublasLt, handle.opDesc,
                alpha2, A2, handle.Adesc, B2 + offset * ldb, Bdesc,
                beta2, C2 + offset * ldc, Cdesc, C2 + offset * ldc, Cdesc,
                &heur->algo, handle.workspace, handle.workspaceSizeInBytes, stream);

            cublasLtMatmul(
                handle.cublasLt, handle.opDesc,
                alpha3, A3, handle.Adesc, B3 + offset * ldb, Bdesc,
                beta3, C3 + offset * ldc, Cdesc, C3 + offset * ldc, Cdesc,
                &heur->algo, handle.workspace, handle.workspaceSizeInBytes, stream);

            offset += nn;
            rem -= nn;
        }
    }
}

//------------------------------
// 3x FP8-GEMM using cuBLASLt
//------------------------------
inline void gemm_low_prec_f8x3(
    const cudaStream_t stream, Handle_t &handle,
    int m, int n, int k,
    const float *alpha1, const float *alpha2, const float *alpha3,
    const __nv_fp8_e4m3 *A1, const __nv_fp8_e4m3 *A2, const __nv_fp8_e4m3 *A3, size_t lda,
    const __nv_fp8_e4m3 *B1, const __nv_fp8_e4m3 *B2, const __nv_fp8_e4m3 *B3, size_t ldb,
    const float *beta1, const float *beta2, const float *beta3,
    float *C1, float *C2, float *C3, size_t ldc //
) {
    const int blk = 8192;
    int rem       = n;
    int offset    = 0;

    while (rem > 0) {
        int nn = (rem <= 12288) ? rem : blk;

        const bool is_blk = (nn == handle.nn_blk);
        auto Bdesc        = is_blk ? handle.Bdesc_blk : handle.Bdesc_tail;
        auto Cdesc        = is_blk ? handle.Cdesc_blk : handle.Cdesc_tail;
        auto heur         = is_blk ? &handle.heur_blk : &handle.heur_tail;

        cublasLtMatmul(
            handle.cublasLt, handle.opDesc,
            alpha1, A1, handle.Adesc, B1 + offset * ldb, Bdesc,
            beta1, C1 + offset * ldc, Cdesc, C1 + offset * ldc, Cdesc,
            &heur->algo, handle.workspace, handle.workspaceSizeInBytes, stream);

        cublasLtMatmul(
            handle.cublasLt, handle.opDesc,
            alpha2, A2, handle.Adesc, B2 + offset * ldb, Bdesc,
            beta2, C2 + offset * ldc, Cdesc, C2 + offset * ldc, Cdesc,
            &heur->algo, handle.workspace, handle.workspaceSizeInBytes, stream);

        cublasLtMatmul(
            handle.cublasLt, handle.opDesc,
            alpha3, A3, handle.Adesc, B3 + offset * ldb, Bdesc,
            beta3, C3 + offset * ldc, Cdesc, C3 + offset * ldc, Cdesc,
            &heur->algo, handle.workspace, handle.workspaceSizeInBytes, stream);

        offset += nn;
        rem -= nn;
    }
}

//------------------------------
// 9x FP8-GEMM using cuBLASLt
//------------------------------
inline void gemm_low_prec_f8x9(
    const cudaStream_t stream, Handle_t &handle,
    int m, int n, int k,
    const float *alpha,
    __nv_fp8_e4m3 **A, size_t lda,
    __nv_fp8_e4m3 **B, size_t ldb,
    const float *beta,
    float **C, size_t ldc //
) {
    const int blk = 8192;
    int rem       = n;
    int offset    = 0;

    while (rem > 0) {
        int nn = (rem <= 12288) ? rem : blk;

        const bool is_blk = (nn == handle.nn_blk);
        auto Bdesc        = is_blk ? handle.Bdesc_blk : handle.Bdesc_tail;
        auto Cdesc        = is_blk ? handle.Cdesc_blk : handle.Cdesc_tail;
        auto heur         = is_blk ? &handle.heur_blk : &handle.heur_tail;

#pragma unroll
        for (int i = 0; i < 9; ++i) {
            cublasLtMatmul(
                handle.cublasLt, handle.opDesc,
                alpha, A[i], handle.Adesc, B[i] + offset * ldb, Bdesc,
                beta, C[i] + offset * ldc, Cdesc, C[i] + offset * ldc, Cdesc,
                &heur->algo, handle.workspace, handle.workspaceSizeInBytes, stream);
        }

        offset += nn;
        rem -= nn;
    }
}
