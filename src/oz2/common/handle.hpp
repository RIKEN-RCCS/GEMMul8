#pragma once
#include "include.hpp"

namespace gemmul8::common {

template <Func FUNC>
inline int block_size_setting(int arch, cublasFillMode_t uplo_A, cublasFillMode_t uplo_B) noexcept {
    if constexpr (FUNC == Func::syrk || FUNC == Func::syrkx || FUNC == Func::herk || FUNC == Func::herkx) {
        switch (arch) {
        case 90: return 2048;
        case 100: return 2048;
        case 103: return 2048;
        default: return 1024;
        }
    } else if constexpr (FUNC == Func::trmm) {
        switch (arch) {
        case 90: return 2048;
        case 100: return 2048;
        case 103: return 2048;
        default: return 1024;
        }
    } else if constexpr (FUNC == Func::trtrmm) {
        switch (arch) {
        case 90: {
            if (uplo_A == CUBLAS_FILL_MODE_UPPER && uplo_B == CUBLAS_FILL_MODE_LOWER) return 4096;
            if (uplo_A == CUBLAS_FILL_MODE_LOWER && uplo_B == CUBLAS_FILL_MODE_UPPER) return 2048;
            return 1024;
        }
        case 100: {
            if (uplo_A == CUBLAS_FILL_MODE_UPPER && uplo_B == CUBLAS_FILL_MODE_LOWER) return 3072;
            if (uplo_A == CUBLAS_FILL_MODE_LOWER && uplo_B == CUBLAS_FILL_MODE_UPPER) return 3072;
            return 2048;
        }
        case 103: {
            if (uplo_A == CUBLAS_FILL_MODE_UPPER && uplo_B == CUBLAS_FILL_MODE_LOWER) return 3072;
            if (uplo_A == CUBLAS_FILL_MODE_LOWER && uplo_B == CUBLAS_FILL_MODE_UPPER) return 3072;
            return 2048;
        }
        default: return 1024;
        }
    } else {
        // gemm, symm, syr2k, her2k, hemm, her2k, trsm
        switch (arch) {
        case 90: return 8192;
        default: return 32768;
        }
    }
}

template <Func FUNC>
inline int block_size(int &arch, cublasFillMode_t uplo_A, cublasFillMode_t uplo_B) noexcept {
    if (arch == 0) {
#if defined(__CUDACC__)
    #if defined(GPU_ARCH)
        arch = GPU_ARCH;
    #else
        int dev = 0;
        if (cudaGetDevice(&dev) == cudaSuccess) {
            int major = 0, minor = 0;
            if (cudaDeviceGetAttribute(&major, cudaDevAttrComputeCapabilityMajor, dev) == cudaSuccess &&
                cudaDeviceGetAttribute(&minor, cudaDevAttrComputeCapabilityMinor, dev) == cudaSuccess) {
                arch = major * 10 + minor;
            }
        }
    #endif
#endif
    }
    return block_size_setting<FUNC>(arch, uplo_A, uplo_B);
}

enum class HandleKind : unsigned char {
    cuBLAS,
    cuBLASLt
};
struct CublasTag {};
struct CublasLtTag {};

enum LtBatchKind {
    None         = 0,
    Strided      = 1,
    PointerArray = 2
};

struct LtMatmulKey {
    int m, n, k;
    int batchCount;
    int64_t lda, ldb, ldc;
    int64_t strideA, strideB, strideC;
    LtBatchKind batchKind;

    bool operator==(const LtMatmulKey &o) const noexcept {
        return m == o.m &&
               n == o.n &&
               k == o.k &&
               batchCount == o.batchCount &&
               lda == o.lda &&
               ldb == o.ldb &&
               ldc == o.ldc &&
               strideA == o.strideA &&
               strideB == o.strideB &&
               strideC == o.strideC &&
               batchKind == o.batchKind;
    }
};

struct LtMatmulKeyHash {
    static inline void combine(size_t &h, const size_t v) noexcept {
        h ^= v + size_t(0x9e3779b97f4a7c15ull) + (h << 6) + (h >> 2);
    }
    size_t operator()(const LtMatmulKey &x) const noexcept {
        size_t h = 0;
        combine(h, static_cast<size_t>(x.m));
        combine(h, static_cast<size_t>(x.n));
        combine(h, static_cast<size_t>(x.k));
        combine(h, static_cast<size_t>(x.batchCount));

        combine(h, static_cast<size_t>(x.lda));
        combine(h, static_cast<size_t>(x.ldb));
        combine(h, static_cast<size_t>(x.ldc));

        combine(h, static_cast<size_t>(x.strideA));
        combine(h, static_cast<size_t>(x.strideB));
        combine(h, static_cast<size_t>(x.strideC));

        combine(h, static_cast<size_t>(x.batchKind));
        return h;
    }
};

struct LtMatmulPlan {
    cublasLtMatmulHeuristicResult_t heur{};
    size_t ws = 0;
};

struct Handle_t {
    HandleKind kind;
    cublasHandle_t cublas{};
    cublasLtHandle_t cublasLt{};

    cublasPointerMode_t ptrMode = CUBLAS_POINTER_MODE_HOST;

    void *workspace             = nullptr;
    size_t workspaceSizeInBytes = 0;

    void **Aarray = nullptr;
    void **Barray = nullptr;
    void **Carray = nullptr;

    cublasLtMatmulDesc_t opDesc     = nullptr;
    cublasLtMatmulPreference_t pref = nullptr;
    cublasLtMatrixLayout_t Adesc    = nullptr;
    cublasLtMatrixLayout_t Bdesc    = nullptr;
    cublasLtMatrixLayout_t Cdesc    = nullptr;

    int arch = 0;
    int nB   = 0;
    std::unordered_map<LtMatmulKey, LtMatmulPlan, LtMatmulKeyHash> plan_cache;

    bool fp8_k_blocking       = false;
    int fp8_k_block_first     = 0;
    int fp8_k_block_next      = 0;
    float fp8_modulus         = 0.0f;
    float fp8_neg_inv_modulus = 0.0f;
    float fp8_half_modulus    = 0.0f;

    Handle_t(CublasTag, cublasHandle_t h) : kind(HandleKind::cuBLAS), cublas(h) {}
    Handle_t(CublasLtTag, cublasLtHandle_t h) : kind(HandleKind::cuBLASLt), cublasLt(h) {}
};

template <Backend BACKEND, Func FUNC>
inline void set_handle(
    const cudaStream_t stream, Handle_t &h,
    int m, int n, int k,
    size_t lda, size_t ldb, size_t ldc,
    size_t workspaceSizeInBytes,
    cublasFillMode_t uplo_A = CUBLAS_FILL_MODE_FULL,
    cublasFillMode_t uplo_B = CUBLAS_FILL_MODE_FULL //
) {
    constexpr auto CUDA_R_LOW          = (BACKEND == Backend::INT8) ? CUDA_R_8I : CUDA_R_8F_E4M3;
    constexpr auto CUDA_R_HIGH         = (BACKEND == Backend::INT8) ? CUDA_R_32I : CUDA_R_32F;
    constexpr auto CUBLAS_COMPUTE_TYPE = (BACKEND == Backend::INT8) ? CUBLAS_COMPUTE_32I : CUBLAS_COMPUTE_32F;
    constexpr cublasOperation_t transa = CUBLAS_OP_T;
    constexpr cublasOperation_t transb = CUBLAS_OP_N;

    h.workspaceSizeInBytes = workspaceSizeInBytes;

    if (h.kind == HandleKind::cuBLAS) {
        cublasGetPointerMode(h.cublas, &h.ptrMode);
        cublasSetPointerMode(h.cublas, CUBLAS_POINTER_MODE_HOST);
    }

    h.nB = block_size<FUNC>(h.arch, uplo_A, uplo_B);

    if (h.kind == HandleKind::cuBLASLt) {
        cublasLtMatmulDescCreate(&h.opDesc, CUBLAS_COMPUTE_TYPE, CUDA_R_HIGH);
        cublasLtMatmulDescSetAttribute(h.opDesc, CUBLASLT_MATMUL_DESC_TRANSA, &transa, sizeof(transa));
        cublasLtMatmulDescSetAttribute(h.opDesc, CUBLASLT_MATMUL_DESC_TRANSB, &transb, sizeof(transb));

        if constexpr (BACKEND == Backend::FP8) {
            int32_t scale_mode = CUBLASLT_MATMUL_MATRIX_SCALE_SCALAR_32F;
            cublasLtMatmulDescSetAttribute(h.opDesc, CUBLASLT_MATMUL_DESC_A_SCALE_MODE, &scale_mode, sizeof(scale_mode));
            cublasLtMatmulDescSetAttribute(h.opDesc, CUBLASLT_MATMUL_DESC_B_SCALE_MODE, &scale_mode, sizeof(scale_mode));

#if defined(__CUDACC__)
            int8_t fast_accum = 0;
            cublasLtMatmulDescSetAttribute(h.opDesc, CUBLASLT_MATMUL_DESC_FAST_ACCUM, &fast_accum, sizeof(fast_accum));
            cublasLtMatmulDescSetAttribute(h.opDesc, CUBLASLT_MATMUL_DESC_SCALE_TYPE, &CUDA_R_HIGH, sizeof(CUDA_R_HIGH));
            cublasLtMatmulDescSetAttribute(h.opDesc, CUBLASLT_MATMUL_DESC_COMPUTE_TYPE, &CUBLAS_COMPUTE_TYPE, sizeof(CUBLAS_COMPUTE_TYPE));
#endif
        }

        cublasLtMatmulPreferenceCreate(&h.pref);
        cublasLtMatmulPreferenceSetAttribute(h.pref, CUBLASLT_MATMUL_PREF_MAX_WORKSPACE_BYTES, &workspaceSizeInBytes, sizeof(workspaceSizeInBytes));

        cublasLtMatrixLayoutCreate(&h.Adesc, CUDA_R_LOW, k, m, (int64_t)lda);
        cublasLtMatrixLayoutCreate(&h.Bdesc, CUDA_R_LOW, k, n, (int64_t)ldb);
        cublasLtMatrixLayoutCreate(&h.Cdesc, CUDA_R_HIGH, m, n, (int64_t)ldc);
    }
}

inline void cleanup_handle(Handle_t &h) {
    h.plan_cache.clear();
    if (h.kind == HandleKind::cuBLASLt) {
        if (h.Adesc) cublasLtMatrixLayoutDestroy(h.Adesc), h.Adesc = nullptr;
        if (h.Bdesc) cublasLtMatrixLayoutDestroy(h.Bdesc), h.Bdesc = nullptr;
        if (h.Cdesc) cublasLtMatrixLayoutDestroy(h.Cdesc), h.Cdesc = nullptr;
        if (h.pref) cublasLtMatmulPreferenceDestroy(h.pref), h.pref = nullptr;
        if (h.opDesc) cublasLtMatmulDescDestroy(h.opDesc), h.opDesc = nullptr;
    }
    if (h.kind == HandleKind::cuBLAS) {
        cublasSetPointerMode(h.cublas, h.ptrMode);
    }
}

inline void update_lt_layouts(
    Handle_t &h,
    int m, int n, int k,
    size_t lda, size_t ldb, size_t ldc,
    int batchCount  = 1,
    int64_t strideA = 0,
    int64_t strideB = 0,
    int64_t strideC = 0 //
) {
    const int64_t rowsA = k, colsA = m, ldA = (int64_t)lda;
    const int64_t rowsB = k, colsB = n, ldB = (int64_t)ldb;
    const int64_t rowsC = m, colsC = n, ldC = (int64_t)ldc;

    cublasLtMatrixLayoutSetAttribute(h.Adesc, CUBLASLT_MATRIX_LAYOUT_ROWS, &rowsA, sizeof(rowsA));
    cublasLtMatrixLayoutSetAttribute(h.Adesc, CUBLASLT_MATRIX_LAYOUT_COLS, &colsA, sizeof(colsA));
    cublasLtMatrixLayoutSetAttribute(h.Adesc, CUBLASLT_MATRIX_LAYOUT_LD, &ldA, sizeof(ldA));
    cublasLtMatrixLayoutSetAttribute(h.Adesc, CUBLASLT_MATRIX_LAYOUT_BATCH_COUNT, &batchCount, sizeof(batchCount));
    cublasLtMatrixLayoutSetAttribute(h.Adesc, CUBLASLT_MATRIX_LAYOUT_STRIDED_BATCH_OFFSET, &strideA, sizeof(strideA));

    cublasLtMatrixLayoutSetAttribute(h.Bdesc, CUBLASLT_MATRIX_LAYOUT_ROWS, &rowsB, sizeof(rowsB));
    cublasLtMatrixLayoutSetAttribute(h.Bdesc, CUBLASLT_MATRIX_LAYOUT_COLS, &colsB, sizeof(colsB));
    cublasLtMatrixLayoutSetAttribute(h.Bdesc, CUBLASLT_MATRIX_LAYOUT_LD, &ldB, sizeof(ldB));
    cublasLtMatrixLayoutSetAttribute(h.Bdesc, CUBLASLT_MATRIX_LAYOUT_BATCH_COUNT, &batchCount, sizeof(batchCount));
    cublasLtMatrixLayoutSetAttribute(h.Bdesc, CUBLASLT_MATRIX_LAYOUT_STRIDED_BATCH_OFFSET, &strideB, sizeof(strideB));

    cublasLtMatrixLayoutSetAttribute(h.Cdesc, CUBLASLT_MATRIX_LAYOUT_ROWS, &rowsC, sizeof(rowsC));
    cublasLtMatrixLayoutSetAttribute(h.Cdesc, CUBLASLT_MATRIX_LAYOUT_COLS, &colsC, sizeof(colsC));
    cublasLtMatrixLayoutSetAttribute(h.Cdesc, CUBLASLT_MATRIX_LAYOUT_LD, &ldC, sizeof(ldC));
    cublasLtMatrixLayoutSetAttribute(h.Cdesc, CUBLASLT_MATRIX_LAYOUT_BATCH_COUNT, &batchCount, sizeof(batchCount));
    cublasLtMatrixLayoutSetAttribute(h.Cdesc, CUBLASLT_MATRIX_LAYOUT_STRIDED_BATCH_OFFSET, &strideC, sizeof(strideC));

#if defined(__CUDACC__)
    const cublasLtBatchMode_t mode = CUBLASLT_BATCH_MODE_STRIDED;
    cublasLtMatrixLayoutSetAttribute(h.Adesc, CUBLASLT_MATRIX_LAYOUT_BATCH_MODE, &mode, sizeof(mode));
    cublasLtMatrixLayoutSetAttribute(h.Bdesc, CUBLASLT_MATRIX_LAYOUT_BATCH_MODE, &mode, sizeof(mode));
    cublasLtMatrixLayoutSetAttribute(h.Cdesc, CUBLASLT_MATRIX_LAYOUT_BATCH_MODE, &mode, sizeof(mode));
#endif
}

#if defined(__CUDACC__)
inline void update_lt_layouts_pointer_array(
    Handle_t &h,
    int m, int n, int k,
    size_t lda, size_t ldb, size_t ldc,
    int batchCount //
) {
    const int64_t rowsA = k, colsA = m, ldA = static_cast<int64_t>(lda);
    const int64_t rowsB = k, colsB = n, ldB = static_cast<int64_t>(ldb);
    const int64_t rowsC = m, colsC = n, ldC = static_cast<int64_t>(ldc);

    cublasLtMatrixLayoutSetAttribute(h.Adesc, CUBLASLT_MATRIX_LAYOUT_ROWS, &rowsA, sizeof(rowsA));
    cublasLtMatrixLayoutSetAttribute(h.Adesc, CUBLASLT_MATRIX_LAYOUT_COLS, &colsA, sizeof(colsA));
    cublasLtMatrixLayoutSetAttribute(h.Adesc, CUBLASLT_MATRIX_LAYOUT_LD, &ldA, sizeof(ldA));
    cublasLtMatrixLayoutSetAttribute(h.Adesc, CUBLASLT_MATRIX_LAYOUT_BATCH_COUNT, &batchCount, sizeof(batchCount));

    cublasLtMatrixLayoutSetAttribute(h.Bdesc, CUBLASLT_MATRIX_LAYOUT_ROWS, &rowsB, sizeof(rowsB));
    cublasLtMatrixLayoutSetAttribute(h.Bdesc, CUBLASLT_MATRIX_LAYOUT_COLS, &colsB, sizeof(colsB));
    cublasLtMatrixLayoutSetAttribute(h.Bdesc, CUBLASLT_MATRIX_LAYOUT_LD, &ldB, sizeof(ldB));
    cublasLtMatrixLayoutSetAttribute(h.Bdesc, CUBLASLT_MATRIX_LAYOUT_BATCH_COUNT, &batchCount, sizeof(batchCount));

    cublasLtMatrixLayoutSetAttribute(h.Cdesc, CUBLASLT_MATRIX_LAYOUT_ROWS, &rowsC, sizeof(rowsC));
    cublasLtMatrixLayoutSetAttribute(h.Cdesc, CUBLASLT_MATRIX_LAYOUT_COLS, &colsC, sizeof(colsC));
    cublasLtMatrixLayoutSetAttribute(h.Cdesc, CUBLASLT_MATRIX_LAYOUT_LD, &ldC, sizeof(ldC));
    cublasLtMatrixLayoutSetAttribute(h.Cdesc, CUBLASLT_MATRIX_LAYOUT_BATCH_COUNT, &batchCount, sizeof(batchCount));

    const cublasLtBatchMode_t mode = CUBLASLT_BATCH_MODE_POINTER_ARRAY;
    cublasLtMatrixLayoutSetAttribute(h.Adesc, CUBLASLT_MATRIX_LAYOUT_BATCH_MODE, &mode, sizeof(mode));
    cublasLtMatrixLayoutSetAttribute(h.Bdesc, CUBLASLT_MATRIX_LAYOUT_BATCH_MODE, &mode, sizeof(mode));
    cublasLtMatrixLayoutSetAttribute(h.Cdesc, CUBLASLT_MATRIX_LAYOUT_BATCH_MODE, &mode, sizeof(mode));
}
#endif

template <Backend BACKEND>
inline LtMatmulPlan &get_or_create_plan(
    Handle_t &h,
    int m, int n, int k,
    int batchCount,
    size_t lda, size_t ldb, size_t ldc,
    int64_t strideA,
    int64_t strideB,
    int64_t strideC,
    LtBatchKind batchKind //
) {
    LtMatmulKey key{
        m, n, k,
        batchCount,
        static_cast<int64_t>(lda),
        static_cast<int64_t>(ldb),
        static_cast<int64_t>(ldc),
        strideA,
        strideB,
        strideC,
        batchKind};

    auto it = h.plan_cache.find(key);
    if (it != h.plan_cache.end()) return it->second;

    LtMatmulPlan plan{};
    int returned = 0;

    cublasLtMatmulAlgoGetHeuristic(
        h.cublasLt,
        h.opDesc,
        h.Adesc, h.Bdesc, h.Cdesc, h.Cdesc,
        h.pref,
        1,
        &plan.heur,
        &returned);

    assert(returned > 0);

    plan.ws = plan.heur.workspaceSize;

    auto [pos, _] = h.plan_cache.emplace(key, plan);
    return pos->second;
}

} // namespace gemmul8::common
