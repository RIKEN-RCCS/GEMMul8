#pragma once
#include "include.hpp"
#include <cstdio>
#include "../../config/config.hpp"

#ifndef GEMMUL8_LT_TRACE
    #define GEMMUL8_LT_TRACE 0
#endif

namespace gemmul8::common {

inline void lt_check_status(cublasStatus_t status, const char *operation) {
    if (status != CUBLAS_STATUS_SUCCESS) {
        std::fprintf(stderr, "GEMMul8: %s failed (status=%d).\n", operation, int(status));
        std::abort();
    }
}

template <Func FUNC>
inline int block_size_setting(int arch, cublasFillMode_t uplo_A, cublasFillMode_t uplo_B) noexcept {
    if constexpr (FUNC == Func::syrk || FUNC == Func::syrkx || FUNC == Func::herk || FUNC == Func::herkx) {
        switch (arch) {
        case 90: return 2048;
        case 100: return 2048;
        case 103: return 2048;
        case 107: return 2048;
        default: return 1024;
        }
    } else if constexpr (FUNC == Func::trmm) {
        switch (arch) {
        case 90: return 2048;
        case 100: return 2048;
        case 103: return 2048;
        case 107: return 2048;
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
        case 107: {
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
        default: return 16384;
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
    size_t workspaceLimit;

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
               batchKind == o.batchKind &&
               workspaceLimit == o.workspaceLimit;
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
        combine(h, x.workspaceLimit);
        return h;
    }
};

struct LtMatmulPlan {
    cublasLtMatmulHeuristicResult_t heur{};
    size_t ws                   = 0;
    cublasStatus_t query_status = CUBLAS_STATUS_NOT_SUPPORTED;
    int returned                = 0;
    bool valid                  = false;
    LtMatmulKey key{};
};

struct Handle_t {
    HandleKind kind;
    config::ConfigSnapshot config{};
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

    LtMatmulKey active_layout{};
    bool active_layout_valid = false;

    int arch = 0;
    int nB   = 0;
    std::unordered_map<LtMatmulKey, LtMatmulPlan, LtMatmulKeyHash> plan_cache;

    bool modulus_complex      = false;
    bool matprod_k_blocking   = false;
    int matprod_k_block_first = 0;
    int matprod_k_block_next  = 0;
    unsigned modulus_idx      = 0;
    unsigned fp8_num_moduli   = 0;

    Handle_t(CublasTag, cublasHandle_t h)
        : kind(HandleKind::cuBLAS), config(config::get_config(h)), cublas(h) {}
    Handle_t(CublasLtTag, cublasLtHandle_t h)
        : kind(HandleKind::cuBLASLt), config(config::get_configLt(h)), cublasLt(h) {}
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

    h.plan_cache.clear();
    h.active_layout_valid  = false;
    h.modulus_complex      = false;
    h.matprod_k_blocking   = false;
    h.workspaceSizeInBytes = workspaceSizeInBytes;

    if (h.kind == HandleKind::cuBLAS) {
        cublasGetPointerMode(h.cublas, &h.ptrMode);
        cublasSetPointerMode(h.cublas, CUBLAS_POINTER_MODE_HOST);
    }

    h.nB = block_size<FUNC>(h.arch, uplo_A, uplo_B);

    if (h.kind == HandleKind::cuBLASLt) {
        lt_check_status(cublasLtMatmulDescCreate(&h.opDesc, CUBLAS_COMPUTE_TYPE, CUDA_R_HIGH), "cublasLtMatmulDescCreate");
        lt_check_status(cublasLtMatmulDescSetAttribute(h.opDesc, CUBLASLT_MATMUL_DESC_TRANSA, &transa, sizeof(transa)), "cublasLtMatmulDescSetAttribute");
        lt_check_status(cublasLtMatmulDescSetAttribute(h.opDesc, CUBLASLT_MATMUL_DESC_TRANSB, &transb, sizeof(transb)), "cublasLtMatmulDescSetAttribute");

        if constexpr (BACKEND == Backend::FP8) {
            int32_t scale_mode = CUBLASLT_MATMUL_MATRIX_SCALE_SCALAR_32F;
            lt_check_status(cublasLtMatmulDescSetAttribute(h.opDesc, CUBLASLT_MATMUL_DESC_A_SCALE_MODE, &scale_mode, sizeof(scale_mode)), "cublasLtMatmulDescSetAttribute");
            lt_check_status(cublasLtMatmulDescSetAttribute(h.opDesc, CUBLASLT_MATMUL_DESC_B_SCALE_MODE, &scale_mode, sizeof(scale_mode)), "cublasLtMatmulDescSetAttribute");

#if defined(__CUDACC__)
            int8_t fast_accum = 0;
            lt_check_status(cublasLtMatmulDescSetAttribute(h.opDesc, CUBLASLT_MATMUL_DESC_FAST_ACCUM, &fast_accum, sizeof(fast_accum)), "cublasLtMatmulDescSetAttribute");
            lt_check_status(cublasLtMatmulDescSetAttribute(h.opDesc, CUBLASLT_MATMUL_DESC_SCALE_TYPE, &CUDA_R_HIGH, sizeof(CUDA_R_HIGH)), "cublasLtMatmulDescSetAttribute");
            lt_check_status(cublasLtMatmulDescSetAttribute(h.opDesc, CUBLASLT_MATMUL_DESC_COMPUTE_TYPE, &CUBLAS_COMPUTE_TYPE, sizeof(CUBLAS_COMPUTE_TYPE)), "cublasLtMatmulDescSetAttribute");
#endif
        }

        lt_check_status(cublasLtMatmulPreferenceCreate(&h.pref), "cublasLtMatmulPreferenceCreate");
        lt_check_status(cublasLtMatmulPreferenceSetAttribute(h.pref, CUBLASLT_MATMUL_PREF_MAX_WORKSPACE_BYTES, &workspaceSizeInBytes, sizeof(workspaceSizeInBytes)), "cublasLtMatmulPreferenceSetAttribute");

        lt_check_status(cublasLtMatrixLayoutCreate(&h.Adesc, CUDA_R_LOW, k, m, (int64_t)lda), "cublasLtMatrixLayoutCreate");
        lt_check_status(cublasLtMatrixLayoutCreate(&h.Bdesc, CUDA_R_LOW, k, n, (int64_t)ldb), "cublasLtMatrixLayoutCreate");
        lt_check_status(cublasLtMatrixLayoutCreate(&h.Cdesc, CUDA_R_HIGH, m, n, (int64_t)ldc), "cublasLtMatrixLayoutCreate");
    }
}

inline void cleanup_handle(Handle_t &h) {
    h.plan_cache.clear();
    h.active_layout_valid = false;
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

    lt_check_status(cublasLtMatrixLayoutSetAttribute(h.Adesc, CUBLASLT_MATRIX_LAYOUT_ROWS, &rowsA, sizeof(rowsA)), "cublasLtMatrixLayoutSetAttribute");
    lt_check_status(cublasLtMatrixLayoutSetAttribute(h.Adesc, CUBLASLT_MATRIX_LAYOUT_COLS, &colsA, sizeof(colsA)), "cublasLtMatrixLayoutSetAttribute");
    lt_check_status(cublasLtMatrixLayoutSetAttribute(h.Adesc, CUBLASLT_MATRIX_LAYOUT_LD, &ldA, sizeof(ldA)), "cublasLtMatrixLayoutSetAttribute");
    lt_check_status(cublasLtMatrixLayoutSetAttribute(h.Adesc, CUBLASLT_MATRIX_LAYOUT_BATCH_COUNT, &batchCount, sizeof(batchCount)), "cublasLtMatrixLayoutSetAttribute");
    lt_check_status(cublasLtMatrixLayoutSetAttribute(h.Adesc, CUBLASLT_MATRIX_LAYOUT_STRIDED_BATCH_OFFSET, &strideA, sizeof(strideA)), "cublasLtMatrixLayoutSetAttribute");

    lt_check_status(cublasLtMatrixLayoutSetAttribute(h.Bdesc, CUBLASLT_MATRIX_LAYOUT_ROWS, &rowsB, sizeof(rowsB)), "cublasLtMatrixLayoutSetAttribute");
    lt_check_status(cublasLtMatrixLayoutSetAttribute(h.Bdesc, CUBLASLT_MATRIX_LAYOUT_COLS, &colsB, sizeof(colsB)), "cublasLtMatrixLayoutSetAttribute");
    lt_check_status(cublasLtMatrixLayoutSetAttribute(h.Bdesc, CUBLASLT_MATRIX_LAYOUT_LD, &ldB, sizeof(ldB)), "cublasLtMatrixLayoutSetAttribute");
    lt_check_status(cublasLtMatrixLayoutSetAttribute(h.Bdesc, CUBLASLT_MATRIX_LAYOUT_BATCH_COUNT, &batchCount, sizeof(batchCount)), "cublasLtMatrixLayoutSetAttribute");
    lt_check_status(cublasLtMatrixLayoutSetAttribute(h.Bdesc, CUBLASLT_MATRIX_LAYOUT_STRIDED_BATCH_OFFSET, &strideB, sizeof(strideB)), "cublasLtMatrixLayoutSetAttribute");

    lt_check_status(cublasLtMatrixLayoutSetAttribute(h.Cdesc, CUBLASLT_MATRIX_LAYOUT_ROWS, &rowsC, sizeof(rowsC)), "cublasLtMatrixLayoutSetAttribute");
    lt_check_status(cublasLtMatrixLayoutSetAttribute(h.Cdesc, CUBLASLT_MATRIX_LAYOUT_COLS, &colsC, sizeof(colsC)), "cublasLtMatrixLayoutSetAttribute");
    lt_check_status(cublasLtMatrixLayoutSetAttribute(h.Cdesc, CUBLASLT_MATRIX_LAYOUT_LD, &ldC, sizeof(ldC)), "cublasLtMatrixLayoutSetAttribute");
    lt_check_status(cublasLtMatrixLayoutSetAttribute(h.Cdesc, CUBLASLT_MATRIX_LAYOUT_BATCH_COUNT, &batchCount, sizeof(batchCount)), "cublasLtMatrixLayoutSetAttribute");
    lt_check_status(cublasLtMatrixLayoutSetAttribute(h.Cdesc, CUBLASLT_MATRIX_LAYOUT_STRIDED_BATCH_OFFSET, &strideC, sizeof(strideC)), "cublasLtMatrixLayoutSetAttribute");

#if defined(__CUDACC__)
    const cublasLtBatchMode_t mode = CUBLASLT_BATCH_MODE_STRIDED;
    lt_check_status(cublasLtMatrixLayoutSetAttribute(h.Adesc, CUBLASLT_MATRIX_LAYOUT_BATCH_MODE, &mode, sizeof(mode)), "cublasLtMatrixLayoutSetAttribute");
    lt_check_status(cublasLtMatrixLayoutSetAttribute(h.Bdesc, CUBLASLT_MATRIX_LAYOUT_BATCH_MODE, &mode, sizeof(mode)), "cublasLtMatrixLayoutSetAttribute");
    lt_check_status(cublasLtMatrixLayoutSetAttribute(h.Cdesc, CUBLASLT_MATRIX_LAYOUT_BATCH_MODE, &mode, sizeof(mode)), "cublasLtMatrixLayoutSetAttribute");
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

    lt_check_status(cublasLtMatrixLayoutSetAttribute(h.Adesc, CUBLASLT_MATRIX_LAYOUT_ROWS, &rowsA, sizeof(rowsA)), "cublasLtMatrixLayoutSetAttribute");
    lt_check_status(cublasLtMatrixLayoutSetAttribute(h.Adesc, CUBLASLT_MATRIX_LAYOUT_COLS, &colsA, sizeof(colsA)), "cublasLtMatrixLayoutSetAttribute");
    lt_check_status(cublasLtMatrixLayoutSetAttribute(h.Adesc, CUBLASLT_MATRIX_LAYOUT_LD, &ldA, sizeof(ldA)), "cublasLtMatrixLayoutSetAttribute");
    lt_check_status(cublasLtMatrixLayoutSetAttribute(h.Adesc, CUBLASLT_MATRIX_LAYOUT_BATCH_COUNT, &batchCount, sizeof(batchCount)), "cublasLtMatrixLayoutSetAttribute");

    lt_check_status(cublasLtMatrixLayoutSetAttribute(h.Bdesc, CUBLASLT_MATRIX_LAYOUT_ROWS, &rowsB, sizeof(rowsB)), "cublasLtMatrixLayoutSetAttribute");
    lt_check_status(cublasLtMatrixLayoutSetAttribute(h.Bdesc, CUBLASLT_MATRIX_LAYOUT_COLS, &colsB, sizeof(colsB)), "cublasLtMatrixLayoutSetAttribute");
    lt_check_status(cublasLtMatrixLayoutSetAttribute(h.Bdesc, CUBLASLT_MATRIX_LAYOUT_LD, &ldB, sizeof(ldB)), "cublasLtMatrixLayoutSetAttribute");
    lt_check_status(cublasLtMatrixLayoutSetAttribute(h.Bdesc, CUBLASLT_MATRIX_LAYOUT_BATCH_COUNT, &batchCount, sizeof(batchCount)), "cublasLtMatrixLayoutSetAttribute");

    lt_check_status(cublasLtMatrixLayoutSetAttribute(h.Cdesc, CUBLASLT_MATRIX_LAYOUT_ROWS, &rowsC, sizeof(rowsC)), "cublasLtMatrixLayoutSetAttribute");
    lt_check_status(cublasLtMatrixLayoutSetAttribute(h.Cdesc, CUBLASLT_MATRIX_LAYOUT_COLS, &colsC, sizeof(colsC)), "cublasLtMatrixLayoutSetAttribute");
    lt_check_status(cublasLtMatrixLayoutSetAttribute(h.Cdesc, CUBLASLT_MATRIX_LAYOUT_LD, &ldC, sizeof(ldC)), "cublasLtMatrixLayoutSetAttribute");
    lt_check_status(cublasLtMatrixLayoutSetAttribute(h.Cdesc, CUBLASLT_MATRIX_LAYOUT_BATCH_COUNT, &batchCount, sizeof(batchCount)), "cublasLtMatrixLayoutSetAttribute");

    const cublasLtBatchMode_t mode = CUBLASLT_BATCH_MODE_POINTER_ARRAY;
    lt_check_status(cublasLtMatrixLayoutSetAttribute(h.Adesc, CUBLASLT_MATRIX_LAYOUT_BATCH_MODE, &mode, sizeof(mode)), "cublasLtMatrixLayoutSetAttribute");
    lt_check_status(cublasLtMatrixLayoutSetAttribute(h.Bdesc, CUBLASLT_MATRIX_LAYOUT_BATCH_MODE, &mode, sizeof(mode)), "cublasLtMatrixLayoutSetAttribute");
    lt_check_status(cublasLtMatrixLayoutSetAttribute(h.Cdesc, CUBLASLT_MATRIX_LAYOUT_BATCH_MODE, &mode, sizeof(mode)), "cublasLtMatrixLayoutSetAttribute");
}
#endif

template <Backend BACKEND>
inline void report_lt_plan(const Handle_t &h, const LtMatmulKey &key,
                           const LtMatmulPlan &plan, const char *reason,
                           cublasStatus_t launch_status = CUBLAS_STATUS_SUCCESS) {
    std::fprintf(stderr,
                 "GEMMul8 Lt: %s; backend=%s arch=%d m=%d n=%d k=%d batch=%d kind=%d "
                 "lda=%lld ldb=%lld ldc=%lld strideA=%lld strideB=%lld strideC=%lld "
                 "workspace=%zu query_status=%d returned=%d candidate_state=%d launch_status=%d\n",
                 reason, BACKEND == Backend::FP8 ? "FP8" : "INT8", h.arch,
                 key.m, key.n, key.k, key.batchCount, int(key.batchKind),
                 (long long)key.lda, (long long)key.ldb, (long long)key.ldc,
                 (long long)key.strideA, (long long)key.strideB, (long long)key.strideC,
                 key.workspaceLimit, int(plan.query_status), plan.returned,
                 int(plan.heur.state), int(launch_status));
}

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
        m, n, k, batchCount,
        static_cast<int64_t>(lda), static_cast<int64_t>(ldb), static_cast<int64_t>(ldc),
        strideA, strideB, strideC, batchKind, h.workspaceSizeInBytes};

    auto it = h.plan_cache.find(key);
    // Negative entries mean an unsupported configuration, not an initialized
    // algorithm. Avoid repeatedly switching descriptors to a rejected mode.
    if (it != h.plan_cache.end() && !it->second.valid) return it->second;

    if (!h.active_layout_valid || !(h.active_layout == key)) {
#if defined(__CUDACC__)
        if (batchKind == LtBatchKind::PointerArray) {
            update_lt_layouts_pointer_array(h, m, n, k, lda, ldb, ldc, batchCount);
        } else
#endif
        {
            update_lt_layouts(h, m, n, k, lda, ldb, ldc, batchCount, strideA, strideB, strideC);
        }
        h.active_layout       = key;
        h.active_layout_valid = true;
    }
    if (it != h.plan_cache.end()) return it->second;

    // A previous set_handle may have used a different workspace budget.
    lt_check_status(cublasLtMatmulPreferenceSetAttribute(
                        h.pref, CUBLASLT_MATMUL_PREF_MAX_WORKSPACE_BYTES,
                        &h.workspaceSizeInBytes, sizeof(h.workspaceSizeInBytes)),
                    "cublasLtMatmulPreferenceSetAttribute(workspace)");

    LtMatmulPlan plan{};
    plan.key            = key;
    plan.heur.state     = CUBLAS_STATUS_NOT_SUPPORTED;
    constexpr int count = 8;
    cublasLtMatmulHeuristicResult_t candidates[count]{};
    plan.query_status = cublasLtMatmulAlgoGetHeuristic(
        h.cublasLt, h.opDesc, h.Adesc, h.Bdesc, h.Cdesc, h.Cdesc,
        h.pref, count, candidates, &plan.returned);

    if (plan.query_status != CUBLAS_STATUS_SUCCESS &&
        plan.query_status != CUBLAS_STATUS_NOT_SUPPORTED) {
        report_lt_plan<BACKEND>(h, key, plan, "heuristic API failed");
        std::abort();
    }
    if (plan.query_status == CUBLAS_STATUS_SUCCESS) {
        for (int i = 0; i < std::min(plan.returned, count); ++i) {
            plan.heur = candidates[i];
            if (candidates[i].state == CUBLAS_STATUS_SUCCESS &&
                candidates[i].workspaceSize <= h.workspaceSizeInBytes) {
                plan.ws    = candidates[i].workspaceSize;
                plan.valid = true;
                break;
            }
        }
    }
    if (!plan.valid && GEMMUL8_LT_TRACE) {
        report_lt_plan<BACKEND>(h, key, plan, "no usable algorithm; caller may fall back");
    }
    auto [pos, inserted] = h.plan_cache.emplace(key, plan);
    return pos->second;
}

template <Backend BACKEND>
inline bool supports_lt_pointer_batch(Handle_t &h, int m, int n, int k, int batch,
                                      size_t lda, size_t ldb, size_t ldc) {
    const LtMatmulKey key{m, n, k, batch, int64_t(lda), int64_t(ldb), int64_t(ldc),
                          0, 0, 0, LtBatchKind::PointerArray, h.workspaceSizeInBytes};
    const auto it = h.plan_cache.find(key);
    if (it != h.plan_cache.end()) return it->second.valid;
    return get_or_create_plan<BACKEND>(h, m, n, k, batch, lda, ldb, ldc, 0, 0, 0,
                                       LtBatchKind::PointerArray)
        .valid;
}

template <Backend BACKEND>
inline void launch_lt_plan(
    const cudaStream_t stream, Handle_t &h, const LtMatmulPlan &plan,
    const void *alpha, const void *A, const void *B,
    const void *beta, const void *C, void *D) {
    if (!plan.valid || plan.ws > h.workspaceSizeInBytes || (plan.ws && !h.workspace)) {
        report_lt_plan<BACKEND>(h, plan.key, plan, "no executable plan/workspace");
        std::abort();
    }
    const auto status = cublasLtMatmul(h.cublasLt, h.opDesc,
                                       alpha, A, h.Adesc, B, h.Bdesc, beta, C, h.Cdesc, D, h.Cdesc,
                                       &plan.heur.algo, h.workspace, plan.ws, stream);
    if (status != CUBLAS_STATUS_SUCCESS) {
        report_lt_plan<BACKEND>(h, plan.key, plan, "cublasLtMatmul failed", status);
        std::abort();
    }
}

} // namespace gemmul8::common
