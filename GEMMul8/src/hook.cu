#include "../include/gemmul8.hpp"
#include "self_hipify.hpp"
#include <algorithm>
#include <cstdlib>
#include <dlfcn.h>
#include <iostream>
#include <string>
#include <unordered_map>

// NOTE: skip_scalA/B rely on pointer identity; actual A/B contents are not verified.
//       Use GEMMUL8_SKIP_SCALE_* only when A/B data are unchanged.
/*
| Variable                      | Default | Applies to | Description                                                                          |
| :---------------------------- | :------ | :--------- | :----------------------------------------------------------------------------------- |
| `GEMMUL8_NUM_MOD_D`           | `0`     | DGEMM      | Number of moduli (`unsigned num_moduli`) used in DGEMM emulation. Controls accuracy. |
| `GEMMUL8_NUM_MOD_S`           | `0`     | SGEMM      | Number of moduli (`unsigned num_moduli`) used in SGEMM emulation. Controls accuracy. |
| `GEMMUL8_NUM_MOD_Z`           | `0`     | ZGEMM      | Number of moduli (`unsigned num_moduli`) used in ZGEMM emulation. Controls accuracy. |
| `GEMMUL8_NUM_MOD_C`           | `0`     | CGEMM      | Number of moduli (`unsigned num_moduli`) used in CGEMM emulation. Controls accuracy. |
| `GEMMUL8_FASTMODE_D`          | `0`     | DGEMM      | Enables fast mode (`1` = fast mode, `0` = accurate mode).                            |
| `GEMMUL8_FASTMODE_S`          | `0`     | SGEMM      | Enables fast mode (`1` = fast mode, `0` = accurate mode).                            |
| `GEMMUL8_FASTMODE_Z`          | `0`     | ZGEMM      | Enables fast mode (`1` = fast mode, `0` = accurate mode).                            |
| `GEMMUL8_FASTMODE_C`          | `0`     | CGEMM      | Enables fast mode (`1` = fast mode, `0` = accurate mode).                            |
| `GEMMUL8_MAX_M`               | `0`     | all        | Maximum value of `M` used to preallocate workspace memory.                           |
| `GEMMUL8_MAX_N`               | `0`     | all        | Maximum value of `N` used to preallocate workspace memory.                           |
| `GEMMUL8_MAX_K`               | `0`     | all        | Maximum value of `K` used to preallocate workspace memory.                           |
| `GEMMUL8_MAX_NUM_MOD`         | `2`     | all        | Maximum number of moduli used when computing the size of the preallocated workspace. |
| `GEMMUL8_SKIP_SCALE_A`        | `0`     | all        | Enables skipping redundant preprocessing for `A` (`1` = enable, `0` = disable).      |
| `GEMMUL8_SKIP_SCALE_B`        | `0`     | all        | Enables skipping redundant preprocessing for `B` (`1` = enable, `0` = disable).      |
| `GEMMUL8_USE_EXTRA_WORKSPACE` | `1`     | all        | Enables extra workspace for intermediate buffers (`1` = enable, `0` = disable).      |
*/

namespace {

// ---- Default initialization values ----
namespace initval {

inline constexpr size_t MAX_M         = 0u;    // default M size
inline constexpr size_t MAX_N         = 0u;    // default N size
inline constexpr size_t MAX_K         = 0u;    // default K size
inline constexpr unsigned MAX_NUM_MOD = 0u;    // default modulus count
inline constexpr unsigned NUM_MOD_D   = 0u;    // default double moduli
inline constexpr unsigned NUM_MOD_S   = 0u;    // default float moduli
inline constexpr unsigned NUM_MOD_Z   = 0u;    // default double-complex moduli
inline constexpr unsigned NUM_MOD_C   = 0u;    // default float-complex moduli
inline constexpr bool FASTMODE_D      = false; // default double fastmode
inline constexpr bool FASTMODE_S      = false; // default float fastmode
inline constexpr bool FASTMODE_Z      = false; // default double-complex fastmode
inline constexpr bool FASTMODE_C      = false; // default float-complex fastmode
inline constexpr bool SCALE_A         = false; // default skip_scalA_switch
inline constexpr bool SCALE_B         = false; // default skip_scalB_switch
inline constexpr bool EXTRA_WORKSPACE = true;  // default UseExtraWorkspace

} // namespace initval

// ---- Structure holding GEMM configuration ----
struct Info_t {
    unsigned num_moduli    = 2;
    cublasOperation_t op_A = CUBLAS_OP_N;
    cublasOperation_t op_B = CUBLAS_OP_N;
    size_t m               = 0;
    size_t n               = 0;
    size_t k               = 0;
    size_t lda             = 0;
    size_t ldb             = 0;
    const void *A          = nullptr;
    const void *B          = nullptr;
    int8_t *workA          = nullptr;
    int8_t *workB          = nullptr;
    char Type              = 'N';
    bool fastmode          = false;
};

// ---- Global caches (per-handle) ----
static std::unordered_map<cublasHandle_t, Info_t> last_info;  // previous GEMM info
static std::unordered_map<cublasHandle_t, void *> work_cache; // workspace memory cache
static std::unordered_map<cublasHandle_t, size_t> work_size;  // cached size per handle
static size_t max_workSize    = 0;                            // global workspace limit in byte
static size_t max_workSizeA   = 0;                            // global workspace limit in byte
static size_t max_workSizeB   = 0;                            // global workspace limit in byte
static bool skip_scalA_switch = false;                        // env flag: skip A scaling
static bool skip_scalB_switch = false;                        // env flag: skip B scaling

// ---- Initialize maximum workspace size ----
static void init_max_workspace() {
    if (max_workSize != 0) return; // already initialized

    size_t max_m        = initval::MAX_M;
    size_t max_n        = initval::MAX_N;
    size_t max_k        = initval::MAX_K;
    unsigned max_moduli = initval::MAX_NUM_MOD;

    // Read environment overrides (optional)
    const char *sm   = getenv("GEMMUL8_MAX_M");
    const char *sn   = getenv("GEMMUL8_MAX_N");
    const char *sk   = getenv("GEMMUL8_MAX_K");
    const char *smod = getenv("GEMMUL8_MAX_NUM_MOD");

    // Parse env values if present
    if (sm) {
        try {
            max_m = std::stoull(sm);
        } catch (...) {}
    }
    if (sn) {
        try {
            max_n = std::stoull(sn);
        } catch (...) {}
    }
    if (sk) {
        try {
            max_k = std::stoull(sk);
        } catch (...) {}
    }
    if (smod) {
        try {
            max_moduli = std::stoul(smod);
        } catch (...) {}
    }

    const char *nmz = getenv("GEMMUL8_NUM_MOD_Z");
    const char *nmc = getenv("GEMMUL8_NUM_MOD_C");
    const char *uew = getenv("GEMMUL8_USE_EXTRA_WORKSPACE");
    unsigned NMOD_Z = 0u;
    if (nmz) {
        try {
            NMOD_Z = std::stoul(nmz);
        } catch (...) {}
    }
    unsigned NMOD_C = 0u;
    if (nmc) {
        try {
            NMOD_C = std::stoul(nmc);
        } catch (...) {}
    }
    bool USE_EXTRA_WORKSPACE = (uew ? std::string(uew) == std::string("1") : initval::EXTRA_WORKSPACE);
    if (USE_EXTRA_WORKSPACE) {
        if (NMOD_Z > 0 || NMOD_C > 0)
            max_workSize = gemmul8::workSize<true, true>(max_m, max_n, max_k, max_moduli, true, true, &max_workSizeA, &max_workSizeB);
        else
            max_workSize = gemmul8::workSize<false, true>(max_m, max_n, max_k, max_moduli, true, true, &max_workSizeA, &max_workSizeB);
    } else {
        if (NMOD_Z > 0 || NMOD_C > 0)
            max_workSize = gemmul8::workSize<true, false>(max_m, max_n, max_k, max_moduli, true, true, &max_workSizeA, &max_workSizeB);
        else
            max_workSize = gemmul8::workSize<false, false>(max_m, max_n, max_k, max_moduli, true, true, &max_workSizeA, &max_workSizeB);
    }
}

static inline void get_env_d(unsigned &num_moduli, bool &fastmode, bool &extra_workspace) {
    const char *skipA = getenv("GEMMUL8_SKIP_SCALE_A");
    const char *skipB = getenv("GEMMUL8_SKIP_SCALE_B");
    const char *nm    = getenv("GEMMUL8_NUM_MOD_D");
    const char *fm    = getenv("GEMMUL8_FASTMODE_D");
    const char *uew   = getenv("GEMMUL8_USE_EXTRA_WORKSPACE");
    num_moduli        = initval::NUM_MOD_D;
    if (nm) {
        try {
            num_moduli = std::stoul(nm);
        } catch (...) {}
    }
    fastmode          = (fm ? std::string(fm) == std::string("1") : initval::FASTMODE_D);
    extra_workspace   = (uew ? std::string(uew) == std::string("1") : initval::EXTRA_WORKSPACE);
    skip_scalA_switch = (skipA ? std::string(skipA) == std::string("1") : initval::SCALE_A);
    skip_scalB_switch = (skipB ? std::string(skipB) == std::string("1") : initval::SCALE_B);
}

static inline void get_env_s(unsigned &num_moduli, bool &fastmode, bool &extra_workspace) {
    const char *skipA = getenv("GEMMUL8_SKIP_SCALE_A");
    const char *skipB = getenv("GEMMUL8_SKIP_SCALE_B");
    const char *nm    = getenv("GEMMUL8_NUM_MOD_S");
    const char *fm    = getenv("GEMMUL8_FASTMODE_S");
    const char *uew   = getenv("GEMMUL8_USE_EXTRA_WORKSPACE");
    num_moduli        = initval::NUM_MOD_S;
    if (nm) {
        try {
            num_moduli = std::stoul(nm);
        } catch (...) {}
    }
    fastmode          = (fm ? std::string(fm) == std::string("1") : initval::FASTMODE_S);
    extra_workspace   = (uew ? std::string(uew) == std::string("1") : initval::EXTRA_WORKSPACE);
    skip_scalA_switch = (skipA ? std::string(skipA) == std::string("1") : initval::SCALE_A);
    skip_scalB_switch = (skipB ? std::string(skipB) == std::string("1") : initval::SCALE_B);
}

static inline void get_env_z(unsigned &num_moduli, bool &fastmode, bool &extra_workspace) {
    const char *skipA = getenv("GEMMUL8_SKIP_SCALE_A");
    const char *skipB = getenv("GEMMUL8_SKIP_SCALE_B");
    const char *nm    = getenv("GEMMUL8_NUM_MOD_Z");
    const char *fm    = getenv("GEMMUL8_FASTMODE_Z");
    const char *uew   = getenv("GEMMUL8_USE_EXTRA_WORKSPACE");
    num_moduli        = initval::NUM_MOD_Z;
    if (nm) {
        try {
            num_moduli = std::stoul(nm);
        } catch (...) {}
    }
    fastmode          = (fm ? std::string(fm) == std::string("1") : initval::FASTMODE_Z);
    extra_workspace   = (uew ? std::string(uew) == std::string("1") : initval::EXTRA_WORKSPACE);
    skip_scalA_switch = (skipA ? std::string(skipA) == std::string("1") : initval::SCALE_A);
    skip_scalB_switch = (skipB ? std::string(skipB) == std::string("1") : initval::SCALE_B);
}

static inline void get_env_c(unsigned &num_moduli, bool &fastmode, bool &extra_workspace) {
    const char *skipA = getenv("GEMMUL8_SKIP_SCALE_A");
    const char *skipB = getenv("GEMMUL8_SKIP_SCALE_B");
    const char *nm    = getenv("GEMMUL8_NUM_MOD_C");
    const char *fm    = getenv("GEMMUL8_FASTMODE_C");
    const char *uew   = getenv("GEMMUL8_USE_EXTRA_WORKSPACE");
    num_moduli        = initval::NUM_MOD_C;
    if (nm) {
        try {
            num_moduli = std::stoul(nm);
        } catch (...) {}
    }
    fastmode          = (fm ? std::string(fm) == std::string("1") : initval::FASTMODE_C);
    extra_workspace   = (uew ? std::string(uew) == std::string("1") : initval::EXTRA_WORKSPACE);
    skip_scalA_switch = (skipA ? std::string(skipA) == std::string("1") : initval::SCALE_A);
    skip_scalB_switch = (skipB ? std::string(skipB) == std::string("1") : initval::SCALE_B);
}

static cublasStatus_t get_work(cublasHandle_t handle, size_t req_size, void **ptr) {
    if (!ptr) return CUBLAS_STATUS_INVALID_VALUE;
    init_max_workspace();

    if (req_size == 0) {
        *ptr = nullptr;
        return CUBLAS_STATUS_SUCCESS;
    }

    void *&current_cache = work_cache[handle];
    size_t &current_size = work_size[handle];
    req_size             = max(req_size, max_workSize);

    if (current_cache && current_size >= req_size) {
        *ptr = current_cache;
        return CUBLAS_STATUS_SUCCESS;
    }

    if (current_cache) {
        cudaError_t free_err = cudaFree(current_cache);
        if (free_err != cudaSuccess) {
            std::cerr << "[GEMMUL8 HOOK] Warning: cudaFree failed ("
                      << cudaGetErrorString(free_err) << ")" << std::endl;
        }
        current_cache = nullptr;
        current_size  = 0;
    }

    cudaError_t err = cudaMalloc(&current_cache, req_size);
    if (err != cudaSuccess) {
        current_cache = nullptr;
        current_size  = 0;
        *ptr          = nullptr;
        std::cerr << "[GEMMUL8 HOOK] Malloc failed for size " << req_size
                  << " bytes. Error: " << cudaGetErrorString(err) << std::endl;
        return CUBLAS_STATUS_ALLOC_FAILED;
    }

    current_size = req_size;
    *ptr         = current_cache;
    return CUBLAS_STATUS_SUCCESS;
}

static void cleanup_work(cublasHandle_t handle) {
    auto it = work_cache.find(handle);
    if (it != work_cache.end()) {
        cudaError_t free_err = cudaFree(it->second);
        if (free_err != cudaSuccess) {
            std::cerr << "[GEMMUL8 HOOK] cublasDestroy: Warning: cudaFree failed ("
                      << cudaGetErrorString(free_err) << ")" << std::endl;
        }
        work_cache.erase(it);
        work_size.erase(handle);
        last_info.erase(handle);
    }
}
} // namespace

// =======================
// Hook: cublasDestroy
// =======================
extern "C" cublasStatus_t cublasDestroy_v2(cublasHandle_t handle) {
#ifdef __CUDA_ARCH__
    return CUBLAS_STATUS_NOT_SUPPORTED;
#else
    cleanup_work(handle);
    using cublasDestroy_t                     = cublasStatus_t (*)(cublasHandle_t);
    static cublasDestroy_t real_cublasDestroy = (cublasDestroy_t)dlsym(RTLD_NEXT, STR(cublasDestroy_v2));
    if (!real_cublasDestroy) return CUBLAS_STATUS_NOT_INITIALIZED;
    return real_cublasDestroy(handle);
#endif
}

// =======================
// Hook: cublasSgemm_v2
// =======================
extern "C" cublasStatus_t cublasSgemm_v2(
    cublasHandle_t handle,
    cublasOperation_t transa, cublasOperation_t transb,
    int m, int n, int k,
    const float *alpha,
    const float *A, int lda,
    const float *B, int ldb,
    const float *beta,
    float *C, int ldc //
) {
#ifdef __CUDA_ARCH__
    return CUBLAS_STATUS_NOT_SUPPORTED;

#else
    if (m <= 0 || n <= 0 || k <= 0) return CUBLAS_STATUS_SUCCESS;
    if (!A || !B || !C) return CUBLAS_STATUS_INVALID_VALUE;

    unsigned num_moduli;
    bool fastmode;
    bool extra_workspace;
    get_env_s(num_moduli, fastmode, extra_workspace);

    if (num_moduli < 2u || 20u < num_moduli) {
        using gemm_t            = cublasStatus_t (*)(cublasHandle_t, cublasOperation_t, cublasOperation_t, int, int, int, const float *, const float *, int, const float *, int, const float *, float *, int);
        static gemm_t real_gemm = (gemm_t)dlsym(RTLD_NEXT, STR(cublasSgemm_v2));
        if (!real_gemm) return CUBLAS_STATUS_NOT_INITIALIZED;
        return real_gemm(handle, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
    }

    size_t wsizeA, wsizeB;
    size_t wsize;
    if (extra_workspace)
        wsize = gemmul8::workSize<false, true>(m, n, k, num_moduli, skip_scalA_switch, skip_scalB_switch, &wsizeA, &wsizeB);
    else
        wsize = gemmul8::workSize<false, false>(m, n, k, num_moduli, skip_scalA_switch, skip_scalB_switch, &wsizeA, &wsizeB);
    void *work            = nullptr;
    cublasStatus_t status = get_work(handle, wsize, &work);
    if (status != CUBLAS_STATUS_SUCCESS) return status;

    int8_t *workA        = reinterpret_cast<int8_t *>(work);
    const size_t offsetA = max(max_workSizeA, wsizeA);
    const size_t offsetB = max(max_workSizeB, wsizeB);
    int8_t *workB        = workA + offsetA;
    int8_t *workC        = workB + offsetB;

    bool skip_scalA  = false; // false (unskip scaling_A) or true (skip scaling_A)
    bool skip_scalB  = false; // false (unskip scaling_B) or true (skip scaling_B)
    Info_t &info_pre = last_info[handle];
    if (info_pre.num_moduli == num_moduli && info_pre.k == k && info_pre.Type == 'S' && info_pre.fastmode == fastmode) {
        if (skip_scalA_switch && info_pre.workA == workA && info_pre.A == A && info_pre.m == m && info_pre.lda == lda && info_pre.op_A == transa) {
            skip_scalA = skip_scalA_switch && true;
        }
        if (skip_scalB_switch && info_pre.workB == workB && info_pre.B == B && info_pre.n == n && info_pre.ldb == ldb && info_pre.op_B == transb) {
            skip_scalB = skip_scalB_switch && true;
        }
    }

    if (extra_workspace) {
        gemmul8::gemm<float, true>(handle,
                                   transa, transb, m, n, k,
                                   alpha, A, lda, B, ldb,
                                   beta, C, ldc,
                                   num_moduli, fastmode,
                                   reinterpret_cast<void *>(workC),
                                   reinterpret_cast<void *>(workA),
                                   reinterpret_cast<void *>(workB),
                                   skip_scalA_switch, skip_scalB_switch,
                                   skip_scalA, skip_scalB);
    } else {
        gemmul8::gemm<float, false>(handle,
                                    transa, transb, m, n, k,
                                    alpha, A, lda, B, ldb,
                                    beta, C, ldc,
                                    num_moduli, fastmode,
                                    reinterpret_cast<void *>(workC),
                                    reinterpret_cast<void *>(workA),
                                    reinterpret_cast<void *>(workB),
                                    skip_scalA_switch, skip_scalB_switch,
                                    skip_scalA, skip_scalB);
    }

    info_pre.num_moduli = num_moduli;
    info_pre.op_A       = transa;
    info_pre.op_B       = transb;
    info_pre.m          = m;
    info_pre.n          = n;
    info_pre.k          = k;
    info_pre.lda        = lda;
    info_pre.ldb        = ldb;
    info_pre.A          = A;
    info_pre.B          = B;
    info_pre.workA      = workA;
    info_pre.workB      = workB;
    info_pre.Type       = 'S';
    info_pre.fastmode   = fastmode;

    return CUBLAS_STATUS_SUCCESS;

#endif
}

// =======================
// Hook: cublasDgemm_v2
// =======================
extern "C" cublasStatus_t cublasDgemm_v2(
    cublasHandle_t handle,
    cublasOperation_t transa, cublasOperation_t transb,
    int m, int n, int k,
    const double *alpha,
    const double *A, int lda,
    const double *B, int ldb,
    const double *beta,
    double *C, int ldc //
) {
#ifdef __CUDA_ARCH__
    return CUBLAS_STATUS_NOT_SUPPORTED;

#else
    if (m <= 0 || n <= 0 || k <= 0) return CUBLAS_STATUS_SUCCESS;
    if (!A || !B || !C) return CUBLAS_STATUS_INVALID_VALUE;

    unsigned num_moduli;
    bool fastmode;
    bool extra_workspace;
    get_env_d(num_moduli, fastmode, extra_workspace);

    if (num_moduli < 2u || 20u < num_moduli) {
        using gemm_t            = cublasStatus_t (*)(cublasHandle_t, cublasOperation_t, cublasOperation_t, int, int, int, const double *, const double *, int, const double *, int, const double *, double *, int);
        static gemm_t real_gemm = (gemm_t)dlsym(RTLD_NEXT, STR(cublasDgemm_v2));
        if (!real_gemm) return CUBLAS_STATUS_NOT_INITIALIZED;
        return real_gemm(handle, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
    }

    size_t wsizeA, wsizeB;
    size_t wsize;
    if (extra_workspace)
        wsize = gemmul8::workSize<false, true>(m, n, k, num_moduli, skip_scalA_switch, skip_scalB_switch, &wsizeA, &wsizeB);
    else
        wsize = gemmul8::workSize<false, false>(m, n, k, num_moduli, skip_scalA_switch, skip_scalB_switch, &wsizeA, &wsizeB);
    void *work            = nullptr;
    cublasStatus_t status = get_work(handle, wsize, &work);
    if (status != CUBLAS_STATUS_SUCCESS) return status;

    int8_t *workA        = reinterpret_cast<int8_t *>(work);
    const size_t offsetA = max(max_workSizeA, wsizeA);
    const size_t offsetB = max(max_workSizeB, wsizeB);
    int8_t *workB        = workA + offsetA;
    int8_t *workC        = workB + offsetB;

    bool skip_scalA  = false; // false (unskip scaling_A) or true (skip scaling_A)
    bool skip_scalB  = false; // false (unskip scaling_B) or true (skip scaling_B)
    Info_t &info_pre = last_info[handle];
    if (info_pre.num_moduli == num_moduli && info_pre.k == k && info_pre.Type == 'D' && info_pre.fastmode == fastmode) {
        if (skip_scalA_switch && info_pre.workA == workA && info_pre.A == A && info_pre.m == m && info_pre.lda == lda && info_pre.op_A == transa) {
            skip_scalA = true;
        }
        if (skip_scalB_switch && info_pre.workB == workB && info_pre.B == B && info_pre.n == n && info_pre.ldb == ldb && info_pre.op_B == transb) {
            skip_scalB = true;
        }
    }

    if (extra_workspace) {
        gemmul8::gemm<double, true>(handle,
                                    transa, transb, m, n, k,
                                    alpha, A, lda, B, ldb,
                                    beta, C, ldc,
                                    num_moduli, fastmode,
                                    reinterpret_cast<void *>(workC),
                                    reinterpret_cast<void *>(workA),
                                    reinterpret_cast<void *>(workB),
                                    skip_scalA_switch, skip_scalB_switch,
                                    skip_scalA, skip_scalB);
    } else {
        gemmul8::gemm<double, false>(handle,
                                     transa, transb, m, n, k,
                                     alpha, A, lda, B, ldb,
                                     beta, C, ldc,
                                     num_moduli, fastmode,
                                     reinterpret_cast<void *>(workC),
                                     reinterpret_cast<void *>(workA),
                                     reinterpret_cast<void *>(workB),
                                     skip_scalA_switch, skip_scalB_switch,
                                     skip_scalA, skip_scalB);
    }

    info_pre.num_moduli = num_moduli;
    info_pre.op_A       = transa;
    info_pre.op_B       = transb;
    info_pre.m          = m;
    info_pre.n          = n;
    info_pre.k          = k;
    info_pre.lda        = lda;
    info_pre.ldb        = ldb;
    info_pre.A          = A;
    info_pre.B          = B;
    info_pre.workA      = workA;
    info_pre.workB      = workB;
    info_pre.Type       = 'D';
    info_pre.fastmode   = fastmode;

    return CUBLAS_STATUS_SUCCESS;

#endif
}

// =======================
// Hook: cublasCgemm_v2
// =======================
extern "C" cublasStatus_t cublasCgemm_v2(
    cublasHandle_t handle,
    cublasOperation_t transa, cublasOperation_t transb,
    int m, int n, int k,
    const cuComplex *alpha,
    const cuComplex *A, int lda,
    const cuComplex *B, int ldb,
    const cuComplex *beta,
    cuComplex *C, int ldc //
) {
#ifdef __CUDA_ARCH__
    return CUBLAS_STATUS_NOT_SUPPORTED;

#else
    if (m <= 0 || n <= 0 || k <= 0) return CUBLAS_STATUS_SUCCESS;
    if (!A || !B || !C) return CUBLAS_STATUS_INVALID_VALUE;

    unsigned num_moduli;
    bool fastmode;
    bool extra_workspace;
    get_env_c(num_moduli, fastmode, extra_workspace);

    if (num_moduli < 2u || 20u < num_moduli) {
        using gemm_t            = cublasStatus_t (*)(cublasHandle_t, cublasOperation_t, cublasOperation_t, int, int, int, const cuComplex *, const cuComplex *, int, const cuComplex *, int, const cuComplex *, cuComplex *, int);
        static gemm_t real_gemm = (gemm_t)dlsym(RTLD_NEXT, STR(cublasCgemm_v2));
        if (!real_gemm) return CUBLAS_STATUS_NOT_INITIALIZED;
        return real_gemm(handle, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
    }

    size_t wsizeA, wsizeB;
    size_t wsize;
    if (extra_workspace)
        wsize = gemmul8::workSize<true, true>(m, n, k, num_moduli, skip_scalA_switch, skip_scalB_switch, &wsizeA, &wsizeB);
    else
        wsize = gemmul8::workSize<true, false>(m, n, k, num_moduli, skip_scalA_switch, skip_scalB_switch, &wsizeA, &wsizeB);
    void *work            = nullptr;
    cublasStatus_t status = get_work(handle, wsize, &work);
    if (status != CUBLAS_STATUS_SUCCESS) return status;

    int8_t *workA        = reinterpret_cast<int8_t *>(work);
    const size_t offsetA = max(max_workSizeA, wsizeA);
    const size_t offsetB = max(max_workSizeB, wsizeB);
    int8_t *workB        = workA + offsetA;
    int8_t *workC        = workB + offsetB;

    bool skip_scalA  = false; // false (unskip scaling_A) or true (skip scaling_A)
    bool skip_scalB  = false; // false (unskip scaling_B) or true (skip scaling_B)
    Info_t &info_pre = last_info[handle];
    if (info_pre.num_moduli == num_moduli && info_pre.k == k && info_pre.Type == 'C' && info_pre.fastmode == fastmode) {
        if (skip_scalA_switch && info_pre.workA == workA && info_pre.A == A && info_pre.m == m && info_pre.lda == lda && info_pre.op_A == transa) {
            skip_scalA = skip_scalA_switch && true;
        }
        if (skip_scalB_switch && info_pre.workB == workB && info_pre.B == B && info_pre.n == n && info_pre.ldb == ldb && info_pre.op_B == transb) {
            skip_scalB = skip_scalB_switch && true;
        }
    }

    if (extra_workspace) {
        gemmul8::gemm<cuComplex, true>(handle,
                                       transa, transb, m, n, k,
                                       alpha, A, lda, B, ldb,
                                       beta, C, ldc,
                                       num_moduli, fastmode,
                                       reinterpret_cast<void *>(workC),
                                       reinterpret_cast<void *>(workA),
                                       reinterpret_cast<void *>(workB),
                                       skip_scalA_switch, skip_scalB_switch,
                                       skip_scalA, skip_scalB);
    } else {
        gemmul8::gemm<cuComplex, false>(handle,
                                        transa, transb, m, n, k,
                                        alpha, A, lda, B, ldb,
                                        beta, C, ldc,
                                        num_moduli, fastmode,
                                        reinterpret_cast<void *>(workC),
                                        reinterpret_cast<void *>(workA),
                                        reinterpret_cast<void *>(workB),
                                        skip_scalA_switch, skip_scalB_switch,
                                        skip_scalA, skip_scalB);
    }

    info_pre.num_moduli = num_moduli;
    info_pre.op_A       = transa;
    info_pre.op_B       = transb;
    info_pre.m          = m;
    info_pre.n          = n;
    info_pre.k          = k;
    info_pre.lda        = lda;
    info_pre.ldb        = ldb;
    info_pre.A          = A;
    info_pre.B          = B;
    info_pre.workA      = workA;
    info_pre.workB      = workB;
    info_pre.Type       = 'C';
    info_pre.fastmode   = fastmode;

    return CUBLAS_STATUS_SUCCESS;

#endif
}

// =======================
// Hook: cublasZgemm_v2
// =======================
extern "C" cublasStatus_t cublasZgemm_v2(
    cublasHandle_t handle,
    cublasOperation_t transa, cublasOperation_t transb,
    int m, int n, int k,
    const cuDoubleComplex *alpha,
    const cuDoubleComplex *A, int lda,
    const cuDoubleComplex *B, int ldb,
    const cuDoubleComplex *beta,
    cuDoubleComplex *C, int ldc //
) {
#ifdef __CUDA_ARCH__
    return CUBLAS_STATUS_NOT_SUPPORTED;

#else
    if (m <= 0 || n <= 0 || k <= 0) return CUBLAS_STATUS_SUCCESS;
    if (!A || !B || !C) return CUBLAS_STATUS_INVALID_VALUE;

    unsigned num_moduli;
    bool fastmode;
    bool extra_workspace;
    get_env_z(num_moduli, fastmode, extra_workspace);

    if (num_moduli < 2u || 20u < num_moduli) {
        using gemm_t            = cublasStatus_t (*)(cublasHandle_t, cublasOperation_t, cublasOperation_t, int, int, int, const cuDoubleComplex *, const cuDoubleComplex *, int, const cuDoubleComplex *, int, const cuDoubleComplex *, cuDoubleComplex *, int);
        static gemm_t real_gemm = (gemm_t)dlsym(RTLD_NEXT, STR(cublasZgemm_v2));
        if (!real_gemm) return CUBLAS_STATUS_NOT_INITIALIZED;
        return real_gemm(handle, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
    }

    size_t wsizeA, wsizeB;
    size_t wsize;
    if (extra_workspace)
        wsize = gemmul8::workSize<true, true>(m, n, k, num_moduli, skip_scalA_switch, skip_scalB_switch, &wsizeA, &wsizeB);
    else
        wsize = gemmul8::workSize<true, false>(m, n, k, num_moduli, skip_scalA_switch, skip_scalB_switch, &wsizeA, &wsizeB);
    void *work            = nullptr;
    cublasStatus_t status = get_work(handle, wsize, &work);
    if (status != CUBLAS_STATUS_SUCCESS) return status;

    int8_t *workA        = reinterpret_cast<int8_t *>(work);
    const size_t offsetA = max(max_workSizeA, wsizeA);
    const size_t offsetB = max(max_workSizeB, wsizeB);
    int8_t *workB        = workA + offsetA;
    int8_t *workC        = workB + offsetB;

    bool skip_scalA  = false; // false (unskip scaling_A) or true (skip scaling_A)
    bool skip_scalB  = false; // false (unskip scaling_B) or true (skip scaling_B)
    Info_t &info_pre = last_info[handle];
    if (info_pre.num_moduli == num_moduli && info_pre.k == k && info_pre.Type == 'Z' && info_pre.fastmode == fastmode) {
        if (skip_scalA_switch && info_pre.workA == workA && info_pre.A == A && info_pre.m == m && info_pre.lda == lda && info_pre.op_A == transa) {
            skip_scalA = true;
        }
        if (skip_scalB_switch && info_pre.workB == workB && info_pre.B == B && info_pre.n == n && info_pre.ldb == ldb && info_pre.op_B == transb) {
            skip_scalB = true;
        }
    }

    if (extra_workspace) {
        gemmul8::gemm<cuDoubleComplex, true>(handle,
                                             transa, transb, m, n, k,
                                             alpha, A, lda, B, ldb,
                                             beta, C, ldc,
                                             num_moduli, fastmode,
                                             reinterpret_cast<void *>(workC),
                                             reinterpret_cast<void *>(workA),
                                             reinterpret_cast<void *>(workB),
                                             skip_scalA_switch, skip_scalB_switch,
                                             skip_scalA, skip_scalB);
    } else {
        gemmul8::gemm<cuDoubleComplex, false>(handle,
                                              transa, transb, m, n, k,
                                              alpha, A, lda, B, ldb,
                                              beta, C, ldc,
                                              num_moduli, fastmode,
                                              reinterpret_cast<void *>(workC),
                                              reinterpret_cast<void *>(workA),
                                              reinterpret_cast<void *>(workB),
                                              skip_scalA_switch, skip_scalB_switch,
                                              skip_scalA, skip_scalB);
    }

    info_pre.num_moduli = num_moduli;
    info_pre.op_A       = transa;
    info_pre.op_B       = transb;
    info_pre.m          = m;
    info_pre.n          = n;
    info_pre.k          = k;
    info_pre.lda        = lda;
    info_pre.ldb        = ldb;
    info_pre.A          = A;
    info_pre.B          = B;
    info_pre.workA      = workA;
    info_pre.workB      = workB;
    info_pre.Type       = 'Z';
    info_pre.fastmode   = fastmode;

    return CUBLAS_STATUS_SUCCESS;

#endif
}

// =======================
// Hook: cublasGemmEx
// =======================
extern "C" cublasStatus_t cublasGemmEx(
    cublasHandle_t handle,
    cublasOperation_t transa, cublasOperation_t transb,
    int m, int n, int k,
    const void *alpha,
    const void *A, cudaDataType Atype, int lda,
    const void *B, cudaDataType Btype, int ldb,
    const void *beta,
    void *C, cudaDataType Ctype, int ldc,
    cublasComputeType_t computeType, cublasGemmAlgo_t algo //
) {
#ifdef __CUDA_ARCH__
    return CUBLAS_STATUS_NOT_SUPPORTED;

#else
    if (m <= 0 || n <= 0 || k <= 0) return CUBLAS_STATUS_SUCCESS;
    if (!A || !B || !C) return CUBLAS_STATUS_INVALID_VALUE;

    // SGEMM
    if (computeType == CUBLAS_COMPUTE_32F &&
        Atype == CUDA_R_32F &&
        Btype == CUDA_R_32F &&
        Ctype == CUDA_R_32F) {
        unsigned num_moduli;
        bool fastmode;
        bool extra_workspace;
        get_env_s(num_moduli, fastmode, extra_workspace);

        if (2u <= num_moduli && num_moduli <= 20u) {
            size_t wsizeA, wsizeB;
            size_t wsize;
            if (extra_workspace)
                wsize = gemmul8::workSize<false, true>(m, n, k, num_moduli, skip_scalA_switch, skip_scalB_switch, &wsizeA, &wsizeB);
            else
                wsize = gemmul8::workSize<false, false>(m, n, k, num_moduli, skip_scalA_switch, skip_scalB_switch, &wsizeA, &wsizeB);
            void *work            = nullptr;
            cublasStatus_t status = get_work(handle, wsize, &work);
            if (status != CUBLAS_STATUS_SUCCESS) return status;

            int8_t *workA        = reinterpret_cast<int8_t *>(work);
            const size_t offsetA = max(max_workSizeA, wsizeA);
            const size_t offsetB = max(max_workSizeB, wsizeB);
            int8_t *workB        = workA + offsetA;
            int8_t *workC        = workB + offsetB;

            bool skip_scalA  = false; // false (unskip scaling_A) or true (skip scaling_A)
            bool skip_scalB  = false; // false (unskip scaling_B) or true (skip scaling_B)
            Info_t &info_pre = last_info[handle];
            if (info_pre.num_moduli == num_moduli && info_pre.k == k && info_pre.Type == 'S' && info_pre.fastmode == fastmode) {
                if (skip_scalA_switch && info_pre.workA == workA && info_pre.A == A && info_pre.m == m && info_pre.lda == lda && info_pre.op_A == transa) {
                    skip_scalA = true;
                }
                if (skip_scalB_switch && info_pre.workB == workB && info_pre.B == B && info_pre.n == n && info_pre.ldb == ldb && info_pre.op_B == transb) {
                    skip_scalB = true;
                }
            }

            if (extra_workspace) {
                gemmul8::gemm<float, true>(handle,
                                           transa, transb, m, n, k,
                                           reinterpret_cast<const float *>(alpha),
                                           reinterpret_cast<const float *>(A), lda,
                                           reinterpret_cast<const float *>(B), ldb,
                                           reinterpret_cast<const float *>(beta),
                                           reinterpret_cast<float *>(C), ldc,
                                           num_moduli, fastmode,
                                           reinterpret_cast<void *>(workC),
                                           reinterpret_cast<void *>(workA),
                                           reinterpret_cast<void *>(workB),
                                           skip_scalA_switch, skip_scalB_switch,
                                           skip_scalA, skip_scalB);
            } else {
                gemmul8::gemm<float, false>(handle,
                                            transa, transb, m, n, k,
                                            reinterpret_cast<const float *>(alpha),
                                            reinterpret_cast<const float *>(A), lda,
                                            reinterpret_cast<const float *>(B), ldb,
                                            reinterpret_cast<const float *>(beta),
                                            reinterpret_cast<float *>(C), ldc,
                                            num_moduli, fastmode,
                                            reinterpret_cast<void *>(workC),
                                            reinterpret_cast<void *>(workA),
                                            reinterpret_cast<void *>(workB),
                                            skip_scalA_switch, skip_scalB_switch,
                                            skip_scalA, skip_scalB);
            }

            info_pre.num_moduli = num_moduli;
            info_pre.op_A       = transa;
            info_pre.op_B       = transb;
            info_pre.m          = m;
            info_pre.n          = n;
            info_pre.k          = k;
            info_pre.lda        = lda;
            info_pre.ldb        = ldb;
            info_pre.A          = A;
            info_pre.B          = B;
            info_pre.workA      = workA;
            info_pre.workB      = workB;
            info_pre.Type       = 'S';
            info_pre.fastmode   = fastmode;

            return CUBLAS_STATUS_SUCCESS;
        }
    }

    // DGEMM
    if (computeType == CUBLAS_COMPUTE_64F &&
        Atype == CUDA_R_64F &&
        Btype == CUDA_R_64F &&
        Ctype == CUDA_R_64F) {
        unsigned num_moduli;
        bool fastmode;
        bool extra_workspace;
        get_env_d(num_moduli, fastmode, extra_workspace);

        if (2u <= num_moduli && num_moduli <= 20u) {
            size_t wsizeA, wsizeB;
            size_t wsize;
            if (extra_workspace)
                wsize = gemmul8::workSize<false, true>(m, n, k, num_moduli, skip_scalA_switch, skip_scalB_switch, &wsizeA, &wsizeB);
            else
                wsize = gemmul8::workSize<false, false>(m, n, k, num_moduli, skip_scalA_switch, skip_scalB_switch, &wsizeA, &wsizeB);
            void *work            = nullptr;
            cublasStatus_t status = get_work(handle, wsize, &work);
            if (status != CUBLAS_STATUS_SUCCESS) return status;

            int8_t *workA        = reinterpret_cast<int8_t *>(work);
            const size_t offsetA = max(max_workSizeA, wsizeA);
            const size_t offsetB = max(max_workSizeB, wsizeB);
            int8_t *workB        = workA + offsetA;
            int8_t *workC        = workB + offsetB;

            bool skip_scalA  = false; // false (unskip scaling_A) or true (skip scaling_A)
            bool skip_scalB  = false; // false (unskip scaling_B) or true (skip scaling_B)
            Info_t &info_pre = last_info[handle];
            if (info_pre.num_moduli == num_moduli && info_pre.k == k && info_pre.Type == 'D' && info_pre.fastmode == fastmode) {
                if (skip_scalA_switch && info_pre.workA == workA && info_pre.A == A && info_pre.m == m && info_pre.lda == lda && info_pre.op_A == transa) {
                    skip_scalA = true;
                }
                if (skip_scalB_switch && info_pre.workB == workB && info_pre.B == B && info_pre.n == n && info_pre.ldb == ldb && info_pre.op_B == transb) {
                    skip_scalB = true;
                }
            }

            if (extra_workspace) {
                gemmul8::gemm<double, true>(handle,
                                            transa, transb, m, n, k,
                                            reinterpret_cast<const double *>(alpha),
                                            reinterpret_cast<const double *>(A), lda,
                                            reinterpret_cast<const double *>(B), ldb,
                                            reinterpret_cast<const double *>(beta),
                                            reinterpret_cast<double *>(C), ldc,
                                            num_moduli, fastmode,
                                            reinterpret_cast<void *>(workC),
                                            reinterpret_cast<void *>(workA),
                                            reinterpret_cast<void *>(workB),
                                            skip_scalA_switch, skip_scalB_switch,
                                            skip_scalA, skip_scalB);
            } else {
                gemmul8::gemm<double, false>(handle,
                                             transa, transb, m, n, k,
                                             reinterpret_cast<const double *>(alpha),
                                             reinterpret_cast<const double *>(A), lda,
                                             reinterpret_cast<const double *>(B), ldb,
                                             reinterpret_cast<const double *>(beta),
                                             reinterpret_cast<double *>(C), ldc,
                                             num_moduli, fastmode,
                                             reinterpret_cast<void *>(workC),
                                             reinterpret_cast<void *>(workA),
                                             reinterpret_cast<void *>(workB),
                                             skip_scalA_switch, skip_scalB_switch,
                                             skip_scalA, skip_scalB);
            }

            info_pre.num_moduli = num_moduli;
            info_pre.op_A       = transa;
            info_pre.op_B       = transb;
            info_pre.m          = m;
            info_pre.n          = n;
            info_pre.k          = k;
            info_pre.lda        = lda;
            info_pre.ldb        = ldb;
            info_pre.A          = A;
            info_pre.B          = B;
            info_pre.workA      = workA;
            info_pre.workB      = workB;
            info_pre.Type       = 'D';
            info_pre.fastmode   = fastmode;

            return CUBLAS_STATUS_SUCCESS;
        }
    }

    // CGEMM
    if (computeType == CUBLAS_COMPUTE_32F &&
        Atype == CUDA_C_32F &&
        Btype == CUDA_C_32F &&
        Ctype == CUDA_C_32F) {
        unsigned num_moduli;
        bool fastmode;
        bool extra_workspace;
        get_env_c(num_moduli, fastmode, extra_workspace);

        if (2u <= num_moduli && num_moduli <= 20u) {
            size_t wsizeA, wsizeB;
            size_t wsize;
            if (extra_workspace)
                wsize = gemmul8::workSize<true, true>(m, n, k, num_moduli, skip_scalA_switch, skip_scalB_switch, &wsizeA, &wsizeB);
            else
                wsize = gemmul8::workSize<true, false>(m, n, k, num_moduli, skip_scalA_switch, skip_scalB_switch, &wsizeA, &wsizeB);
            void *work            = nullptr;
            cublasStatus_t status = get_work(handle, wsize, &work);
            if (status != CUBLAS_STATUS_SUCCESS) return status;

            int8_t *workA        = reinterpret_cast<int8_t *>(work);
            const size_t offsetA = max(max_workSizeA, wsizeA);
            const size_t offsetB = max(max_workSizeB, wsizeB);
            int8_t *workB        = workA + offsetA;
            int8_t *workC        = workB + offsetB;

            bool skip_scalA  = false; // false (unskip scaling_A) or true (skip scaling_A)
            bool skip_scalB  = false; // false (unskip scaling_B) or true (skip scaling_B)
            Info_t &info_pre = last_info[handle];
            if (info_pre.num_moduli == num_moduli && info_pre.k == k && info_pre.Type == 'C' && info_pre.fastmode == fastmode) {
                if (skip_scalA_switch && info_pre.workA == workA && info_pre.A == A && info_pre.m == m && info_pre.lda == lda && info_pre.op_A == transa) {
                    skip_scalA = true;
                }
                if (skip_scalB_switch && info_pre.workB == workB && info_pre.B == B && info_pre.n == n && info_pre.ldb == ldb && info_pre.op_B == transb) {
                    skip_scalB = true;
                }
            }

            if (extra_workspace) {
                gemmul8::gemm<cuComplex, true>(handle,
                                               transa, transb, m, n, k,
                                               reinterpret_cast<const cuComplex *>(alpha),
                                               reinterpret_cast<const cuComplex *>(A), lda,
                                               reinterpret_cast<const cuComplex *>(B), ldb,
                                               reinterpret_cast<const cuComplex *>(beta),
                                               reinterpret_cast<cuComplex *>(C), ldc,
                                               num_moduli, fastmode,
                                               reinterpret_cast<void *>(workC),
                                               reinterpret_cast<void *>(workA),
                                               reinterpret_cast<void *>(workB),
                                               skip_scalA_switch, skip_scalB_switch,
                                               skip_scalA, skip_scalB);
            } else {
                gemmul8::gemm<cuComplex, false>(handle,
                                                transa, transb, m, n, k,
                                                reinterpret_cast<const cuComplex *>(alpha),
                                                reinterpret_cast<const cuComplex *>(A), lda,
                                                reinterpret_cast<const cuComplex *>(B), ldb,
                                                reinterpret_cast<const cuComplex *>(beta),
                                                reinterpret_cast<cuComplex *>(C), ldc,
                                                num_moduli, fastmode,
                                                reinterpret_cast<void *>(workC),
                                                reinterpret_cast<void *>(workA),
                                                reinterpret_cast<void *>(workB),
                                                skip_scalA_switch, skip_scalB_switch,
                                                skip_scalA, skip_scalB);
            }

            info_pre.num_moduli = num_moduli;
            info_pre.op_A       = transa;
            info_pre.op_B       = transb;
            info_pre.m          = m;
            info_pre.n          = n;
            info_pre.k          = k;
            info_pre.lda        = lda;
            info_pre.ldb        = ldb;
            info_pre.A          = A;
            info_pre.B          = B;
            info_pre.workA      = workA;
            info_pre.workB      = workB;
            info_pre.Type       = 'C';
            info_pre.fastmode   = fastmode;

            return CUBLAS_STATUS_SUCCESS;
        }
    }

    // ZGEMM
    if (computeType == CUBLAS_COMPUTE_64F &&
        Atype == CUDA_C_64F &&
        Btype == CUDA_C_64F &&
        Ctype == CUDA_C_64F) {
        unsigned num_moduli;
        bool fastmode;
        bool extra_workspace;
        get_env_z(num_moduli, fastmode, extra_workspace);

        if (2u <= num_moduli && num_moduli <= 20u) {
            size_t wsizeA, wsizeB;
            size_t wsize;
            if (extra_workspace)
                wsize = gemmul8::workSize<true, true>(m, n, k, num_moduli, skip_scalA_switch, skip_scalB_switch, &wsizeA, &wsizeB);
            else
                wsize = gemmul8::workSize<true, false>(m, n, k, num_moduli, skip_scalA_switch, skip_scalB_switch, &wsizeA, &wsizeB);
            void *work            = nullptr;
            cublasStatus_t status = get_work(handle, wsize, &work);
            if (status != CUBLAS_STATUS_SUCCESS) return status;

            int8_t *workA        = reinterpret_cast<int8_t *>(work);
            const size_t offsetA = max(max_workSizeA, wsizeA);
            const size_t offsetB = max(max_workSizeB, wsizeB);
            int8_t *workB        = workA + offsetA;
            int8_t *workC        = workB + offsetB;

            bool skip_scalA  = false; // false (unskip scaling_A) or true (skip scaling_A)
            bool skip_scalB  = false; // false (unskip scaling_B) or true (skip scaling_B)
            Info_t &info_pre = last_info[handle];
            if (info_pre.num_moduli == num_moduli && info_pre.k == k && info_pre.Type == 'Z' && info_pre.fastmode == fastmode) {
                if (skip_scalA_switch && info_pre.workA == workA && info_pre.A == A && info_pre.m == m && info_pre.lda == lda && info_pre.op_A == transa) {
                    skip_scalA = true;
                }
                if (skip_scalB_switch && info_pre.workB == workB && info_pre.B == B && info_pre.n == n && info_pre.ldb == ldb && info_pre.op_B == transb) {
                    skip_scalB = true;
                }
            }

            if (extra_workspace) {
                gemmul8::gemm<cuDoubleComplex, true>(handle,
                                                     transa, transb, m, n, k,
                                                     reinterpret_cast<const cuDoubleComplex *>(alpha),
                                                     reinterpret_cast<const cuDoubleComplex *>(A), lda,
                                                     reinterpret_cast<const cuDoubleComplex *>(B), ldb,
                                                     reinterpret_cast<const cuDoubleComplex *>(beta),
                                                     reinterpret_cast<cuDoubleComplex *>(C), ldc,
                                                     num_moduli, fastmode,
                                                     reinterpret_cast<void *>(workC),
                                                     reinterpret_cast<void *>(workA),
                                                     reinterpret_cast<void *>(workB),
                                                     skip_scalA_switch, skip_scalB_switch,
                                                     skip_scalA, skip_scalB);
            } else {
                gemmul8::gemm<cuDoubleComplex, false>(handle,
                                                      transa, transb, m, n, k,
                                                      reinterpret_cast<const cuDoubleComplex *>(alpha),
                                                      reinterpret_cast<const cuDoubleComplex *>(A), lda,
                                                      reinterpret_cast<const cuDoubleComplex *>(B), ldb,
                                                      reinterpret_cast<const cuDoubleComplex *>(beta),
                                                      reinterpret_cast<cuDoubleComplex *>(C), ldc,
                                                      num_moduli, fastmode,
                                                      reinterpret_cast<void *>(workC),
                                                      reinterpret_cast<void *>(workA),
                                                      reinterpret_cast<void *>(workB),
                                                      skip_scalA_switch, skip_scalB_switch,
                                                      skip_scalA, skip_scalB);
            }

            info_pre.num_moduli = num_moduli;
            info_pre.op_A       = transa;
            info_pre.op_B       = transb;
            info_pre.m          = m;
            info_pre.n          = n;
            info_pre.k          = k;
            info_pre.lda        = lda;
            info_pre.ldb        = ldb;
            info_pre.A          = A;
            info_pre.B          = B;
            info_pre.workA      = workA;
            info_pre.workB      = workB;
            info_pre.Type       = 'Z';
            info_pre.fastmode   = fastmode;

            return CUBLAS_STATUS_SUCCESS;
        }
    }

    // otherwise
    using gemm_t = cublasStatus_t (*)(cublasHandle_t,
                                      cublasOperation_t, cublasOperation_t,
                                      int, int, int,
                                      const void *,
                                      const void *, cudaDataType, int,
                                      const void *, cudaDataType, int,
                                      const void *,
                                      void *, cudaDataType, int,
                                      cublasComputeType_t, cublasGemmAlgo_t);

    static gemm_t real_gemm = (gemm_t)dlsym(RTLD_NEXT, STR(cublasGemmEx));
    if (!real_gemm) return CUBLAS_STATUS_NOT_INITIALIZED;
    return real_gemm(handle, transa, transb, m, n, k, alpha, A, Atype, lda, B, Btype, ldb, beta, C, Ctype, ldc, computeType, algo);

#endif
}
