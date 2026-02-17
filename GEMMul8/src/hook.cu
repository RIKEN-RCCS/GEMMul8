#include "../include/gemmul8.hpp"
#include "self_hipify.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <memory>
#include <mutex>
#include <string>
#include <unordered_map>
#include <dlfcn.h>

// NOTE:
// - skip_scalA/B rely on pointer identity; actual A/B contents are not verified.
// - Use GEMMUL8_SKIP_SCALE_* only when A/B data are unchanged.
// - When GEMMUL8_NUM_MOD_X < 2 or 20 < GEMMUL8_NUM_MOD_X, native GEMM routines are used.
/*
| Variable                | Default | Applies to | Description                                                                                               |
| :---------------------- | :------ | :--------- | :-------------------------------------------------------------------------------------------------------- |
| `GEMMUL8_BACKEND`       | `0`     | all        | Selects the emulation backend (`0` or `INT8` = INT8-based emulation, `1` or `FP8` = FP8-based emulation). |
| `GEMMUL8_NUM_MOD_D`     | `0`     | DGEMM      | Number of moduli used in DGEMM emulation. When num_moduli < 2 or 20 < num_moduli, native DGEMM is used.   |
| `GEMMUL8_NUM_MOD_S`     | `0`     | SGEMM      | Number of moduli used in SGEMM emulation. When num_moduli < 2 or 14 < num_moduli, native SGEMM is used.   |
| `GEMMUL8_NUM_MOD_Z`     | `0`     | ZGEMM      | Number of moduli used in ZGEMM emulation. When num_moduli < 2 or 20 < num_moduli, native ZGEMM is used.   |
| `GEMMUL8_NUM_MOD_C`     | `0`     | CGEMM      | Number of moduli used in CGEMM emulation. When num_moduli < 2 or 14 < num_moduli, native CGEMM is used.   |
| `GEMMUL8_FASTMODE_D`    | `0`     | DGEMM      | Enables fast mode (`1` = fast mode, `0` = accurate mode).                                                 |
| `GEMMUL8_FASTMODE_S`    | `0`     | SGEMM      | Enables fast mode (`1` = fast mode, `0` = accurate mode).                                                 |
| `GEMMUL8_FASTMODE_Z`    | `0`     | ZGEMM      | Enables fast mode (`1` = fast mode, `0` = accurate mode).                                                 |
| `GEMMUL8_FASTMODE_C`    | `0`     | CGEMM      | Enables fast mode (`1` = fast mode, `0` = accurate mode).                                                 |
| `GEMMUL8_MAXWS_BACKEND` | `0`     | all        | Max workspace calc target (`0` or `INT8`, `1` or `FP8`, `2` or `BOTH`). Default is INT8.                  |
| `GEMMUL8_MAX_M`         | `0`     | all        | Maximum value of `M` used to preallocate workspace memory.                                                |
| `GEMMUL8_MAX_N`         | `0`     | all        | Maximum value of `N` used to preallocate workspace memory.                                                |
| `GEMMUL8_MAX_K`         | `0`     | all        | Maximum value of `K` used to preallocate workspace memory.                                                |
| `GEMMUL8_MAX_NUM_MOD`   | `2`     | all        | Maximum number of moduli used when computing the size of the preallocated workspace.                      |
| `GEMMUL8_SKIP_SCALE_A`  | `0`     | all        | Enables skipping redundant preprocessing for `A` (`1` = enable, `0` = disable).                           |
| `GEMMUL8_SKIP_SCALE_B`  | `0`     | all        | Enables skipping redundant preprocessing for `B` (`1` = enable, `0` = disable).                           |
*/

namespace {

template <typename T> inline constexpr unsigned num_moduli_threshold  = 20u;
template <> inline constexpr unsigned num_moduli_threshold<float>     = 14u;
template <> inline constexpr unsigned num_moduli_threshold<cuComplex> = 14u;

// ---- Default initialization values ----
namespace initial_vals {

inline constexpr size_t MAX_M = 0u; // default M size
inline constexpr size_t MAX_N = 0u; // default N size
inline constexpr size_t MAX_K = 0u; // default K size

inline constexpr unsigned MAX_NUM_MOD = 2u; // default modulus count
inline constexpr unsigned NUM_MOD_D   = 0u; // default double moduli
inline constexpr unsigned NUM_MOD_S   = 0u; // default float moduli
inline constexpr unsigned NUM_MOD_Z   = 0u; // default double-complex moduli
inline constexpr unsigned NUM_MOD_C   = 0u; // default float-complex moduli

inline constexpr bool FASTMODE_D = false; // default double fastmode
inline constexpr bool FASTMODE_S = false; // default float fastmode
inline constexpr bool FASTMODE_Z = false; // default double-complex fastmode
inline constexpr bool FASTMODE_C = false; // default float-complex fastmode

inline constexpr bool SCALE_A = false; // default skip_scalA_switch
inline constexpr bool SCALE_B = false; // default skip_scalB_switch

} // namespace initial_vals

// ---- Structure holding GEMM configuration ----
struct Info_t {
    unsigned num_moduli      = 2;
    cublasOperation_t op_A   = CUBLAS_OP_N;
    cublasOperation_t op_B   = CUBLAS_OP_N;
    size_t m                 = 0;
    size_t n                 = 0;
    size_t k                 = 0;
    size_t lda               = 0;
    size_t ldb               = 0;
    const void *A            = nullptr;
    const void *B            = nullptr;
    int8_t *workA            = nullptr;
    int8_t *workB            = nullptr;
    char Type                = 'N';
    bool fastmode            = false;
    gemmul8::Backend backend = gemmul8::Backend::INT8;

    __host__ __forceinline__ bool match_core(unsigned nm, size_t kk, char t, bool fm, gemmul8::Backend b) const {
        return (num_moduli == nm) && (k == kk) && (Type == t) && (fastmode == fm) && (backend == b);
    }

    __host__ __forceinline__ bool can_skip_A(
        bool enable_switch, int8_t *cur_workA, const void *cur_A,
        size_t cur_m, size_t cur_lda, cublasOperation_t cur_opA //
    ) const {
        return enable_switch &&
               (workA == cur_workA) && (A == cur_A) &&
               (m == cur_m) && (lda == cur_lda) && (op_A == cur_opA);
    }

    __host__ __forceinline__ bool can_skip_B(
        bool enable_switch, int8_t *cur_workB, const void *cur_B,
        size_t cur_n, size_t cur_ldb, cublasOperation_t cur_opB //
    ) const {
        return enable_switch &&
               (workB == cur_workB) && (B == cur_B) &&
               (n == cur_n) && (ldb == cur_ldb) && (op_B == cur_opB);
    }
};

// HandleState + global map
struct HandleState {
    std::mutex mtx; // per-handle lock
    void *workA         = nullptr;
    size_t workA_size   = 0;
    void *workB         = nullptr;
    size_t workB_size   = 0;
    void *workC         = nullptr;
    size_t workC_size   = 0;
    cublasLtHandle_t lt = nullptr;
    Info_t last{};

    cudaStream_t last_stream = nullptr;
    cudaEvent_t last_event   = nullptr;
    bool last_stream_valid   = false;
};
static std::mutex g_state_map_mtx;
static std::unordered_map<cublasHandle_t, std::shared_ptr<HandleState>> g_state_map;

static inline std::shared_ptr<HandleState> get_state(cublasHandle_t h) {
    std::lock_guard<std::mutex> g(g_state_map_mtx);
    auto &p = g_state_map[h];
    if (!p) p = std::make_shared<HandleState>();
    return p;
}

static inline void erase_state(cublasHandle_t h) {
    std::lock_guard<std::mutex> g(g_state_map_mtx);
    g_state_map.erase(h);
}

static inline cublasStatus_t ensure_stream_ordered_locked(HandleState &st, cudaStream_t cur_stream) {
    if (!st.last_stream_valid) {
        st.last_stream       = cur_stream;
        st.last_stream_valid = true;
        return CUBLAS_STATUS_SUCCESS;
    }
    if (st.last_stream == cur_stream) return CUBLAS_STATUS_SUCCESS;

    if (!st.last_event) {
        cudaError_t e = cudaEventCreateWithFlags(&st.last_event, cudaEventDisableTiming);
        if (e != cudaSuccess) return CUBLAS_STATUS_INTERNAL_ERROR;
    }

    cudaError_t e1 = cudaEventRecord(st.last_event, st.last_stream);
    if (e1 != cudaSuccess) return CUBLAS_STATUS_INTERNAL_ERROR;

    cudaError_t e2 = cudaStreamWaitEvent(cur_stream, st.last_event, 0);
    if (e2 != cudaSuccess) return CUBLAS_STATUS_INTERNAL_ERROR;

    st.last_stream = cur_stream;
    return CUBLAS_STATUS_SUCCESS;
}

// ---- Global max workspace (process-wide) ----
static size_t max_workSizeA   = 0;     // global workspace limit in byte
static size_t max_workSizeB   = 0;     // global workspace limit in byte
static size_t max_workSizeC   = 0;     // global workspace limit in byte

// ---- Small helpers ----
static inline bool env_is_one(const char *s, bool def) {
    return s ? (std::strcmp(s, "1") == 0) : def;
}

static inline unsigned env_u32(const char *s, unsigned def) {
    if (!s) return def;
    try {
        return std::stoul(s);
    } catch (...) {
        return def;
    }
}

static inline size_t env_u64(const char *s, size_t def) {
    if (!s) return def;
    try {
        return static_cast<size_t>(std::stoull(s));
    } catch (...) {
        return def;
    }
}

enum class MaxWSBackend : unsigned {
    INT8 = 0,
    FP8  = 1,
    BOTH = 2,
};

static inline MaxWSBackend env_maxws_backend(const char *s, MaxWSBackend def) {
    if (!s) return def;

    if (std::strcmp(s, "0") == 0) return MaxWSBackend::INT8;
    if (std::strcmp(s, "1") == 0) return MaxWSBackend::FP8;
    if (std::strcmp(s, "2") == 0) return MaxWSBackend::BOTH;

    if (std::strcmp(s, "INT8") == 0) return MaxWSBackend::INT8;
    if (std::strcmp(s, "FP8") == 0) return MaxWSBackend::FP8;
    if (std::strcmp(s, "BOTH") == 0) return MaxWSBackend::BOTH;

    return def;
}

static inline MaxWSBackend requested_maxws_backend() {
    return env_maxws_backend(getenv("GEMMUL8_MAXWS_BACKEND"), MaxWSBackend::INT8);
}

static inline gemmul8::Backend env_backend(const char *s, gemmul8::Backend def) {
    if (!s) return def;
    if (std::strcmp(s, "0") == 0) return gemmul8::Backend::INT8;
    if (std::strcmp(s, "1") == 0) return gemmul8::Backend::FP8;
    if (std::strcmp(s, "INT8") == 0) return gemmul8::Backend::INT8;
    if (std::strcmp(s, "FP8") == 0) return gemmul8::Backend::FP8;
    return def;
}

static inline gemmul8::Backend requested_backend() {
    return env_backend(getenv("GEMMUL8_BACKEND"), gemmul8::Backend::INT8);
}

// ---- Initialize maximum workspace size ----
static std::once_flag g_maxws_once;

static void init_max_workspace() {
    std::call_once(g_maxws_once, [] {
        size_t max_m        = initial_vals::MAX_M;
        size_t max_n        = initial_vals::MAX_N;
        size_t max_k        = initial_vals::MAX_K;
        unsigned max_moduli = initial_vals::MAX_NUM_MOD;

        max_m      = env_u64(getenv("GEMMUL8_MAX_M"), max_m);
        max_n      = env_u64(getenv("GEMMUL8_MAX_N"), max_n);
        max_k      = env_u64(getenv("GEMMUL8_MAX_K"), max_k);
        max_moduli = env_u32(getenv("GEMMUL8_MAX_NUM_MOD"), max_moduli);

        unsigned NMOD_Z             = env_u32(getenv("GEMMUL8_NUM_MOD_Z"), 0u);
        unsigned NMOD_C             = env_u32(getenv("GEMMUL8_NUM_MOD_C"), 0u);
        const bool want_complex_max = (NMOD_Z > 0u) || (NMOD_C > 0u);

        const MaxWSBackend mws = requested_maxws_backend();
        const bool do_int8     = (mws == MaxWSBackend::INT8) || (mws == MaxWSBackend::BOTH);
        const bool do_fp8      = (mws == MaxWSBackend::FP8) || (mws == MaxWSBackend::BOTH);

        size_t wA_i = 0, wB_i = 0, w_i = 0, wC_i = 0;
        size_t wA_f = 0, wB_f = 0, w_f = 0, wC_f = 0;

        if (do_int8) {
            if (want_complex_max) {
                w_i = gemmul8::workSize<true, gemmul8::Backend::INT8>(
                    max_m, max_n, max_k, max_moduli, true, true, &wA_i, &wB_i);
            } else {
                w_i = gemmul8::workSize<false, gemmul8::Backend::INT8>(
                    max_m, max_n, max_k, max_moduli, true, true, &wA_i, &wB_i);
            }
            wC_i = (w_i > (wA_i + wB_i)) ? (w_i - wA_i - wB_i) : 0;
        }

        if (do_fp8) {
            if (want_complex_max) {
                w_f = gemmul8::workSize<true, gemmul8::Backend::FP8>(
                    max_m, max_n, max_k, max_moduli, true, true, &wA_f, &wB_f);
            } else {
                w_f = gemmul8::workSize<false, gemmul8::Backend::FP8>(
                    max_m, max_n, max_k, max_moduli, true, true, &wA_f, &wB_f);
            }
            wC_f = (w_f > (wA_f + wB_f)) ? (w_f - wA_f - wB_f) : 0;
        }

        max_workSizeA = std::max(wA_i, wA_f);
        max_workSizeB = std::max(wB_i, wB_f);
        max_workSizeC = std::max(wC_i, wC_f);
    });
}

// ---- Environment per type ----
static inline void get_env_d(unsigned &num_moduli, bool &fastmode, bool &enable_skipA, bool &enable_skipB) {
    num_moduli   = env_u32(getenv("GEMMUL8_NUM_MOD_D"), initial_vals::NUM_MOD_D);
    fastmode     = env_is_one(getenv("GEMMUL8_FASTMODE_D"), initial_vals::FASTMODE_D);
    enable_skipA = env_is_one(getenv("GEMMUL8_SKIP_SCALE_A"), initial_vals::SCALE_A);
    enable_skipB = env_is_one(getenv("GEMMUL8_SKIP_SCALE_B"), initial_vals::SCALE_B);
}

static inline void get_env_s(unsigned &num_moduli, bool &fastmode, bool &enable_skipA, bool &enable_skipB) {
    num_moduli   = env_u32(getenv("GEMMUL8_NUM_MOD_S"), initial_vals::NUM_MOD_S);
    fastmode     = env_is_one(getenv("GEMMUL8_FASTMODE_S"), initial_vals::FASTMODE_S);
    enable_skipA = env_is_one(getenv("GEMMUL8_SKIP_SCALE_A"), initial_vals::SCALE_A);
    enable_skipB = env_is_one(getenv("GEMMUL8_SKIP_SCALE_B"), initial_vals::SCALE_B);
}

static inline void get_env_z(unsigned &num_moduli, bool &fastmode, bool &enable_skipA, bool &enable_skipB) {
    num_moduli   = env_u32(getenv("GEMMUL8_NUM_MOD_Z"), initial_vals::NUM_MOD_Z);
    fastmode     = env_is_one(getenv("GEMMUL8_FASTMODE_Z"), initial_vals::FASTMODE_Z);
    enable_skipA = env_is_one(getenv("GEMMUL8_SKIP_SCALE_A"), initial_vals::SCALE_A);
    enable_skipB = env_is_one(getenv("GEMMUL8_SKIP_SCALE_B"), initial_vals::SCALE_B);
}

static inline void get_env_c(unsigned &num_moduli, bool &fastmode, bool &enable_skipA, bool &enable_skipB) {
    num_moduli   = env_u32(getenv("GEMMUL8_NUM_MOD_C"), initial_vals::NUM_MOD_C);
    fastmode     = env_is_one(getenv("GEMMUL8_FASTMODE_C"), initial_vals::FASTMODE_C);
    enable_skipA = env_is_one(getenv("GEMMUL8_SKIP_SCALE_A"), initial_vals::SCALE_A);
    enable_skipB = env_is_one(getenv("GEMMUL8_SKIP_SCALE_B"), initial_vals::SCALE_B);
}

// ---- Lt handle management ----
static cublasStatus_t ensure_lt_handle_locked(HandleState &st, cublasLtHandle_t *out) {
    if (!out) return CUBLAS_STATUS_INVALID_VALUE;
    if (st.lt) {
        *out = st.lt;
        return CUBLAS_STATUS_SUCCESS;
    }
    cublasLtHandle_t lt = nullptr;
    cublasStatus_t s    = cublasLtCreate(&lt);
    if (s != CUBLAS_STATUS_SUCCESS) {
        *out = nullptr;
        return s;
    }
    st.lt = lt;
    *out  = lt;
    return CUBLAS_STATUS_SUCCESS;
}

// ---- Workspace management ----
static cublasStatus_t get_work_locked(void *&buf, size_t &buf_size, size_t req_size, void **out, const char *tag, cudaStream_t stream) {

    if (!out) return CUBLAS_STATUS_INVALID_VALUE;
    init_max_workspace(); // call_once

    if (req_size == 0) {
        *out = nullptr;
        return CUBLAS_STATUS_SUCCESS;
    }

    // grow-only policy (never shrink)
    if (buf && buf_size >= req_size) {
        *out = buf;
        return CUBLAS_STATUS_SUCCESS;
    }

    void *newp      = nullptr;
    cudaError_t err = cudaMallocAsync(&newp, req_size, stream);
    if (err != cudaSuccess) {
        *out = nullptr;
        std::cerr << "[GEMMUL8 HOOK] cudaMallocAsync failed for " << (tag ? tag : "workspace")
                  << " size " << req_size << " bytes. Error: " << cudaGetErrorString(err) << "\n";
        return CUBLAS_STATUS_ALLOC_FAILED;
    }

    if (buf) {
        cudaError_t free_err = cudaFreeAsync(buf, stream);
        if (free_err != cudaSuccess) {
            std::cerr << "[GEMMUL8 HOOK] Warning: cudaFreeAsync failed for " << (tag ? tag : "workspace")
                      << " (" << cudaGetErrorString(free_err) << ")\n";

            cudaError_t free_new_err = cudaFreeAsync(newp, stream);
            if (free_new_err != cudaSuccess) {
                std::cerr << "[GEMMUL8 HOOK] Warning: cudaFreeAsync failed for newly-allocated "
                          << (tag ? tag : "workspace") << " (" << cudaGetErrorString(free_new_err) << ")\n";
            }

            *out = nullptr;
            return CUBLAS_STATUS_INTERNAL_ERROR;
        }
    }

    buf      = newp;
    buf_size = req_size;
    *out     = buf;
    return CUBLAS_STATUS_SUCCESS;
}

static void cleanup_work(cublasHandle_t handle) {
    auto sp = get_state(handle);

    cudaStream_t stream = nullptr;
    bool have_stream    = false;

    {
        std::lock_guard<std::mutex> lk(sp->mtx);
        // 可能なら「最後に使った stream」を優先
        if (sp->last_stream_valid) {
            stream      = sp->last_stream;
            have_stream = true;
        }
    }

    if (!have_stream) {
        cudaStream_t s = 0;
        if (cublasGetStream(handle, &s) == CUBLAS_STATUS_SUCCESS) {
            stream      = s;
            have_stream = true;
        }
    }

    if (!have_stream) {
        cudaError_t sync_err = cudaDeviceSynchronize();
        if (sync_err != cudaSuccess) {
            std::cerr << "[GEMMUL8 HOOK] cleanup_work: cudaDeviceSynchronize failed ("
                      << cudaGetErrorString(sync_err) << ")\n";
        }
    }

    {
        std::lock_guard<std::mutex> lk(sp->mtx);

        auto free_one = [&](void *&ptr, size_t &sz, const char *name) {
            if (!ptr) return;

            if (have_stream) {
                cudaError_t e = cudaFreeAsync(ptr, stream);
                if (e != cudaSuccess) {
                    std::cerr << "[GEMMUL8 HOOK] cublasDestroy: cudaFreeAsync " << (name ? name : "workspace")
                              << " failed (" << cudaGetErrorString(e) << ")\n";
                    cudaError_t se = cudaStreamSynchronize(stream);
                    if (se == cudaSuccess) {
                        cudaError_t fe = cudaFree(ptr);
                        if (fe != cudaSuccess) {
                            std::cerr << "[GEMMUL8 HOOK] cublasDestroy: cudaFree fallback failed for "
                                      << (name ? name : "workspace") << " (" << cudaGetErrorString(fe) << ")\n";
                        }
                    } else {
                        std::cerr << "[GEMMUL8 HOOK] cublasDestroy: cudaStreamSynchronize failed ("
                                  << cudaGetErrorString(se) << ")\n";
                    }
                }
            } else {
                cudaError_t fe = cudaFree(ptr);
                if (fe != cudaSuccess) {
                    std::cerr << "[GEMMUL8 HOOK] cublasDestroy: cudaFree (no-stream) failed for "
                              << (name ? name : "workspace") << " (" << cudaGetErrorString(fe) << ")\n";
                }
            }

            ptr = nullptr;
            sz  = 0;
        };

        free_one(sp->workA, sp->workA_size, "workA");
        free_one(sp->workB, sp->workB_size, "workB");
        free_one(sp->workC, sp->workC_size, "workC");

        if (sp->lt) {
            (void)cublasLtDestroy(sp->lt);
            sp->lt = nullptr;
        }

        if (sp->last_event) {
            cudaEventDestroy(sp->last_event);
            sp->last_event = nullptr;
        }

        sp->last_stream       = nullptr;
        sp->last_stream_valid = false;

        sp->last = Info_t{};
    }

    erase_state(handle);
}

// ---- dlsym loader ----
template <typename Fn>
static inline Fn load_real(const char *sym) {
    return reinterpret_cast<Fn>(dlsym(RTLD_NEXT, sym));
}

// ---- Traits to unify duplicated GEMM bodies ----
template <typename T> struct GemmTraits;

template <> struct GemmTraits<float> {
    static constexpr bool isComplex = false;
    static constexpr char typeChar  = 'S';
    static inline void get_env(unsigned &nm, bool &fm, bool &enA, bool &enB) { get_env_s(nm, fm, enA, enB); }

    using GemmV2Fn = cublasStatus_t (*)(cublasHandle_t, cublasOperation_t, cublasOperation_t,
                                        int, int, int,
                                        const float *, const float *, int,
                                        const float *, int,
                                        const float *, float *, int);

    static constexpr const char *gemm_v2_sym = STR(cublasSgemm_v2);
};

template <> struct GemmTraits<double> {
    static constexpr bool isComplex = false;
    static constexpr char typeChar  = 'D';
    static inline void get_env(unsigned &nm, bool &fm, bool &enA, bool &enB) { get_env_d(nm, fm, enA, enB); }

    using GemmV2Fn = cublasStatus_t (*)(cublasHandle_t, cublasOperation_t, cublasOperation_t,
                                        int, int, int,
                                        const double *, const double *, int,
                                        const double *, int,
                                        const double *, double *, int);

    static constexpr const char *gemm_v2_sym = STR(cublasDgemm_v2);
};

template <> struct GemmTraits<cuComplex> {
    static constexpr bool isComplex = true;
    static constexpr char typeChar  = 'C';
    static inline void get_env(unsigned &nm, bool &fm, bool &enA, bool &enB) { get_env_c(nm, fm, enA, enB); }

    using GemmV2Fn = cublasStatus_t (*)(cublasHandle_t, cublasOperation_t, cublasOperation_t,
                                        int, int, int,
                                        const cuComplex *, const cuComplex *, int,
                                        const cuComplex *, int,
                                        const cuComplex *, cuComplex *, int);

    static constexpr const char *gemm_v2_sym = STR(cublasCgemm_v2);
};

template <> struct GemmTraits<cuDoubleComplex> {
    static constexpr bool isComplex = true;
    static constexpr char typeChar  = 'Z';
    static inline void get_env(unsigned &nm, bool &fm, bool &enA, bool &enB) { get_env_z(nm, fm, enA, enB); }

    using GemmV2Fn = cublasStatus_t (*)(cublasHandle_t, cublasOperation_t, cublasOperation_t,
                                        int, int, int,
                                        const cuDoubleComplex *, const cuDoubleComplex *, int,
                                        const cuDoubleComplex *, int,
                                        const cuDoubleComplex *, cuDoubleComplex *, int);

    static constexpr const char *gemm_v2_sym = STR(cublasZgemm_v2);
};

template <typename T> static inline typename GemmTraits<T>::GemmV2Fn real_gemm_v2() {
    using Fn     = typename GemmTraits<T>::GemmV2Fn;
    static Fn fn = load_real<Fn>(GemmTraits<T>::gemm_v2_sym);
    return fn;
}

// ---- GEMM call wrapper (INT8 or FP8) ----
template <typename T> static inline cublasStatus_t call_gemmul8_gemm(
    gemmul8::Backend backend,
    cublasHandle_t handle,
    HandleState &hst,
    cublasOperation_t transa, cublasOperation_t transb,
    int m, int n, int k,
    const T *alpha, const T *A, int lda, const T *B, int ldb,
    const T *beta, T *C, int ldc,
    unsigned num_moduli, bool fastmode,
    void *workC, void *workA, void *workB,
    bool enable_skip_A, bool enable_skip_B,
    bool skip_A, bool skip_B,
    cudaStream_t stream //
) {
    if (backend == gemmul8::Backend::INT8) {
        (void)gemmul8::gemm<T, gemmul8::Backend::INT8>(
            handle,
            transa, transb,
            static_cast<size_t>(m), static_cast<size_t>(n), static_cast<size_t>(k),
            alpha,
            A, static_cast<size_t>(lda),
            B, static_cast<size_t>(ldb),
            beta,
            C, static_cast<size_t>(ldc),
            num_moduli, fastmode,
            workC, workA, workB,
            enable_skip_A, enable_skip_B,
            skip_A, skip_B);
        return CUBLAS_STATUS_SUCCESS;
    }

    cublasLtHandle_t lt  = nullptr;
    cublasStatus_t st_lt = ensure_lt_handle_locked(hst, &lt);
    if (st_lt != CUBLAS_STATUS_SUCCESS) return st_lt;

    (void)gemmul8::gemm<T, gemmul8::Backend::FP8>(
        lt,
        transa, transb,
        static_cast<size_t>(m), static_cast<size_t>(n), static_cast<size_t>(k),
        alpha,
        A, static_cast<size_t>(lda),
        B, static_cast<size_t>(ldb),
        beta,
        C, static_cast<size_t>(ldc),
        num_moduli, fastmode,
        workC, workA, workB,
        enable_skip_A, enable_skip_B,
        skip_A, skip_B,
        stream);

    return CUBLAS_STATUS_SUCCESS;
}

// ---- workSize wrapper (INT8 or FP8) ----
template <typename T> static inline size_t call_gemmul8_workSize(
    gemmul8::Backend backend,
    int m, int n, int k,
    unsigned num_moduli,
    bool enable_skip_A, bool enable_skip_B,
    size_t *wA, size_t *wB //
) {
    if (backend == gemmul8::Backend::INT8) {
        return gemmul8::workSize<GemmTraits<T>::isComplex, gemmul8::Backend::INT8>(
            static_cast<size_t>(m), static_cast<size_t>(n), static_cast<size_t>(k),
            num_moduli, enable_skip_A, enable_skip_B, wA, wB);
    }

    return gemmul8::workSize<GemmTraits<T>::isComplex, gemmul8::Backend::FP8>(
        static_cast<size_t>(m), static_cast<size_t>(n), static_cast<size_t>(k),
        num_moduli, enable_skip_A, enable_skip_B, wA, wB);
}

// ---- Shared body for *gemm_v2 hooks ----
template <typename T> static cublasStatus_t gemm_v2_impl(
    cublasHandle_t handle,
    cublasOperation_t transa, cublasOperation_t transb,
    int m, int n, int k,
    const T *alpha, const T *A, int lda, const T *B, int ldb,
    const T *beta, T *C, int ldc //
) {
    if (m <= 0 || n <= 0 || k <= 0) return CUBLAS_STATUS_SUCCESS;
    if (!A || !B || !C) return CUBLAS_STATUS_INVALID_VALUE;

    unsigned num_moduli = 0;
    bool fastmode       = false;
    bool enable_skipA   = false;
    bool enable_skipB   = false;
    GemmTraits<T>::get_env(num_moduli, fastmode, enable_skipA, enable_skipB);

    if (num_moduli < 2u || num_moduli_threshold<T> < num_moduli) {
        auto fn = real_gemm_v2<T>();
        if (!fn) return CUBLAS_STATUS_NOT_INITIALIZED;
        return fn(handle, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
    }

    const gemmul8::Backend backend = requested_backend();

    auto sp = get_state(handle);
    std::lock_guard<std::mutex> lk(sp->mtx);

    init_max_workspace();

    cudaStream_t stream      = 0;
    cublasStatus_t st_stream = cublasGetStream(handle, &stream);
    if (st_stream != CUBLAS_STATUS_SUCCESS) return st_stream;

    cublasStatus_t st_ord = ensure_stream_ordered_locked(*sp, stream);
    if (st_ord != CUBLAS_STATUS_SUCCESS) return st_ord;

    size_t wsizeA = 0, wsizeB = 0;
    const size_t wsize = call_gemmul8_workSize<T>(
        backend, m, n, k, num_moduli,
        enable_skipA, enable_skipB,
        &wsizeA, &wsizeB);

    // wsize is assumed to cover A+B+rest; if not, layout is inconsistent.
    if (wsize < wsizeA + wsizeB) return CUBLAS_STATUS_INVALID_VALUE;
    const size_t needA = wsizeA;
    const size_t needB = wsizeB;
    const size_t needC = wsize - needA - needB;

    // Allocate 3 independent workspaces (grow-only).
    size_t reqA = needA;
    size_t reqB = needB;
    size_t reqC = needC;

    const bool enforce_maxws = (enable_skipA || enable_skipB);
    if (enforce_maxws) {
        if (enable_skipA) reqA = std::max(reqA, max_workSizeA);
        if (enable_skipB) reqB = std::max(reqB, max_workSizeB);
        reqC = std::max(reqC, max_workSizeC);
    }

    void *workA_raw = nullptr;
    void *workB_raw = nullptr;
    void *workC_raw = nullptr;

    cublasStatus_t st = get_work_locked(sp->workA, sp->workA_size, reqA, &workA_raw, "workA", stream);
    if (st != CUBLAS_STATUS_SUCCESS) return st;
    st = get_work_locked(sp->workB, sp->workB_size, reqB, &workB_raw, "workB", stream);
    if (st != CUBLAS_STATUS_SUCCESS) return st;
    st = get_work_locked(sp->workC, sp->workC_size, reqC, &workC_raw, "workC", stream);
    if (st != CUBLAS_STATUS_SUCCESS) return st;

    int8_t *workA = reinterpret_cast<int8_t *>(workA_raw);
    int8_t *workB = reinterpret_cast<int8_t *>(workB_raw);
    int8_t *workC = reinterpret_cast<int8_t *>(workC_raw);

    bool skipA = false;
    bool skipB = false;

    Info_t &info_pre = sp->last;
    if (info_pre.match_core(num_moduli, static_cast<size_t>(k), GemmTraits<T>::typeChar, fastmode, backend)) {
        if (info_pre.can_skip_A(enable_skipA, workA, A, static_cast<size_t>(m), static_cast<size_t>(lda), transa)) skipA = true;
        if (info_pre.can_skip_B(enable_skipB, workB, B, static_cast<size_t>(n), static_cast<size_t>(ldb), transb)) skipB = true;
    }

    st = call_gemmul8_gemm<T>(
        backend,
        handle,
        *sp,
        transa, transb,
        m, n, k,
        alpha, A, lda,
        B, ldb,
        beta,
        C, ldc,
        num_moduli, fastmode,
        reinterpret_cast<void *>(workC),
        reinterpret_cast<void *>(workA),
        reinterpret_cast<void *>(workB),
        enable_skipA, enable_skipB,
        skipA, skipB, stream);

    if (st != CUBLAS_STATUS_SUCCESS) return st;

    // update cache
    info_pre.num_moduli = num_moduli;
    info_pre.op_A       = transa;
    info_pre.op_B       = transb;
    info_pre.m          = static_cast<size_t>(m);
    info_pre.n          = static_cast<size_t>(n);
    info_pre.k          = static_cast<size_t>(k);
    info_pre.lda        = static_cast<size_t>(lda);
    info_pre.ldb        = static_cast<size_t>(ldb);
    info_pre.A          = A;
    info_pre.B          = B;
    info_pre.workA      = workA;
    info_pre.workB      = workB;
    info_pre.Type       = GemmTraits<T>::typeChar;
    info_pre.fastmode   = fastmode;
    info_pre.backend    = backend;

    return CUBLAS_STATUS_SUCCESS;
}

// ---- Shared body for GemmEx “emulation path” ----
template <typename T> static cublasStatus_t gemm_ex_impl(
    cublasHandle_t handle,
    cublasOperation_t transa, cublasOperation_t transb,
    int m, int n, int k,
    const void *alpha, const void *A, int lda, const void *B, int ldb,
    const void *beta, void *C, int ldc,
    unsigned num_moduli, bool fastmode,
    gemmul8::Backend backend,
    bool enable_skipA, bool enable_skipB //
) {
    auto sp = get_state(handle);
    std::lock_guard<std::mutex> lk(sp->mtx);

    init_max_workspace();

    cudaStream_t stream      = 0;
    cublasStatus_t st_stream = cublasGetStream(handle, &stream);
    if (st_stream != CUBLAS_STATUS_SUCCESS) return st_stream;

    cublasStatus_t st_ord = ensure_stream_ordered_locked(*sp, stream);
    if (st_ord != CUBLAS_STATUS_SUCCESS) return st_ord;

    size_t wsizeA = 0, wsizeB = 0;
    const size_t wsize = call_gemmul8_workSize<T>(
        backend, m, n, k, num_moduli,
        enable_skipA, enable_skipB,
        &wsizeA, &wsizeB);

    if (wsize < wsizeA + wsizeB) return CUBLAS_STATUS_INVALID_VALUE;
    const size_t needA = wsizeA;
    const size_t needB = wsizeB;
    const size_t needC = wsize - needA - needB;

    size_t reqA = needA;
    size_t reqB = needB;
    size_t reqC = needC;

    const bool enforce_maxws = (enable_skipA || enable_skipB);
    if (enforce_maxws) {
        if (enable_skipA) reqA = std::max(reqA, max_workSizeA);
        if (enable_skipB) reqB = std::max(reqB, max_workSizeB);
        reqC = std::max(reqC, max_workSizeC);
    }

    void *workA_raw = nullptr;
    void *workB_raw = nullptr;
    void *workC_raw = nullptr;

    cublasStatus_t st = get_work_locked(sp->workA, sp->workA_size, reqA, &workA_raw, "workA", stream);
    if (st != CUBLAS_STATUS_SUCCESS) return st;
    st = get_work_locked(sp->workB, sp->workB_size, reqB, &workB_raw, "workB", stream);
    if (st != CUBLAS_STATUS_SUCCESS) return st;
    st = get_work_locked(sp->workC, sp->workC_size, reqC, &workC_raw, "workC", stream);
    if (st != CUBLAS_STATUS_SUCCESS) return st;

    int8_t *workA = reinterpret_cast<int8_t *>(workA_raw);
    int8_t *workB = reinterpret_cast<int8_t *>(workB_raw);
    int8_t *workC = reinterpret_cast<int8_t *>(workC_raw);

    bool skipA = false;
    bool skipB = false;

    Info_t &info_pre = sp->last;
    if (info_pre.match_core(num_moduli, static_cast<size_t>(k), GemmTraits<T>::typeChar, fastmode, backend)) {
        if (info_pre.can_skip_A(enable_skipA, workA, A, static_cast<size_t>(m), static_cast<size_t>(lda), transa)) skipA = true;
        if (info_pre.can_skip_B(enable_skipB, workB, B, static_cast<size_t>(n), static_cast<size_t>(ldb), transb)) skipB = true;
    }

    st = call_gemmul8_gemm<T>(
        backend,
        handle,
        *sp,
        transa, transb,
        m, n, k,
        reinterpret_cast<const T *>(alpha),
        reinterpret_cast<const T *>(A), lda,
        reinterpret_cast<const T *>(B), ldb,
        reinterpret_cast<const T *>(beta),
        reinterpret_cast<T *>(C), ldc,
        num_moduli, fastmode,
        reinterpret_cast<void *>(workC),
        reinterpret_cast<void *>(workA),
        reinterpret_cast<void *>(workB),
        enable_skipA, enable_skipB,
        skipA, skipB, stream);

    if (st != CUBLAS_STATUS_SUCCESS) return st;

    // update cache
    info_pre.num_moduli = num_moduli;
    info_pre.op_A       = transa;
    info_pre.op_B       = transb;
    info_pre.m          = static_cast<size_t>(m);
    info_pre.n          = static_cast<size_t>(n);
    info_pre.k          = static_cast<size_t>(k);
    info_pre.lda        = static_cast<size_t>(lda);
    info_pre.ldb        = static_cast<size_t>(ldb);
    info_pre.A          = A;
    info_pre.B          = B;
    info_pre.workA      = workA;
    info_pre.workB      = workB;
    info_pre.Type       = GemmTraits<T>::typeChar;
    info_pre.fastmode   = fastmode;
    info_pre.backend    = backend;

    return CUBLAS_STATUS_SUCCESS;
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
    using Fn       = cublasStatus_t (*)(cublasHandle_t);
    static Fn real = load_real<Fn>(STR(cublasDestroy_v2));
    if (!real) return CUBLAS_STATUS_NOT_INITIALIZED;
    return real(handle);
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
    return gemm_v2_impl<float>(handle, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
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
    return gemm_v2_impl<double>(handle, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
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
    return gemm_v2_impl<cuComplex>(handle, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
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
    return gemm_v2_impl<cuDoubleComplex>(handle, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
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

    const gemmul8::Backend backend = requested_backend(); // same env for all T

    // SGEMM case
    if (computeType == CUBLAS_COMPUTE_32F &&
        Atype == CUDA_R_32F && Btype == CUDA_R_32F && Ctype == CUDA_R_32F) {

        unsigned num_moduli = 0;
        bool fastmode       = false;
        bool enable_skipA   = false;
        bool enable_skipB   = false;
        get_env_s(num_moduli, fastmode, enable_skipA, enable_skipB);

        if (2u <= num_moduli && num_moduli <= num_moduli_threshold<float>) {
            return gemm_ex_impl<float>(
                handle, transa, transb, m, n, k,
                alpha, A, lda, B, ldb, beta, C, ldc,
                num_moduli, fastmode, backend, enable_skipA, enable_skipB);
        }
    }

    // DGEMM case
    if (computeType == CUBLAS_COMPUTE_64F &&
        Atype == CUDA_R_64F && Btype == CUDA_R_64F && Ctype == CUDA_R_64F) {

        unsigned num_moduli = 0;
        bool fastmode       = false;
        bool enable_skipA   = false;
        bool enable_skipB   = false;
        get_env_d(num_moduli, fastmode, enable_skipA, enable_skipB);

        if (2u <= num_moduli && num_moduli <= num_moduli_threshold<double>) {
            return gemm_ex_impl<double>(
                handle, transa, transb, m, n, k,
                alpha, A, lda, B, ldb, beta, C, ldc,
                num_moduli, fastmode, backend, enable_skipA, enable_skipB);
        }
    }

    // CGEMM case
    if (computeType == CUBLAS_COMPUTE_32F &&
        Atype == CUDA_C_32F && Btype == CUDA_C_32F && Ctype == CUDA_C_32F) {

        unsigned num_moduli = 0;
        bool fastmode       = false;
        bool enable_skipA   = false;
        bool enable_skipB   = false;
        get_env_c(num_moduli, fastmode, enable_skipA, enable_skipB);

        if (2u <= num_moduli && num_moduli <= num_moduli_threshold<cuComplex>) {
            return gemm_ex_impl<cuComplex>(
                handle, transa, transb, m, n, k,
                alpha, A, lda, B, ldb, beta, C, ldc,
                num_moduli, fastmode, backend, enable_skipA, enable_skipB);
        }
    }

    // ZGEMM case
    if (computeType == CUBLAS_COMPUTE_64F &&
        Atype == CUDA_C_64F && Btype == CUDA_C_64F && Ctype == CUDA_C_64F) {

        unsigned num_moduli = 0;
        bool fastmode       = false;
        bool enable_skipA   = false;
        bool enable_skipB   = false;
        get_env_z(num_moduli, fastmode, enable_skipA, enable_skipB);

        if (2u <= num_moduli && num_moduli <= num_moduli_threshold<cuDoubleComplex>) {
            return gemm_ex_impl<cuDoubleComplex>(
                handle, transa, transb, m, n, k,
                alpha, A, lda, B, ldb, beta, C, ldc,
                num_moduli, fastmode, backend, enable_skipA, enable_skipB);
        }
    }

    // otherwise: call native GemmEx
    using Fn = cublasStatus_t (*)(cublasHandle_t,
                                  cublasOperation_t, cublasOperation_t,
                                  int, int, int,
                                  const void *,
                                  const void *, cudaDataType, int,
                                  const void *, cudaDataType, int,
                                  const void *,
                                  void *, cudaDataType, int,
                                  cublasComputeType_t, cublasGemmAlgo_t);

    static Fn real = load_real<Fn>(STR(cublasGemmEx));
    if (!real) return CUBLAS_STATUS_NOT_INITIALIZED;

    return real(handle, transa, transb,
                m, n, k,
                alpha,
                A, Atype, lda,
                B, Btype, ldb,
                beta,
                C, Ctype, ldc,
                computeType, algo);
#endif
}
