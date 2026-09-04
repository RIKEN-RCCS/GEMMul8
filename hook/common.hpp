/**
 * Environment variables for GEMMul8 hook mode
 * -------------------------------------------
 *
 * GEMMul8 hook mode can intercept standard cuBLAS/hipBLAS routines and route
 * selected calls to GEMMul8 emulation.  Runtime control is specified by
 * operation-specific environment variables.
 *
 * Runtime emulation variables:
 *
 *   | Variable pattern                         | Default | Description                                                                         |
 *   | :--------------------------------------- | :------ | :---------------------------------------------------------------------------------- |
 *   | GEMMUL8_BACKEND_<OP>                     | INT8    | Emulation backend. 0 or INT8 = INT8 backend; 1 or FP8 = FP8 backend.                |
 *   | GEMMUL8_NUM_MOD_S_<OP>                   | 0       | Number of moduli for FP32 real routines. Native BLAS is used if outside [2, 13].    |
 *   | GEMMUL8_NUM_MOD_D_<OP>                   | 0       | Number of moduli for FP64 real routines. Native BLAS is used if outside [2, 20].    |
 *   | GEMMUL8_NUM_MOD_C_<OP>                   | 0       | Number of moduli for FP32 complex routines. Native BLAS is used if outside [2, 13]. |
 *   | GEMMUL8_NUM_MOD_Z_<OP>                   | 0       | Number of moduli for FP64 complex routines. Native BLAS is used if outside [2, 20]. |
 *   | GEMMUL8_FASTMODE_S_<OP>                  | 1       | Fast mode switch for FP32 real routines. 1 = fast mode; 0 = accurate mode.          |
 *   | GEMMUL8_FASTMODE_D_<OP>                  | 1       | Fast mode switch for FP64 real routines. 1 = fast mode; 0 = accurate mode.          |
 *   | GEMMUL8_FASTMODE_C_<OP>                  | 1       | Fast mode switch for FP32 complex routines. 1 = fast mode; 0 = accurate mode.       |
 *   | GEMMUL8_FASTMODE_Z_<OP>                  | 1       | Fast mode switch for FP64 complex routines. 1 = fast mode; 0 = accurate mode.       |
 *   | GEMMUL8_MEMORY_SAVING                    | unset   | Memory-saving override. 1 = enable; 0 = disable.                                    |
 *   | GEMMUL8_MAX_WORKSIZE                     | unset   | Workspace-size limit in bytes used when memory saving is enabled.                   |
 *   | GEMMUL8_BLK_SIZE_TRSM                    | unset   | TRSM block-size override. >0 = use the specified size; <=0 = automatic selection.   |
 *   | GEMMUL8_SKIP_SCALE_A                     | 0       | Enables reuse of preprocessed/scaled A when the operand cache key matches.          |
 *   | GEMMUL8_SKIP_SCALE_B                     | 0       | Enables reuse of preprocessed/scaled B when the operand cache key matches.          |
 *
 * Operation suffixes:
 *
 *   GEMM
 *   SYMM_LEFT, SYMM_RIGHT
 *   SYRK
 *   SYR2K
 *   SYRKX
 *   HEMM_LEFT, HEMM_RIGHT
 *   HERK
 *   HER2K
 *   HERKX
 *   TRMM_LEFT, TRMM_RIGHT
 *   TRSM_LEFT, TRSM_RIGHT
 *
 * Max-workspace preallocation variables:
 *
 *   GEMMul8 normally grows workspace buffers on demand.
 *   To stabilize workspace addresses, avoid reallocating workspace, and improve skip-scaling reuse,
 *   define the maximum BLAS size arguments for the operations that will be used.
 *   An operation is included in max-workspace calculation only
 *   when its required size variables are defined.
 *
 *   | Operation suffix | Required size variables                                    | Internal workspace query       |
 *   | :--------------- | :--------------------------------------------------------- | :----------------------------- |
 *   | GEMM             | GEMMUL8_MAX_M_GEMM, GEMMUL8_MAX_N_GEMM, GEMMUL8_MAX_K_GEMM | workSize(m, n, k, ...)         |
 *   | SYMM_LEFT        | GEMMUL8_MAX_M_SYMM_LEFT, GEMMUL8_MAX_N_SYMM_LEFT           | workSize(m, n, m, ...)         |
 *   | SYMM_RIGHT       | GEMMUL8_MAX_M_SYMM_RIGHT, GEMMUL8_MAX_N_SYMM_RIGHT         | workSize(m, n, n, ...)         |
 *   | SYRK             | GEMMUL8_MAX_N_SYRK, GEMMUL8_MAX_K_SYRK                     | workSize(n, n, k, ...)         |
 *   | SYR2K            | GEMMUL8_MAX_N_SYR2K, GEMMUL8_MAX_K_SYR2K                   | workSize(n, n, k, ...)         |
 *   | SYRKX            | GEMMUL8_MAX_N_SYRKX, GEMMUL8_MAX_K_SYRKX                   | workSize(n, n, k, ...)         |
 *   | HEMM_LEFT        | GEMMUL8_MAX_M_HEMM_LEFT, GEMMUL8_MAX_N_HEMM_LEFT           | workSize(m, n, m, ...)         |
 *   | HEMM_RIGHT       | GEMMUL8_MAX_M_HEMM_RIGHT, GEMMUL8_MAX_N_HEMM_RIGHT         | workSize(m, n, n, ...)         |
 *   | HERK             | GEMMUL8_MAX_N_HERK, GEMMUL8_MAX_K_HERK                     | workSize(n, n, k, ...)         |
 *   | HER2K            | GEMMUL8_MAX_N_HER2K, GEMMUL8_MAX_K_HER2K                   | workSize(n, n, k, ...)         |
 *   | HERKX            | GEMMUL8_MAX_N_HERKX, GEMMUL8_MAX_K_HERKX                   | workSize(n, n, k, ...)         |
 *   | TRMM_LEFT        | GEMMUL8_MAX_M_TRMM_LEFT, GEMMUL8_MAX_N_TRMM_LEFT           | workSize(m, n, m, ...)         |
 *   | TRMM_RIGHT       | GEMMUL8_MAX_M_TRMM_RIGHT, GEMMUL8_MAX_N_TRMM_RIGHT         | workSize(m, n, n, ...)         |
 *   | TRSM_LEFT        | GEMMUL8_MAX_M_TRSM_LEFT, GEMMUL8_MAX_N_TRSM_LEFT           | workSizeTrsm(LEFT, m, n, ...)  |
 *   | TRSM_RIGHT       | GEMMUL8_MAX_M_TRSM_RIGHT, GEMMUL8_MAX_N_TRSM_RIGHT         | workSizeTrsm(RIGHT, m, n, ...) |
 *
 *   Additional max-workspace variables:
 *
 *   | Variable pattern                         | Default | Description                                                           |
 *   | :--------------------------------------- | :------ | :-------------------------------------------------------------------- |
 *   | GEMMUL8_MAXWS_BACKEND_<OP>               | INT8    | Backend used for max-workspace calculation. 0/INT8, 1/FP8, or 2/BOTH. |
 *   | GEMMUL8_MAX_NUM_MOD_<OP>                 | 2       | Number of moduli used for max-workspace calculation.                  |
 *
 *   For GEMM only, the following old names are also accepted when the
 *   corresponding "_GEMM" variables are not defined:
 *
 *     GEMMUL8_BACKEND
 *     GEMMUL8_NUM_MOD_S, GEMMUL8_NUM_MOD_D, GEMMUL8_NUM_MOD_C, GEMMUL8_NUM_MOD_Z
 *     GEMMUL8_FASTMODE_S, GEMMUL8_FASTMODE_D, GEMMUL8_FASTMODE_C, GEMMUL8_FASTMODE_Z
 *     GEMMUL8_MAXWS_BACKEND
 *     GEMMUL8_MAX_M, GEMMUL8_MAX_N, GEMMUL8_MAX_K, GEMMUL8_MAX_NUM_MOD
 *
 * EXAMPLE:
 *
 *   # GEMM
 *   export GEMMUL8_BACKEND_GEMM=INT8
 *   export GEMMUL8_NUM_MOD_D_GEMM=15
 *   export GEMMUL8_FASTMODE_D_GEMM=1
 *   export GEMMUL8_MAXWS_BACKEND_GEMM=BOTH
 *   export GEMMUL8_MAX_M_GEMM=32768
 *   export GEMMUL8_MAX_N_GEMM=32768
 *   export GEMMUL8_MAX_K_GEMM=32768
 *   export GEMMUL8_MAX_NUM_MOD_GEMM=15
 *
 *   # SYMM with side == LEFT
 *   export GEMMUL8_BACKEND_SYMM_LEFT=INT8
 *   export GEMMUL8_NUM_MOD_D_SYMM_LEFT=15
 *   export GEMMUL8_FASTMODE_D_SYMM_LEFT=1
 *   export GEMMUL8_MAX_M_SYMM_LEFT=32768
 *   export GEMMUL8_MAX_N_SYMM_LEFT=32768
 *   export GEMMUL8_MAX_NUM_MOD_SYMM_LEFT=15
 *
 *   # SYRK
 *   export GEMMUL8_BACKEND_SYRK=INT8
 *   export GEMMUL8_NUM_MOD_D_SYRK=15
 *   export GEMMUL8_FASTMODE_D_SYRK=1
 *   export GEMMUL8_MAX_N_SYRK=32768
 *   export GEMMUL8_MAX_K_SYRK=32768
 *   export GEMMUL8_MAX_NUM_MOD_SYRK=15
 *
 *   # TRMM with side == RIGHT
 *   export GEMMUL8_BACKEND_TRMM_RIGHT=FP8
 *   export GEMMUL8_NUM_MOD_D_TRMM_RIGHT=12
 *   export GEMMUL8_FASTMODE_D_TRMM_RIGHT=1
 *   export GEMMUL8_MAX_M_TRMM_RIGHT=32768
 *   export GEMMUL8_MAX_N_TRMM_RIGHT=32768
 *   export GEMMUL8_MAX_NUM_MOD_TRMM_RIGHT=12
 *
 *   # TRSM with side == LEFT
 *   export GEMMUL8_BACKEND_TRSM_LEFT=INT8
 *   export GEMMUL8_NUM_MOD_D_TRSM_LEFT=10
 *   export GEMMUL8_FASTMODE_D_TRSM_LEFT=0
 *   export GEMMUL8_MAX_M_TRSM_LEFT=32768
 *   export GEMMUL8_MAX_N_TRSM_LEFT=32768
 *   export GEMMUL8_MAX_NUM_MOD_TRSM_LEFT=10
 *
 *   # Global memory-saving settings
 *   export GEMMUL8_MEMORY_SAVING=1
 *   export GEMMUL8_MAX_WORKSIZE=4294967296   # 4 GiB, in bytes
 *
 *   # Global skip-scaling switches
 *   export GEMMUL8_SKIP_SCALE_A=1
 *   export GEMMUL8_SKIP_SCALE_B=1
 *
 * NOTE:
 *
 *   - GEMMUL8_SKIP_SCALE_A and GEMMUL8_SKIP_SCALE_B are global switches.
 *     They are not operation-specific because the skip-scaling cache is
 *     keyed by operand metadata rather than by routine family.
 *
 *   - Skip scaling relies on pointer identity and operand metadata;
 *     matrix contents are not checked.
 *
 *   - Enable GEMMUL8_SKIP_SCALE_A/B only when the corresponding input data are
 *     unchanged between calls.
 *
 *   - GEMMUL8_MEMORY_SAVING and GEMMUL8_MAX_WORKSIZE override the corresponding
 *     handle-local settings only when the environment variables are explicitly defined.
 *
 *   - Handle-local memory saving disables skip-scaling/reuse even when
 *     GEMMUL8_SKIP_SCALE_A/B is enabled.
 *
 *   - GEMMUL8_BLK_SIZE_TRSM overrides the handle-local TRSM block size only
 *     when the environment variable is explicitly defined.
 *
 *   - GEMMUL8_MEMORY_SAVING, GEMMUL8_MAX_WORKSIZE, and GEMMUL8_BLK_SIZE_TRSM are read
 *     on every intercepted GEMMul8 call (GEMMUL8_BLK_SIZE_TRSM on every intercepted TRSM call).
 *     Changing their values during execution therefore takes effect from the next matching call.
 *     The value is copied to the handle-local setting; unsetting the environment variable
 *     does not restore an earlier value. Use 0 or the public setter to reset the setting.
 *
 *   - If an operand workspace is reallocated, the corresponding cache entry is
 *     invalidated automatically.
 */
#pragma once
#include "../include/gemmul8.hpp"
#include "../src/config/config.hpp"
#include "../src/oz2/common/include.hpp"
#include "../src/oz2/trsm/worksize.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#if !defined(_WIN32)
    #include <dlfcn.h>
#endif
#include <iostream>
#include <memory>
#include <mutex>
#include <string>
#include <unordered_map>

namespace gemmul8::hook {

inline size_t max_workSizeA = 0;
inline size_t max_workSizeB = 0;
inline size_t max_workSizeC = 0;
inline std::once_flag g_maxws_once;

template <typename T> inline constexpr int num_moduli_threshold       = 20;
template <> inline constexpr int num_moduli_threshold<float>          = 13;
template <> inline constexpr int num_moduli_threshold<cuFloatComplex> = 13;

enum class MaxWSBackend : unsigned {
    INT8 = 0,
    FP8  = 1,
    BOTH = 2,
};

enum class HookOp : unsigned {
    GEMM,
    SYMM_LEFT,
    SYMM_RIGHT,
    SYRK,
    SYR2K,
    SYRKX,
    HEMM_LEFT,
    HEMM_RIGHT,
    HERK,
    HER2K,
    HERKX,
    TRMM_LEFT,
    TRMM_RIGHT,
    TRSM_LEFT,
    TRSM_RIGHT,
};

static inline const char *hook_op_suffix(const HookOp op) {
    switch (op) {
    case HookOp::GEMM: return "GEMM";
    case HookOp::SYMM_LEFT: return "SYMM_LEFT";
    case HookOp::SYMM_RIGHT: return "SYMM_RIGHT";
    case HookOp::SYRK: return "SYRK";
    case HookOp::SYR2K: return "SYR2K";
    case HookOp::SYRKX: return "SYRKX";
    case HookOp::HEMM_LEFT: return "HEMM_LEFT";
    case HookOp::HEMM_RIGHT: return "HEMM_RIGHT";
    case HookOp::HERK: return "HERK";
    case HookOp::HER2K: return "HER2K";
    case HookOp::HERKX: return "HERKX";
    case HookOp::TRMM_LEFT: return "TRMM_LEFT";
    case HookOp::TRMM_RIGHT: return "TRMM_RIGHT";
    case HookOp::TRSM_LEFT: return "TRSM_LEFT";
    case HookOp::TRSM_RIGHT: return "TRSM_RIGHT";
    }

    return "GEMM";
}

static inline HookOp hook_op_from_symm_side(const cublasSideMode_t side) {
    return (side == CUBLAS_SIDE_LEFT) ? HookOp::SYMM_LEFT : HookOp::SYMM_RIGHT;
}

static inline HookOp hook_op_from_hemm_side(const cublasSideMode_t side) {
    return (side == CUBLAS_SIDE_LEFT) ? HookOp::HEMM_LEFT : HookOp::HEMM_RIGHT;
}

static inline HookOp hook_op_from_trmm_side(const cublasSideMode_t side) {
    return (side == CUBLAS_SIDE_LEFT) ? HookOp::TRMM_LEFT : HookOp::TRMM_RIGHT;
}

static inline HookOp hook_op_from_trsm_side(const cublasSideMode_t side) {
    return (side == CUBLAS_SIDE_LEFT) ? HookOp::TRSM_LEFT : HookOp::TRSM_RIGHT;
}

// ---- Default initialization values ----
namespace initial_vals {

inline constexpr Backend BACKEND = Backend::INT8;

inline constexpr size_t MAX_M               = 0u; // default M size
inline constexpr size_t MAX_N               = 0u; // default N size
inline constexpr size_t MAX_K               = 0u; // default K size
inline constexpr MaxWSBackend MaxWS_BACKEND = MaxWSBackend::INT8;

inline constexpr int MAX_NUM_MOD = 2; // default modulus count
inline constexpr int NUM_MOD_D   = 0; // default double moduli
inline constexpr int NUM_MOD_S   = 0; // default float moduli
inline constexpr int NUM_MOD_Z   = 0; // default double-complex moduli
inline constexpr int NUM_MOD_C   = 0; // default float-complex moduli

inline constexpr bool FASTMODE_D = true; // default double fastmode
inline constexpr bool FASTMODE_S = true; // default float fastmode
inline constexpr bool FASTMODE_Z = true; // default double-complex fastmode
inline constexpr bool FASTMODE_C = true; // default float-complex fastmode

inline constexpr int BLK_SIZE_TRSM = 0; // default TRSM block-size override

inline constexpr bool SCALE_A = false; // default skip_scalA_switch
inline constexpr bool SCALE_B = false; // default skip_scalB_switch

} // namespace initial_vals

inline thread_local unsigned hook_depth = 0;

struct HookGuard {
    HookGuard() { ++hook_depth; }
    ~HookGuard() { --hook_depth; }

    HookGuard(const HookGuard &)            = delete;
    HookGuard &operator=(const HookGuard &) = delete;
};

inline bool inside_hook() {
    return hook_depth != 0;
}

enum class ScalarKind : char {
    S = 'S',
    D = 'D',
    C = 'C',
    Z = 'Z',
};

enum class MatrixKind : unsigned char {
    General,
    Triangular,
    Symmetric,
    Hermitian,
};

template <typename T> inline constexpr ScalarKind scalar_kind_v        = ScalarKind::D;
template <> inline constexpr ScalarKind scalar_kind_v<float>           = ScalarKind::S;
template <> inline constexpr ScalarKind scalar_kind_v<double>          = ScalarKind::D;
template <> inline constexpr ScalarKind scalar_kind_v<cuFloatComplex>  = ScalarKind::C;
template <> inline constexpr ScalarKind scalar_kind_v<cuDoubleComplex> = ScalarKind::Z;

struct ScaledOperandKey {
    const void *ptr = nullptr;

    size_t rows = 0;
    size_t cols = 0;
    size_t ld   = 0;

    cublasOperation_t op  = CUBLAS_OP_N;
    cublasFillMode_t uplo = CUBLAS_FILL_MODE_FULL;
    cublasDiagType_t diag = CUBLAS_DIAG_NON_UNIT;

    MatrixKind matrix_kind = MatrixKind::General;
    ScalarKind scalar_kind = ScalarKind::D;

    Backend backend = Backend::INT8;
    int num_moduli  = 2;
    bool fastmode   = false;

    friend bool operator==(const ScaledOperandKey &a, const ScaledOperandKey &b) {
        return a.ptr == b.ptr &&
               a.rows == b.rows &&
               a.cols == b.cols &&
               a.ld == b.ld &&
               a.op == b.op &&
               a.uplo == b.uplo &&
               a.diag == b.diag &&
               a.matrix_kind == b.matrix_kind &&
               a.scalar_kind == b.scalar_kind &&
               a.backend == b.backend &&
               a.num_moduli == b.num_moduli &&
               a.fastmode == b.fastmode;
    }
};

struct ScaledOperandEntry {
    bool valid = false;
    ScaledOperandKey key{};
    int8_t *work = nullptr;

    bool matches(const ScaledOperandKey &cur_key, const int8_t *cur_work) const {
        return valid && cur_work != nullptr && work == cur_work && key == cur_key;
    }
};

struct HandleState {
    std::mutex mtx;

    cublasHandle_t parent = nullptr;

    void *workA       = nullptr;
    size_t workA_size = 0;

    void *workB       = nullptr;
    size_t workB_size = 0;

    void *workC       = nullptr;
    size_t workC_size = 0;

    cublasLtHandle_t lt = nullptr;

    std::array<ScaledOperandEntry, 2> scaled_operands{};
    unsigned next_scaled_slot = 0;

    cudaStream_t last_stream = nullptr;
    cudaEvent_t last_event   = nullptr;
    bool last_stream_valid   = false;
};

inline std::mutex g_state_map_mtx;
inline std::unordered_map<cublasHandle_t, std::shared_ptr<HandleState>> g_state_map;

inline bool can_skip_scaled_operand_locked(
    const HandleState &st,
    const ScaledOperandKey &key,
    const int8_t *work,
    bool enable_skip //
) {
    if (!enable_skip) return false;

    for (const auto &entry : st.scaled_operands) {
        if (entry.matches(key, work)) return true;
    }

    return false;
}

inline void update_scaled_operand_locked(
    HandleState &st,
    const ScaledOperandKey &key,
    int8_t *work //
) {
    if (!work) return;

    for (auto &entry : st.scaled_operands) {
        if (entry.valid && entry.work == work) {
            entry.key = key;
            return;
        }
    }

    auto &entry = st.scaled_operands[st.next_scaled_slot % st.scaled_operands.size()];
    entry.valid = true;
    entry.key   = key;
    entry.work  = work;

    ++st.next_scaled_slot;
}

inline void invalidate_scaled_operand_work_locked(
    HandleState &st,
    const void *work //
) {
    for (auto &entry : st.scaled_operands) {
        if (entry.valid && entry.work == work) {
            entry = ScaledOperandEntry{};
        }
    }
}

inline void clear_scaled_operands_locked(HandleState &st) {
    for (auto &entry : st.scaled_operands) {
        entry = ScaledOperandEntry{};
    }
    st.next_scaled_slot = 0;
}

static inline std::shared_ptr<HandleState> get_state(cublasHandle_t h) {
    std::lock_guard<std::mutex> g(g_state_map_mtx);
    auto &p = g_state_map[h];
    if (!p) {
        p         = std::make_shared<HandleState>();
        p->parent = h;
    }
    return p;
}

static inline std::shared_ptr<HandleState> find_state(cublasHandle_t h) {
    std::lock_guard<std::mutex> g(g_state_map_mtx);
    const auto it = g_state_map.find(h);
    return (it == g_state_map.end()) ? nullptr : it->second;
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

// ---- Small helpers ----
static inline bool env_is_one(const char *s, bool def) {
    return s ? (std::strcmp(s, "1") == 0) : def;
}

static inline bool env_bool01(const char *s, bool def) {
    if (!s) return def;
    if (std::strcmp(s, "1") == 0) return true;
    if (std::strcmp(s, "0") == 0) return false;
    return def;
}

static inline int env_i32(const char *s, int def) {
    if (!s) return def;
    try {
        return std::stoi(s);
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

static inline const char *getenv_op(const char *base, const HookOp op) {
    const std::string name = std::string(base) + "_" + hook_op_suffix(op);

    if (const char *v = std::getenv(name.c_str())) {
        return v;
    }

    // Backward compatibility for GEMM only.
    if (op == HookOp::GEMM) {
        return std::getenv(base);
    }

    return nullptr;
}

static inline bool has_env_op(const char *base, const HookOp op) {
    const std::string name = std::string(base) + "_" + hook_op_suffix(op);

    if (std::getenv(name.c_str())) {
        return true;
    }

    // Backward compatibility for GEMM only.
    if (op == HookOp::GEMM) {
        return std::getenv(base) != nullptr;
    }

    return false;
}

static inline bool has_maxws_mn(const HookOp op) {
    return has_env_op("GEMMUL8_MAX_M", op) &&
           has_env_op("GEMMUL8_MAX_N", op);
}

static inline bool has_maxws_nk(const HookOp op) {
    return has_env_op("GEMMUL8_MAX_N", op) &&
           has_env_op("GEMMUL8_MAX_K", op);
}

static inline bool has_maxws_mnk(const HookOp op) {
    return has_env_op("GEMMUL8_MAX_M", op) &&
           has_env_op("GEMMUL8_MAX_N", op) &&
           has_env_op("GEMMUL8_MAX_K", op);
}

static inline bool has_maxws_config(const HookOp op) {
    switch (op) {
    case HookOp::GEMM:
        return has_maxws_mnk(op);

    case HookOp::SYMM_LEFT:
    case HookOp::SYMM_RIGHT:
    case HookOp::HEMM_LEFT:
    case HookOp::HEMM_RIGHT:
    case HookOp::TRMM_LEFT:
    case HookOp::TRMM_RIGHT:
    case HookOp::TRSM_LEFT:
    case HookOp::TRSM_RIGHT:
        return has_maxws_mn(op);

    case HookOp::SYRK:
    case HookOp::SYR2K:
    case HookOp::SYRKX:
    case HookOp::HERK:
    case HookOp::HER2K:
    case HookOp::HERKX:
        return has_maxws_nk(op);
    }

    return false;
}

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

static inline MaxWSBackend requested_maxws_backend(const HookOp op) {
    return env_maxws_backend(getenv_op("GEMMUL8_MAXWS_BACKEND", op), initial_vals::MaxWS_BACKEND);
}

static inline size_t requested_max_m(const HookOp op) {
    return env_u64(getenv_op("GEMMUL8_MAX_M", op), initial_vals::MAX_M);
}

static inline size_t requested_max_n(const HookOp op) {
    return env_u64(getenv_op("GEMMUL8_MAX_N", op), initial_vals::MAX_N);
}

static inline size_t requested_max_k(const HookOp op) {
    return env_u64(getenv_op("GEMMUL8_MAX_K", op), initial_vals::MAX_K);
}

static inline int requested_max_num_mod(const HookOp op) {
    return env_i32(getenv_op("GEMMUL8_MAX_NUM_MOD", op), initial_vals::MAX_NUM_MOD);
}

static inline Backend env_backend(const char *s, Backend def) {
    if (!s) return def;
    if (std::strcmp(s, "0") == 0) return Backend::INT8;
    if (std::strcmp(s, "1") == 0) return Backend::FP8;
    if (std::strcmp(s, "INT8") == 0) return Backend::INT8;
    if (std::strcmp(s, "FP8") == 0) return Backend::FP8;
    return def;
}

static inline Backend requested_backend(const HookOp op) {
    return env_backend(getenv_op("GEMMUL8_BACKEND", op), initial_vals::BACKEND);
}

static inline bool has_requested_block_size_trsm() {
    return std::getenv("GEMMUL8_BLK_SIZE_TRSM") != nullptr;
}

static inline int requested_block_size_trsm() {
    const int nB = env_i32(std::getenv("GEMMUL8_BLK_SIZE_TRSM"), initial_vals::BLK_SIZE_TRSM);
    return (nB > 0) ? nB : 0;
}

// Apply the TRSM environment override to the handle
static inline void apply_block_size_trsm_environment(cublasHandle_t handle) {
    if (!has_requested_block_size_trsm()) return;
    gemmul8::set_block_size_trsm(handle, requested_block_size_trsm());
}

// Apply workspace-related environment overrides to the handle
static inline gemmul8::config::ConfigSnapshot apply_workspace_environment(cublasHandle_t handle) {
    auto config = gemmul8::config::get_config(handle);

    if (const char *s = std::getenv("GEMMUL8_MEMORY_SAVING")) {
        const bool value = env_bool01(s, config.memory_saving);
        if (value != config.memory_saving) {
            gemmul8::set_memory_saving(handle, value);
            config.memory_saving = value;
        }
    }

    if (const char *s = std::getenv("GEMMUL8_MAX_WORKSIZE")) {
        const size_t value = env_u64(s, config.max_worksize);
        if (value != config.max_worksize) {
            gemmul8::set_max_worksize(handle, value);
            config.max_worksize = value;
        }
    }

    return config;
}

struct WorkspaceConfig {
    bool memory_saving_mode = false;
    bool memory_saving      = false;
    size_t limit            = 0;
    bool enable_skipA       = false;
    bool enable_skipB       = false;
};

static inline WorkspaceConfig workspace_config(
    cublasHandle_t handle,
    bool requested_skipA,
    bool requested_skipB //
) {
    const auto config = apply_workspace_environment(handle);

    WorkspaceConfig out;
    out.memory_saving_mode = config.memory_saving;
    out.memory_saving      = config.memory_saving && config.max_worksize > 0;
    out.limit              = config.max_worksize;
    out.enable_skipA       = config.memory_saving ? false : requested_skipA;
    out.enable_skipB       = config.memory_saving ? false : requested_skipB;
    return out;
}

struct WorkspaceRequest {
    size_t workA = 0;
    size_t workB = 0;
    size_t workC = 0;
};

static inline WorkspaceRequest workspace_request(
    const WorkspaceConfig &config,
    size_t total,
    size_t workA,
    size_t workB //
) {
    WorkspaceRequest out;

    if (config.memory_saving && total > config.limit) {
        out.workC = config.limit;
        return out;
    }

    out.workA = workA;
    out.workB = workB;
    out.workC = (total >= workA + workB) ? (total - workA - workB) : 0;

    if (!config.memory_saving_mode) {
        if (config.enable_skipA) out.workA = std::max(out.workA, max_workSizeA);
        if (config.enable_skipB) out.workB = std::max(out.workB, max_workSizeB);
        if (config.enable_skipA || config.enable_skipB) {
            out.workC = std::max(out.workC, max_workSizeC);
        }
    }

    return out;
}

template <bool COMPLEX, Backend BACKEND, Func FUNC>
static inline void update_maxws_regular(
    size_t m, size_t n, size_t k, int num_moduli,
    bool enable_skipA, bool enable_skipB,
    size_t &maxA, size_t &maxB, size_t &maxC //
) {
    size_t wA = 0;
    size_t wB = 0;

    const size_t w = workSize<COMPLEX, BACKEND, FUNC>(
        m, n, k, num_moduli, enable_skipA, enable_skipB, &wA, &wB);

    const size_t wC = (w > (wA + wB)) ? (w - wA - wB) : 0;

    maxA = std::max(maxA, wA);
    maxB = std::max(maxB, wB);
    maxC = std::max(maxC, wC);
}

template <typename T, Backend BACKEND>
static inline void update_maxws_trsm(
    cublasSideMode_t side, size_t m, size_t n, int num_moduli,
    size_t &maxC //
) {
    const int nB_override = requested_block_size_trsm();

    const size_t w = (side == CUBLAS_SIDE_LEFT)
                         ? gemmul8::oz2::trsm::workSize_left<T, BACKEND>(
                               m, n, static_cast<unsigned>(num_moduli), nB_override)
                         : gemmul8::oz2::trsm::workSize_right<T, BACKEND>(
                               m, n, static_cast<unsigned>(num_moduli), nB_override);

    maxC = std::max(maxC, w);
}

template <bool COMPLEX, Backend BACKEND>
static inline void update_maxws_for_hook_op_typed(
    HookOp op,
    size_t &maxA,
    size_t &maxB,
    size_t &maxC //
) {
    const int num_moduli = requested_max_num_mod(op);

    switch (op) {
    case HookOp::GEMM: {
        const size_t m = requested_max_m(op);
        const size_t n = requested_max_n(op);
        const size_t k = requested_max_k(op);

        update_maxws_regular<COMPLEX, BACKEND, Func::gemm>(
            m, n, k, num_moduli,
            true, true,
            maxA, maxB, maxC);
        break;
    }

    case HookOp::SYMM_LEFT: {
        const size_t m = requested_max_m(op);
        const size_t n = requested_max_n(op);

        update_maxws_regular<COMPLEX, BACKEND, Func::symm>(
            m, n, m, num_moduli,
            true, true,
            maxA, maxB, maxC);
        break;
    }

    case HookOp::SYMM_RIGHT: {
        const size_t m = requested_max_m(op);
        const size_t n = requested_max_n(op);

        update_maxws_regular<COMPLEX, BACKEND, Func::symm>(
            m, n, n, num_moduli,
            true, true,
            maxA, maxB, maxC);
        break;
    }

    case HookOp::SYRK: {
        const size_t n = requested_max_n(op);
        const size_t k = requested_max_k(op);

        update_maxws_regular<COMPLEX, BACKEND, Func::syrk>(
            n, n, k, num_moduli,
            true, false,
            maxA, maxB, maxC);
        break;
    }

    case HookOp::SYR2K: {
        const size_t n = requested_max_n(op);
        const size_t k = requested_max_k(op);

        update_maxws_regular<COMPLEX, BACKEND, Func::syr2k>(
            n, n, k, num_moduli,
            true, true,
            maxA, maxB, maxC);
        break;
    }

    case HookOp::SYRKX: {
        const size_t n = requested_max_n(op);
        const size_t k = requested_max_k(op);

        update_maxws_regular<COMPLEX, BACKEND, Func::syrkx>(
            n, n, k, num_moduli,
            true, true,
            maxA, maxB, maxC);
        break;
    }

    case HookOp::HEMM_LEFT: {
        if constexpr (COMPLEX) {
            const size_t m = requested_max_m(op);
            const size_t n = requested_max_n(op);

            update_maxws_regular<true, BACKEND, Func::hemm>(
                m, n, m, num_moduli,
                true, true,
                maxA, maxB, maxC);
        }
        break;
    }

    case HookOp::HEMM_RIGHT: {
        if constexpr (COMPLEX) {
            const size_t m = requested_max_m(op);
            const size_t n = requested_max_n(op);

            update_maxws_regular<true, BACKEND, Func::hemm>(
                m, n, n, num_moduli,
                true, true,
                maxA, maxB, maxC);
        }
        break;
    }

    case HookOp::HERK: {
        if constexpr (COMPLEX) {
            const size_t n = requested_max_n(op);
            const size_t k = requested_max_k(op);

            update_maxws_regular<true, BACKEND, Func::herk>(
                n, n, k, num_moduli,
                true, false,
                maxA, maxB, maxC);
        }
        break;
    }

    case HookOp::HER2K: {
        if constexpr (COMPLEX) {
            const size_t n = requested_max_n(op);
            const size_t k = requested_max_k(op);

            update_maxws_regular<true, BACKEND, Func::her2k>(
                n, n, k, num_moduli,
                true, true,
                maxA, maxB, maxC);
        }
        break;
    }

    case HookOp::HERKX: {
        if constexpr (COMPLEX) {
            const size_t n = requested_max_n(op);
            const size_t k = requested_max_k(op);

            update_maxws_regular<true, BACKEND, Func::herkx>(
                n, n, k, num_moduli,
                true, true,
                maxA, maxB, maxC);
        }
        break;
    }

    case HookOp::TRMM_LEFT: {
        const size_t m = requested_max_m(op);
        const size_t n = requested_max_n(op);

        update_maxws_regular<COMPLEX, BACKEND, Func::trmm>(
            m, n, m, num_moduli,
            true, true,
            maxA, maxB, maxC);
        break;
    }

    case HookOp::TRMM_RIGHT: {
        const size_t m = requested_max_m(op);
        const size_t n = requested_max_n(op);

        update_maxws_regular<COMPLEX, BACKEND, Func::trmm>(
            m, n, n, num_moduli,
            true, true,
            maxA, maxB, maxC);
        break;
    }

    case HookOp::TRSM_LEFT: {
        const size_t m = requested_max_m(op);
        const size_t n = requested_max_n(op);

        if constexpr (COMPLEX) {
            update_maxws_trsm<cuDoubleComplex, BACKEND>(
                CUBLAS_SIDE_LEFT, m, n, num_moduli, maxC);
        } else {
            update_maxws_trsm<double, BACKEND>(
                CUBLAS_SIDE_LEFT, m, n, num_moduli, maxC);
        }
        break;
    }

    case HookOp::TRSM_RIGHT: {
        const size_t m = requested_max_m(op);
        const size_t n = requested_max_n(op);

        if constexpr (COMPLEX) {
            update_maxws_trsm<cuDoubleComplex, BACKEND>(
                CUBLAS_SIDE_RIGHT, m, n, num_moduli, maxC);
        } else {
            update_maxws_trsm<double, BACKEND>(
                CUBLAS_SIDE_RIGHT, m, n, num_moduli, maxC);
        }
        break;
    }
    }
}

static inline bool hook_op_is_hermitian_family(const HookOp op) {
    return op == HookOp::HEMM_LEFT ||
           op == HookOp::HEMM_RIGHT ||
           op == HookOp::HERK ||
           op == HookOp::HER2K ||
           op == HookOp::HERKX;
}

static inline void update_maxws_for_one_op(HookOp op, size_t &maxA, size_t &maxB, size_t &maxC) {
    if (!has_maxws_config(op)) {
        return;
    }

    const MaxWSBackend mws = requested_maxws_backend(op);

    const bool do_int8 = (mws == MaxWSBackend::INT8) || (mws == MaxWSBackend::BOTH);
    const bool do_fp8  = (mws == MaxWSBackend::FP8) || (mws == MaxWSBackend::BOTH);

    const int nmod_c = env_i32(
        getenv_op("GEMMUL8_NUM_MOD_C", op),
        initial_vals::NUM_MOD_C);

    const int nmod_z = env_i32(
        getenv_op("GEMMUL8_NUM_MOD_Z", op),
        initial_vals::NUM_MOD_Z);

    const bool want_complex =
        hook_op_is_hermitian_family(op) || (nmod_c > 0) || (nmod_z > 0);

    if (do_int8) {
        update_maxws_for_hook_op_typed<false, Backend::INT8>(
            op, maxA, maxB, maxC);

        if (want_complex) {
            update_maxws_for_hook_op_typed<true, Backend::INT8>(
                op, maxA, maxB, maxC);
        }
    }

    if (do_fp8) {
        update_maxws_for_hook_op_typed<false, Backend::FP8>(
            op, maxA, maxB, maxC);

        if (want_complex) {
            update_maxws_for_hook_op_typed<true, Backend::FP8>(
                op, maxA, maxB, maxC);
        }
    }
}

static void init_max_workspace() {
    std::call_once(g_maxws_once, [] {
        size_t maxA = 0;
        size_t maxB = 0;
        size_t maxC = 0;

        update_maxws_for_one_op(HookOp::GEMM, maxA, maxB, maxC);
        update_maxws_for_one_op(HookOp::SYMM_LEFT, maxA, maxB, maxC);
        update_maxws_for_one_op(HookOp::SYMM_RIGHT, maxA, maxB, maxC);
        update_maxws_for_one_op(HookOp::SYRK, maxA, maxB, maxC);
        update_maxws_for_one_op(HookOp::SYR2K, maxA, maxB, maxC);
        update_maxws_for_one_op(HookOp::SYRKX, maxA, maxB, maxC);
        update_maxws_for_one_op(HookOp::HEMM_LEFT, maxA, maxB, maxC);
        update_maxws_for_one_op(HookOp::HEMM_RIGHT, maxA, maxB, maxC);
        update_maxws_for_one_op(HookOp::HERK, maxA, maxB, maxC);
        update_maxws_for_one_op(HookOp::HER2K, maxA, maxB, maxC);
        update_maxws_for_one_op(HookOp::HERKX, maxA, maxB, maxC);
        update_maxws_for_one_op(HookOp::TRMM_LEFT, maxA, maxB, maxC);
        update_maxws_for_one_op(HookOp::TRMM_RIGHT, maxA, maxB, maxC);
        update_maxws_for_one_op(HookOp::TRSM_LEFT, maxA, maxB, maxC);
        update_maxws_for_one_op(HookOp::TRSM_RIGHT, maxA, maxB, maxC);

        max_workSizeA = maxA;
        max_workSizeB = maxB;
        max_workSizeC = maxC;
    });
}

// ---- Environment per type ----
static inline void get_env_d(
    const HookOp op,
    int &num_moduli, bool &fastmode,
    bool &enable_skipA, bool &enable_skipB //
) {
    num_moduli   = env_i32(getenv_op("GEMMUL8_NUM_MOD_D", op), initial_vals::NUM_MOD_D);
    fastmode     = env_is_one(getenv_op("GEMMUL8_FASTMODE_D", op), initial_vals::FASTMODE_D);
    enable_skipA = env_is_one(std::getenv("GEMMUL8_SKIP_SCALE_A"), initial_vals::SCALE_A);
    enable_skipB = env_is_one(std::getenv("GEMMUL8_SKIP_SCALE_B"), initial_vals::SCALE_B);
}

static inline void get_env_s(
    const HookOp op,
    int &num_moduli, bool &fastmode,
    bool &enable_skipA, bool &enable_skipB //
) {
    num_moduli   = env_i32(getenv_op("GEMMUL8_NUM_MOD_S", op), initial_vals::NUM_MOD_S);
    fastmode     = env_is_one(getenv_op("GEMMUL8_FASTMODE_S", op), initial_vals::FASTMODE_S);
    enable_skipA = env_is_one(std::getenv("GEMMUL8_SKIP_SCALE_A"), initial_vals::SCALE_A);
    enable_skipB = env_is_one(std::getenv("GEMMUL8_SKIP_SCALE_B"), initial_vals::SCALE_B);
}

static inline void get_env_z(
    const HookOp op,
    int &num_moduli, bool &fastmode,
    bool &enable_skipA, bool &enable_skipB //
) {
    num_moduli   = env_i32(getenv_op("GEMMUL8_NUM_MOD_Z", op), initial_vals::NUM_MOD_Z);
    fastmode     = env_is_one(getenv_op("GEMMUL8_FASTMODE_Z", op), initial_vals::FASTMODE_Z);
    enable_skipA = env_is_one(std::getenv("GEMMUL8_SKIP_SCALE_A"), initial_vals::SCALE_A);
    enable_skipB = env_is_one(std::getenv("GEMMUL8_SKIP_SCALE_B"), initial_vals::SCALE_B);
}

static inline void get_env_c(
    const HookOp op,
    int &num_moduli, bool &fastmode,
    bool &enable_skipA, bool &enable_skipB //
) {
    num_moduli   = env_i32(getenv_op("GEMMUL8_NUM_MOD_C", op), initial_vals::NUM_MOD_C);
    fastmode     = env_is_one(getenv_op("GEMMUL8_FASTMODE_C", op), initial_vals::FASTMODE_C);
    enable_skipA = env_is_one(std::getenv("GEMMUL8_SKIP_SCALE_A"), initial_vals::SCALE_A);
    enable_skipB = env_is_one(std::getenv("GEMMUL8_SKIP_SCALE_B"), initial_vals::SCALE_B);
}

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
    gemmul8::config::bind_config(st.parent, lt);
    *out = lt;
    return CUBLAS_STATUS_SUCCESS;
}

static inline cublasStatus_t release_work_locked(
    HandleState &st,
    void *&buf,
    size_t &buf_size,
    const char *tag,
    cudaStream_t stream //
) {
    if (!buf) {
        buf_size = 0;
        return CUBLAS_STATUS_SUCCESS;
    }

    void *old                  = buf;
    const cudaError_t free_err = cudaFreeAsync(old, stream);
    if (free_err != cudaSuccess) {
        std::cerr << "[GEMMUL8 HOOK] cudaFreeAsync failed for "
                  << (tag ? tag : "workspace")
                  << " (" << cudaGetErrorString(free_err) << ")\n";
        return CUBLAS_STATUS_INTERNAL_ERROR;
    }

    invalidate_scaled_operand_work_locked(st, old);
    buf      = nullptr;
    buf_size = 0;

    const cudaError_t sync_err = cudaStreamSynchronize(stream);
    if (sync_err != cudaSuccess) {
        std::cerr << "[GEMMUL8 HOOK] cudaStreamSynchronize failed after freeing "
                  << (tag ? tag : "workspace")
                  << " (" << cudaGetErrorString(sync_err) << ")\n";
        return CUBLAS_STATUS_INTERNAL_ERROR;
    }

    return CUBLAS_STATUS_SUCCESS;
}

static inline cublasStatus_t prepare_memory_saving_workspaces_locked(
    HandleState &st,
    const WorkspaceRequest &req,
    size_t limit,
    cudaStream_t stream //
) {
    if (limit == 0) return CUBLAS_STATUS_SUCCESS;

    const size_t projected =
        std::max(st.workA_size, req.workA) +
        std::max(st.workB_size, req.workB) +
        std::max(st.workC_size, req.workC);

    if (projected <= limit) return CUBLAS_STATUS_SUCCESS;

    cublasStatus_t st_release = CUBLAS_STATUS_SUCCESS;

    if (st.workA_size > req.workA) {
        st_release = release_work_locked(st, st.workA, st.workA_size, "workA", stream);
        if (st_release != CUBLAS_STATUS_SUCCESS) return st_release;
    }
    if (st.workB_size > req.workB) {
        st_release = release_work_locked(st, st.workB, st.workB_size, "workB", stream);
        if (st_release != CUBLAS_STATUS_SUCCESS) return st_release;
    }
    if (st.workC_size > req.workC) {
        st_release = release_work_locked(st, st.workC, st.workC_size, "workC", stream);
        if (st_release != CUBLAS_STATUS_SUCCESS) return st_release;
    }

    return CUBLAS_STATUS_SUCCESS;
}

static inline cublasStatus_t get_work_locked(
    HandleState &st,
    void *&buf,
    size_t &buf_size,
    size_t req_size,
    void **out,
    const char *tag,
    cudaStream_t stream //
) {

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

    if (buf) {
        const cublasStatus_t st_release = release_work_locked(st, buf, buf_size, tag, stream);
        if (st_release != CUBLAS_STATUS_SUCCESS) {
            *out = nullptr;
            return st_release;
        }
    }

    void *newp      = nullptr;
    cudaError_t err = cudaMallocAsync(&newp, req_size, stream);
    if (err != cudaSuccess) {
        *out = nullptr;
        std::cerr << "[GEMMUL8 HOOK] cudaMallocAsync failed for " << (tag ? tag : "workspace")
                  << " size " << req_size << " bytes. Error: " << cudaGetErrorString(err) << "\n";
        return CUBLAS_STATUS_ALLOC_FAILED;
    }

    buf      = newp;
    buf_size = req_size;
    *out     = buf;
    return CUBLAS_STATUS_SUCCESS;
}

static inline void cleanup_work(cublasHandle_t handle) {
    auto sp = find_state(handle);
    if (!sp) return;

    cudaStream_t stream = nullptr;
    bool have_stream    = false;

    {
        std::lock_guard<std::mutex> lk(sp->mtx);
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
            gemmul8::config::unbind_configLt(sp->lt);
            (void)cublasLtDestroy(sp->lt);
            sp->lt = nullptr;
        }

        if (sp->last_event) {
            cudaEventDestroy(sp->last_event);
            sp->last_event = nullptr;
        }

        sp->last_stream       = nullptr;
        sp->last_stream_valid = false;

        clear_scaled_operands_locked(*sp);
    }

    erase_state(handle);
}

// ---- dlsym loader ----
template <typename Fn>
static inline Fn load_real(const char *sym) {
    return reinterpret_cast<Fn>(dlsym(RTLD_NEXT, sym));
}

} // namespace gemmul8::hook
