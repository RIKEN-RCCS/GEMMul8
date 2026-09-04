/**
 * GEMM hook
 * ---------
 * Hook targets:
 *   - cublas{S,D,C,Z}gemm_v2
 *   - cublas{S,D,C,Z}gemm_v2_64
 *   - cublasGemmEx
 *   - cublasGemmEx_64
 *
 * CUDA-only hook targets:
 *   - cublas{C,Z}gemm3m
 *   - cublas{C,Z}gemm3m_64
 *   - cublasSgemmEx
 *   - cublasSgemmEx_64
 *   - cublasCgemmEx
 *   - cublasCgemmEx_64
 *   - cublasCgemm3mEx
 *   - cublasCgemm3mEx_64
 *
 * Notes:
 *   - cublas{S,D,C,Z}gemm and cublas{S,D,C,Z}gemm_64 are not defined here
 *     because the v2 variants are used when v2 exists.
 */
#include "common.hpp"

namespace {

template <typename T> struct HookScalarIsComplex {
    static constexpr bool value = false;
};
template <> struct HookScalarIsComplex<cuFloatComplex> {
    static constexpr bool value = true;
};
template <> struct HookScalarIsComplex<cuDoubleComplex> {
    static constexpr bool value = true;
};

static inline bool is_valid_gemm_op(const cublasOperation_t op) {
    return op == CUBLAS_OP_N || op == CUBLAS_OP_T || op == CUBLAS_OP_C;
}

static inline int64_t gemm_lda_min(const cublasOperation_t transa, const int64_t m, const int64_t k) {
    return std::max<int64_t>(1, (transa == CUBLAS_OP_N) ? m : k);
}

static inline int64_t gemm_ldb_min(const cublasOperation_t transb, const int64_t n, const int64_t k) {
    return std::max<int64_t>(1, (transb == CUBLAS_OP_N) ? k : n);
}

static inline int64_t gemm_ldc_min(const int64_t m) {
    return std::max<int64_t>(1, m);
}

static inline cublasStatus_t validate_gemm_params(
    const cublasHandle_t handle,
    const cublasOperation_t transa,
    const cublasOperation_t transb,
    const int64_t m, const int64_t n, const int64_t k,
    const void *alpha,
    const void *A, const int64_t lda,
    const void *B, const int64_t ldb,
    const void *beta,
    const void *C, const int64_t ldc //
) {
    if (!handle) return CUBLAS_STATUS_NOT_INITIALIZED;
    if (!is_valid_gemm_op(transa) || !is_valid_gemm_op(transb)) return CUBLAS_STATUS_INVALID_VALUE;
    if (!alpha || !beta || !A || !B || !C) return CUBLAS_STATUS_INVALID_VALUE;
    if (m < 0 || n < 0 || k < 0) return CUBLAS_STATUS_INVALID_VALUE;
    if (lda < gemm_lda_min(transa, m, k)) return CUBLAS_STATUS_INVALID_VALUE;
    if (ldb < gemm_ldb_min(transb, n, k)) return CUBLAS_STATUS_INVALID_VALUE;
    if (ldc < gemm_ldc_min(m)) return CUBLAS_STATUS_INVALID_VALUE;

    return CUBLAS_STATUS_SUCCESS;
}

template <typename T>
static inline gemmul8::hook::ScaledOperandKey make_gemm_A_key(
    const void *A,
    int64_t m,
    int64_t k,
    int64_t lda,
    cublasOperation_t transa,
    int num_moduli,
    bool fastmode,
    gemmul8::Backend backend) {
    return gemmul8::hook::ScaledOperandKey{
        .ptr         = A,
        .rows        = static_cast<size_t>(m),
        .cols        = static_cast<size_t>(k),
        .ld          = static_cast<size_t>(lda),
        .op          = transa,
        .uplo        = CUBLAS_FILL_MODE_FULL,
        .diag        = CUBLAS_DIAG_NON_UNIT,
        .matrix_kind = gemmul8::hook::MatrixKind::General,
        .scalar_kind = gemmul8::hook::scalar_kind_v<T>,
        .backend     = backend,
        .num_moduli  = num_moduli,
        .fastmode    = fastmode,
    };
}

template <typename T>
static inline gemmul8::hook::ScaledOperandKey make_gemm_B_key(
    const void *B,
    int64_t n,
    int64_t k,
    int64_t ldb,
    cublasOperation_t transb,
    int num_moduli,
    bool fastmode,
    gemmul8::Backend backend) {
    return gemmul8::hook::ScaledOperandKey{
        .ptr         = B,
        .rows        = static_cast<size_t>(k),
        .cols        = static_cast<size_t>(n),
        .ld          = static_cast<size_t>(ldb),
        .op          = transb,
        .uplo        = CUBLAS_FILL_MODE_FULL,
        .diag        = CUBLAS_DIAG_NON_UNIT,
        .matrix_kind = gemmul8::hook::MatrixKind::General,
        .scalar_kind = gemmul8::hook::scalar_kind_v<T>,
        .backend     = backend,
        .num_moduli  = num_moduli,
        .fastmode    = fastmode,
    };
}

// ---- Traits to unify duplicated GEMM bodies ----
template <typename T> struct GemmTraits;

#define DEFINE_GEMM_TRAITS(TYPE, GETENV_FN, SYM, SYM64)                                                 \
    template <> struct GemmTraits<TYPE> {                                                               \
        static constexpr bool isComplex = HookScalarIsComplex<TYPE>::value;                             \
        static inline void get_env(gemmul8::hook::HookOp op, int &nm, bool &fm, bool &enA, bool &enB) { \
            gemmul8::hook::GETENV_FN(op, nm, fm, enA, enB);                                             \
        }                                                                                               \
        using GemmV2Fn = cublasStatus_t (*)(                                                            \
            cublasHandle_t, cublasOperation_t, cublasOperation_t,                                       \
            int, int, int,                                                                              \
            const TYPE *, const TYPE *, int,                                                            \
            const TYPE *, int,                                                                          \
            const TYPE *, TYPE *, int);                                                                 \
        using GemmV2Fn64 = cublasStatus_t (*)(                                                          \
            cublasHandle_t, cublasOperation_t, cublasOperation_t,                                       \
            int64_t, int64_t, int64_t,                                                                  \
            const TYPE *, const TYPE *, int64_t,                                                        \
            const TYPE *, int64_t,                                                                      \
            const TYPE *, TYPE *, int64_t);                                                             \
        static constexpr const char *gemm_v2_sym    = STR(SYM);                                         \
        static constexpr const char *gemm_v2_64_sym = STR(SYM64);                                       \
    };

DEFINE_GEMM_TRAITS(float, get_env_s, cublasSgemm_v2, cublasSgemm_v2_64)
DEFINE_GEMM_TRAITS(double, get_env_d, cublasDgemm_v2, cublasDgemm_v2_64)
DEFINE_GEMM_TRAITS(cuFloatComplex, get_env_c, cublasCgemm_v2, cublasCgemm_v2_64)
DEFINE_GEMM_TRAITS(cuDoubleComplex, get_env_z, cublasZgemm_v2, cublasZgemm_v2_64)

#undef DEFINE_GEMM_TRAITS

template <typename T>
static inline typename GemmTraits<T>::GemmV2Fn real_gemm_v2() {
    using Fn     = typename GemmTraits<T>::GemmV2Fn;
    static Fn fn = gemmul8::hook::load_real<Fn>(GemmTraits<T>::gemm_v2_sym);
    return fn;
}

template <typename T>
static inline typename GemmTraits<T>::GemmV2Fn64 real_gemm_v2_64() {
    using Fn     = typename GemmTraits<T>::GemmV2Fn64;
    static Fn fn = gemmul8::hook::load_real<Fn>(GemmTraits<T>::gemm_v2_64_sym);
    return fn;
}

using GemmExFn = cublasStatus_t (*)(
    cublasHandle_t, cublasOperation_t, cublasOperation_t,
    int, int, int,
    const void *,
    const void *, cudaDataType, int,
    const void *, cudaDataType, int,
    const void *,
    void *, cudaDataType, int,
    cublasComputeType_t, cublasGemmAlgo_t);

using GemmExFn64 = cublasStatus_t (*)(
    cublasHandle_t, cublasOperation_t, cublasOperation_t,
    int64_t, int64_t, int64_t,
    const void *,
    const void *, cudaDataType, int64_t,
    const void *, cudaDataType, int64_t,
    const void *,
    void *, cudaDataType, int64_t,
    cublasComputeType_t, cublasGemmAlgo_t);

static inline GemmExFn real_gemm_ex() {
    static GemmExFn fn = gemmul8::hook::load_real<GemmExFn>(STR(cublasGemmEx));
    return fn;
}

static inline GemmExFn64 real_gemm_ex_64() {
    static GemmExFn64 fn = gemmul8::hook::load_real<GemmExFn64>(STR(cublasGemmEx_64));
    return fn;
}

// ---- GEMM call wrapper (INT8 or FP8) ----
template <typename T>
static inline cublasStatus_t call_gemmul8_gemm(
    gemmul8::Backend backend,
    cublasHandle_t handle,
    gemmul8::hook::HandleState &hst,
    cublasOperation_t transa, cublasOperation_t transb,
    int64_t m, int64_t n, int64_t k,
    const T *alpha, const T *A, int64_t lda, const T *B, int64_t ldb,
    const T *beta, T *C, int64_t ldc,
    int num_moduli, bool fastmode,
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
    cublasStatus_t st_lt = gemmul8::hook::ensure_lt_handle_locked(hst, &lt);
    if (st_lt != CUBLAS_STATUS_SUCCESS) return st_lt;

    (void)gemmul8::gemmLt<T, gemmul8::Backend::FP8>(
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
template <typename T>
static inline size_t call_gemmul8_workSize(
    gemmul8::Backend backend,
    int64_t m, int64_t n, int64_t k,
    int num_moduli,
    bool enable_skip_A, bool enable_skip_B,
    size_t *wA, size_t *wB //
) {
    constexpr bool COMPLEX       = GemmTraits<T>::isComplex;
    constexpr gemmul8::Func FUNC = gemmul8::Func::gemm;

    if (backend == gemmul8::Backend::INT8) {
        constexpr gemmul8::Backend BACKEND = gemmul8::Backend::INT8;
        return gemmul8::workSize<COMPLEX, BACKEND, FUNC>(
            static_cast<size_t>(m), static_cast<size_t>(n), static_cast<size_t>(k),
            num_moduli, enable_skip_A, enable_skip_B, wA, wB);
    } else {
        constexpr gemmul8::Backend BACKEND = gemmul8::Backend::FP8;
        return gemmul8::workSize<COMPLEX, BACKEND, FUNC>(
            static_cast<size_t>(m), static_cast<size_t>(n), static_cast<size_t>(k),
            num_moduli, enable_skip_A, enable_skip_B, wA, wB);
    }
}

template <typename T>
struct GemmArgs {
    cublasHandle_t handle    = nullptr;
    cublasOperation_t transa = CUBLAS_OP_N;
    cublasOperation_t transb = CUBLAS_OP_N;
    int64_t m                = 0;
    int64_t n                = 0;
    int64_t k                = 0;
    const T *alpha           = nullptr;
    const T *A               = nullptr;
    int64_t lda              = 0;
    const T *B               = nullptr;
    int64_t ldb              = 0;
    const T *beta            = nullptr;
    T *C                     = nullptr;
    int64_t ldc              = 0;
};

struct HookGemmEnv {
    int num_moduli           = 0;
    bool fastmode            = false;
    bool enable_skipA        = false;
    bool enable_skipB        = false;
    gemmul8::Backend backend = gemmul8::Backend::INT8;
};

template <typename T>
static inline cublasStatus_t run_gemmul8_gemm_emulation(const GemmArgs<T> &a, const HookGemmEnv &env) {
    auto sp = gemmul8::hook::get_state(a.handle);
    std::lock_guard<std::mutex> lk(sp->mtx);

    gemmul8::hook::init_max_workspace();

    cudaStream_t stream      = 0;
    cublasStatus_t st_stream = cublasGetStream(a.handle, &stream);
    if (st_stream != CUBLAS_STATUS_SUCCESS) return st_stream;

    cublasStatus_t st_ord = gemmul8::hook::ensure_stream_ordered_locked(*sp, stream);
    if (st_ord != CUBLAS_STATUS_SUCCESS) return st_ord;

    const auto ws_config = gemmul8::hook::workspace_config(
        a.handle, env.enable_skipA, env.enable_skipB);

    size_t wsizeA      = 0;
    size_t wsizeB      = 0;
    const size_t wsize = call_gemmul8_workSize<T>(
        env.backend,
        a.m, a.n, a.k,
        env.num_moduli,
        ws_config.enable_skipA,
        ws_config.enable_skipB,
        &wsizeA, &wsizeB);

    if (wsize < wsizeA + wsizeB) return CUBLAS_STATUS_INVALID_VALUE;

    const auto request = gemmul8::hook::workspace_request(ws_config, wsize, wsizeA, wsizeB);

    cublasStatus_t st = CUBLAS_STATUS_SUCCESS;

    if (ws_config.memory_saving) {
        st = gemmul8::hook::prepare_memory_saving_workspaces_locked(*sp, request, ws_config.limit, stream);
        if (st != CUBLAS_STATUS_SUCCESS) return st;
    }

    void *workA_raw = nullptr;
    void *workB_raw = nullptr;
    void *workC_raw = nullptr;

    st = gemmul8::hook::get_work_locked(*sp, sp->workA, sp->workA_size, request.workA, &workA_raw, "workA", stream);
    if (st != CUBLAS_STATUS_SUCCESS) return st;

    st = gemmul8::hook::get_work_locked(*sp, sp->workB, sp->workB_size, request.workB, &workB_raw, "workB", stream);
    if (st != CUBLAS_STATUS_SUCCESS) return st;

    st = gemmul8::hook::get_work_locked(*sp, sp->workC, sp->workC_size, request.workC, &workC_raw, "workC", stream);
    if (st != CUBLAS_STATUS_SUCCESS) return st;

    int8_t *workA = reinterpret_cast<int8_t *>(workA_raw);
    int8_t *workB = reinterpret_cast<int8_t *>(workB_raw);
    int8_t *workC = reinterpret_cast<int8_t *>(workC_raw);

    const auto keyA = make_gemm_A_key<T>(
        a.A, a.m, a.k, a.lda, a.transa,
        env.num_moduli, env.fastmode, env.backend);

    const auto keyB = make_gemm_B_key<T>(
        a.B, a.n, a.k, a.ldb, a.transb,
        env.num_moduli, env.fastmode, env.backend);

    const bool skipA = gemmul8::hook::can_skip_scaled_operand_locked(*sp, keyA, workA, ws_config.enable_skipA);
    const bool skipB = gemmul8::hook::can_skip_scaled_operand_locked(*sp, keyB, workB, ws_config.enable_skipB);

    st = call_gemmul8_gemm<T>(
        env.backend,
        a.handle,
        *sp,
        a.transa, a.transb,
        a.m, a.n, a.k,
        a.alpha, a.A, a.lda, a.B, a.ldb,
        a.beta, a.C, a.ldc,
        env.num_moduli, env.fastmode,
        reinterpret_cast<void *>(workC),
        reinterpret_cast<void *>(workA),
        reinterpret_cast<void *>(workB),
        ws_config.enable_skipA, ws_config.enable_skipB,
        skipA, skipB,
        stream);

    if (st != CUBLAS_STATUS_SUCCESS) return st;

    gemmul8::hook::update_scaled_operand_locked(*sp, keyA, workA);
    gemmul8::hook::update_scaled_operand_locked(*sp, keyB, workB);

    return CUBLAS_STATUS_SUCCESS;
}

template <typename T, typename NativeCall>
static inline cublasStatus_t gemm_common_impl(
    cublasHandle_t handle,
    cublasOperation_t transa,
    cublasOperation_t transb,
    int64_t m, int64_t n, int64_t k,
    const T *alpha,
    const T *A, int64_t lda,
    const T *B, int64_t ldb,
    const T *beta,
    T *C, int64_t ldc,
    NativeCall call_native //
) {
    constexpr gemmul8::hook::HookOp OP = gemmul8::hook::HookOp::GEMM;

    HookGemmEnv env{};
    GemmTraits<T>::get_env(OP, env.num_moduli, env.fastmode, env.enable_skipA, env.enable_skipB);
    env.backend = gemmul8::hook::requested_backend(OP);

    constexpr int num_moduli_min = 2;
    constexpr int num_moduli_max = gemmul8::hook::num_moduli_threshold<T>;
    if (env.num_moduli < num_moduli_min || num_moduli_max < env.num_moduli) {
        return call_native();
    }

    if (m < 0 || n < 0 || k < 0) return CUBLAS_STATUS_INVALID_VALUE;
    if (m == 0 || n == 0) return CUBLAS_STATUS_SUCCESS;
    if (k == 0) return call_native();

    const cublasStatus_t st_param = validate_gemm_params(
        handle, transa, transb, m, n, k,
        alpha, A, lda, B, ldb, beta, C, ldc);
    if (st_param != CUBLAS_STATUS_SUCCESS) return st_param;

    GemmArgs<T> args{};
    args.handle = handle;
    args.transa = transa;
    args.transb = transb;
    args.m      = m;
    args.n      = n;
    args.k      = k;
    args.alpha  = alpha;
    args.A      = A;
    args.lda    = lda;
    args.B      = B;
    args.ldb    = ldb;
    args.beta   = beta;
    args.C      = C;
    args.ldc    = ldc;

    return run_gemmul8_gemm_emulation<T>(args, env);
}

template <typename T>
static inline cublasStatus_t gemm_v2_impl(
    cublasHandle_t handle,
    cublasOperation_t transa,
    cublasOperation_t transb,
    int m, int n, int k,
    const T *alpha,
    const T *A, int lda,
    const T *B, int ldb,
    const T *beta,
    T *C, int ldc //
) {
    auto call_native = [&]() -> cublasStatus_t {
        auto fn = real_gemm_v2<T>();
        if (!fn) return CUBLAS_STATUS_NOT_INITIALIZED;
        return fn(handle, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
    };

    return gemm_common_impl<T>(
        handle, transa, transb, m, n, k,
        alpha, A, lda, B, ldb, beta, C, ldc,
        call_native);
}

template <typename T>
static inline cublasStatus_t gemm_v2_64_impl(
    cublasHandle_t handle,
    cublasOperation_t transa,
    cublasOperation_t transb,
    int64_t m, int64_t n, int64_t k,
    const T *alpha,
    const T *A, int64_t lda,
    const T *B, int64_t ldb,
    const T *beta,
    T *C, int64_t ldc //
) {
    auto call_native = [&]() -> cublasStatus_t {
        auto fn = real_gemm_v2_64<T>();
        if (!fn) return CUBLAS_STATUS_NOT_INITIALIZED;
        return fn(handle, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
    };

    return gemm_common_impl<T>(
        handle, transa, transb, m, n, k,
        alpha, A, lda, B, ldb, beta, C, ldc,
        call_native);
}

template <typename T>
static inline cublasStatus_t gemm_ex_impl(
    cublasHandle_t handle,
    cublasOperation_t transa,
    cublasOperation_t transb,
    int64_t m, int64_t n, int64_t k,
    const void *alpha,
    const void *A, int64_t lda,
    const void *B, int64_t ldb,
    const void *beta,
    void *C, int64_t ldc,
    int num_moduli, bool fastmode,
    gemmul8::Backend backend,
    bool enable_skipA, bool enable_skipB //
) {
    GemmArgs<T> args{};
    args.handle = handle;
    args.transa = transa;
    args.transb = transb;
    args.m      = m;
    args.n      = n;
    args.k      = k;
    args.alpha  = reinterpret_cast<const T *>(alpha);
    args.A      = reinterpret_cast<const T *>(A);
    args.lda    = lda;
    args.B      = reinterpret_cast<const T *>(B);
    args.ldb    = ldb;
    args.beta   = reinterpret_cast<const T *>(beta);
    args.C      = reinterpret_cast<T *>(C);
    args.ldc    = ldc;

    const HookGemmEnv env{num_moduli, fastmode, enable_skipA, enable_skipB, backend};
    return run_gemmul8_gemm_emulation<T>(args, env);
}

template <typename T, typename NativeCall>
static inline cublasStatus_t gemm_ex_typed_dispatch(
    cublasHandle_t handle,
    cublasOperation_t transa,
    cublasOperation_t transb,
    int64_t m, int64_t n, int64_t k,
    const void *alpha,
    const void *A, int64_t lda,
    const void *B, int64_t ldb,
    const void *beta,
    void *C, int64_t ldc,
    int num_moduli, bool fastmode,
    gemmul8::Backend backend,
    bool enable_skipA, bool enable_skipB,
    NativeCall call_native //
) {
    constexpr int num_moduli_min = 2;
    constexpr int num_moduli_max = gemmul8::hook::num_moduli_threshold<T>;
    if (num_moduli < num_moduli_min || num_moduli_max < num_moduli) {
        return call_native();
    }

    if (m < 0 || n < 0 || k < 0) return CUBLAS_STATUS_INVALID_VALUE;
    if (m == 0 || n == 0) return CUBLAS_STATUS_SUCCESS;
    if (k == 0) return call_native();

    const cublasStatus_t st_param = validate_gemm_params(
        handle, transa, transb, m, n, k,
        alpha, A, lda, B, ldb, beta, C, ldc);
    if (st_param != CUBLAS_STATUS_SUCCESS) return st_param;

    return gemm_ex_impl<T>(
        handle, transa, transb, m, n, k,
        alpha, A, lda, B, ldb, beta, C, ldc,
        num_moduli, fastmode, backend, enable_skipA, enable_skipB);
}

} // namespace

#define DEFINE_GEMM_V2_HOOK(NAME, TYPE)                                                                  \
    extern "C" cublasStatus_t NAME(                                                                      \
        cublasHandle_t handle,                                                                           \
        cublasOperation_t transa, cublasOperation_t transb,                                              \
        int m, int n, int k,                                                                             \
        const TYPE *alpha,                                                                               \
        const TYPE *A, int lda,                                                                          \
        const TYPE *B, int ldb,                                                                          \
        const TYPE *beta,                                                                                \
        TYPE *C, int ldc) {                                                                              \
        if (gemmul8::hook::inside_hook()) {                                                              \
            auto fn = real_gemm_v2<TYPE>();                                                              \
            if (!fn) return CUBLAS_STATUS_NOT_INITIALIZED;                                               \
            return fn(handle, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);             \
        }                                                                                                \
        gemmul8::hook::HookGuard guard;                                                                  \
        return gemm_v2_impl<TYPE>(handle, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc); \
    }

#define DEFINE_GEMM_V2_64_HOOK(NAME, TYPE)                                                                  \
    extern "C" cublasStatus_t NAME(                                                                         \
        cublasHandle_t handle,                                                                              \
        cublasOperation_t transa, cublasOperation_t transb,                                                 \
        int64_t m, int64_t n, int64_t k,                                                                    \
        const TYPE *alpha,                                                                                  \
        const TYPE *A, int64_t lda,                                                                         \
        const TYPE *B, int64_t ldb,                                                                         \
        const TYPE *beta,                                                                                   \
        TYPE *C, int64_t ldc) {                                                                             \
        if (gemmul8::hook::inside_hook()) {                                                                 \
            auto fn = real_gemm_v2_64<TYPE>();                                                              \
            if (!fn) return CUBLAS_STATUS_NOT_INITIALIZED;                                                  \
            return fn(handle, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);                \
        }                                                                                                   \
        gemmul8::hook::HookGuard guard;                                                                     \
        return gemm_v2_64_impl<TYPE>(handle, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc); \
    }

DEFINE_GEMM_V2_HOOK(cublasSgemm_v2, float)
DEFINE_GEMM_V2_HOOK(cublasDgemm_v2, double)
DEFINE_GEMM_V2_HOOK(cublasCgemm_v2, cuFloatComplex)
DEFINE_GEMM_V2_HOOK(cublasZgemm_v2, cuDoubleComplex)

DEFINE_GEMM_V2_64_HOOK(cublasSgemm_v2_64, float)
DEFINE_GEMM_V2_64_HOOK(cublasDgemm_v2_64, double)
DEFINE_GEMM_V2_64_HOOK(cublasCgemm_v2_64, cuFloatComplex)
DEFINE_GEMM_V2_64_HOOK(cublasZgemm_v2_64, cuDoubleComplex)

#undef DEFINE_GEMM_V2_HOOK
#undef DEFINE_GEMM_V2_64_HOOK

#if defined(__CUDACC__) && !defined(__HIPCC__)
namespace {

template <typename T> struct Gemm3mTraits;

template <> struct Gemm3mTraits<cuFloatComplex> {
    using Fn                            = typename GemmTraits<cuFloatComplex>::GemmV2Fn;
    using Fn64                          = typename GemmTraits<cuFloatComplex>::GemmV2Fn64;
    static constexpr const char *sym    = STR(cublasCgemm3m);
    static constexpr const char *sym_64 = STR(cublasCgemm3m_64);
};

template <> struct Gemm3mTraits<cuDoubleComplex> {
    using Fn                            = typename GemmTraits<cuDoubleComplex>::GemmV2Fn;
    using Fn64                          = typename GemmTraits<cuDoubleComplex>::GemmV2Fn64;
    static constexpr const char *sym    = STR(cublasZgemm3m);
    static constexpr const char *sym_64 = STR(cublasZgemm3m_64);
};

template <typename T>
static inline typename Gemm3mTraits<T>::Fn real_gemm3m() {
    using Fn     = typename Gemm3mTraits<T>::Fn;
    static Fn fn = gemmul8::hook::load_real<Fn>(Gemm3mTraits<T>::sym);
    return fn;
}

template <typename T>
static inline typename Gemm3mTraits<T>::Fn64 real_gemm3m_64() {
    using Fn     = typename Gemm3mTraits<T>::Fn64;
    static Fn fn = gemmul8::hook::load_real<Fn>(Gemm3mTraits<T>::sym_64);
    return fn;
}

template <typename T, typename NativeCall>
static inline cublasStatus_t gemm3m_common_impl(
    cublasHandle_t handle,
    cublasOperation_t transa,
    cublasOperation_t transb,
    int64_t m, int64_t n, int64_t k,
    const T *alpha,
    const T *A, int64_t lda,
    const T *B, int64_t ldb,
    const T *beta,
    T *C, int64_t ldc,
    NativeCall call_native //
) {
    return gemm_common_impl<T>(
        handle, transa, transb, m, n, k,
        alpha, A, lda, B, ldb, beta, C, ldc,
        call_native);
}

} // namespace

    #define DEFINE_GEMM3M_HOOK(NAME, TYPE)                                                                                      \
        extern "C" cublasStatus_t NAME(                                                                                         \
            cublasHandle_t handle,                                                                                              \
            cublasOperation_t transa, cublasOperation_t transb,                                                                 \
            int m, int n, int k,                                                                                                \
            const TYPE *alpha,                                                                                                  \
            const TYPE *A, int lda,                                                                                             \
            const TYPE *B, int ldb,                                                                                             \
            const TYPE *beta,                                                                                                   \
            TYPE *C, int ldc) {                                                                                                 \
            if (gemmul8::hook::inside_hook()) {                                                                                 \
                auto fn = real_gemm3m<TYPE>();                                                                                  \
                if (!fn) return CUBLAS_STATUS_NOT_INITIALIZED;                                                                  \
                return fn(handle, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);                                \
            }                                                                                                                   \
            gemmul8::hook::HookGuard guard;                                                                                     \
            auto call_native = [&]() -> cublasStatus_t {                                                                        \
                auto fn = real_gemm3m<TYPE>();                                                                                  \
                if (!fn) return CUBLAS_STATUS_NOT_INITIALIZED;                                                                  \
                return fn(handle, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);                                \
            };                                                                                                                  \
            return gemm3m_common_impl<TYPE>(handle, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc, call_native); \
        }

    #define DEFINE_GEMM3M_64_HOOK(NAME, TYPE)                                                                                   \
        extern "C" cublasStatus_t NAME(                                                                                         \
            cublasHandle_t handle,                                                                                              \
            cublasOperation_t transa, cublasOperation_t transb,                                                                 \
            int64_t m, int64_t n, int64_t k,                                                                                    \
            const TYPE *alpha,                                                                                                  \
            const TYPE *A, int64_t lda,                                                                                         \
            const TYPE *B, int64_t ldb,                                                                                         \
            const TYPE *beta,                                                                                                   \
            TYPE *C, int64_t ldc) {                                                                                             \
            if (gemmul8::hook::inside_hook()) {                                                                                 \
                auto fn = real_gemm3m_64<TYPE>();                                                                               \
                if (!fn) return CUBLAS_STATUS_NOT_INITIALIZED;                                                                  \
                return fn(handle, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);                                \
            }                                                                                                                   \
            gemmul8::hook::HookGuard guard;                                                                                     \
            auto call_native = [&]() -> cublasStatus_t {                                                                        \
                auto fn = real_gemm3m_64<TYPE>();                                                                               \
                if (!fn) return CUBLAS_STATUS_NOT_INITIALIZED;                                                                  \
                return fn(handle, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);                                \
            };                                                                                                                  \
            return gemm3m_common_impl<TYPE>(handle, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc, call_native); \
        }

DEFINE_GEMM3M_HOOK(cublasCgemm3m, cuFloatComplex)
DEFINE_GEMM3M_HOOK(cublasZgemm3m, cuDoubleComplex)
DEFINE_GEMM3M_64_HOOK(cublasCgemm3m_64, cuFloatComplex)
DEFINE_GEMM3M_64_HOOK(cublasZgemm3m_64, cuDoubleComplex)

    #undef DEFINE_GEMM3M_HOOK
    #undef DEFINE_GEMM3M_64_HOOK
#endif

namespace {

template <typename NativeCall>
static inline cublasStatus_t gemm_ex_common_impl(
    cublasHandle_t handle,
    cublasOperation_t transa,
    cublasOperation_t transb,
    int64_t m, int64_t n, int64_t k,
    const void *alpha,
    const void *A, cudaDataType Atype, int64_t lda,
    const void *B, cudaDataType Btype, int64_t ldb,
    const void *beta,
    void *C, cudaDataType Ctype, int64_t ldc,
    cublasComputeType_t computeType,
    NativeCall call_native //
) {
    if (gemmul8::hook::inside_hook()) return call_native();

    gemmul8::hook::HookGuard guard;

    constexpr gemmul8::hook::HookOp OP = gemmul8::hook::HookOp::GEMM;
    const gemmul8::Backend backend     = gemmul8::hook::requested_backend(OP);

    if (computeType == CUBLAS_COMPUTE_32F &&
        Atype == CUDA_R_32F && Btype == CUDA_R_32F && Ctype == CUDA_R_32F) {
        int num_moduli    = 0;
        bool fastmode     = false;
        bool enable_skipA = false;
        bool enable_skipB = false;
        gemmul8::hook::get_env_s(OP, num_moduli, fastmode, enable_skipA, enable_skipB);
        return gemm_ex_typed_dispatch<float>(
            handle, transa, transb, m, n, k,
            alpha, A, lda, B, ldb, beta, C, ldc,
            num_moduli, fastmode, backend, enable_skipA, enable_skipB,
            call_native);
    }

    if (computeType == CUBLAS_COMPUTE_64F &&
        Atype == CUDA_R_64F && Btype == CUDA_R_64F && Ctype == CUDA_R_64F) {
        int num_moduli    = 0;
        bool fastmode     = false;
        bool enable_skipA = false;
        bool enable_skipB = false;
        gemmul8::hook::get_env_d(OP, num_moduli, fastmode, enable_skipA, enable_skipB);
        return gemm_ex_typed_dispatch<double>(
            handle, transa, transb, m, n, k,
            alpha, A, lda, B, ldb, beta, C, ldc,
            num_moduli, fastmode, backend, enable_skipA, enable_skipB,
            call_native);
    }

    if (computeType == CUBLAS_COMPUTE_32F &&
        Atype == CUDA_C_32F && Btype == CUDA_C_32F && Ctype == CUDA_C_32F) {
        int num_moduli    = 0;
        bool fastmode     = false;
        bool enable_skipA = false;
        bool enable_skipB = false;
        gemmul8::hook::get_env_c(OP, num_moduli, fastmode, enable_skipA, enable_skipB);
        return gemm_ex_typed_dispatch<cuFloatComplex>(
            handle, transa, transb, m, n, k,
            alpha, A, lda, B, ldb, beta, C, ldc,
            num_moduli, fastmode, backend, enable_skipA, enable_skipB,
            call_native);
    }

    if (computeType == CUBLAS_COMPUTE_64F &&
        Atype == CUDA_C_64F && Btype == CUDA_C_64F && Ctype == CUDA_C_64F) {
        int num_moduli    = 0;
        bool fastmode     = false;
        bool enable_skipA = false;
        bool enable_skipB = false;
        gemmul8::hook::get_env_z(OP, num_moduli, fastmode, enable_skipA, enable_skipB);
        return gemm_ex_typed_dispatch<cuDoubleComplex>(
            handle, transa, transb, m, n, k,
            alpha, A, lda, B, ldb, beta, C, ldc,
            num_moduli, fastmode, backend, enable_skipA, enable_skipB,
            call_native);
    }

    return call_native();
}

} // namespace

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
    auto call_native = [&]() -> cublasStatus_t {
        auto fn = real_gemm_ex();
        if (!fn) return CUBLAS_STATUS_NOT_INITIALIZED;
        return fn(handle, transa, transb, m, n, k,
                  alpha, A, Atype, lda, B, Btype, ldb,
                  beta, C, Ctype, ldc,
                  computeType, algo);
    };

    return gemm_ex_common_impl(
        handle, transa, transb, m, n, k,
        alpha, A, Atype, lda, B, Btype, ldb,
        beta, C, Ctype, ldc,
        computeType,
        call_native);
#endif
}

extern "C" cublasStatus_t cublasGemmEx_64(
    cublasHandle_t handle,
    cublasOperation_t transa, cublasOperation_t transb,
    int64_t m, int64_t n, int64_t k,
    const void *alpha,
    const void *A, cudaDataType Atype, int64_t lda,
    const void *B, cudaDataType Btype, int64_t ldb,
    const void *beta,
    void *C, cudaDataType Ctype, int64_t ldc,
    cublasComputeType_t computeType, cublasGemmAlgo_t algo //
) {
#ifdef __CUDA_ARCH__
    return CUBLAS_STATUS_NOT_SUPPORTED;
#else
    auto call_native = [&]() -> cublasStatus_t {
        auto fn = real_gemm_ex_64();
        if (!fn) return CUBLAS_STATUS_NOT_INITIALIZED;
        return fn(handle, transa, transb, m, n, k,
                  alpha, A, Atype, lda, B, Btype, ldb,
                  beta, C, Ctype, ldc,
                  computeType, algo);
    };

    return gemm_ex_common_impl(
        handle, transa, transb, m, n, k,
        alpha, A, Atype, lda, B, Btype, ldb,
        beta, C, Ctype, ldc,
        computeType,
        call_native);
#endif
}

#if defined(__CUDACC__) && !defined(__HIPCC__)
namespace {

using SgemmExFn = cublasStatus_t (*)(
    cublasHandle_t,
    cublasOperation_t, cublasOperation_t,
    int, int, int,
    const float *,
    const void *, cudaDataType, int,
    const void *, cudaDataType, int,
    const float *,
    void *, cudaDataType, int);

using SgemmExFn64 = cublasStatus_t (*)(
    cublasHandle_t,
    cublasOperation_t, cublasOperation_t,
    int64_t, int64_t, int64_t,
    const float *,
    const void *, cudaDataType, int64_t,
    const void *, cudaDataType, int64_t,
    const float *,
    void *, cudaDataType, int64_t);

using CgemmExFn = cublasStatus_t (*)(
    cublasHandle_t,
    cublasOperation_t, cublasOperation_t,
    int, int, int,
    const cuFloatComplex *,
    const void *, cudaDataType, int,
    const void *, cudaDataType, int,
    const cuFloatComplex *,
    void *, cudaDataType, int);

using CgemmExFn64 = cublasStatus_t (*)(
    cublasHandle_t,
    cublasOperation_t, cublasOperation_t,
    int64_t, int64_t, int64_t,
    const cuFloatComplex *,
    const void *, cudaDataType, int64_t,
    const void *, cudaDataType, int64_t,
    const cuFloatComplex *,
    void *, cudaDataType, int64_t);

static inline SgemmExFn real_sgemm_ex() {
    static SgemmExFn fn = gemmul8::hook::load_real<SgemmExFn>(STR(cublasSgemmEx));
    return fn;
}

static inline SgemmExFn64 real_sgemm_ex_64() {
    static SgemmExFn64 fn = gemmul8::hook::load_real<SgemmExFn64>(STR(cublasSgemmEx_64));
    return fn;
}

static inline CgemmExFn real_cgemm_ex() {
    static CgemmExFn fn = gemmul8::hook::load_real<CgemmExFn>(STR(cublasCgemmEx));
    return fn;
}

static inline CgemmExFn64 real_cgemm_ex_64() {
    static CgemmExFn64 fn = gemmul8::hook::load_real<CgemmExFn64>(STR(cublasCgemmEx_64));
    return fn;
}

static inline CgemmExFn real_cgemm3m_ex() {
    static CgemmExFn fn = gemmul8::hook::load_real<CgemmExFn>(STR(cublasCgemm3mEx));
    return fn;
}

static inline CgemmExFn64 real_cgemm3m_ex_64() {
    static CgemmExFn64 fn = gemmul8::hook::load_real<CgemmExFn64>(STR(cublasCgemm3mEx_64));
    return fn;
}

} // namespace

extern "C" cublasStatus_t cublasSgemmEx(
    cublasHandle_t handle,
    cublasOperation_t transa, cublasOperation_t transb,
    int m, int n, int k,
    const float *alpha,
    const void *A, cudaDataType Atype, int lda,
    const void *B, cudaDataType Btype, int ldb,
    const float *beta,
    void *C, cudaDataType Ctype, int ldc //
) {
    #ifdef __CUDA_ARCH__
    return CUBLAS_STATUS_NOT_SUPPORTED;
    #else
    auto call_native = [&]() -> cublasStatus_t {
        auto fn = real_sgemm_ex();
        if (!fn) return CUBLAS_STATUS_NOT_INITIALIZED;

        return fn(handle, transa, transb, m, n, k,
                  alpha, A, Atype, lda, B, Btype, ldb,
                  beta, C, Ctype, ldc);
    };

    return gemm_ex_common_impl(
        handle, transa, transb,
        static_cast<int64_t>(m),
        static_cast<int64_t>(n),
        static_cast<int64_t>(k),
        alpha,
        A, Atype, static_cast<int64_t>(lda),
        B, Btype, static_cast<int64_t>(ldb),
        beta,
        C, Ctype, static_cast<int64_t>(ldc),
        CUBLAS_COMPUTE_32F,
        call_native);
    #endif
}

extern "C" cublasStatus_t cublasSgemmEx_64(
    cublasHandle_t handle,
    cublasOperation_t transa, cublasOperation_t transb,
    int64_t m, int64_t n, int64_t k,
    const float *alpha,
    const void *A, cudaDataType Atype, int64_t lda,
    const void *B, cudaDataType Btype, int64_t ldb,
    const float *beta,
    void *C, cudaDataType Ctype, int64_t ldc //
) {
    #ifdef __CUDA_ARCH__
    return CUBLAS_STATUS_NOT_SUPPORTED;
    #else
    auto call_native = [&]() -> cublasStatus_t {
        auto fn = real_sgemm_ex_64();
        if (!fn) return CUBLAS_STATUS_NOT_INITIALIZED;

        return fn(handle, transa, transb, m, n, k,
                  alpha, A, Atype, lda, B, Btype, ldb,
                  beta, C, Ctype, ldc);
    };

    return gemm_ex_common_impl(
        handle, transa, transb, m, n, k,
        alpha, A, Atype, lda, B, Btype, ldb,
        beta, C, Ctype, ldc, CUBLAS_COMPUTE_32F,
        call_native);
    #endif
}

extern "C" cublasStatus_t cublasCgemmEx(
    cublasHandle_t handle,
    cublasOperation_t transa, cublasOperation_t transb,
    int m, int n, int k,
    const cuFloatComplex *alpha,
    const void *A, cudaDataType Atype, int lda,
    const void *B, cudaDataType Btype, int ldb,
    const cuFloatComplex *beta,
    void *C, cudaDataType Ctype, int ldc //
) {
    #ifdef __CUDA_ARCH__
    return CUBLAS_STATUS_NOT_SUPPORTED;
    #else
    auto call_native = [&]() -> cublasStatus_t {
        auto fn = real_cgemm_ex();
        if (!fn) return CUBLAS_STATUS_NOT_INITIALIZED;

        return fn(handle, transa, transb, m, n, k,
                  alpha, A, Atype, lda, B, Btype, ldb,
                  beta, C, Ctype, ldc);
    };

    return gemm_ex_common_impl(
        handle, transa, transb,
        static_cast<int64_t>(m),
        static_cast<int64_t>(n),
        static_cast<int64_t>(k),
        alpha,
        A, Atype, static_cast<int64_t>(lda),
        B, Btype, static_cast<int64_t>(ldb),
        beta,
        C, Ctype, static_cast<int64_t>(ldc),
        CUBLAS_COMPUTE_32F,
        call_native);
    #endif
}

extern "C" cublasStatus_t cublasCgemmEx_64(
    cublasHandle_t handle,
    cublasOperation_t transa, cublasOperation_t transb,
    int64_t m, int64_t n, int64_t k,
    const cuFloatComplex *alpha,
    const void *A, cudaDataType Atype, int64_t lda,
    const void *B, cudaDataType Btype, int64_t ldb,
    const cuFloatComplex *beta,
    void *C, cudaDataType Ctype, int64_t ldc //
) {
    #ifdef __CUDA_ARCH__
    return CUBLAS_STATUS_NOT_SUPPORTED;
    #else
    auto call_native = [&]() -> cublasStatus_t {
        auto fn = real_cgemm_ex_64();
        if (!fn) return CUBLAS_STATUS_NOT_INITIALIZED;

        return fn(handle, transa, transb, m, n, k,
                  alpha, A, Atype, lda, B, Btype, ldb,
                  beta, C, Ctype, ldc);
    };

    return gemm_ex_common_impl(
        handle, transa, transb, m, n, k,
        alpha, A, Atype, lda, B, Btype, ldb,
        beta, C, Ctype, ldc,
        CUBLAS_COMPUTE_32F,
        call_native);
    #endif
}

extern "C" cublasStatus_t cublasCgemm3mEx(
    cublasHandle_t handle,
    cublasOperation_t transa, cublasOperation_t transb,
    int m, int n, int k,
    const cuFloatComplex *alpha,
    const void *A, cudaDataType Atype, int lda,
    const void *B, cudaDataType Btype, int ldb,
    const cuFloatComplex *beta,
    void *C, cudaDataType Ctype, int ldc //
) {
    #ifdef __CUDA_ARCH__
    return CUBLAS_STATUS_NOT_SUPPORTED;
    #else
    auto call_native = [&]() -> cublasStatus_t {
        auto fn = real_cgemm3m_ex();
        if (!fn) return CUBLAS_STATUS_NOT_INITIALIZED;

        return fn(handle, transa, transb, m, n, k,
                  alpha, A, Atype, lda, B, Btype, ldb,
                  beta, C, Ctype, ldc);
    };

    return gemm_ex_common_impl(
        handle, transa, transb,
        static_cast<int64_t>(m),
        static_cast<int64_t>(n),
        static_cast<int64_t>(k),
        alpha,
        A, Atype, static_cast<int64_t>(lda),
        B, Btype, static_cast<int64_t>(ldb),
        beta,
        C, Ctype, static_cast<int64_t>(ldc),
        CUBLAS_COMPUTE_32F,
        call_native);
    #endif
}

extern "C" cublasStatus_t cublasCgemm3mEx_64(
    cublasHandle_t handle,
    cublasOperation_t transa, cublasOperation_t transb,
    int64_t m, int64_t n, int64_t k,
    const cuFloatComplex *alpha,
    const void *A, cudaDataType Atype, int64_t lda,
    const void *B, cudaDataType Btype, int64_t ldb,
    const cuFloatComplex *beta,
    void *C, cudaDataType Ctype, int64_t ldc //
) {
    #ifdef __CUDA_ARCH__
    return CUBLAS_STATUS_NOT_SUPPORTED;
    #else
    auto call_native = [&]() -> cublasStatus_t {
        auto fn = real_cgemm3m_ex_64();
        if (!fn) return CUBLAS_STATUS_NOT_INITIALIZED;

        return fn(handle, transa, transb, m, n, k,
                  alpha, A, Atype, lda, B, Btype, ldb,
                  beta, C, Ctype, ldc);
    };

    return gemm_ex_common_impl(
        handle, transa, transb, m, n, k,
        alpha, A, Atype, lda, B, Btype, ldb,
        beta, C, Ctype, ldc, CUBLAS_COMPUTE_32F,
        call_native);
    #endif
}

#endif
