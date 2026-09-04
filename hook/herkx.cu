/**
 * HERKX hook
 * ----------
 * Hook targets:
 *   - cublas{C,Z}herkx
 *   - cublas{C,Z}herkx_64
 *
 * Notes:
 *   - HERKX is a complex Hermitian routine; there are no S/D hooks.
 *   - If a v2 routine exists in cuBLAS, only the v2/v2_64 symbols are defined.
 */
#include "common.hpp"

namespace {

template <typename T> struct HerkxTraits;
template <typename T> struct HerkxTraitsScalar;
template <> struct HerkxTraitsScalar<cuFloatComplex> {
    using type = float;
};
template <> struct HerkxTraitsScalar<cuDoubleComplex> {
    using type = double;
};
template <typename T> using HerkxTraitsScalarT = typename HerkxTraitsScalar<T>::type;

static inline bool is_valid_uplo(const cublasFillMode_t uplo) {
    return uplo == CUBLAS_FILL_MODE_UPPER || uplo == CUBLAS_FILL_MODE_LOWER;
}

static inline bool is_valid_herkx_op(const cublasOperation_t trans) {
    return trans == CUBLAS_OP_N || trans == CUBLAS_OP_C;
}

static inline int64_t herkx_lda_min(const cublasOperation_t trans, const int64_t n, const int64_t k) {
    return std::max<int64_t>(1, (trans == CUBLAS_OP_N) ? n : k);
}

static inline int64_t herkx_ldc_min(const int64_t n) {
    return std::max<int64_t>(1, n);
}

static inline cublasStatus_t validate_herkx_params(
    const cublasHandle_t handle,
    const cublasFillMode_t uplo,
    const cublasOperation_t trans,
    const int64_t n, const int64_t k,
    const void *alpha,
    const void *A, const int64_t lda,
    const void *B, const int64_t ldb,
    const void *beta,
    const void *C, const int64_t ldc //
) {
    if (!handle) return CUBLAS_STATUS_NOT_INITIALIZED;
    if (!is_valid_uplo(uplo) || !is_valid_herkx_op(trans)) return CUBLAS_STATUS_INVALID_VALUE;
    if (!alpha || !beta || !A || !B || !C) return CUBLAS_STATUS_INVALID_VALUE;
    if (n < 0 || k < 0) return CUBLAS_STATUS_INVALID_VALUE;
    if (lda < herkx_lda_min(trans, n, k)) return CUBLAS_STATUS_INVALID_VALUE;
    if (ldb < herkx_lda_min(trans, n, k)) return CUBLAS_STATUS_INVALID_VALUE;
    if (ldc < herkx_ldc_min(n)) return CUBLAS_STATUS_INVALID_VALUE;
    return CUBLAS_STATUS_SUCCESS;
}

template <typename T>
static inline gemmul8::hook::ScaledOperandKey make_herkx_operand_key(
    const void *ptr,
    int64_t n,
    int64_t k,
    int64_t ld,
    cublasOperation_t trans,
    int num_moduli,
    bool fastmode,
    gemmul8::Backend backend) {
    return gemmul8::hook::ScaledOperandKey{
        .ptr         = ptr,
        .rows        = static_cast<size_t>(n),
        .cols        = static_cast<size_t>(k),
        .ld          = static_cast<size_t>(ld),
        .op          = trans,
        .uplo        = CUBLAS_FILL_MODE_FULL,
        .diag        = CUBLAS_DIAG_NON_UNIT,
        .matrix_kind = gemmul8::hook::MatrixKind::General,
        .scalar_kind = gemmul8::hook::scalar_kind_v<T>,
        .backend     = backend,
        .num_moduli  = num_moduli,
        .fastmode    = fastmode,
    };
}

#define DEFINE_HERKX_TRAITS(TYPE, GETENV_FN, SYM, SYM64)                                                \
    template <> struct HerkxTraits<TYPE> {                                                              \
        static constexpr bool isComplex = true;                                                         \
        static inline void get_env(gemmul8::hook::HookOp op, int &nm, bool &fm, bool &enA, bool &enB) { \
            gemmul8::hook::GETENV_FN(op, nm, fm, enA, enB);                                             \
        }                                                                                               \
        using Scalar = HerkxTraitsScalarT<TYPE>;                                                        \
        using Fn     = cublasStatus_t (*)(                                                              \
            cublasHandle_t, cublasFillMode_t, cublasOperation_t,                                        \
            int, int,                                                                                   \
            const TYPE *, const TYPE *, int,                                                            \
            const TYPE *, int,                                                                          \
            const Scalar *, TYPE *, int);                                                               \
        using Fn64 = cublasStatus_t (*)(                                                                \
            cublasHandle_t, cublasFillMode_t, cublasOperation_t,                                        \
            int64_t, int64_t,                                                                           \
            const TYPE *, const TYPE *, int64_t,                                                        \
            const TYPE *, int64_t,                                                                      \
            const Scalar *, TYPE *, int64_t);                                                           \
        static constexpr const char *sym    = STR(SYM);                                                 \
        static constexpr const char *sym_64 = STR(SYM64);                                               \
    };

DEFINE_HERKX_TRAITS(cuFloatComplex, get_env_c, cublasCherkx, cublasCherkx_64)
DEFINE_HERKX_TRAITS(cuDoubleComplex, get_env_z, cublasZherkx, cublasZherkx_64)

#undef DEFINE_HERKX_TRAITS

template <typename T>
static inline typename HerkxTraits<T>::Fn real_herkx() {
    using Fn     = typename HerkxTraits<T>::Fn;
    static Fn fn = gemmul8::hook::load_real<Fn>(HerkxTraits<T>::sym);
    return fn;
}

template <typename T>
static inline typename HerkxTraits<T>::Fn64 real_herkx_64() {
    using Fn     = typename HerkxTraits<T>::Fn64;
    static Fn fn = gemmul8::hook::load_real<Fn>(HerkxTraits<T>::sym_64);
    return fn;
}

template <typename T>
static inline cublasStatus_t call_gemmul8_herkx(
    gemmul8::Backend backend,
    cublasHandle_t handle,
    gemmul8::hook::HandleState &hst,
    cublasFillMode_t uplo, cublasOperation_t trans,
    int64_t n, int64_t k,
    const T *alpha,
    const T *A, int64_t lda,
    const T *B, int64_t ldb,
    const HerkxTraitsScalarT<T> *beta,
    T *C, int64_t ldc,
    int num_moduli, bool fastmode,
    void *workC, void *workA, void *workB,
    bool enable_skip_A, bool enable_skip_B,
    bool skip_A, bool skip_B,
    cudaStream_t stream //
) {
    if (backend == gemmul8::Backend::INT8) {

        (void)gemmul8::herkx<T, gemmul8::Backend::INT8>(
            handle,
            uplo, trans,
            static_cast<size_t>(n), static_cast<size_t>(k),
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

    (void)gemmul8::herkxLt<T, gemmul8::Backend::FP8>(
        lt,
        uplo, trans,
        static_cast<size_t>(n), static_cast<size_t>(k),
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

template <typename T>
static inline size_t call_gemmul8_herkx_workSize(
    gemmul8::Backend backend,
    int64_t n, int64_t k,
    int num_moduli,
    bool enable_skip_A,
    bool enable_skip_B,
    size_t *wA,
    size_t *wB //
) {
    constexpr bool COMPLEX       = true;
    constexpr gemmul8::Func FUNC = gemmul8::Func::herkx;
    if (backend == gemmul8::Backend::INT8) {
        constexpr gemmul8::Backend BACKEND = gemmul8::Backend::INT8;
        return gemmul8::workSize<COMPLEX, BACKEND, FUNC>(
            static_cast<size_t>(n), static_cast<size_t>(n), static_cast<size_t>(k),
            num_moduli, enable_skip_A, enable_skip_B, wA, wB);
    } else {
        constexpr gemmul8::Backend BACKEND = gemmul8::Backend::FP8;
        return gemmul8::workSize<COMPLEX, BACKEND, FUNC>(
            static_cast<size_t>(n), static_cast<size_t>(n), static_cast<size_t>(k),
            num_moduli, enable_skip_A, enable_skip_B, wA, wB);
    }
}

template <typename T>
struct HerkxArgs {
    cublasHandle_t handle             = nullptr;
    cublasFillMode_t uplo             = CUBLAS_FILL_MODE_UPPER;
    cublasOperation_t trans           = CUBLAS_OP_N;
    int64_t n                         = 0;
    int64_t k                         = 0;
    const T *alpha                    = nullptr;
    const T *A                        = nullptr;
    int64_t lda                       = 0;
    const T *B                        = nullptr;
    int64_t ldb                       = 0;
    const HerkxTraitsScalarT<T> *beta = nullptr;
    T *C                              = nullptr;
    int64_t ldc                       = 0;
};

struct HookHerkxEnv {
    int num_moduli           = 0;
    bool fastmode            = false;
    bool enable_skipA        = false;
    bool enable_skipB        = false;
    gemmul8::Backend backend = gemmul8::Backend::INT8;
};

template <typename T>
static inline cublasStatus_t run_gemmul8_herkx_emulation(const HerkxArgs<T> &a, const HookHerkxEnv &env) {
    auto sp = gemmul8::hook::get_state(a.handle);
    std::lock_guard<std::mutex> lk(sp->mtx);

    gemmul8::hook::init_max_workspace();

    cudaStream_t stream      = 0;
    cublasStatus_t st_stream = cublasGetStream(a.handle, &stream);
    if (st_stream != CUBLAS_STATUS_SUCCESS) return st_stream;

    cublasStatus_t st_ord = gemmul8::hook::ensure_stream_ordered_locked(*sp, stream);
    if (st_ord != CUBLAS_STATUS_SUCCESS) return st_ord;

    const auto ws_config = gemmul8::hook::workspace_config(a.handle, env.enable_skipA, env.enable_skipB);

    size_t wsizeA      = 0;
    size_t wsizeB      = 0;
    const size_t wsize = call_gemmul8_herkx_workSize<T>(
        env.backend,
        a.n, a.k,
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

    const auto keyA = make_herkx_operand_key<T>(
        a.A, a.n, a.k, a.lda, a.trans,
        env.num_moduli, env.fastmode, env.backend);

    const auto keyB = make_herkx_operand_key<T>(
        a.B, a.n, a.k, a.ldb, a.trans,
        env.num_moduli, env.fastmode, env.backend);

    const bool skipA = gemmul8::hook::can_skip_scaled_operand_locked(*sp, keyA, workA, ws_config.enable_skipA);
    const bool skipB = gemmul8::hook::can_skip_scaled_operand_locked(*sp, keyB, workB, ws_config.enable_skipB);

    st = call_gemmul8_herkx<T>(
        env.backend,
        a.handle,
        *sp,
        a.uplo, a.trans,
        a.n, a.k,
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
static inline cublasStatus_t herkx_common_impl(
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    int64_t n, int64_t k,
    const T *alpha,
    const T *A, int64_t lda,
    const T *B, int64_t ldb,
    const HerkxTraitsScalarT<T> *beta,
    T *C, int64_t ldc,
    NativeCall call_native //
) {
    constexpr gemmul8::hook::HookOp OP = gemmul8::hook::HookOp::HERKX;

    HookHerkxEnv env{};
    HerkxTraits<T>::get_env(OP, env.num_moduli, env.fastmode, env.enable_skipA, env.enable_skipB);
    env.backend = gemmul8::hook::requested_backend(OP);

    constexpr int num_moduli_min = 2;
    constexpr int num_moduli_max = gemmul8::hook::num_moduli_threshold<T>;
    if (env.num_moduli < num_moduli_min || num_moduli_max < env.num_moduli) return call_native();

    if (n < 0 || k < 0) return CUBLAS_STATUS_INVALID_VALUE;
    if (n == 0) return CUBLAS_STATUS_SUCCESS;
    if (k == 0) return call_native();

    const cublasStatus_t st_param = validate_herkx_params(
        handle, uplo, trans, n, k,
        alpha, A, lda, B, ldb, beta, C, ldc);
    if (st_param != CUBLAS_STATUS_SUCCESS) return st_param;

    HerkxArgs<T> args{};
    args.handle = handle;
    args.uplo   = uplo;
    args.trans  = trans;
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

    return run_gemmul8_herkx_emulation<T>(args, env);
}

template <typename T>
static inline cublasStatus_t herkx_impl(
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    int n, int k,
    const T *alpha,
    const T *A, int lda,
    const T *B, int ldb,
    const HerkxTraitsScalarT<T> *beta,
    T *C, int ldc //
) {
    auto call_native = [&]() -> cublasStatus_t {
        auto fn = real_herkx<T>();
        if (!fn) return CUBLAS_STATUS_NOT_INITIALIZED;
        return fn(handle, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
    };

    return herkx_common_impl<T>(
        handle, uplo, trans,
        n, k,
        alpha, A, lda,
        B, ldb,
        beta, C, ldc,
        call_native);
}

template <typename T>
static inline cublasStatus_t herkx_64_impl(
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    int64_t n, int64_t k,
    const T *alpha,
    const T *A, int64_t lda,
    const T *B, int64_t ldb,
    const HerkxTraitsScalarT<T> *beta,
    T *C, int64_t ldc //
) {
    auto call_native = [&]() -> cublasStatus_t {
        auto fn = real_herkx_64<T>();
        if (!fn) return CUBLAS_STATUS_NOT_INITIALIZED;
        return fn(handle, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
    };

    return herkx_common_impl<T>(
        handle, uplo, trans,
        n, k,
        alpha, A, lda,
        B, ldb,
        beta, C, ldc,
        call_native);
}

} // namespace

#define DEFINE_HERKX_HOOK(NAME, TYPE)                                                            \
    extern "C" cublasStatus_t NAME(                                                              \
        cublasHandle_t handle,                                                                   \
        cublasFillMode_t uplo,                                                                   \
        cublasOperation_t trans,                                                                 \
        int n, int k,                                                                            \
        const TYPE *alpha,                                                                       \
        const TYPE *A, int lda,                                                                  \
        const TYPE *B, int ldb,                                                                  \
        const HerkxTraitsScalarT<TYPE> *beta,                                                    \
        TYPE *C, int ldc) {                                                                      \
        if (gemmul8::hook::inside_hook()) {                                                      \
            auto fn = real_herkx<TYPE>();                                                        \
            if (!fn) return CUBLAS_STATUS_NOT_INITIALIZED;                                       \
            return fn(handle, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc);           \
        }                                                                                        \
        gemmul8::hook::HookGuard guard;                                                          \
        return herkx_impl<TYPE>(handle, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc); \
    }

#define DEFINE_HERKX_64_HOOK(NAME, TYPE)                                                            \
    extern "C" cublasStatus_t NAME(                                                                 \
        cublasHandle_t handle,                                                                      \
        cublasFillMode_t uplo,                                                                      \
        cublasOperation_t trans,                                                                    \
        int64_t n, int64_t k,                                                                       \
        const TYPE *alpha,                                                                          \
        const TYPE *A, int64_t lda,                                                                 \
        const TYPE *B, int64_t ldb,                                                                 \
        const HerkxTraitsScalarT<TYPE> *beta,                                                       \
        TYPE *C, int64_t ldc) {                                                                     \
        if (gemmul8::hook::inside_hook()) {                                                         \
            auto fn = real_herkx_64<TYPE>();                                                        \
            if (!fn) return CUBLAS_STATUS_NOT_INITIALIZED;                                          \
            return fn(handle, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc);              \
        }                                                                                           \
        gemmul8::hook::HookGuard guard;                                                             \
        return herkx_64_impl<TYPE>(handle, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc); \
    }

DEFINE_HERKX_HOOK(cublasCherkx, cuFloatComplex)
DEFINE_HERKX_HOOK(cublasZherkx, cuDoubleComplex)

DEFINE_HERKX_64_HOOK(cublasCherkx_64, cuFloatComplex)
DEFINE_HERKX_64_HOOK(cublasZherkx_64, cuDoubleComplex)

#undef DEFINE_HERKX_HOOK
#undef DEFINE_HERKX_64_HOOK
