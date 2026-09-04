/**
 * HEMM hook
 * ---------
 * Hook targets:
 *   - cublas{C,Z}hemm_v2
 *   - cublas{C,Z}hemm_v2_64
 *
 * Notes:
 *   - HEMM is a complex Hermitian routine; there are no S/D HEMM hooks.
 *   - cublas{C,Z}hemm and cublas{C,Z}hemm_64 are not defined here
 *     because the v2 variants are used when v2 exists.
 */
#include "common.hpp"

namespace {

// Local scalar-kind helper for compile-time complex dispatch.
template <typename T> struct HookScalarIsComplex {
    static constexpr bool value = false;
};
template <> struct HookScalarIsComplex<cuFloatComplex> {
    static constexpr bool value = true;
};
template <> struct HookScalarIsComplex<cuDoubleComplex> {
    static constexpr bool value = true;
};

static inline bool is_valid_side(const cublasSideMode_t side) {
    return side == CUBLAS_SIDE_LEFT || side == CUBLAS_SIDE_RIGHT;
}

static inline bool is_valid_uplo(const cublasFillMode_t uplo) {
    return uplo == CUBLAS_FILL_MODE_UPPER || uplo == CUBLAS_FILL_MODE_LOWER;
}

static inline int64_t hemm_order_A(const cublasSideMode_t side, const int64_t m, const int64_t n) {
    return (side == CUBLAS_SIDE_LEFT) ? m : n;
}

static inline int64_t hemm_lda_min(const cublasSideMode_t side, const int64_t m, const int64_t n) {
    return std::max<int64_t>(1, hemm_order_A(side, m, n));
}

static inline int64_t hemm_ldb_min(const int64_t m) {
    return std::max<int64_t>(1, m);
}

static inline int64_t hemm_ldc_min(const int64_t m) {
    return std::max<int64_t>(1, m);
}

static inline cublasStatus_t validate_hemm_params(
    const cublasHandle_t handle,
    const cublasSideMode_t side,
    const cublasFillMode_t uplo,
    const int64_t m, const int64_t n,
    const void *alpha,
    const void *A, const int64_t lda,
    const void *B, const int64_t ldb,
    const void *beta,
    const void *C, const int64_t ldc //
) {
    if (!handle) return CUBLAS_STATUS_NOT_INITIALIZED;
    if (!is_valid_side(side) || !is_valid_uplo(uplo)) return CUBLAS_STATUS_INVALID_VALUE;
    if (!alpha || !beta || !A || !B || !C) return CUBLAS_STATUS_INVALID_VALUE;
    if (m < 0 || n < 0) return CUBLAS_STATUS_INVALID_VALUE;
    if (lda < hemm_lda_min(side, m, n)) return CUBLAS_STATUS_INVALID_VALUE;
    if (ldb < hemm_ldb_min(m)) return CUBLAS_STATUS_INVALID_VALUE;
    if (ldc < hemm_ldc_min(m)) return CUBLAS_STATUS_INVALID_VALUE;
    return CUBLAS_STATUS_SUCCESS;
}

template <typename T>
static inline gemmul8::hook::ScaledOperandKey make_hemm_A_key(
    const void *A,
    int64_t order,
    int64_t lda,
    cublasFillMode_t uplo,
    int num_moduli,
    bool fastmode,
    gemmul8::Backend backend //
) {
    return gemmul8::hook::ScaledOperandKey{
        .ptr         = A,
        .rows        = static_cast<size_t>(order),
        .cols        = static_cast<size_t>(order),
        .ld          = static_cast<size_t>(lda),
        .op          = CUBLAS_OP_N,
        .uplo        = uplo,
        .diag        = CUBLAS_DIAG_NON_UNIT,
        .matrix_kind = gemmul8::hook::MatrixKind::Hermitian,
        .scalar_kind = gemmul8::hook::scalar_kind_v<T>,
        .backend     = backend,
        .num_moduli  = num_moduli,
        .fastmode    = fastmode,
    };
}

template <typename T>
static inline gemmul8::hook::ScaledOperandKey make_hemm_B_key(
    const void *B,
    int64_t m,
    int64_t n,
    int64_t ldb,
    int num_moduli,
    bool fastmode,
    gemmul8::Backend backend //
) {
    return gemmul8::hook::ScaledOperandKey{
        .ptr         = B,
        .rows        = static_cast<size_t>(m),
        .cols        = static_cast<size_t>(n),
        .ld          = static_cast<size_t>(ldb),
        .op          = CUBLAS_OP_N,
        .uplo        = CUBLAS_FILL_MODE_FULL,
        .diag        = CUBLAS_DIAG_NON_UNIT,
        .matrix_kind = gemmul8::hook::MatrixKind::General,
        .scalar_kind = gemmul8::hook::scalar_kind_v<T>,
        .backend     = backend,
        .num_moduli  = num_moduli,
        .fastmode    = fastmode,
    };
}

template <typename T> struct HemmTraits;

template <> struct HemmTraits<cuFloatComplex> {
    static constexpr bool isComplex = true;
    static inline void get_env(gemmul8::hook::HookOp op, int &nm, bool &fm, bool &enA, bool &enB) {
        gemmul8::hook::get_env_c(op, nm, fm, enA, enB);
    }
    using HemmV2Fn = cublasStatus_t (*)(
        cublasHandle_t, cublasSideMode_t, cublasFillMode_t,
        int, int,
        const cuFloatComplex *, const cuFloatComplex *, int,
        const cuFloatComplex *, int,
        const cuFloatComplex *, cuFloatComplex *, int);
    using HemmV2Fn64 = cublasStatus_t (*)(
        cublasHandle_t, cublasSideMode_t, cublasFillMode_t,
        int64_t, int64_t,
        const cuFloatComplex *, const cuFloatComplex *, int64_t,
        const cuFloatComplex *, int64_t,
        const cuFloatComplex *, cuFloatComplex *, int64_t);
    static constexpr const char *hemm_v2_sym    = STR(cublasChemm_v2);
    static constexpr const char *hemm_v2_64_sym = STR(cublasChemm_v2_64);
};
template <> struct HemmTraits<cuDoubleComplex> {
    static constexpr bool isComplex = true;
    static inline void get_env(gemmul8::hook::HookOp op, int &nm, bool &fm, bool &enA, bool &enB) {
        gemmul8::hook::get_env_z(op, nm, fm, enA, enB);
    }
    using HemmV2Fn = cublasStatus_t (*)(
        cublasHandle_t, cublasSideMode_t, cublasFillMode_t,
        int, int,
        const cuDoubleComplex *, const cuDoubleComplex *, int,
        const cuDoubleComplex *, int,
        const cuDoubleComplex *, cuDoubleComplex *, int);
    using HemmV2Fn64 = cublasStatus_t (*)(
        cublasHandle_t, cublasSideMode_t, cublasFillMode_t,
        int64_t, int64_t,
        const cuDoubleComplex *, const cuDoubleComplex *, int64_t,
        const cuDoubleComplex *, int64_t,
        const cuDoubleComplex *, cuDoubleComplex *, int64_t);
    static constexpr const char *hemm_v2_sym    = STR(cublasZhemm_v2);
    static constexpr const char *hemm_v2_64_sym = STR(cublasZhemm_v2_64);
};

template <typename T>
static inline typename HemmTraits<T>::HemmV2Fn real_hemm_v2() {
    using Fn     = typename HemmTraits<T>::HemmV2Fn;
    static Fn fn = gemmul8::hook::load_real<Fn>(HemmTraits<T>::hemm_v2_sym);
    return fn;
}

template <typename T>
static inline typename HemmTraits<T>::HemmV2Fn64 real_hemm_v2_64() {
    using Fn     = typename HemmTraits<T>::HemmV2Fn64;
    static Fn fn = gemmul8::hook::load_real<Fn>(HemmTraits<T>::hemm_v2_64_sym);
    return fn;
}

// ---- HEMM call wrapper (INT8 or FP8) ----
template <typename T>
static inline cublasStatus_t call_gemmul8_hemm(
    gemmul8::Backend backend,
    cublasHandle_t handle,
    gemmul8::hook::HandleState &hst,
    cublasSideMode_t side, cublasFillMode_t uplo,
    int64_t m, int64_t n,
    const T *alpha, const T *A, int64_t lda,
    const T *B, int64_t ldb,
    const T *beta, T *C, int64_t ldc,
    int num_moduli, bool fastmode,
    void *workC, void *workA, void *workB,
    bool enable_skip_A, bool enable_skip_B,
    bool skip_A, bool skip_B,
    cudaStream_t stream //
) {
    if (backend == gemmul8::Backend::INT8) {

        (void)gemmul8::hemm<T, gemmul8::Backend::INT8>(
            handle,
            side, uplo,
            static_cast<size_t>(m), static_cast<size_t>(n),
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

    (void)gemmul8::hemmLt<T, gemmul8::Backend::FP8>(
        lt,
        side, uplo,
        static_cast<size_t>(m), static_cast<size_t>(n),
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
static inline size_t call_gemmul8_hemm_workSize(
    gemmul8::Backend backend,
    cublasSideMode_t side,
    int64_t m, int64_t n,
    int num_moduli,
    bool enable_skip_A, bool enable_skip_B,
    size_t *wA, size_t *wB //
) {
    constexpr bool COMPLEX       = HemmTraits<T>::isComplex;
    constexpr gemmul8::Func FUNC = gemmul8::Func::hemm;

    const size_t mm = static_cast<size_t>(m);
    const size_t nn = static_cast<size_t>(n);
    const size_t kk = static_cast<size_t>((side == CUBLAS_SIDE_LEFT) ? m : n);

    if (backend == gemmul8::Backend::INT8) {
        constexpr gemmul8::Backend BACKEND = gemmul8::Backend::INT8;
        return gemmul8::workSize<COMPLEX, BACKEND, FUNC>(
            mm, nn, kk,
            num_moduli, enable_skip_A, enable_skip_B, wA, wB);
    } else {
        constexpr gemmul8::Backend BACKEND = gemmul8::Backend::FP8;
        return gemmul8::workSize<COMPLEX, BACKEND, FUNC>(
            mm, nn, kk,
            num_moduli, enable_skip_A, enable_skip_B, wA, wB);
    }
}

template <typename T>
struct HemmArgs {
    cublasHandle_t handle = nullptr;
    cublasSideMode_t side = CUBLAS_SIDE_LEFT;
    cublasFillMode_t uplo = CUBLAS_FILL_MODE_UPPER;
    int64_t m             = 0;
    int64_t n             = 0;
    const T *alpha        = nullptr;
    const T *A            = nullptr;
    int64_t lda           = 0;
    const T *B            = nullptr;
    int64_t ldb           = 0;
    const T *beta         = nullptr;
    T *C                  = nullptr;
    int64_t ldc           = 0;
};

struct HookHemmEnv {
    int num_moduli           = 0;
    bool fastmode            = false;
    bool enable_skipA        = false;
    bool enable_skipB        = false;
    gemmul8::Backend backend = gemmul8::Backend::INT8;
};

template <typename T>
static inline cublasStatus_t run_gemmul8_hemm_emulation(const HemmArgs<T> &a, const HookHemmEnv &env) {
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
    const size_t wsize = call_gemmul8_hemm_workSize<T>(
        env.backend,
        a.side,
        a.m, a.n,
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

    const int64_t orderA = hemm_order_A(a.side, a.m, a.n);

    const auto keyA = make_hemm_A_key<T>(
        a.A, orderA, a.lda, a.uplo,
        env.num_moduli, env.fastmode, env.backend);

    const auto keyB = make_hemm_B_key<T>(
        a.B, a.m, a.n, a.ldb,
        env.num_moduli, env.fastmode, env.backend);

    const bool skipA = gemmul8::hook::can_skip_scaled_operand_locked(*sp, keyA, workA, ws_config.enable_skipA);
    const bool skipB = gemmul8::hook::can_skip_scaled_operand_locked(*sp, keyB, workB, ws_config.enable_skipB);

    st = call_gemmul8_hemm<T>(
        env.backend,
        a.handle,
        *sp,
        a.side, a.uplo,
        a.m, a.n,
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
static inline cublasStatus_t hemm_common_impl(
    cublasHandle_t handle,
    cublasSideMode_t side,
    cublasFillMode_t uplo,
    int64_t m, int64_t n,
    const T *alpha,
    const T *A, int64_t lda,
    const T *B, int64_t ldb,
    const T *beta,
    T *C, int64_t ldc,
    NativeCall call_native //
) {
    if (!is_valid_side(side)) return CUBLAS_STATUS_INVALID_VALUE;

    const gemmul8::hook::HookOp OP = gemmul8::hook::hook_op_from_hemm_side(side);

    HookHemmEnv env{};
    HemmTraits<T>::get_env(OP, env.num_moduli, env.fastmode, env.enable_skipA, env.enable_skipB);
    env.backend = gemmul8::hook::requested_backend(OP);

    constexpr int num_moduli_min = 2;
    constexpr int num_moduli_max = gemmul8::hook::num_moduli_threshold<T>;
    if (env.num_moduli < num_moduli_min || num_moduli_max < env.num_moduli) return call_native();

    if (m < 0 || n < 0) return CUBLAS_STATUS_INVALID_VALUE;
    if (m == 0 || n == 0) return CUBLAS_STATUS_SUCCESS;

    const cublasStatus_t st_param = validate_hemm_params(
        handle, side, uplo, m, n,
        alpha, A, lda, B, ldb, beta, C, ldc);
    if (st_param != CUBLAS_STATUS_SUCCESS) return st_param;

    HemmArgs<T> args{};
    args.handle = handle;
    args.side   = side;
    args.uplo   = uplo;
    args.m      = m;
    args.n      = n;
    args.alpha  = alpha;
    args.A      = A;
    args.lda    = lda;
    args.B      = B;
    args.ldb    = ldb;
    args.beta   = beta;
    args.C      = C;
    args.ldc    = ldc;

    return run_gemmul8_hemm_emulation<T>(args, env);
}

template <typename T>
static inline cublasStatus_t hemm_v2_impl(
    cublasHandle_t handle,
    cublasSideMode_t side,
    cublasFillMode_t uplo,
    int m, int n,
    const T *alpha,
    const T *A, int lda,
    const T *B, int ldb,
    const T *beta,
    T *C, int ldc //
) {
    auto call_native = [&]() -> cublasStatus_t {
        auto fn = real_hemm_v2<T>();
        if (!fn) return CUBLAS_STATUS_NOT_INITIALIZED;
        return fn(handle, side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc);
    };

    return hemm_common_impl<T>(
        handle, side, uplo,
        m, n,
        alpha, A, lda,
        B, ldb,
        beta, C, ldc,
        call_native);
}

template <typename T>
static inline cublasStatus_t hemm_v2_64_impl(
    cublasHandle_t handle,
    cublasSideMode_t side,
    cublasFillMode_t uplo,
    int64_t m, int64_t n,
    const T *alpha,
    const T *A, int64_t lda,
    const T *B, int64_t ldb,
    const T *beta,
    T *C, int64_t ldc //
) {
    auto call_native = [&]() -> cublasStatus_t {
        auto fn = real_hemm_v2_64<T>();
        if (!fn) return CUBLAS_STATUS_NOT_INITIALIZED;
        return fn(handle, side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc);
    };

    return hemm_common_impl<T>(
        handle, side, uplo,
        m, n,
        alpha, A, lda,
        B, ldb,
        beta, C, ldc,
        call_native);
}

} // namespace

#define DEFINE_HEMM_V2_HOOK(NAME, TYPE)                                                           \
    extern "C" cublasStatus_t NAME(                                                               \
        cublasHandle_t handle,                                                                    \
        cublasSideMode_t side,                                                                    \
        cublasFillMode_t uplo,                                                                    \
        int m, int n,                                                                             \
        const TYPE *alpha,                                                                        \
        const TYPE *A, int lda,                                                                   \
        const TYPE *B, int ldb,                                                                   \
        const TYPE *beta,                                                                         \
        TYPE *C, int ldc) {                                                                       \
        if (gemmul8::hook::inside_hook()) {                                                       \
            auto fn = real_hemm_v2<TYPE>();                                                       \
            if (!fn) return CUBLAS_STATUS_NOT_INITIALIZED;                                        \
            return fn(handle, side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc);             \
        }                                                                                         \
        gemmul8::hook::HookGuard guard;                                                           \
        return hemm_v2_impl<TYPE>(handle, side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc); \
    }

#define DEFINE_HEMM_V2_64_HOOK(NAME, TYPE)                                                           \
    extern "C" cublasStatus_t NAME(                                                                  \
        cublasHandle_t handle,                                                                       \
        cublasSideMode_t side,                                                                       \
        cublasFillMode_t uplo,                                                                       \
        int64_t m, int64_t n,                                                                        \
        const TYPE *alpha,                                                                           \
        const TYPE *A, int64_t lda,                                                                  \
        const TYPE *B, int64_t ldb,                                                                  \
        const TYPE *beta,                                                                            \
        TYPE *C, int64_t ldc) {                                                                      \
        if (gemmul8::hook::inside_hook()) {                                                          \
            auto fn = real_hemm_v2_64<TYPE>();                                                       \
            if (!fn) return CUBLAS_STATUS_NOT_INITIALIZED;                                           \
            return fn(handle, side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc);                \
        }                                                                                            \
        gemmul8::hook::HookGuard guard;                                                              \
        return hemm_v2_64_impl<TYPE>(handle, side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc); \
    }

DEFINE_HEMM_V2_HOOK(cublasChemm_v2, cuFloatComplex)
DEFINE_HEMM_V2_HOOK(cublasZhemm_v2, cuDoubleComplex)

DEFINE_HEMM_V2_64_HOOK(cublasChemm_v2_64, cuFloatComplex)
DEFINE_HEMM_V2_64_HOOK(cublasZhemm_v2_64, cuDoubleComplex)

#undef DEFINE_HEMM_V2_HOOK
#undef DEFINE_HEMM_V2_64_HOOK
