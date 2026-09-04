/**
 * TRSM hook
 * ---------
 * Hook targets:
 *   - cublas{S,D,C,Z}trsm_v2
 *   - cublas{S,D,C,Z}trsm_v2_64
 *
 * Notes:
 *   - cublas{S,D,C,Z}trsm and cublas{S,D,C,Z}trsm_64 are not defined here
 *     because the v2 variants are used when v2 exists.
 *   - TRSM overwrites B in place, following the BLAS convention.
 *   - The GEMMul8 TRSM direct interface uses a single workspace and does not
 *     expose workA/workB or skip-scaling arguments.
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

static inline bool is_valid_side(const cublasSideMode_t side) {
    return side == CUBLAS_SIDE_LEFT || side == CUBLAS_SIDE_RIGHT;
}

static inline bool is_valid_uplo(const cublasFillMode_t uplo) {
    return uplo == CUBLAS_FILL_MODE_UPPER || uplo == CUBLAS_FILL_MODE_LOWER;
}

static inline bool is_valid_trsm_op(const cublasOperation_t trans) {
    return trans == CUBLAS_OP_N || trans == CUBLAS_OP_T || trans == CUBLAS_OP_C;
}

static inline bool is_valid_diag(const cublasDiagType_t diag) {
    return diag == CUBLAS_DIAG_NON_UNIT || diag == CUBLAS_DIAG_UNIT;
}

static inline int64_t trsm_order_A(const cublasSideMode_t side, const int64_t m, const int64_t n) {
    return (side == CUBLAS_SIDE_LEFT) ? m : n;
}

static inline int64_t trsm_lda_min(const cublasSideMode_t side, const int64_t m, const int64_t n) {
    return std::max<int64_t>(1, trsm_order_A(side, m, n));
}

static inline int64_t trsm_ldb_min(const int64_t m) {
    return std::max<int64_t>(1, m);
}

static inline cublasStatus_t validate_trsm_params(
    const cublasHandle_t handle,
    const cublasSideMode_t side,
    const cublasFillMode_t uplo,
    const cublasOperation_t trans,
    const cublasDiagType_t diag,
    const int64_t m, const int64_t n,
    const void *alpha,
    const void *A, const int64_t lda,
    const void *B, const int64_t ldb //
) {
    if (!handle) return CUBLAS_STATUS_NOT_INITIALIZED;
    if (!is_valid_side(side) || !is_valid_uplo(uplo) || !is_valid_trsm_op(trans) || !is_valid_diag(diag)) {
        return CUBLAS_STATUS_INVALID_VALUE;
    }
    if (!alpha || !A || !B) return CUBLAS_STATUS_INVALID_VALUE;
    if (m < 0 || n < 0) return CUBLAS_STATUS_INVALID_VALUE;
    if (lda < trsm_lda_min(side, m, n)) return CUBLAS_STATUS_INVALID_VALUE;
    if (ldb < trsm_ldb_min(m)) return CUBLAS_STATUS_INVALID_VALUE;
    return CUBLAS_STATUS_SUCCESS;
}

template <typename T> struct TrsmTraits;

#define DEFINE_TRSM_TRAITS(TYPE, GETENV_FN, SYM, SYM64)                                                 \
    template <> struct TrsmTraits<TYPE> {                                                               \
        static constexpr bool isComplex = HookScalarIsComplex<TYPE>::value;                             \
        static inline void get_env(gemmul8::hook::HookOp op, int &nm, bool &fm, bool &enA, bool &enB) { \
            gemmul8::hook::GETENV_FN(op, nm, fm, enA, enB);                                             \
        }                                                                                               \
        using Fn = cublasStatus_t (*)(                                                                  \
            cublasHandle_t, cublasSideMode_t, cublasFillMode_t, cublasOperation_t, cublasDiagType_t,    \
            int, int,                                                                                   \
            const TYPE *, const TYPE *, int,                                                            \
            TYPE *, int);                                                                               \
        using Fn64 = cublasStatus_t (*)(                                                                \
            cublasHandle_t, cublasSideMode_t, cublasFillMode_t, cublasOperation_t, cublasDiagType_t,    \
            int64_t, int64_t,                                                                           \
            const TYPE *, const TYPE *, int64_t,                                                        \
            TYPE *, int64_t);                                                                           \
        static constexpr const char *sym    = STR(SYM);                                                 \
        static constexpr const char *sym_64 = STR(SYM64);                                               \
    };

DEFINE_TRSM_TRAITS(float, get_env_s, cublasStrsm_v2, cublasStrsm_v2_64)
DEFINE_TRSM_TRAITS(double, get_env_d, cublasDtrsm_v2, cublasDtrsm_v2_64)
DEFINE_TRSM_TRAITS(cuFloatComplex, get_env_c, cublasCtrsm_v2, cublasCtrsm_v2_64)
DEFINE_TRSM_TRAITS(cuDoubleComplex, get_env_z, cublasZtrsm_v2, cublasZtrsm_v2_64)

#undef DEFINE_TRSM_TRAITS

template <typename T>
static inline typename TrsmTraits<T>::Fn real_trsm_v2() {
    using Fn     = typename TrsmTraits<T>::Fn;
    static Fn fn = gemmul8::hook::load_real<Fn>(TrsmTraits<T>::sym);
    return fn;
}

template <typename T>
static inline typename TrsmTraits<T>::Fn64 real_trsm_v2_64() {
    using Fn     = typename TrsmTraits<T>::Fn64;
    static Fn fn = gemmul8::hook::load_real<Fn>(TrsmTraits<T>::sym_64);
    return fn;
}

template <typename T>
static inline cublasStatus_t call_gemmul8_trsm(
    gemmul8::Backend backend,
    cublasHandle_t handle,
    gemmul8::hook::HandleState &hst,
    cublasSideMode_t side, cublasFillMode_t uplo,
    cublasOperation_t trans, cublasDiagType_t diag,
    int64_t m, int64_t n,
    const T *alpha,
    const T *A, int64_t lda,
    T *B, int64_t ldb,
    int num_moduli, bool fastmode,
    void *work,
    cudaStream_t stream //
) {
    if (backend == gemmul8::Backend::INT8) {

        (void)gemmul8::trsm<T, gemmul8::Backend::INT8>(
            handle,
            side, uplo, trans, diag,
            static_cast<size_t>(m), static_cast<size_t>(n),
            alpha,
            A, static_cast<size_t>(lda),
            B, static_cast<size_t>(ldb),
            num_moduli, fastmode,
            work);

        return CUBLAS_STATUS_SUCCESS;
    }

    cublasLtHandle_t lt  = nullptr;
    cublasStatus_t st_lt = gemmul8::hook::ensure_lt_handle_locked(hst, &lt);
    if (st_lt != CUBLAS_STATUS_SUCCESS) return st_lt;

    (void)gemmul8::trsmLt<T, gemmul8::Backend::FP8>(
        lt,
        side, uplo, trans, diag,
        static_cast<size_t>(m), static_cast<size_t>(n),
        alpha,
        A, static_cast<size_t>(lda),
        B, static_cast<size_t>(ldb),
        num_moduli, fastmode,
        work,
        stream);

    return CUBLAS_STATUS_SUCCESS;
}

template <typename T>
static inline size_t call_gemmul8_trsm_workSize(
    gemmul8::Backend backend,
    cublasHandle_t handle,
    cublasLtHandle_t lt,
    cublasSideMode_t side,
    int64_t m, int64_t n,
    int num_moduli //
) {
    const size_t mm = static_cast<size_t>(m);
    const size_t nn = static_cast<size_t>(n);

    if (backend == gemmul8::Backend::INT8) {
        constexpr gemmul8::Backend BACKEND = gemmul8::Backend::INT8;
        return gemmul8::workSizeTrsm<T, BACKEND>(handle, side, mm, nn, num_moduli);
    } else {
        constexpr gemmul8::Backend BACKEND = gemmul8::Backend::FP8;
        return gemmul8::workSizeTrsmLt<T, BACKEND>(lt, side, mm, nn, num_moduli);
    }
}

template <typename T>
struct TrsmArgs {
    cublasHandle_t handle   = nullptr;
    cublasSideMode_t side   = CUBLAS_SIDE_LEFT;
    cublasFillMode_t uplo   = CUBLAS_FILL_MODE_UPPER;
    cublasOperation_t trans = CUBLAS_OP_N;
    cublasDiagType_t diag   = CUBLAS_DIAG_NON_UNIT;
    int64_t m               = 0;
    int64_t n               = 0;
    const T *alpha          = nullptr;
    const T *A              = nullptr;
    int64_t lda             = 0;
    T *B                    = nullptr;
    int64_t ldb             = 0;
};

struct HookTrsmEnv {
    int num_moduli           = 0;
    bool fastmode            = false;
    gemmul8::Backend backend = gemmul8::Backend::INT8;
};

template <typename T>
static inline cublasStatus_t run_gemmul8_trsm_emulation(const TrsmArgs<T> &a, const HookTrsmEnv &env) {
    auto sp = gemmul8::hook::get_state(a.handle);
    std::lock_guard<std::mutex> lk(sp->mtx);

    gemmul8::hook::init_max_workspace();

    cudaStream_t stream      = 0;
    cublasStatus_t st_stream = cublasGetStream(a.handle, &stream);
    if (st_stream != CUBLAS_STATUS_SUCCESS) return st_stream;

    cublasStatus_t st_ord = gemmul8::hook::ensure_stream_ordered_locked(*sp, stream);
    if (st_ord != CUBLAS_STATUS_SUCCESS) return st_ord;

    cublasLtHandle_t lt = nullptr;
    if (env.backend == gemmul8::Backend::FP8) {
        const cublasStatus_t st_lt = gemmul8::hook::ensure_lt_handle_locked(*sp, &lt);
        if (st_lt != CUBLAS_STATUS_SUCCESS) return st_lt;
    }

    const size_t need = call_gemmul8_trsm_workSize<T>(
        env.backend,
        a.handle, lt,
        a.side,
        a.m, a.n,
        env.num_moduli);

    const auto ws_config = gemmul8::hook::workspace_config(a.handle, false, false);
    auto request         = gemmul8::hook::workspace_request(ws_config, need, 0, 0);

    if (!ws_config.memory_saving_mode) {
        request.workC = std::max(request.workC, gemmul8::hook::max_workSizeC);
    }

    cublasStatus_t st = CUBLAS_STATUS_SUCCESS;
    if (ws_config.memory_saving) {
        st = gemmul8::hook::prepare_memory_saving_workspaces_locked(*sp, request, ws_config.limit, stream);
        if (st != CUBLAS_STATUS_SUCCESS) return st;
    }

    void *work_raw = nullptr;

    st = gemmul8::hook::get_work_locked(
        *sp,
        sp->workC,
        sp->workC_size,
        request.workC,
        &work_raw,
        "workC",
        stream);
    if (st != CUBLAS_STATUS_SUCCESS) return st;

    st = call_gemmul8_trsm<T>(
        env.backend,
        a.handle,
        *sp,
        a.side, a.uplo, a.trans, a.diag,
        a.m, a.n,
        a.alpha, a.A, a.lda,
        a.B, a.ldb,
        env.num_moduli, env.fastmode,
        work_raw,
        stream);

    return st;
}

template <typename T, typename NativeCall>
static inline cublasStatus_t trsm_common_impl(
    cublasHandle_t handle,
    cublasSideMode_t side,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    cublasDiagType_t diag,
    int64_t m, int64_t n,
    const T *alpha,
    const T *A, int64_t lda,
    T *B, int64_t ldb,
    NativeCall call_native //
) {
    if (!is_valid_side(side)) return CUBLAS_STATUS_INVALID_VALUE;

    const gemmul8::hook::HookOp OP = gemmul8::hook::hook_op_from_trsm_side(side);

    int dummy_num_moduli     = 0;
    bool dummy_fastmode      = false;
    bool dummy_enable_skip_A = false;
    bool dummy_enable_skip_B = false;

    HookTrsmEnv env{};
    TrsmTraits<T>::get_env(OP, dummy_num_moduli, dummy_fastmode, dummy_enable_skip_A, dummy_enable_skip_B);
    env.num_moduli = dummy_num_moduli;
    env.fastmode   = dummy_fastmode;
    env.backend    = gemmul8::hook::requested_backend(OP);

    constexpr int num_moduli_min = 2;
    constexpr int num_moduli_max = gemmul8::hook::num_moduli_threshold<T>;
    if (env.num_moduli < num_moduli_min || num_moduli_max < env.num_moduli) return call_native();

    if (m < 0 || n < 0) return CUBLAS_STATUS_INVALID_VALUE;
    if (m == 0 || n == 0) return CUBLAS_STATUS_SUCCESS;

    const cublasStatus_t st_param = validate_trsm_params(
        handle, side, uplo, trans, diag, m, n,
        alpha, A, lda, B, ldb);
    if (st_param != CUBLAS_STATUS_SUCCESS) return st_param;

    gemmul8::hook::apply_block_size_trsm_environment(handle);

    TrsmArgs<T> args{};
    args.handle = handle;
    args.side   = side;
    args.uplo   = uplo;
    args.trans  = trans;
    args.diag   = diag;
    args.m      = m;
    args.n      = n;
    args.alpha  = alpha;
    args.A      = A;
    args.lda    = lda;
    args.B      = B;
    args.ldb    = ldb;

    return run_gemmul8_trsm_emulation<T>(args, env);
}

template <typename T>
static inline cublasStatus_t trsm_v2_impl(
    cublasHandle_t handle,
    cublasSideMode_t side,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    cublasDiagType_t diag,
    int m, int n,
    const T *alpha,
    const T *A, int lda,
    T *B, int ldb //
) {
    auto call_native = [&]() -> cublasStatus_t {
        auto fn = real_trsm_v2<T>();
        if (!fn) return CUBLAS_STATUS_NOT_INITIALIZED;
        return fn(handle, side, uplo, trans, diag, m, n, alpha, A, lda, B, ldb);
    };

    return trsm_common_impl<T>(
        handle, side, uplo, trans, diag,
        m, n,
        alpha, A, lda,
        B, ldb,
        call_native);
}

template <typename T>
static inline cublasStatus_t trsm_v2_64_impl(
    cublasHandle_t handle,
    cublasSideMode_t side,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    cublasDiagType_t diag,
    int64_t m, int64_t n,
    const T *alpha,
    const T *A, int64_t lda,
    T *B, int64_t ldb //
) {
    auto call_native = [&]() -> cublasStatus_t {
        auto fn = real_trsm_v2_64<T>();
        if (!fn) return CUBLAS_STATUS_NOT_INITIALIZED;
        return fn(handle, side, uplo, trans, diag, m, n, alpha, A, lda, B, ldb);
    };

    return trsm_common_impl<T>(
        handle, side, uplo, trans, diag,
        m, n,
        alpha, A, lda,
        B, ldb,
        call_native);
}

} // namespace

#define DEFINE_TRSM_V2_HOOK(NAME, TYPE)                                                               \
    extern "C" cublasStatus_t NAME(                                                                   \
        cublasHandle_t handle,                                                                        \
        cublasSideMode_t side, cublasFillMode_t uplo, cublasOperation_t trans, cublasDiagType_t diag, \
        int m, int n,                                                                                 \
        const TYPE *alpha,                                                                            \
        const TYPE *A, int lda,                                                                       \
        TYPE *B, int ldb) {                                                                           \
        if (gemmul8::hook::inside_hook()) {                                                           \
            auto fn = real_trsm_v2<TYPE>();                                                           \
            if (!fn) return CUBLAS_STATUS_NOT_INITIALIZED;                                            \
            return fn(handle, side, uplo, trans, diag, m, n, alpha, A, lda, B, ldb);                  \
        }                                                                                             \
        gemmul8::hook::HookGuard guard;                                                               \
        return trsm_v2_impl<TYPE>(handle, side, uplo, trans, diag, m, n, alpha, A, lda, B, ldb);      \
    }

#define DEFINE_TRSM_V2_64_HOOK(NAME, TYPE)                                                            \
    extern "C" cublasStatus_t NAME(                                                                   \
        cublasHandle_t handle,                                                                        \
        cublasSideMode_t side, cublasFillMode_t uplo, cublasOperation_t trans, cublasDiagType_t diag, \
        int64_t m, int64_t n,                                                                         \
        const TYPE *alpha,                                                                            \
        const TYPE *A, int64_t lda,                                                                   \
        TYPE *B, int64_t ldb) {                                                                       \
        if (gemmul8::hook::inside_hook()) {                                                           \
            auto fn = real_trsm_v2_64<TYPE>();                                                        \
            if (!fn) return CUBLAS_STATUS_NOT_INITIALIZED;                                            \
            return fn(handle, side, uplo, trans, diag, m, n, alpha, A, lda, B, ldb);                  \
        }                                                                                             \
        gemmul8::hook::HookGuard guard;                                                               \
        return trsm_v2_64_impl<TYPE>(handle, side, uplo, trans, diag, m, n, alpha, A, lda, B, ldb);   \
    }

DEFINE_TRSM_V2_HOOK(cublasStrsm_v2, float)
DEFINE_TRSM_V2_HOOK(cublasDtrsm_v2, double)
DEFINE_TRSM_V2_HOOK(cublasCtrsm_v2, cuFloatComplex)
DEFINE_TRSM_V2_HOOK(cublasZtrsm_v2, cuDoubleComplex)

DEFINE_TRSM_V2_64_HOOK(cublasStrsm_v2_64, float)
DEFINE_TRSM_V2_64_HOOK(cublasDtrsm_v2_64, double)
DEFINE_TRSM_V2_64_HOOK(cublasCtrsm_v2_64, cuFloatComplex)
DEFINE_TRSM_V2_64_HOOK(cublasZtrsm_v2_64, cuDoubleComplex)

#undef DEFINE_TRSM_V2_HOOK
#undef DEFINE_TRSM_V2_64_HOOK
