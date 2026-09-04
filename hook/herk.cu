
/**
 * HERK hook
 * ---------
 * Hook targets:
 *   - cublas{C,Z}herk_v2
 *   - cublas{C,Z}herk_v2_64
 *
 * CUDA-only hook targets:
 *   - cublasCherkEx
 *   - cublasCherkEx_64
 *   - cublasCherk3mEx
 *   - cublasCherk3mEx_64
 *
 * Notes:
 *   - HERK is a complex Hermitian rank-k update; there are no S/D HERK hooks.
 *   - cublas{C,Z}herk and cublas{C,Z}herk_64 are not defined here because
 *     the v2 variants are used when v2 exists.
 */
#include "common.hpp"

namespace {

template <typename T> struct HerkScalar;
template <> struct HerkScalar<cuFloatComplex> {
    using type = float;
};
template <> struct HerkScalar<cuDoubleComplex> {
    using type = double;
};
template <typename T> using HerkScalarT = typename HerkScalar<T>::type;

static inline bool is_valid_uplo(const cublasFillMode_t uplo) {
    return uplo == CUBLAS_FILL_MODE_UPPER || uplo == CUBLAS_FILL_MODE_LOWER;
}

static inline bool is_valid_herk_op(const cublasOperation_t trans) {
    return trans == CUBLAS_OP_N || trans == CUBLAS_OP_C;
}

static inline int64_t herk_lda_min(const cublasOperation_t trans, const int64_t n, const int64_t k) {
    return std::max<int64_t>(1, (trans == CUBLAS_OP_N) ? n : k);
}

static inline int64_t herk_ldc_min(const int64_t n) {
    return std::max<int64_t>(1, n);
}

static inline cublasStatus_t validate_herk_params(
    const cublasHandle_t handle,
    const cublasFillMode_t uplo,
    const cublasOperation_t trans,
    const int64_t n, const int64_t k,
    const void *alpha,
    const void *A, const int64_t lda,
    const void *beta,
    const void *C, const int64_t ldc //
) {
    if (!handle) return CUBLAS_STATUS_NOT_INITIALIZED;
    if (!is_valid_uplo(uplo) || !is_valid_herk_op(trans)) return CUBLAS_STATUS_INVALID_VALUE;
    if (!alpha || !beta || !A || !C) return CUBLAS_STATUS_INVALID_VALUE;
    if (n < 0 || k < 0) return CUBLAS_STATUS_INVALID_VALUE;
    if (lda < herk_lda_min(trans, n, k)) return CUBLAS_STATUS_INVALID_VALUE;
    if (ldc < herk_ldc_min(n)) return CUBLAS_STATUS_INVALID_VALUE;
    return CUBLAS_STATUS_SUCCESS;
}

template <typename T>
static inline gemmul8::hook::ScaledOperandKey make_herk_A_key(
    const void *A,
    int64_t n,
    int64_t k,
    int64_t lda,
    cublasOperation_t trans,
    int num_moduli,
    bool fastmode,
    gemmul8::Backend backend) {
    return gemmul8::hook::ScaledOperandKey{
        .ptr         = A,
        .rows        = static_cast<size_t>(n),
        .cols        = static_cast<size_t>(k),
        .ld          = static_cast<size_t>(lda),
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

template <typename T> struct HerkTraits;

template <> struct HerkTraits<cuFloatComplex> {
    using Scalar                    = float;
    static constexpr bool isComplex = true;
    static inline void get_env(gemmul8::hook::HookOp op, int &nm, bool &fm, bool &enA, bool &enB) {
        gemmul8::hook::get_env_c(op, nm, fm, enA, enB);
    }
    using HerkV2Fn = cublasStatus_t (*)(
        cublasHandle_t, cublasFillMode_t, cublasOperation_t,
        int, int,
        const float *, const cuFloatComplex *, int,
        const float *, cuFloatComplex *, int);
    using HerkV2Fn64 = cublasStatus_t (*)(
        cublasHandle_t, cublasFillMode_t, cublasOperation_t,
        int64_t, int64_t,
        const float *, const cuFloatComplex *, int64_t,
        const float *, cuFloatComplex *, int64_t);
    static constexpr const char *herk_v2_sym    = STR(cublasCherk_v2);
    static constexpr const char *herk_v2_64_sym = STR(cublasCherk_v2_64);
};

template <> struct HerkTraits<cuDoubleComplex> {
    using Scalar                    = double;
    static constexpr bool isComplex = true;
    static inline void get_env(gemmul8::hook::HookOp op, int &nm, bool &fm, bool &enA, bool &enB) {
        gemmul8::hook::get_env_z(op, nm, fm, enA, enB);
    }
    using HerkV2Fn = cublasStatus_t (*)(
        cublasHandle_t, cublasFillMode_t, cublasOperation_t,
        int, int,
        const double *, const cuDoubleComplex *, int,
        const double *, cuDoubleComplex *, int);
    using HerkV2Fn64 = cublasStatus_t (*)(
        cublasHandle_t, cublasFillMode_t, cublasOperation_t,
        int64_t, int64_t,
        const double *, const cuDoubleComplex *, int64_t,
        const double *, cuDoubleComplex *, int64_t);
    static constexpr const char *herk_v2_sym    = STR(cublasZherk_v2);
    static constexpr const char *herk_v2_64_sym = STR(cublasZherk_v2_64);
};

template <typename T>
static inline typename HerkTraits<T>::HerkV2Fn real_herk_v2() {
    using Fn     = typename HerkTraits<T>::HerkV2Fn;
    static Fn fn = gemmul8::hook::load_real<Fn>(HerkTraits<T>::herk_v2_sym);
    return fn;
}

template <typename T>
static inline typename HerkTraits<T>::HerkV2Fn64 real_herk_v2_64() {
    using Fn     = typename HerkTraits<T>::HerkV2Fn64;
    static Fn fn = gemmul8::hook::load_real<Fn>(HerkTraits<T>::herk_v2_64_sym);
    return fn;
}

template <typename T>
static inline cublasStatus_t call_gemmul8_herk(
    gemmul8::Backend backend,
    cublasHandle_t handle,
    gemmul8::hook::HandleState &hst,
    cublasFillMode_t uplo, cublasOperation_t trans,
    int64_t n, int64_t k,
    const HerkScalarT<T> *alpha,
    const T *A, int64_t lda,
    const HerkScalarT<T> *beta,
    T *C, int64_t ldc,
    int num_moduli, bool fastmode,
    void *workC, void *workA,
    bool enable_skip_A,
    bool skip_A,
    cudaStream_t stream //
) {
    if (backend == gemmul8::Backend::INT8) {

        (void)gemmul8::herk<T, gemmul8::Backend::INT8, T>(
            handle,
            uplo, trans,
            static_cast<size_t>(n), static_cast<size_t>(k),
            alpha,
            A, static_cast<size_t>(lda),
            beta,
            C, static_cast<size_t>(ldc),
            num_moduli, fastmode,
            workC, workA,
            enable_skip_A, skip_A);

        return CUBLAS_STATUS_SUCCESS;
    }

    cublasLtHandle_t lt  = nullptr;
    cublasStatus_t st_lt = gemmul8::hook::ensure_lt_handle_locked(hst, &lt);
    if (st_lt != CUBLAS_STATUS_SUCCESS) return st_lt;

    (void)gemmul8::herkLt<T, gemmul8::Backend::FP8, T>(
        lt,
        uplo, trans,
        static_cast<size_t>(n), static_cast<size_t>(k),
        alpha,
        A, static_cast<size_t>(lda),
        beta,
        C, static_cast<size_t>(ldc),
        num_moduli, fastmode,
        workC, workA,
        enable_skip_A, skip_A,
        stream);

    return CUBLAS_STATUS_SUCCESS;
}

template <typename T>
static inline size_t call_gemmul8_herk_workSize(
    gemmul8::Backend backend,
    int64_t n, int64_t k,
    int num_moduli,
    bool enable_skip_A,
    size_t *wA //
) {
    constexpr bool COMPLEX       = true;
    constexpr gemmul8::Func FUNC = gemmul8::Func::herk;

    size_t wB = 0;
    if (backend == gemmul8::Backend::INT8) {
        constexpr gemmul8::Backend BACKEND = gemmul8::Backend::INT8;
        return gemmul8::workSize<COMPLEX, BACKEND, FUNC>(
            static_cast<size_t>(n), static_cast<size_t>(n), static_cast<size_t>(k),
            num_moduli, enable_skip_A, false, wA, &wB);
    } else {
        constexpr gemmul8::Backend BACKEND = gemmul8::Backend::FP8;
        return gemmul8::workSize<COMPLEX, BACKEND, FUNC>(
            static_cast<size_t>(n), static_cast<size_t>(n), static_cast<size_t>(k),
            num_moduli, enable_skip_A, false, wA, &wB);
    }
}

template <typename T>
struct HerkArgs {
    cublasHandle_t handle       = nullptr;
    cublasFillMode_t uplo       = CUBLAS_FILL_MODE_UPPER;
    cublasOperation_t trans     = CUBLAS_OP_N;
    int64_t n                   = 0;
    int64_t k                   = 0;
    const HerkScalarT<T> *alpha = nullptr;
    const T *A                  = nullptr;
    int64_t lda                 = 0;
    const HerkScalarT<T> *beta  = nullptr;
    T *C                        = nullptr;
    int64_t ldc                 = 0;
};

struct HookHerkEnv {
    int num_moduli           = 0;
    bool fastmode            = false;
    bool enable_skipA        = false;
    gemmul8::Backend backend = gemmul8::Backend::INT8;
};

template <typename T>
static inline cublasStatus_t run_gemmul8_herk_emulation(const HerkArgs<T> &a, const HookHerkEnv &env) {
    auto sp = gemmul8::hook::get_state(a.handle);
    std::lock_guard<std::mutex> lk(sp->mtx);

    gemmul8::hook::init_max_workspace();

    cudaStream_t stream      = 0;
    cublasStatus_t st_stream = cublasGetStream(a.handle, &stream);
    if (st_stream != CUBLAS_STATUS_SUCCESS) return st_stream;

    cublasStatus_t st_ord = gemmul8::hook::ensure_stream_ordered_locked(*sp, stream);
    if (st_ord != CUBLAS_STATUS_SUCCESS) return st_ord;

    const auto ws_config = gemmul8::hook::workspace_config(a.handle, env.enable_skipA, false);

    size_t wsizeA      = 0;
    const size_t wsize = call_gemmul8_herk_workSize<T>(
        env.backend,
        a.n, a.k,
        env.num_moduli,
        ws_config.enable_skipA,
        &wsizeA);

    if (wsize < wsizeA) return CUBLAS_STATUS_INVALID_VALUE;

    const auto request = gemmul8::hook::workspace_request(ws_config, wsize, wsizeA, 0);

    cublasStatus_t st = CUBLAS_STATUS_SUCCESS;

    if (ws_config.memory_saving) {
        st = gemmul8::hook::prepare_memory_saving_workspaces_locked(*sp, request, ws_config.limit, stream);
        if (st != CUBLAS_STATUS_SUCCESS) return st;
    }

    void *workA_raw = nullptr;
    void *workC_raw = nullptr;

    st = gemmul8::hook::get_work_locked(*sp, sp->workA, sp->workA_size, request.workA, &workA_raw, "workA", stream);
    if (st != CUBLAS_STATUS_SUCCESS) return st;

    st = gemmul8::hook::get_work_locked(*sp, sp->workC, sp->workC_size, request.workC, &workC_raw, "workC", stream);
    if (st != CUBLAS_STATUS_SUCCESS) return st;

    int8_t *workA = reinterpret_cast<int8_t *>(workA_raw);
    int8_t *workC = reinterpret_cast<int8_t *>(workC_raw);

    const auto keyA = make_herk_A_key<T>(
        a.A, a.n, a.k, a.lda, a.trans,
        env.num_moduli, env.fastmode, env.backend);

    const bool skipA = gemmul8::hook::can_skip_scaled_operand_locked(*sp, keyA, workA, ws_config.enable_skipA);

    st = call_gemmul8_herk<T>(
        env.backend,
        a.handle,
        *sp,
        a.uplo, a.trans,
        a.n, a.k,
        a.alpha, a.A, a.lda,
        a.beta, a.C, a.ldc,
        env.num_moduli, env.fastmode,
        reinterpret_cast<void *>(workC),
        reinterpret_cast<void *>(workA),
        ws_config.enable_skipA,
        skipA,
        stream);

    if (st != CUBLAS_STATUS_SUCCESS) return st;

    gemmul8::hook::update_scaled_operand_locked(*sp, keyA, workA);

    return CUBLAS_STATUS_SUCCESS;
}

template <typename T, typename NativeCall>
static inline cublasStatus_t herk_common_impl(
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    int64_t n, int64_t k,
    const HerkScalarT<T> *alpha,
    const T *A, int64_t lda,
    const HerkScalarT<T> *beta,
    T *C, int64_t ldc,
    NativeCall call_native //
) {
    constexpr gemmul8::hook::HookOp OP = gemmul8::hook::HookOp::HERK;

    HookHerkEnv env{};
    bool dummy_skipB = false;
    HerkTraits<T>::get_env(OP, env.num_moduli, env.fastmode, env.enable_skipA, dummy_skipB);
    env.backend = gemmul8::hook::requested_backend(OP);

    constexpr int num_moduli_min = 2;
    constexpr int num_moduli_max = gemmul8::hook::num_moduli_threshold<T>;
    if (env.num_moduli < num_moduli_min || num_moduli_max < env.num_moduli) return call_native();

    if (n < 0 || k < 0) return CUBLAS_STATUS_INVALID_VALUE;
    if (n == 0) return CUBLAS_STATUS_SUCCESS;
    if (k == 0) return call_native();

    const cublasStatus_t st_param = validate_herk_params(
        handle, uplo, trans, n, k,
        alpha, A, lda, beta, C, ldc);
    if (st_param != CUBLAS_STATUS_SUCCESS) return st_param;

    HerkArgs<T> args{};
    args.handle = handle;
    args.uplo   = uplo;
    args.trans  = trans;
    args.n      = n;
    args.k      = k;
    args.alpha  = alpha;
    args.A      = A;
    args.lda    = lda;
    args.beta   = beta;
    args.C      = C;
    args.ldc    = ldc;

    return run_gemmul8_herk_emulation<T>(args, env);
}

template <typename T>
static inline cublasStatus_t herk_v2_impl(
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    int n, int k,
    const HerkScalarT<T> *alpha,
    const T *A, int lda,
    const HerkScalarT<T> *beta,
    T *C, int ldc //
) {
    auto call_native = [&]() -> cublasStatus_t {
        auto fn = real_herk_v2<T>();
        if (!fn) return CUBLAS_STATUS_NOT_INITIALIZED;
        return fn(handle, uplo, trans, n, k, alpha, A, lda, beta, C, ldc);
    };

    return herk_common_impl<T>(
        handle, uplo, trans,
        n, k,
        alpha, A, lda,
        beta, C, ldc,
        call_native);
}

template <typename T>
static inline cublasStatus_t herk_v2_64_impl(
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    int64_t n, int64_t k,
    const HerkScalarT<T> *alpha,
    const T *A, int64_t lda,
    const HerkScalarT<T> *beta,
    T *C, int64_t ldc //
) {
    auto call_native = [&]() -> cublasStatus_t {
        auto fn = real_herk_v2_64<T>();
        if (!fn) return CUBLAS_STATUS_NOT_INITIALIZED;
        return fn(handle, uplo, trans, n, k, alpha, A, lda, beta, C, ldc);
    };

    return herk_common_impl<T>(
        handle, uplo, trans,
        n, k,
        alpha, A, lda,
        beta, C, ldc,
        call_native);
}

template <typename NativeCall>
static inline cublasStatus_t herk_ex_common_impl(
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    int64_t n, int64_t k,
    const float *alpha,
    const void *A, cudaDataType Atype, int64_t lda,
    const float *beta,
    void *C, cudaDataType Ctype, int64_t ldc,
    NativeCall call_native //
) {
    if (gemmul8::hook::inside_hook()) return call_native();
    gemmul8::hook::HookGuard guard;

    if (Atype != CUDA_C_32F || Ctype != CUDA_C_32F) return call_native();

    return herk_common_impl<cuFloatComplex>(
        handle, uplo, trans,
        n, k,
        alpha, reinterpret_cast<const cuFloatComplex *>(A), lda,
        beta, reinterpret_cast<cuFloatComplex *>(C), ldc,
        call_native);
}

} // namespace

#define DEFINE_HERK_V2_HOOK(NAME, TYPE)                                                    \
    extern "C" cublasStatus_t NAME(                                                        \
        cublasHandle_t handle,                                                             \
        cublasFillMode_t uplo,                                                             \
        cublasOperation_t trans,                                                           \
        int n, int k,                                                                      \
        const HerkScalarT<TYPE> *alpha,                                                    \
        const TYPE *A, int lda,                                                            \
        const HerkScalarT<TYPE> *beta,                                                     \
        TYPE *C, int ldc) {                                                                \
        if (gemmul8::hook::inside_hook()) {                                                \
            auto fn = real_herk_v2<TYPE>();                                                \
            if (!fn) return CUBLAS_STATUS_NOT_INITIALIZED;                                 \
            return fn(handle, uplo, trans, n, k, alpha, A, lda, beta, C, ldc);             \
        }                                                                                  \
        gemmul8::hook::HookGuard guard;                                                    \
        return herk_v2_impl<TYPE>(handle, uplo, trans, n, k, alpha, A, lda, beta, C, ldc); \
    }

#define DEFINE_HERK_V2_64_HOOK(NAME, TYPE)                                                    \
    extern "C" cublasStatus_t NAME(                                                           \
        cublasHandle_t handle,                                                                \
        cublasFillMode_t uplo,                                                                \
        cublasOperation_t trans,                                                              \
        int64_t n, int64_t k,                                                                 \
        const HerkScalarT<TYPE> *alpha,                                                       \
        const TYPE *A, int64_t lda,                                                           \
        const HerkScalarT<TYPE> *beta,                                                        \
        TYPE *C, int64_t ldc) {                                                               \
        if (gemmul8::hook::inside_hook()) {                                                   \
            auto fn = real_herk_v2_64<TYPE>();                                                \
            if (!fn) return CUBLAS_STATUS_NOT_INITIALIZED;                                    \
            return fn(handle, uplo, trans, n, k, alpha, A, lda, beta, C, ldc);                \
        }                                                                                     \
        gemmul8::hook::HookGuard guard;                                                       \
        return herk_v2_64_impl<TYPE>(handle, uplo, trans, n, k, alpha, A, lda, beta, C, ldc); \
    }

DEFINE_HERK_V2_HOOK(cublasCherk_v2, cuFloatComplex)
DEFINE_HERK_V2_HOOK(cublasZherk_v2, cuDoubleComplex)

DEFINE_HERK_V2_64_HOOK(cublasCherk_v2_64, cuFloatComplex)
DEFINE_HERK_V2_64_HOOK(cublasZherk_v2_64, cuDoubleComplex)

#undef DEFINE_HERK_V2_HOOK
#undef DEFINE_HERK_V2_64_HOOK

#if defined(__CUDACC__) && !defined(__HIPCC__)
namespace {

using CherkExFn = cublasStatus_t (*)(
    cublasHandle_t,
    cublasFillMode_t, cublasOperation_t,
    int, int,
    const float *,
    const void *, cudaDataType, int,
    const float *,
    void *, cudaDataType, int);

static inline CherkExFn real_cherk_ex() {
    static CherkExFn fn = gemmul8::hook::load_real<CherkExFn>(STR(cublasCherkEx));
    return fn;
}

} // namespace

extern "C" cublasStatus_t cublasCherkEx(
    cublasHandle_t handle,
    cublasFillMode_t uplo, cublasOperation_t trans,
    int n, int k,
    const float *alpha,
    const void *A, cudaDataType Atype, int lda,
    const float *beta,
    void *C, cudaDataType Ctype, int ldc //
) {
    auto call_native = [&]() -> cublasStatus_t {
        auto fn = real_cherk_ex();
        if (!fn) return CUBLAS_STATUS_NOT_INITIALIZED;
        return fn(handle, uplo, trans, n, k, alpha, A, Atype, lda, beta, C, Ctype, ldc);
    };

    return herk_ex_common_impl(
        handle, uplo, trans,
        static_cast<int64_t>(n), static_cast<int64_t>(k),
        alpha, A, Atype, static_cast<int64_t>(lda),
        beta, C, Ctype, static_cast<int64_t>(ldc),
        call_native);
}
#endif

#if defined(__CUDACC__) && !defined(__HIPCC__)
namespace {

using CherkExFn64 = cublasStatus_t (*)(
    cublasHandle_t,
    cublasFillMode_t, cublasOperation_t,
    int64_t, int64_t,
    const float *,
    const void *, cudaDataType, int64_t,
    const float *,
    void *, cudaDataType, int64_t);

static inline CherkExFn64 real_cherk_ex_64() {
    static CherkExFn64 fn = gemmul8::hook::load_real<CherkExFn64>(STR(cublasCherkEx_64));
    return fn;
}

static inline CherkExFn real_cherk3m_ex() {
    static CherkExFn fn = gemmul8::hook::load_real<CherkExFn>(STR(cublasCherk3mEx));
    return fn;
}

static inline CherkExFn64 real_cherk3m_ex_64() {
    static CherkExFn64 fn = gemmul8::hook::load_real<CherkExFn64>(STR(cublasCherk3mEx_64));
    return fn;
}

} // namespace

extern "C" cublasStatus_t cublasCherkEx_64(
    cublasHandle_t handle,
    cublasFillMode_t uplo, cublasOperation_t trans,
    int64_t n, int64_t k,
    const float *alpha,
    const void *A, cudaDataType Atype, int64_t lda,
    const float *beta,
    void *C, cudaDataType Ctype, int64_t ldc //
) {
    auto call_native = [&]() -> cublasStatus_t {
        auto fn = real_cherk_ex_64();
        if (!fn) return CUBLAS_STATUS_NOT_INITIALIZED;
        return fn(handle, uplo, trans, n, k, alpha, A, Atype, lda, beta, C, Ctype, ldc);
    };

    return herk_ex_common_impl(
        handle, uplo, trans,
        n, k,
        alpha, A, Atype, lda,
        beta, C, Ctype, ldc,
        call_native);
}

extern "C" cublasStatus_t cublasCherk3mEx(
    cublasHandle_t handle,
    cublasFillMode_t uplo, cublasOperation_t trans,
    int n, int k,
    const float *alpha,
    const void *A, cudaDataType Atype, int lda,
    const float *beta,
    void *C, cudaDataType Ctype, int ldc //
) {
    auto call_native = [&]() -> cublasStatus_t {
        auto fn = real_cherk3m_ex();
        if (!fn) return CUBLAS_STATUS_NOT_INITIALIZED;
        return fn(handle, uplo, trans, n, k, alpha, A, Atype, lda, beta, C, Ctype, ldc);
    };

    return herk_ex_common_impl(
        handle, uplo, trans,
        static_cast<int64_t>(n), static_cast<int64_t>(k),
        alpha, A, Atype, static_cast<int64_t>(lda),
        beta, C, Ctype, static_cast<int64_t>(ldc),
        call_native);
}

extern "C" cublasStatus_t cublasCherk3mEx_64(
    cublasHandle_t handle,
    cublasFillMode_t uplo, cublasOperation_t trans,
    int64_t n, int64_t k,
    const float *alpha,
    const void *A, cudaDataType Atype, int64_t lda,
    const float *beta,
    void *C, cudaDataType Ctype, int64_t ldc //
) {
    auto call_native = [&]() -> cublasStatus_t {
        auto fn = real_cherk3m_ex_64();
        if (!fn) return CUBLAS_STATUS_NOT_INITIALIZED;
        return fn(handle, uplo, trans, n, k, alpha, A, Atype, lda, beta, C, Ctype, ldc);
    };

    return herk_ex_common_impl(
        handle, uplo, trans,
        n, k,
        alpha, A, Atype, lda,
        beta, C, Ctype, ldc,
        call_native);
}
#endif
