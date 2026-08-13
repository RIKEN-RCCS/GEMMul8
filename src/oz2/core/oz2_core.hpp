#pragma once
#include "../common/common.hpp"
#include "../common/table.hpp"
#include "../common/timer.hpp"
#include "mod_hi2mid.hpp"
#include "undo_scaling.hpp"
#include "matmult.hpp"
#include "scaling.hpp"
#include "helper_matmult.hpp"

namespace gemmul8::oz2::core {

// gemm, symm, syr2k, syrkx, trmm, hemm, her2k, herkx, trtrmm
template <Func FUNC,
          typename TA, typename TB, typename TC,
          Backend BACKEND, unsigned NUM_MODULI,
          typename TAlpha            = TC,
          typename TBeta             = TC,
          common::MatStruct STRUCT_A = common::MatStruct::Full,
          common::MatStruct STRUCT_B = common::MatStruct::Full,
          cublasFillMode_t UPLO_A    = CUBLAS_FILL_MODE_FULL,
          cublasFillMode_t UPLO_B    = CUBLAS_FILL_MODE_FULL,
          cublasDiagType_t DIAG_A    = CUBLAS_DIAG_NON_UNIT,
          cublasDiagType_t DIAG_B    = CUBLAS_DIAG_NON_UNIT,
          cublasFillMode_t UPLO_C    = CUBLAS_FILL_MODE_FULL>
inline std::vector<double> oz2_core(
    common::Handle_t handle,
    cublasOperation_t op_A, cublasOperation_t op_B,
    size_t m, size_t n, size_t k,
    const TAlpha *alpha,
    const TA *const A, size_t lda,
    const TB *const B, size_t ldb,
    const TBeta *beta,
    TC *const C, size_t ldc,
    bool fastmode,
    void *const work, void *const workA, void *const workB,
    bool enable_skip_scalA, bool enable_skip_scalB,
    bool skip_scalA, bool skip_scalB,
    cudaStream_t stream //
) {
    constexpr bool COMPLEX = common::isComplex<TA>;
    using LowT             = common::low_t<BACKEND>;
    using MidT             = common::mid_t<BACKEND, COMPLEX>;
    using HiT              = common::hi_t<BACKEND>;

    // Set constants
    const size_t m_pad     = common::padding(m);
    const size_t n_pad     = common::padding(n);
    const size_t k_pad     = common::padding(k);
    const size_t n_work    = (FUNC == Func::trtrmm) ? n_pad : n;
    const size_t lda_lo    = k_pad;
    const size_t ldb_lo    = k_pad;
    const size_t ldc_hi    = m_pad;
    const size_t sizeA     = k_pad * m_pad;
    const size_t sizeB     = k_pad * n_work;
    const size_t sizeC     = m_pad * n_work;
    const size_t size_vecA = m_pad;
    const size_t size_vecB = n_pad;
    const bool skipA       = skip_scalA && enable_skip_scalA;
    const bool skipB       = skip_scalB && enable_skip_scalB;

    // Set workspace
    constexpr common::MatMulKind KIND               = matmul_kind<FUNC, STRUCT_A, STRUCT_B, UPLO_C>();
    constexpr bool use_pointer_arrays               = common::isCUDA && (BACKEND == Backend::FP8) && (KIND == common::MatMulKind::Gemm);
    constexpr unsigned matmul_per_moduli            = COMPLEX ? 3u : 1u;
    constexpr unsigned num_C_hi                     = (BACKEND == Backend::INT8) ? 1u : 3u;
    constexpr unsigned num_mat                      = common::table::num_mat_v<BACKEND, NUM_MODULI>;
    constexpr unsigned pointer_products_per_modulus = (use_pointer_arrays) ? ((COMPLEX) ? 9u : 3u) : 0u;
    const size_t pointer_array_bytes                = (use_pointer_arrays)
                                                          ? common::padding(size_t(3 * pointer_products_per_modulus * NUM_MODULI) * sizeof(void *))
                                                          : 0u;

    const size_t offsetA        = sizeA * num_mat * matmul_per_moduli;
    const size_t offsetB        = sizeB * num_mat * matmul_per_moduli;
    constexpr size_t lwork_blas = size_t(32) << 20; // 32 MiB
    const size_t sizeC_Mid      = sizeof(MidT) * sizeC;
    const size_t worksizeC_1    = (NUM_MODULI - 1) * sizeC_Mid;
    const size_t worksizeC_2    = std::max<size_t>(pointer_array_bytes + lwork_blas, sizeC_Mid);
    const size_t incC_hi        = matmul_per_moduli * num_C_hi * sizeC;
    const size_t sizeC_Hi       = sizeof(HiT) * incC_hi;
    const size_t worksizeC      = worksizeC_1 + worksizeC_2 + sizeC_Hi;

    void *const work_aligned  = common::align(work);
    void *const workA_aligned = (workA) ? common::align(workA) : nullptr;
    void *const workB_aligned = (workB) ? common::align(workB) : nullptr;

    const size_t incA_lo      = offsetA + ((enable_skip_scalA) ? (sizeA * matmul_per_moduli) : 0u);
    const size_t incsftA      = size_vecA * ((enable_skip_scalA) ? 2u : 1u);
    const size_t incB_lo      = offsetB + ((enable_skip_scalB) ? (sizeB * matmul_per_moduli) : 0u);
    const size_t incsftB      = size_vecB * ((enable_skip_scalB) ? 2u : 1u);
    LowT *const A_lo_base     = reinterpret_cast<LowT *>((workA) ? workA_aligned : work_aligned);
    int16_t *const sftA       = reinterpret_cast<int16_t *>(A_lo_base + incA_lo);
    LowT *const B_lo_base     = reinterpret_cast<LowT *>((workB) ? workB_aligned : ((workA) ? work_aligned : (sftA + incsftA)));
    int16_t *const sftB       = reinterpret_cast<int16_t *>(B_lo_base + incB_lo);
    MidT *const C_mid         = reinterpret_cast<MidT *>((workB) ? ((workA) ? work_aligned : (sftA + incsftA)) : (sftB + incsftB));
    int8_t *const C_work_base = reinterpret_cast<int8_t *>(C_mid);

    // Set handle
    cublasFillMode_t UPLO_A_handle = UPLO_A;
    cublasFillMode_t UPLO_B_handle = UPLO_B;
    if constexpr (KIND == common::MatMulKind::Trtrmm) {
        if (op_A != CUBLAS_OP_N) {
            UPLO_A_handle = (UPLO_A == CUBLAS_FILL_MODE_UPPER) ? CUBLAS_FILL_MODE_LOWER : CUBLAS_FILL_MODE_UPPER;
        }
        if (op_B != CUBLAS_OP_N) {
            UPLO_B_handle = (UPLO_B == CUBLAS_FILL_MODE_UPPER) ? CUBLAS_FILL_MODE_LOWER : CUBLAS_FILL_MODE_UPPER;
        }
    }
    common::set_handle<BACKEND, FUNC>(stream, handle, ldc_hi, n_work, lda_lo, ldb_lo, ldb_lo, ldc_hi, lwork_blas, UPLO_A_handle, UPLO_B_handle);
    if constexpr (BACKEND == Backend::FP8) {
        handle.fp8_k_blocking = false;
    }

    // set timer
    std::vector<double> timer(4, 0.0);
    common::Timer<NUM_MODULI> phase_timer;

    // A,B --(scaling)--> Aint,Bint --(mod)--> A_lo,B_lo
    phase_timer.record(phase_timer.scaling_begin(), stream);

    auto A_lo = common::make_matptr<LowT, COMPLEX>(A_lo_base, sizeA * num_mat);
    auto B_lo = common::make_matptr<LowT, COMPLEX>(B_lo_base, sizeB * num_mat);

    if (fastmode) {

        scaling_fast<TA, TB, BACKEND, NUM_MODULI, STRUCT_A, STRUCT_B, UPLO_A, UPLO_B, DIAG_A, DIAG_B>(
            stream, op_A, op_B, m, n, k,
            A, lda, A_lo, lda_lo, sizeA, sftA, skipA,
            B, ldb, B_lo, ldb_lo, sizeB, sftB, skipB);

    } else {

        auto A_lo_high   = common::make_matptr<LowT, COMPLEX>(A_lo_base + ((enable_skip_scalA) ? offsetA : 0u), sizeA);
        auto B_lo_high   = common::make_matptr<LowT, COMPLEX>(B_lo_base + ((enable_skip_scalB) ? offsetB : 0u), sizeB);
        auto C_hi        = common::make_matptr<HiT, COMPLEX>(reinterpret_cast<HiT *>(C_work_base), sizeC * num_C_hi);
        handle.workspace = static_cast<void *>(C_work_base + sizeC_Hi);

        scaling_accu<FUNC, TA, TB, BACKEND, NUM_MODULI, STRUCT_A, STRUCT_B, UPLO_A, UPLO_B, DIAG_A, DIAG_B, UPLO_C>(
            stream, handle, op_A, op_B, m, n, k,
            A, lda, A_lo, A_lo_high, lda_lo, sizeA, sftA, skipA, enable_skip_scalA,
            B, ldb, B_lo, B_lo_high, ldb_lo, sizeB, sftB, skipB, enable_skip_scalB,
            C_hi, ldc_hi);
    }

    phase_timer.record(phase_timer.scaling_end(), stream);

    for (unsigned i = 0; i < NUM_MODULI;) {

        unsigned bcnt = batch_count<NUM_MODULI>(
            handle.arch, n_work, i, sizeC_Mid, sizeC_Hi, lwork_blas, worksizeC, pointer_products_per_modulus);

        if constexpr (BACKEND == Backend::FP8) {
            bcnt                       = limit_fp8_batch_for_k<NUM_MODULI>(i, bcnt, k_pad);
            const int KB               = int(common::table::k_block_first_fp8[i]);
            const float p              = float(common::table::moduli_fp8[i]);
            handle.fp8_k_blocking      = (k_pad > size_t(KB));
            handle.fp8_k_block_first   = KB;
            handle.fp8_k_block_next    = common::table::k_block_next_fp8[i];
            handle.fp8_modulus         = p;
            handle.fp8_neg_inv_modulus = float(-1.0 / double(common::table::moduli_fp8[i]));
            handle.fp8_half_modulus    = 0.5f * p;
        }

        int8_t *const base = C_work_base + size_t(i) * sizeC_Mid;

        const size_t pointer_batch_count = use_pointer_arrays ? pointer_products_per_modulus * bcnt : 0u;
        const size_t ptr_bytes           = use_pointer_arrays ? 3 * pointer_batch_count * sizeof(void *) : 0u;
        const size_t ptr_gap             = use_pointer_arrays ? common::padding(ptr_bytes) : 0u;

        void **const Aarray = use_pointer_arrays ? reinterpret_cast<void **>(base) : nullptr;
        void **const Barray = use_pointer_arrays ? reinterpret_cast<void **>(base + pointer_batch_count * sizeof(void *)) : nullptr;
        void **const Carray = use_pointer_arrays ? reinterpret_cast<void **>(base + 2 * pointer_batch_count * sizeof(void *)) : nullptr;

        const size_t mid_bytes = bcnt * sizeC_Mid;
        const size_t gap       = std::max<size_t>(ptr_gap + lwork_blas, mid_bytes);

        HiT *const C_hi  = reinterpret_cast<HiT *>(base + gap);
        auto C_hi_matptr = common::make_matptr<HiT, COMPLEX>(C_hi, sizeC * num_C_hi);

        handle.Aarray    = Aarray;
        handle.Barray    = Barray;
        handle.Carray    = Carray;
        handle.workspace = static_cast<void *>(base + ptr_gap);

        // Error-free matrix multiplication C_hi := A_lo[i]*B_lo[i]
        error_free_matmult<FUNC, BACKEND, COMPLEX, STRUCT_A, STRUCT_B, UPLO_A, UPLO_B, UPLO_C>(
            stream, handle, i, bcnt,
            ldc_hi, n_work, lda_lo, ldb_lo,
            sizeA, sizeB, sizeC, op_A, op_B,
            A_lo, B_lo, C_hi_matptr);

        const unsigned gid = phase_timer.begin_group();
        phase_timer.record(phase_timer.mm_end(gid), stream);

        // C_mid[i] := mod(C_hi, p[i])
        mod_hi2mid<FUNC, BACKEND, COMPLEX, UPLO_A, UPLO_B, UPLO_C>(
            stream, op_A, op_B, i, bcnt, ldc_hi, n_work, C_hi_matptr, C_mid, incC_hi);

        phase_timer.record(phase_timer.mod_hi2mid_end(gid), stream);

        i += bcnt;
    }

    // C_mid --(CRT Reduction)--> Aint*Bint --(undo scaling)--> A*B
    // and C := alpha*AB + beta*C
    undo_scaling<FUNC, TC, BACKEND, NUM_MODULI, TAlpha, TBeta, UPLO_A, UPLO_B, UPLO_C>(
        stream, op_A, op_B, m, n, C_mid, ldc_hi, sizeC, C, ldc, sftA, sftB, alpha, beta);

    phase_timer.record(phase_timer.undo_scaling_end(), stream);

    // Collect phase times.
    phase_timer.collect(timer);

    // Clean up Lt memory
    common::cleanup_handle(handle);

    return timer;
}

// syrk, herk
template <Func FUNC,
          typename TA, typename TC,
          Backend BACKEND, unsigned NUM_MODULI,
          typename TAlpha         = TC,
          typename TBeta          = TC,
          cublasFillMode_t UPLO_C = CUBLAS_FILL_MODE_FULL>
inline std::vector<double> oz2_core_rk(
    common::Handle_t handle,
    cublasOperation_t op_A,
    size_t n, size_t k,
    const TAlpha *alpha,
    const TA *const A, size_t lda,
    const TBeta *beta,
    TC *const C, size_t ldc,
    bool fastmode,
    void *const work, void *const workA,
    bool enable_skip_scalA,
    bool skip_scalA,
    cudaStream_t stream //
) {
    static_assert(FUNC == Func::syrk || FUNC == Func::herk,
                  "oz2_core_rk supports only syrk/herk.");
    if constexpr (FUNC == Func::herk) {
        static_assert(common::isComplex<TA>, "herk requires complex input type.");
    }
    static_assert(UPLO_C == CUBLAS_FILL_MODE_UPPER || UPLO_C == CUBLAS_FILL_MODE_LOWER,
                  "oz2_core_rk requires UPLO_C = UPPER or LOWER.");

    constexpr bool COMPLEX            = common::isComplex<TA>;
    using LowT                        = common::low_t<BACKEND>;
    using MidT                        = common::mid_t<BACKEND, COMPLEX>;
    using HiT                         = common::hi_t<BACKEND>;
    constexpr cublasFillMode_t UPLO_A = CUBLAS_FILL_MODE_FULL;
    constexpr cublasFillMode_t UPLO_B = CUBLAS_FILL_MODE_FULL;

    const size_t n_pad     = common::padding(n);
    const size_t k_pad     = common::padding(k);
    const size_t lda_lo    = k_pad;
    const size_t ldc_hi    = n_pad;
    const size_t sizeA     = k_pad * n_pad;
    const size_t sizeC     = n_pad * n;
    const size_t size_vecA = n_pad;
    const bool skipA       = skip_scalA && enable_skip_scalA;

    // Set workspace
    constexpr unsigned matmul_per_moduli = COMPLEX ? 3u : 1u;
    constexpr unsigned num_C_hi          = (BACKEND == Backend::INT8) ? 1u : 3u;
    constexpr unsigned num_mat           = common::table::num_mat_v<BACKEND, NUM_MODULI>;

    const size_t offsetA        = sizeA * num_mat * matmul_per_moduli;
    constexpr size_t lwork_blas = size_t(32) << 20; // 32 MiB
    const size_t sizeC_Mid      = sizeof(MidT) * sizeC;
    const size_t worksizeC_1    = (NUM_MODULI - 1) * sizeC_Mid;
    const size_t worksizeC_2    = std::max<size_t>(lwork_blas, sizeC_Mid);
    const size_t incC_hi        = matmul_per_moduli * num_C_hi * sizeC;
    const size_t sizeC_Hi       = sizeof(HiT) * incC_hi;
    const size_t worksizeC      = worksizeC_1 + worksizeC_2 + sizeC_Hi;

    void *const work_aligned  = common::align(work);
    void *const workA_aligned = (workA) ? common::align(workA) : nullptr;

    const size_t incA_lo = offsetA;
    const size_t incsftA = size_vecA;

    LowT *const A_lo_base     = reinterpret_cast<LowT *>((workA) ? workA_aligned : work_aligned);
    int16_t *const sftA       = reinterpret_cast<int16_t *>(A_lo_base + incA_lo);
    MidT *const C_mid         = reinterpret_cast<MidT *>((workA) ? work_aligned : (sftA + incsftA));
    int8_t *const C_work_base = reinterpret_cast<int8_t *>(C_mid);

    // Set handle
    common::set_handle<BACKEND, FUNC>(stream, handle, ldc_hi, n, lda_lo, lda_lo, lda_lo, ldc_hi, lwork_blas, UPLO_A, UPLO_B);
    handle.Aarray = nullptr;
    handle.Barray = nullptr;
    handle.Carray = nullptr;
    if constexpr (BACKEND == Backend::FP8) {
        handle.fp8_k_blocking = false;
    }

    // set timer
    std::vector<double> timer(4, 0.0);
    common::Timer<NUM_MODULI> phase_timer;

    // A,B --(scaling)--> Aint,Bint --(mod)--> A_lo,B_lo
    phase_timer.record(phase_timer.scaling_begin(), stream);

    auto A_lo = common::make_matptr<LowT, COMPLEX>(A_lo_base, sizeA * num_mat);

    if (fastmode) {

        scaling_fast_rk<FUNC, TA, BACKEND, NUM_MODULI>(
            stream, op_A, n, k, A, lda, A_lo, lda_lo, sizeA, sftA, skipA);

    } else {

        auto A_lo_high   = common::make_matptr<LowT, COMPLEX>(A_lo_base, sizeA);
        auto C_hi        = common::make_matptr<HiT, COMPLEX>(reinterpret_cast<HiT *>(C_work_base), sizeC * num_C_hi);
        handle.workspace = static_cast<void *>(C_work_base + sizeC_Hi);

        scaling_accu_rk<FUNC, TA, BACKEND, NUM_MODULI, UPLO_C>(
            stream, handle, op_A, n, k, A, lda, A_lo, A_lo_high, lda_lo, sizeA, sftA, skipA, C_hi, ldc_hi);
    }

    phase_timer.record(phase_timer.scaling_end(), stream);

    for (unsigned i = 0; i < NUM_MODULI;) {

        unsigned bcnt = batch_count<NUM_MODULI>(
            handle.arch, n, i, sizeC_Mid, sizeC_Hi, lwork_blas, worksizeC, 0u);

        if constexpr (BACKEND == Backend::FP8) {
            bcnt                       = limit_fp8_batch_for_k<NUM_MODULI>(i, bcnt, k_pad);
            const int KB               = int(common::table::k_block_first_fp8[i]);
            const float p              = float(common::table::moduli_fp8[i]);
            handle.fp8_k_blocking      = (k_pad > size_t(KB));
            handle.fp8_k_block_first   = KB;
            handle.fp8_k_block_next    = common::table::k_block_next_fp8[i];
            handle.fp8_modulus         = p;
            handle.fp8_neg_inv_modulus = float(-1.0 / double(common::table::moduli_fp8[i]));
            handle.fp8_half_modulus    = 0.5f * p;
        }

        const size_t mid_bytes = bcnt * sizeC_Mid;
        const size_t gap       = std::max<size_t>(lwork_blas, mid_bytes);

        int8_t *const base = C_work_base + size_t(i) * sizeC_Mid;
        HiT *const C_hi    = reinterpret_cast<HiT *>(base + gap);
        auto C_hi_matptr   = common::make_matptr<HiT, COMPLEX>(C_hi, sizeC * num_C_hi);
        handle.workspace   = static_cast<void *>(base);

        // Error-free matrix multiplication C_hi := A_lo[i]*B_lo[i]
        error_free_matmult_rk<FUNC, BACKEND, COMPLEX, UPLO_C>(
            stream, handle, i, bcnt, ldc_hi, n, lda_lo, sizeA, sizeC, A_lo, C_hi_matptr);

        const unsigned gid = phase_timer.begin_group();
        phase_timer.record(phase_timer.mm_end(gid), stream);

        // C_mid[i] := mod(C_hi, p[i])
        mod_hi2mid<FUNC, BACKEND, COMPLEX, UPLO_A, UPLO_B, UPLO_C>(
            stream, op_A, CUBLAS_OP_N, i, bcnt, ldc_hi, n, C_hi_matptr, C_mid, incC_hi);

        phase_timer.record(phase_timer.mod_hi2mid_end(gid), stream);

        i += bcnt;
    }

    // C_mid --(CRT Reduction)--> Aint*Bint --(undo scaling)--> A*B
    // and C := alpha*AB + beta*C
    undo_scaling<FUNC, TC, BACKEND, NUM_MODULI, TAlpha, TBeta, UPLO_A, UPLO_B, UPLO_C>(
        stream, op_A, CUBLAS_OP_N, n, n, C_mid, ldc_hi, sizeC, C, ldc, sftA, sftA, alpha, beta);

    phase_timer.record(phase_timer.undo_scaling_end(), stream);

    // Collect phase times.
    phase_timer.collect(timer);

    // Clean up Lt memory
    common::cleanup_handle(handle);

    return timer;
}

} // namespace gemmul8::oz2::core
