#pragma once
#include "../common/common.hpp"
#include "../core/block_size.hpp"
#include "../core/blocking.hpp"
#include "../core/oz2_core.hpp"
#include "../gemm/worksize.hpp"
#include "../trmm/worksize.hpp"
#include "worksize.hpp"

namespace gemmul8::oz2::trtrmm {

namespace impl {

template <cublasFillMode_t UA, cublasFillMode_t UB,
          cublasDiagType_t DA, cublasDiagType_t DB,
          typename TA, typename TB, typename TC,
          Backend BACKEND, unsigned NUM_MODULI>
std::vector<double> run(
    common::Handle_t handle,
    cublasOperation_t op_A, cublasOperation_t op_B,
    size_t n,
    const TC *alpha,
    const TA *A, size_t lda,
    const TB *B, size_t ldb,
    const TC *beta,
    TC *C, size_t ldc,
    bool fastmode,
    void *work, void *workA, void *workB,
    bool enable_skip_scalA, bool enable_skip_scalB,
    bool skip_scalA, bool skip_scalB,
    cudaStream_t stream //
) {
    constexpr common::MatStruct TRI  = common::MatStruct::Triangular;
    constexpr common::MatStruct FULL = common::MatStruct::Full;
    constexpr cublasFillMode_t UF    = CUBLAS_FILL_MODE_FULL;
    constexpr cublasDiagType_t DN    = CUBLAS_DIAG_NON_UNIT;
    constexpr cublasFillMode_t UC    = CUBLAS_FILL_MODE_FULL;

    const bool memory_saving_mode = handle.config.memory_saving;
    const size_t limit            = handle.config.max_worksize;
    const bool memory_saving      = memory_saving_mode && limit > 0;
    const bool use_skipA          = memory_saving_mode ? false : enable_skip_scalA;
    const bool use_skipB          = memory_saving_mode ? false : enable_skip_scalB;
    const bool do_skipA           = memory_saving_mode ? false : skip_scalA;
    const bool do_skipB           = memory_saving_mode ? false : skip_scalB;

    const size_t full_worksize = workSize<common::isComplex<TA>, BACKEND>(
        n, n, n, NUM_MODULI, use_skipA, use_skipB, nullptr, nullptr, fastmode);

    if (!memory_saving || full_worksize <= limit) {
        return core::oz2_core<Func::trtrmm, TA, TB, TC,
                              BACKEND, NUM_MODULI, TC, TC,
                              TRI, TRI, UA, UB, DA, DB, UC>(
            handle, op_A, op_B, n, n, n,
            alpha, A, lda, B, ldb, beta, C, ldc,
            fastmode, work, workA, workB,
            use_skipA, use_skipB, do_skipA, do_skipB, stream);
    }

    const cublasFillMode_t effA = core::blocking::effective_uplo(UA, op_A);
    const cublasFillMode_t effB = core::blocking::effective_uplo(UB, op_B);

    auto fits_no_gemm = [&](size_t b) {
        if (workSize<common::isComplex<TA>, BACKEND>(
                b, b, b, NUM_MODULI, false, false, nullptr, nullptr, fastmode) > limit) {
            return false;
        }

        return trmm::workSize<common::isComplex<TA>, BACKEND>(
                   b, b, b, NUM_MODULI, false, false, nullptr, nullptr, fastmode) <= limit;
    };

    auto fits = [&](size_t b) {
        return fits_no_gemm(b) &&
               gemm::workSize<common::isComplex<TA>, BACKEND>(
                   b, b, b, NUM_MODULI, false, false, nullptr, nullptr, fastmode) <= limit;
    };

    size_t sB = 0;

    if (effA == effB) {
        const size_t b2 = core::block_size_from_count(n, 2);
        if (b2 < n && core::block_count(n, b2) == 2 && fits_no_gemm(b2)) {
            sB = b2;
        }
    }

    if (sB == 0) {
        sB = core::find_block_size_trtrmm(n, fits);
    }

    if (sB == 0) {
        assert(false && "The workspace-size limit is too small.");
        return std::vector<double>(4, 0.0);
    }

    core::blocking::OneScalar<TC> one_storage;
    const TC *one = one_storage.get(beta ? static_cast<const void *>(beta) : static_cast<const void *>(alpha), stream);
    if (!one) {
        assert(false && "Failed to create beta=1 scalar for blocked TRTRMM.");
        return std::vector<double>(4, 0.0);
    }

    std::vector<double> timer(4, 0.0);

    for (size_t i = 0; i < n; i += sB) {
        const size_t mi = std::min(sB, n - i);
        for (size_t j = 0; j < n; j += sB) {
            const size_t nj = std::min(sB, n - j);
            TC *Cij         = C + i + j * ldc;
            bool first      = true;

            const auto p_range = core::blocking::triangular_product_p_range(effA, effB, i, j, n, sB);

            for (size_t p = p_range.begin; p < p_range.end; p += sB) {
                const size_t kp      = std::min(sB, n - p);
                const TA *Aip        = core::blocking::matrix_block_ptr(A, lda, op_A, i, p);
                const TB *Bpj        = core::blocking::matrix_block_ptr(B, ldb, op_B, p, j);
                const TC *beta_block = first ? beta : one;

                std::vector<double> t;
                if (i == p && p == j) {
                    t = core::oz2_core<Func::trtrmm, TA, TB, TC,
                                       BACKEND, NUM_MODULI, TC, TC,
                                       TRI, TRI, UA, UB, DA, DB, UC>(
                        handle, op_A, op_B, mi, nj, kp,
                        alpha, Aip, lda, Bpj, ldb, beta_block, Cij, ldc,
                        fastmode, work, nullptr, nullptr,
                        false, false, false, false, stream);
                } else if (i == p) {
                    t = core::oz2_core<Func::trmm, TA, TB, TC,
                                       BACKEND, NUM_MODULI, TC, TC,
                                       TRI, FULL, UA, UF, DA, DN>(
                        handle, op_A, op_B, mi, nj, kp,
                        alpha, Aip, lda, Bpj, ldb, beta_block, Cij, ldc,
                        fastmode, work, nullptr, nullptr,
                        false, false, false, false, stream);
                } else if (p == j) {
                    t = core::oz2_core<Func::trmm, TA, TB, TC,
                                       BACKEND, NUM_MODULI, TC, TC,
                                       FULL, TRI, UF, UB, DN, DB>(
                        handle, op_A, op_B, mi, nj, kp,
                        alpha, Aip, lda, Bpj, ldb, beta_block, Cij, ldc,
                        fastmode, work, nullptr, nullptr,
                        false, false, false, false, stream);
                } else {
                    t = core::oz2_core<Func::gemm, TA, TB, TC, BACKEND, NUM_MODULI>(
                        handle, op_A, op_B, mi, nj, kp,
                        alpha, Aip, lda, Bpj, ldb, beta_block, Cij, ldc,
                        fastmode, work, nullptr, nullptr,
                        false, false, false, false, stream);
                }

                core::blocking::add_timer(timer, t);
                first = false;
            }

            if (first) {
                core::blocking::scale_block(stream, mi, nj, beta, Cij, ldc);
            }
        }
    }

    one_storage.release(stream);
    return timer;
}

} // namespace impl

template <typename TA, typename TB, typename TC, Backend BACKEND, unsigned NUM_MODULI>
std::vector<double> trtrmm_core(
    common::Handle_t handle,
    cublasFillMode_t uplo_A, cublasFillMode_t uplo_B,
    cublasOperation_t trans_A, cublasOperation_t trans_B,
    cublasDiagType_t diag_A, cublasDiagType_t diag_B,
    size_t n,
    const TC *alpha,
    const TA *const A, size_t lda,
    const TB *const B, size_t ldb,
    const TC *beta,
    TC *const C, size_t ldc,
    bool fastmode,
    void *const work, void *const workA, void *const workB,
    bool enable_skip_scalA, bool enable_skip_scalB,
    bool skip_scalA, bool skip_scalB,
    cudaStream_t stream //
) {
    static_assert(common::isComplex<TA> == common::isComplex<TB> &&
                      common::isComplex<TB> == common::isComplex<TC>,
                  "TA, TB, and TC must be all real or all complex");

    if (uplo_A == CUBLAS_FILL_MODE_FULL || uplo_B == CUBLAS_FILL_MODE_FULL) {
        assert(false && "unsupported");
        return std::vector<double>(4, 0.0);
    }

    if (work == nullptr) {
        size_t wa = 0, wb = 0;
        const size_t total = workSize<common::isComplex<TA>, BACKEND>(
            n, n, n, NUM_MODULI, enable_skip_scalA, enable_skip_scalB, &wa, &wb, fastmode);
        std::vector<double> timer(4, 0.0);
        timer[0] = static_cast<double>(total);
        timer[1] = static_cast<double>(wa);
        timer[2] = static_cast<double>(wb);
        return timer;
    }

#define GEMMUL8_TRTRMM_RUN(UA_, UB_, DA_, DB_)                              \
    return impl::run<UA_, UB_, DA_, DB_, TA, TB, TC, BACKEND, NUM_MODULI>(  \
        handle, trans_A, trans_B, n, alpha, A, lda, B, ldb, beta, C, ldc,   \
        fastmode, work, workA, workB, enable_skip_scalA, enable_skip_scalB, \
        skip_scalA, skip_scalB, stream)

    if (uplo_A == CUBLAS_FILL_MODE_UPPER) {
        if (uplo_B == CUBLAS_FILL_MODE_UPPER) {
            if (diag_A == CUBLAS_DIAG_NON_UNIT) {
                if (diag_B == CUBLAS_DIAG_NON_UNIT) GEMMUL8_TRTRMM_RUN(CUBLAS_FILL_MODE_UPPER, CUBLAS_FILL_MODE_UPPER, CUBLAS_DIAG_NON_UNIT, CUBLAS_DIAG_NON_UNIT);
                GEMMUL8_TRTRMM_RUN(CUBLAS_FILL_MODE_UPPER, CUBLAS_FILL_MODE_UPPER, CUBLAS_DIAG_NON_UNIT, CUBLAS_DIAG_UNIT);
            }
            if (diag_B == CUBLAS_DIAG_NON_UNIT) GEMMUL8_TRTRMM_RUN(CUBLAS_FILL_MODE_UPPER, CUBLAS_FILL_MODE_UPPER, CUBLAS_DIAG_UNIT, CUBLAS_DIAG_NON_UNIT);
            GEMMUL8_TRTRMM_RUN(CUBLAS_FILL_MODE_UPPER, CUBLAS_FILL_MODE_UPPER, CUBLAS_DIAG_UNIT, CUBLAS_DIAG_UNIT);
        }
        if (diag_A == CUBLAS_DIAG_NON_UNIT) {
            if (diag_B == CUBLAS_DIAG_NON_UNIT) GEMMUL8_TRTRMM_RUN(CUBLAS_FILL_MODE_UPPER, CUBLAS_FILL_MODE_LOWER, CUBLAS_DIAG_NON_UNIT, CUBLAS_DIAG_NON_UNIT);
            GEMMUL8_TRTRMM_RUN(CUBLAS_FILL_MODE_UPPER, CUBLAS_FILL_MODE_LOWER, CUBLAS_DIAG_NON_UNIT, CUBLAS_DIAG_UNIT);
        }
        if (diag_B == CUBLAS_DIAG_NON_UNIT) GEMMUL8_TRTRMM_RUN(CUBLAS_FILL_MODE_UPPER, CUBLAS_FILL_MODE_LOWER, CUBLAS_DIAG_UNIT, CUBLAS_DIAG_NON_UNIT);
        GEMMUL8_TRTRMM_RUN(CUBLAS_FILL_MODE_UPPER, CUBLAS_FILL_MODE_LOWER, CUBLAS_DIAG_UNIT, CUBLAS_DIAG_UNIT);
    }

    if (uplo_B == CUBLAS_FILL_MODE_UPPER) {
        if (diag_A == CUBLAS_DIAG_NON_UNIT) {
            if (diag_B == CUBLAS_DIAG_NON_UNIT) GEMMUL8_TRTRMM_RUN(CUBLAS_FILL_MODE_LOWER, CUBLAS_FILL_MODE_UPPER, CUBLAS_DIAG_NON_UNIT, CUBLAS_DIAG_NON_UNIT);
            GEMMUL8_TRTRMM_RUN(CUBLAS_FILL_MODE_LOWER, CUBLAS_FILL_MODE_UPPER, CUBLAS_DIAG_NON_UNIT, CUBLAS_DIAG_UNIT);
        }
        if (diag_B == CUBLAS_DIAG_NON_UNIT) GEMMUL8_TRTRMM_RUN(CUBLAS_FILL_MODE_LOWER, CUBLAS_FILL_MODE_UPPER, CUBLAS_DIAG_UNIT, CUBLAS_DIAG_NON_UNIT);
        GEMMUL8_TRTRMM_RUN(CUBLAS_FILL_MODE_LOWER, CUBLAS_FILL_MODE_UPPER, CUBLAS_DIAG_UNIT, CUBLAS_DIAG_UNIT);
    }

    if (diag_A == CUBLAS_DIAG_NON_UNIT) {
        if (diag_B == CUBLAS_DIAG_NON_UNIT) GEMMUL8_TRTRMM_RUN(CUBLAS_FILL_MODE_LOWER, CUBLAS_FILL_MODE_LOWER, CUBLAS_DIAG_NON_UNIT, CUBLAS_DIAG_NON_UNIT);
        GEMMUL8_TRTRMM_RUN(CUBLAS_FILL_MODE_LOWER, CUBLAS_FILL_MODE_LOWER, CUBLAS_DIAG_NON_UNIT, CUBLAS_DIAG_UNIT);
    }
    if (diag_B == CUBLAS_DIAG_NON_UNIT) GEMMUL8_TRTRMM_RUN(CUBLAS_FILL_MODE_LOWER, CUBLAS_FILL_MODE_LOWER, CUBLAS_DIAG_UNIT, CUBLAS_DIAG_NON_UNIT);
    GEMMUL8_TRTRMM_RUN(CUBLAS_FILL_MODE_LOWER, CUBLAS_FILL_MODE_LOWER, CUBLAS_DIAG_UNIT, CUBLAS_DIAG_UNIT);

#undef GEMMUL8_TRTRMM_RUN
}

} // namespace gemmul8::oz2::trtrmm
