#pragma once
#include "../common/common.hpp"
#include "../core/block_size.hpp"
#include "../core/blocking.hpp"
#include "../core/oz2_core.hpp"
#include "../gemm/worksize.hpp"
#include "worksize.hpp"

namespace gemmul8::oz2::herkx {

namespace impl {

template <cublasFillMode_t UPLO_C,
          typename TA, typename TB, typename TC,
          Backend BACKEND, unsigned NUM_MODULI>
std::vector<double> run(
    common::Handle_t handle,
    cublasOperation_t trans,
    size_t n, size_t k,
    const TC *alpha,
    const TA *A, size_t lda,
    const TB *B, size_t ldb,
    const common::underlying_t<TC> *beta,
    TC *C, size_t ldc,
    bool fastmode,
    void *work, void *workA, void *workB,
    bool enable_skip_scalA, bool enable_skip_scalB,
    bool skip_scalA, bool skip_scalB,
    cudaStream_t stream //
) {
    using U                          = common::underlying_t<TC>;
    constexpr common::MatStruct FULL = common::MatStruct::Full;
    constexpr cublasFillMode_t UF    = CUBLAS_FILL_MODE_FULL;
    constexpr cublasDiagType_t DN    = CUBLAS_DIAG_NON_UNIT;

    const bool memory_saving_mode = handle.config.memory_saving;
    const size_t limit            = handle.config.max_worksize;
    const bool memory_saving      = memory_saving_mode && limit > 0;
    const bool use_skipA          = memory_saving_mode ? false : enable_skip_scalA;
    const bool use_skipB          = memory_saving_mode ? false : enable_skip_scalB;
    const bool do_skipA           = memory_saving_mode ? false : skip_scalA;
    const bool do_skipB           = memory_saving_mode ? false : skip_scalB;
    const cublasOperation_t op_B  = (trans == CUBLAS_OP_N) ? CUBLAS_OP_C : CUBLAS_OP_N;

    const size_t full_worksize = workSize<true, BACKEND>(
        n, n, k, NUM_MODULI, use_skipA, use_skipB, nullptr, nullptr, fastmode);

    if (!memory_saving || full_worksize <= limit) {
        return core::oz2_core<Func::herkx, TA, TB, TC,
                              BACKEND, NUM_MODULI, TC, U,
                              FULL, FULL, UF, UF, DN, DN, UPLO_C>(
            handle, trans, op_B, n, n, k,
            alpha, A, lda, B, ldb, beta, C, ldc,
            fastmode, work, workA, workB,
            use_skipA, use_skipB, do_skipA, do_skipB, stream);
    }

    auto fits = [&](size_t nB, size_t kB, bool need_gemm) {
        if (workSize<true, BACKEND>(
                nB, nB, kB, NUM_MODULI, false, false, nullptr, nullptr, fastmode) > limit) {
            return false;
        }

        return !need_gemm ||
               gemm::workSize<true, BACKEND>(
                   nB, nB, kB, NUM_MODULI, false, false, nullptr, nullptr, fastmode) <= limit;
    };

    const core::BlockSize2D block = core::find_block_size_rankk<BACKEND, NUM_MODULI>(n, k, false, fits);
    if (!block) {
        assert(false && "The workspace-size limit is too small.");
        return std::vector<double>(4, 0.0);
    }

    core::blocking::OneScalar<U> one_real_storage;
    core::blocking::OneScalar<TC> one_complex_storage;
    core::blocking::RealToComplexScalar<TC> beta_complex_storage;

    const U *one_real      = nullptr;
    const TC *one_complex  = nullptr;
    const TC *beta_complex = nullptr;

    if (k > block.yB) {
        one_real = one_real_storage.get(beta ? static_cast<const void *>(beta) : static_cast<const void *>(alpha), stream);
        if (!one_real) {
            assert(false && "Failed to create beta=1 scalar for blocked HERKX.");
            return std::vector<double>(4, 0.0);
        }
    }
    if (n > block.xB) {
        beta_complex = beta_complex_storage.get(beta, stream);
        if (k > block.yB) one_complex = one_complex_storage.get(alpha, stream);
        if ((beta && !beta_complex) || (k > block.yB && !one_complex)) {
            assert(false && "Failed to create complex beta scalars for blocked HERKX.");
            one_real_storage.release(stream);
            one_complex_storage.release(stream);
            beta_complex_storage.release(stream);
            return std::vector<double>(4, 0.0);
        }
    }

    std::vector<double> timer(4, 0.0);
    core::blocking::for_each_rankk_tile(
        UPLO_C, n, k, block.xB, block.yB,
        [&](size_t i, size_t ni, size_t p, size_t kp, bool first) {
            TC *Cii      = C + i + i * ldc;
            const TA *Ai = core::blocking::matrix_block_ptr(A, lda, trans, i, p);
            const TB *Bi = core::blocking::matrix_block_ptr(B, ldb, trans, i, p);
            const auto t = core::oz2_core<Func::herkx, TA, TB, TC, BACKEND, NUM_MODULI, TC, U,
                                          FULL, FULL, UF, UF, DN, DN, UPLO_C>(
                handle, trans, op_B, ni, ni, kp,
                alpha, Ai, lda, Bi, ldb, first ? beta : one_real, Cii, ldc,
                fastmode, work, nullptr, nullptr,
                false, false, false, false, stream);
            core::blocking::add_timer(timer, t);
        },
        [&](size_t i, size_t ni, size_t j, size_t nj, size_t p, size_t kp, bool first) {
            TC *Cij      = C + i + j * ldc;
            const TA *Ai = core::blocking::matrix_block_ptr(A, lda, trans, i, p);
            const TB *Bj = core::blocking::matrix_block_ptr(B, ldb, trans, j, p);
            const auto t = core::oz2_core<Func::gemm, TA, TB, TC, BACKEND, NUM_MODULI>(
                handle, trans, op_B, ni, nj, kp,
                alpha, Ai, lda, Bj, ldb, first ? beta_complex : one_complex, Cij, ldc,
                fastmode, work, nullptr, nullptr,
                false, false, false, false, stream);
            core::blocking::add_timer(timer, t);
        });

    one_real_storage.release(stream);
    one_complex_storage.release(stream);
    beta_complex_storage.release(stream);
    return timer;
}

} // namespace impl

template <typename TA, typename TB, typename TC, Backend BACKEND, unsigned NUM_MODULI>
std::vector<double> herkx_core(
    common::Handle_t handle,
    cublasFillMode_t uplo, cublasOperation_t trans,
    size_t n, size_t k,
    const TC *alpha,
    const TA *const A, size_t lda,
    const TB *const B, size_t ldb,
    const common::underlying_t<TC> *beta,
    TC *const C, size_t ldc,
    bool fastmode,
    void *const work, void *const workA, void *const workB,
    bool enable_skip_scalA, bool enable_skip_scalB,
    bool skip_scalA, bool skip_scalB,
    cudaStream_t stream //
) {
    static_assert(common::isComplex<TA> && common::isComplex<TB> && common::isComplex<TC>,
                  "TA, TB, and TC must be all complex");

    if (uplo == CUBLAS_FILL_MODE_FULL || trans == CUBLAS_OP_T) {
        assert(false && "unsupported");
        return std::vector<double>(4, 0.0);
    }

    if (work == nullptr) {
        size_t wa = 0, wb = 0;
        const size_t total = workSize<true, BACKEND>(
            n, n, k, NUM_MODULI, enable_skip_scalA, enable_skip_scalB, &wa, &wb, fastmode);
        std::vector<double> timer(4, 0.0);
        timer[0] = static_cast<double>(total);
        timer[1] = static_cast<double>(wa);
        timer[2] = static_cast<double>(wb);
        return timer;
    }

    if (uplo == CUBLAS_FILL_MODE_UPPER) {
        return impl::run<CUBLAS_FILL_MODE_UPPER, TA, TB, TC, BACKEND, NUM_MODULI>(
            handle, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc,
            fastmode, work, workA, workB,
            enable_skip_scalA, enable_skip_scalB, skip_scalA, skip_scalB, stream);
    }
    return impl::run<CUBLAS_FILL_MODE_LOWER, TA, TB, TC, BACKEND, NUM_MODULI>(
        handle, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc,
        fastmode, work, workA, workB,
        enable_skip_scalA, enable_skip_scalB, skip_scalA, skip_scalB, stream);
}

} // namespace gemmul8::oz2::herkx
