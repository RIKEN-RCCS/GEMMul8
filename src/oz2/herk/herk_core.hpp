#pragma once
#include "../common/common.hpp"
#include "../core/block_size.hpp"
#include "../core/blocking.hpp"
#include "../core/oz2_core.hpp"
#include "../gemm/worksize.hpp"
#include "worksize.hpp"

namespace gemmul8::oz2::herk {

namespace impl {

template <cublasFillMode_t UPLO_C,
          typename TA, typename TC,
          Backend BACKEND, unsigned NUM_MODULI>
std::vector<double> run(
    common::Handle_t handle,
    cublasOperation_t trans,
    size_t n, size_t k,
    const common::underlying_t<TC> *alpha,
    const TA *A, size_t lda,
    const common::underlying_t<TC> *beta,
    TC *C, size_t ldc,
    bool fastmode,
    void *work, void *workA,
    bool enable_skip_scalA, bool skip_scalA,
    cudaStream_t stream //
) {
    using U = common::underlying_t<TC>;

    const bool memory_saving_mode = handle.config.memory_saving;
    const size_t limit            = handle.config.max_worksize;
    const bool memory_saving      = memory_saving_mode && limit > 0;
    const bool use_skipA          = memory_saving_mode ? false : enable_skip_scalA;
    const bool do_skipA           = memory_saving_mode ? false : skip_scalA;

    const size_t full_worksize = workSize<true, BACKEND>(
        n, n, k, NUM_MODULI, use_skipA, false, nullptr, nullptr, fastmode);

    if (!memory_saving || full_worksize <= limit) {
        return core::oz2_core_rk<Func::herk, TA, TC, BACKEND, NUM_MODULI, U, U, UPLO_C>(
            handle, trans, n, k,
            alpha, A, lda, beta, C, ldc,
            fastmode, work, workA,
            use_skipA, do_skipA, stream);
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

    const size_t nB              = block.xB;
    const size_t kB              = block.yB;
    const cublasOperation_t op_B = (trans == CUBLAS_OP_N) ? CUBLAS_OP_C : CUBLAS_OP_N;

    core::blocking::OneScalar<U> one_real_storage;
    core::blocking::OneScalar<TC> one_complex_storage;
    core::blocking::RealToComplexScalar<TC> alpha_complex_storage;
    core::blocking::RealToComplexScalar<TC> beta_complex_storage;

    const U *one_real       = nullptr;
    const TC *one_complex   = nullptr;
    const TC *alpha_complex = nullptr;
    const TC *beta_complex  = nullptr;

    if (k > kB) {
        one_real = one_real_storage.get(beta ? static_cast<const void *>(beta) : static_cast<const void *>(alpha), stream);
        if (!one_real) {
            assert(false && "Failed to create beta=1 scalar for blocked HERK.");
            return std::vector<double>(4, 0.0);
        }
    }

    if (n > nB) {
        alpha_complex = alpha_complex_storage.get(alpha, stream);
        beta_complex  = beta_complex_storage.get(beta, stream);
        if (k > kB) one_complex = one_complex_storage.get(alpha, stream);
        if (!alpha_complex || (beta && !beta_complex) || (k > kB && !one_complex)) {
            assert(false && "Failed to create complex scalars for blocked HERK.");
            one_real_storage.release(stream);
            alpha_complex_storage.release(stream);
            beta_complex_storage.release(stream);
            one_complex_storage.release(stream);
            return std::vector<double>(4, 0.0);
        }
    }

    std::vector<double> timer(4, 0.0);

    core::blocking::for_each_rankk_tile(
        UPLO_C, n, k, nB, kB,
        [&](size_t i, size_t ni, size_t p, size_t kp, bool first) {
            TC *Cii             = C + i + i * ldc;
            const TA *Ai        = core::blocking::matrix_block_ptr(A, lda, trans, i, p);
            const U *beta_block = first ? beta : one_real;
            const auto t        = core::oz2_core_rk<Func::herk, TA, TC, BACKEND, NUM_MODULI, U, U, UPLO_C>(
                handle, trans, ni, kp,
                alpha, Ai, lda, beta_block, Cii, ldc,
                fastmode, work, nullptr,
                false, false, stream);
            core::blocking::add_timer(timer, t);
        },
        [&](size_t i, size_t ni, size_t j, size_t nj, size_t p, size_t kp, bool first) {
            TC *Cij              = C + i + j * ldc;
            const TA *Ai         = core::blocking::matrix_block_ptr(A, lda, trans, i, p);
            const TA *Aj         = core::blocking::matrix_block_ptr(A, lda, trans, j, p);
            const TC *beta_block = first ? beta_complex : one_complex;
            const auto t         = core::oz2_core<Func::gemm, TA, TA, TC, BACKEND, NUM_MODULI>(
                handle, trans, op_B, ni, nj, kp,
                alpha_complex, Ai, lda, Aj, lda, beta_block, Cij, ldc,
                fastmode, work, nullptr, nullptr,
                false, false, false, false, stream);
            core::blocking::add_timer(timer, t);
        });

    one_real_storage.release(stream);
    one_complex_storage.release(stream);
    alpha_complex_storage.release(stream);
    beta_complex_storage.release(stream);
    return timer;
}

} // namespace impl

template <typename TA, typename TC, Backend BACKEND, unsigned NUM_MODULI>
std::vector<double> herk_core(
    common::Handle_t handle,
    cublasFillMode_t uplo, cublasOperation_t trans,
    size_t n, size_t k,
    const common::underlying_t<TC> *alpha,
    const TA *const A, size_t lda,
    const common::underlying_t<TC> *beta,
    TC *const C, size_t ldc,
    bool fastmode,
    void *const work, void *const workA,
    bool enable_skip_scalA,
    bool skip_scalA,
    cudaStream_t stream //
) {
    static_assert(common::isComplex<TA> && common::isComplex<TC>,
                  "TA and TC must be all complex");

    if (uplo == CUBLAS_FILL_MODE_FULL || trans == CUBLAS_OP_T) {
        assert(false && "unsupported");
        return std::vector<double>(4, 0.0);
    }

    if (work == nullptr) {
        size_t workSize_A           = 0;
        const size_t workSize_total = workSize<true, BACKEND>(
            n, n, k, NUM_MODULI, enable_skip_scalA, false, &workSize_A, nullptr, fastmode);
        std::vector<double> timer(4, 0.0);
        timer[0] = static_cast<double>(workSize_total);
        timer[1] = static_cast<double>(workSize_A);
        return timer;
    }

    if (uplo == CUBLAS_FILL_MODE_UPPER) {
        return impl::run<CUBLAS_FILL_MODE_UPPER, TA, TC, BACKEND, NUM_MODULI>(
            handle, trans, n, k, alpha, A, lda, beta, C, ldc,
            fastmode, work, workA, enable_skip_scalA, skip_scalA, stream);
    }
    return impl::run<CUBLAS_FILL_MODE_LOWER, TA, TC, BACKEND, NUM_MODULI>(
        handle, trans, n, k, alpha, A, lda, beta, C, ldc,
        fastmode, work, workA, enable_skip_scalA, skip_scalA, stream);
}

} // namespace gemmul8::oz2::herk
