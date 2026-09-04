#pragma once
#include "../common/common.hpp"
#include "../core/oz2_core.hpp"
#include "../core/block_size.hpp"
#include "../core/blocking.hpp"
#include "worksize.hpp"

namespace gemmul8::oz2::gemm {

template <typename TA, typename TB, typename TC, Backend BACKEND, unsigned NUM_MODULI>
std::vector<double> gemm_core(
    common::Handle_t handle,
    cublasOperation_t op_A, cublasOperation_t op_B,
    size_t m, size_t n, size_t k,
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

    // Return workspace size
    if (work == nullptr) {
        size_t workSize_A, workSize_B;
        size_t workSize_total = workSize<common::isComplex<TA>, BACKEND>(
            m, n, k, NUM_MODULI, enable_skip_scalA, enable_skip_scalB, &workSize_A, &workSize_B, fastmode);
        std::vector<double> timer(4, 0.0);
        timer[0] = static_cast<double>(workSize_total);
        timer[1] = static_cast<double>(workSize_A);
        timer[2] = static_cast<double>(workSize_B);
        return timer;
    }

    const bool memory_saving_mode = handle.config.memory_saving;
    const size_t limit            = handle.config.max_worksize;
    const bool memory_saving      = memory_saving_mode && limit > 0;
    const bool use_skipA          = memory_saving_mode ? false : enable_skip_scalA;
    const bool use_skipB          = memory_saving_mode ? false : enable_skip_scalB;
    const bool do_skipA           = memory_saving_mode ? false : skip_scalA;
    const bool do_skipB           = memory_saving_mode ? false : skip_scalB;
    void *const use_workA         = workA;
    void *const use_workB         = workB;

    const size_t full_worksize = workSize<common::isComplex<TA>, BACKEND>(
        m, n, k, NUM_MODULI, use_skipA, use_skipB, nullptr, nullptr, fastmode);

    if (!memory_saving || full_worksize <= limit) {
        return core::oz2_core<Func::gemm, TA, TB, TC, BACKEND, NUM_MODULI>(
            handle, op_A, op_B, m, n, k,
            alpha, A, lda, B, ldb, beta, C, ldc,
            fastmode, work, use_workA, use_workB,
            use_skipA, use_skipB,
            do_skipA, do_skipB, stream);
    }

    const core::BlockSize3D block = core::find_block_size_gemm<BACKEND, NUM_MODULI>(
        m, n, k, limit,
        [&](size_t mB, size_t nB, size_t kB) {
            return workSize<common::isComplex<TA>, BACKEND>(
                mB, nB, kB, NUM_MODULI, false, false, nullptr, nullptr, fastmode);
        });

    if (!block) {
        assert(false && "The workspace-size limit is too small.");
        return std::vector<double>(4, 0.0);
    }

    std::vector<double> timer(4, 0.0);

    core::blocking::OneScalar<TC> one_storage;
    const TC *one = nullptr;
    if (k > block.kB) {
        const TC *reference = beta ? beta : alpha;
        one                 = one_storage.get(reference, stream);
        if (!one) {
            assert(false && "Failed to create beta=1 scalar for K-blocked GEMM.");
            return std::vector<double>(4, 0.0);
        }
    }

    for (size_t i = 0; i < m; i += block.mB) {
        const size_t mi = std::min<size_t>(block.mB, m - i);

        for (size_t j = 0; j < n; j += block.nB) {
            const size_t nj = std::min<size_t>(block.nB, n - j);
            TC *const Cij   = C + i + j * ldc;

            bool p_first = true;
            for (size_t p = 0; p < k; p += block.kB) {
                const size_t kp = std::min<size_t>(block.kB, k - p);

                const TA *const Aip        = core::blocking::matrix_block_ptr(A, lda, op_A, i, p);
                const TB *const Bpj        = core::blocking::matrix_block_ptr(B, ldb, op_B, p, j);
                const TC *const beta_block = (p_first) ? beta : one;

                const auto t = core::oz2_core<Func::gemm, TA, TB, TC, BACKEND, NUM_MODULI>(
                    handle, op_A, op_B, mi, nj, kp,
                    alpha, Aip, lda, Bpj, ldb,
                    beta_block, Cij, ldc,
                    fastmode, work, nullptr, nullptr,
                    false, false, false, false,
                    stream);

                core::blocking::add_timer(timer, t);
                p_first = false;
            }
        }
    }

    one_storage.release(stream);
    return timer;
}

} // namespace gemmul8::oz2::gemm
