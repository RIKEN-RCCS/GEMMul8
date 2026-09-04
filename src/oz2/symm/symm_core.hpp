#pragma once
#include "../common/common.hpp"
#include "../core/block_size.hpp"
#include "../core/blocking.hpp"
#include "../core/oz2_core.hpp"
#include "../gemm/worksize.hpp"
#include "worksize.hpp"

namespace gemmul8::oz2::symm {

namespace impl {

template <cublasFillMode_t UPLO_S,
          typename TSym, typename TFull, typename TC,
          Backend BACKEND, unsigned NUM_MODULI>
std::vector<double> run(
    common::Handle_t handle,
    cublasSideMode_t side,
    size_t rowsC, size_t colsC,
    const TC *alpha,
    const TSym *Sym, size_t ldSym,
    const TFull *Full, size_t ldFull,
    const TC *beta,
    TC *C, size_t ldc,
    bool fastmode,
    void *work, void *workSym, void *workFull,
    bool enable_skip_scalSym, bool enable_skip_scalFull,
    bool skip_scalSym, bool skip_scalFull,
    cudaStream_t stream //
) {
    constexpr cublasFillMode_t UF  = CUBLAS_FILL_MODE_FULL;
    constexpr cublasDiagType_t DN  = CUBLAS_DIAG_NON_UNIT;
    constexpr common::MatStruct SF = common::MatStruct::Full;
    constexpr common::MatStruct SS = common::MatStruct::Symmetric;

    const bool memory_saving_mode = handle.config.memory_saving;
    const size_t limit            = handle.config.max_worksize;
    const bool memory_saving      = memory_saving_mode && limit > 0;
    const bool use_skipSym        = memory_saving_mode ? false : enable_skip_scalSym;
    const bool use_skipFull       = memory_saving_mode ? false : enable_skip_scalFull;
    const bool do_skipSym         = memory_saving_mode ? false : skip_scalSym;
    const bool do_skipFull        = memory_saving_mode ? false : skip_scalFull;

    const bool left  = side == CUBLAS_SIDE_LEFT;
    const size_t s   = left ? rowsC : colsC;
    const size_t f   = left ? colsC : rowsC;
    const size_t m   = rowsC;
    const size_t n   = colsC;
    const size_t k   = s;
    const bool skipA = left ? use_skipSym : use_skipFull;
    const bool skipB = left ? use_skipFull : use_skipSym;

    const size_t full_worksize = workSize<common::isComplex<TSym>, BACKEND>(
        m, n, k, NUM_MODULI, skipA, skipB, nullptr, nullptr, fastmode);

    if (!memory_saving || full_worksize <= limit) {
        if (left) {
            return core::oz2_core<Func::symm, TSym, TFull, TC,
                                  BACKEND, NUM_MODULI, TC, TC,
                                  SS, SF, UPLO_S, UF, DN, DN>(
                handle, CUBLAS_OP_N, CUBLAS_OP_N, m, n, k,
                alpha, Sym, ldSym, Full, ldFull, beta, C, ldc,
                fastmode, work, workSym, workFull,
                use_skipSym, use_skipFull, do_skipSym, do_skipFull, stream);
        }
        return core::oz2_core<Func::symm, TFull, TSym, TC,
                              BACKEND, NUM_MODULI, TC, TC,
                              SF, SS, UF, UPLO_S, DN, DN>(
            handle, CUBLAS_OP_N, CUBLAS_OP_N, m, n, k,
            alpha, Full, ldFull, Sym, ldSym, beta, C, ldc,
            fastmode, work, workFull, workSym,
            use_skipFull, use_skipSym, do_skipFull, do_skipSym, stream);
    }

    auto fits = [&](size_t sB, size_t fB, bool need_gemm) {
        const size_t sm = left ? sB : fB;
        const size_t sn = left ? fB : sB;

        if (workSize<common::isComplex<TSym>, BACKEND>(
                sm, sn, sB, NUM_MODULI, false, false, nullptr, nullptr, fastmode) > limit) {
            return false;
        }

        return !need_gemm ||
               gemm::workSize<common::isComplex<TSym>, BACKEND>(
                   sm, sn, sB, NUM_MODULI, false, false, nullptr, nullptr, fastmode) <= limit;
    };

    const core::BlockSize2D block = core::find_block_size_structured(s, f, false, fits);
    if (!block) {
        assert(false && "The workspace-size limit is too small.");
        return std::vector<double>(4, 0.0);
    }

    core::blocking::OneScalar<TC> one_storage;
    const TC *one = nullptr;
    if (s > block.xB) {
        one = one_storage.get(beta ? static_cast<const void *>(beta) : static_cast<const void *>(alpha), stream);
        if (!one) {
            assert(false && "Failed to create beta=1 scalar for blocked SYMM.");
            return std::vector<double>(4, 0.0);
        }
    }

    std::vector<double> timer(4, 0.0);
    core::blocking::for_each_structured_tile<false>(
        side, UPLO_S, s, f, block.xB, block.yB,
        [&](size_t q, size_t sq, size_t r, size_t fr) {
            const TSym *Sqq = Sym + q + q * ldSym;
            if (left) {
                const auto t = core::oz2_core<Func::symm, TSym, TFull, TC,
                                              BACKEND, NUM_MODULI, TC, TC,
                                              SS, SF, UPLO_S, UF, DN, DN>(
                    handle, CUBLAS_OP_N, CUBLAS_OP_N, sq, fr, sq,
                    alpha, Sqq, ldSym, Full + q + r * ldFull, ldFull,
                    beta, C + q + r * ldc, ldc,
                    fastmode, work, nullptr, nullptr,
                    false, false, false, false, stream);
                core::blocking::add_timer(timer, t);
            } else {
                const auto t = core::oz2_core<Func::symm, TFull, TSym, TC,
                                              BACKEND, NUM_MODULI, TC, TC,
                                              SF, SS, UF, UPLO_S, DN, DN>(
                    handle, CUBLAS_OP_N, CUBLAS_OP_N, fr, sq, sq,
                    alpha, Full + r + q * ldFull, ldFull, Sqq, ldSym,
                    beta, C + r + q * ldc, ldc,
                    fastmode, work, nullptr, nullptr,
                    false, false, false, false, stream);
                core::blocking::add_timer(timer, t);
            }
        },
        [&](size_t q, size_t sq, size_t p, size_t sp, size_t r, size_t fr) {
            if (left) {
                const auto Sqp = core::blocking::stored_structured_block<false>(Sym, ldSym, UPLO_S, q, p);
                const auto t   = core::oz2_core<Func::gemm, TSym, TFull, TC, BACKEND, NUM_MODULI>(
                    handle, Sqp.op, CUBLAS_OP_N, sq, fr, sp,
                    alpha, Sqp.ptr, ldSym, Full + p + r * ldFull, ldFull,
                    one, C + q + r * ldc, ldc,
                    fastmode, work, nullptr, nullptr,
                    false, false, false, false, stream);
                core::blocking::add_timer(timer, t);
            } else {
                const auto Spq = core::blocking::stored_structured_block<false>(Sym, ldSym, UPLO_S, p, q);
                const auto t   = core::oz2_core<Func::gemm, TFull, TSym, TC, BACKEND, NUM_MODULI>(
                    handle, CUBLAS_OP_N, Spq.op, fr, sq, sp,
                    alpha, Full + r + p * ldFull, ldFull, Spq.ptr, ldSym,
                    one, C + r + q * ldc, ldc,
                    fastmode, work, nullptr, nullptr,
                    false, false, false, false, stream);
                core::blocking::add_timer(timer, t);
            }
        });

    one_storage.release(stream);
    return timer;
}

} // namespace impl

template <typename TSym, typename TFull, typename TC, Backend BACKEND, unsigned NUM_MODULI>
std::vector<double> symm_core(
    common::Handle_t handle,
    cublasSideMode_t side, cublasFillMode_t uplo,
    size_t rowsC, size_t colsC,
    const TC *alpha,
    const TSym *const Sym, size_t ldSym,
    const TFull *const Full, size_t ldFull,
    const TC *beta,
    TC *const C, size_t ldc,
    bool fastmode,
    void *const work, void *const workSym, void *const workFull,
    bool enable_skip_scalSym, bool enable_skip_scalFull,
    bool skip_scalSym, bool skip_scalFull,
    cudaStream_t stream //
) {
    static_assert(common::isComplex<TSym> == common::isComplex<TFull> &&
                      common::isComplex<TFull> == common::isComplex<TC>,
                  "TA, TB, and TC must be all real or all complex");

    if (uplo == CUBLAS_FILL_MODE_FULL) {
        assert(false && "unsupported");
        return std::vector<double>(4, 0.0);
    }

    if (work == nullptr) {
        const size_t m = rowsC;
        const size_t n = colsC;
        const size_t k = (side == CUBLAS_SIDE_LEFT) ? rowsC : colsC;
        size_t wa = 0, wb = 0;
        const size_t total = workSize<common::isComplex<TSym>, BACKEND>(
            m, n, k, NUM_MODULI,
            side == CUBLAS_SIDE_LEFT ? enable_skip_scalSym : enable_skip_scalFull,
            side == CUBLAS_SIDE_LEFT ? enable_skip_scalFull : enable_skip_scalSym,
            &wa, &wb, fastmode);
        std::vector<double> timer(4, 0.0);
        timer[0] = static_cast<double>(total);
        timer[1] = static_cast<double>(side == CUBLAS_SIDE_LEFT ? wa : wb);
        timer[2] = static_cast<double>(side == CUBLAS_SIDE_LEFT ? wb : wa);
        return timer;
    }

    if (uplo == CUBLAS_FILL_MODE_UPPER) {
        return impl::run<CUBLAS_FILL_MODE_UPPER, TSym, TFull, TC, BACKEND, NUM_MODULI>(
            handle, side, rowsC, colsC, alpha, Sym, ldSym, Full, ldFull, beta, C, ldc,
            fastmode, work, workSym, workFull,
            enable_skip_scalSym, enable_skip_scalFull, skip_scalSym, skip_scalFull, stream);
    }
    return impl::run<CUBLAS_FILL_MODE_LOWER, TSym, TFull, TC, BACKEND, NUM_MODULI>(
        handle, side, rowsC, colsC, alpha, Sym, ldSym, Full, ldFull, beta, C, ldc,
        fastmode, work, workSym, workFull,
        enable_skip_scalSym, enable_skip_scalFull, skip_scalSym, skip_scalFull, stream);
}

} // namespace gemmul8::oz2::symm
