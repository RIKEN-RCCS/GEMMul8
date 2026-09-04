#pragma once
#include "../common/common.hpp"
#include "../core/block_size.hpp"
#include "../core/blocking.hpp"
#include "../core/oz2_core.hpp"
#include "../gemm/worksize.hpp"
#include "worksize.hpp"

namespace gemmul8::oz2::hemm {

namespace impl {

template <cublasFillMode_t UPLO_S,
          typename THerm, typename TFull, typename TC,
          Backend BACKEND, unsigned NUM_MODULI>
std::vector<double> run(
    common::Handle_t handle,
    cublasSideMode_t side,
    size_t rowsC, size_t colsC,
    const TC *alpha,
    const THerm *Herm, size_t ldHerm,
    const TFull *Full, size_t ldFull,
    const TC *beta,
    TC *C, size_t ldc,
    bool fastmode,
    void *work, void *workHerm, void *workFull,
    bool enable_skip_scalHerm, bool enable_skip_scalFull,
    bool skip_scalHerm, bool skip_scalFull,
    cudaStream_t stream //
) {
    constexpr cublasFillMode_t UF  = CUBLAS_FILL_MODE_FULL;
    constexpr cublasDiagType_t DN  = CUBLAS_DIAG_NON_UNIT;
    constexpr common::MatStruct SF = common::MatStruct::Full;
    constexpr common::MatStruct SS = common::MatStruct::Hermitian;

    const bool memory_saving_mode = handle.config.memory_saving;
    const size_t limit            = handle.config.max_worksize;
    const bool memory_saving      = memory_saving_mode && limit > 0;
    const bool use_skipHerm       = memory_saving_mode ? false : enable_skip_scalHerm;
    const bool use_skipFull       = memory_saving_mode ? false : enable_skip_scalFull;
    const bool do_skipHerm        = memory_saving_mode ? false : skip_scalHerm;
    const bool do_skipFull        = memory_saving_mode ? false : skip_scalFull;

    const size_t s = (side == CUBLAS_SIDE_LEFT) ? rowsC : colsC;
    const size_t f = (side == CUBLAS_SIDE_LEFT) ? colsC : rowsC;
    const size_t m = rowsC, n = colsC, k = s;

    const bool left            = side == CUBLAS_SIDE_LEFT;
    const bool skipA           = left ? use_skipHerm : use_skipFull;
    const bool skipB           = left ? use_skipFull : use_skipHerm;
    const size_t full_worksize = workSize<common::isComplex<THerm>, BACKEND>(
        m, n, k, NUM_MODULI, skipA, skipB, nullptr, nullptr, fastmode);

    if (!memory_saving || full_worksize <= limit) {
        if (left) {
            return core::oz2_core<Func::hemm, THerm, TFull, TC,
                                  BACKEND, NUM_MODULI, TC, TC,
                                  SS, SF, UPLO_S, UF, DN, DN>(
                handle, CUBLAS_OP_N, CUBLAS_OP_N, m, n, k,
                alpha, Herm, ldHerm, Full, ldFull, beta, C, ldc,
                fastmode, work, workHerm, workFull,
                use_skipHerm, use_skipFull, do_skipHerm, do_skipFull, stream);
        }
        return core::oz2_core<Func::hemm, TFull, THerm, TC,
                              BACKEND, NUM_MODULI, TC, TC,
                              SF, SS, UF, UPLO_S, DN, DN>(
            handle, CUBLAS_OP_N, CUBLAS_OP_N, m, n, k,
            alpha, Full, ldFull, Herm, ldHerm, beta, C, ldc,
            fastmode, work, workFull, workHerm,
            use_skipFull, use_skipHerm, do_skipFull, do_skipHerm, stream);
    }

    auto fits = [&](size_t sB, size_t fB, bool need_gemm) {
        const size_t sm = left ? sB : fB;
        const size_t sn = left ? fB : sB;

        if (workSize<common::isComplex<THerm>, BACKEND>(
                sm, sn, sB, NUM_MODULI, false, false, nullptr, nullptr, fastmode) > limit) {
            return false;
        }

        return !need_gemm ||
               gemm::workSize<common::isComplex<THerm>, BACKEND>(
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
            assert(false && "Failed to create beta=1 scalar for blocked HEMM.");
            return std::vector<double>(4, 0.0);
        }
    }

    std::vector<double> timer(4, 0.0);
    core::blocking::for_each_structured_tile<false>(
        side, UPLO_S, s, f, block.xB, block.yB,
        [&](size_t q, size_t sq, size_t r, size_t fr) {
            const THerm *Sqq = Herm + q + q * ldHerm;
            if (left) {
                const auto t = core::oz2_core<Func::hemm, THerm, TFull, TC,
                                              BACKEND, NUM_MODULI, TC, TC,
                                              SS, SF, UPLO_S, UF, DN, DN>(
                    handle, CUBLAS_OP_N, CUBLAS_OP_N, sq, fr, sq,
                    alpha, Sqq, ldHerm, Full + q + r * ldFull, ldFull,
                    beta, C + q + r * ldc, ldc,
                    fastmode, work, nullptr, nullptr,
                    false, false, false, false, stream);
                core::blocking::add_timer(timer, t);
            } else {
                const auto t = core::oz2_core<Func::hemm, TFull, THerm, TC,
                                              BACKEND, NUM_MODULI, TC, TC,
                                              SF, SS, UF, UPLO_S, DN, DN>(
                    handle, CUBLAS_OP_N, CUBLAS_OP_N, fr, sq, sq,
                    alpha, Full + r + q * ldFull, ldFull, Sqq, ldHerm,
                    beta, C + r + q * ldc, ldc,
                    fastmode, work, nullptr, nullptr,
                    false, false, false, false, stream);
                core::blocking::add_timer(timer, t);
            }
        },
        [&](size_t q, size_t sq, size_t p, size_t sp, size_t r, size_t fr) {
            if (left) {
                const auto Sqp = core::blocking::stored_structured_block<true>(Herm, ldHerm, UPLO_S, q, p);
                const auto t   = core::oz2_core<Func::gemm, THerm, TFull, TC, BACKEND, NUM_MODULI>(
                    handle, Sqp.op, CUBLAS_OP_N, sq, fr, sp,
                    alpha, Sqp.ptr, ldHerm, Full + p + r * ldFull, ldFull,
                    one, C + q + r * ldc, ldc,
                    fastmode, work, nullptr, nullptr,
                    false, false, false, false, stream);
                core::blocking::add_timer(timer, t);
            } else {
                const auto Spq = core::blocking::stored_structured_block<true>(Herm, ldHerm, UPLO_S, p, q);
                const auto t   = core::oz2_core<Func::gemm, TFull, THerm, TC, BACKEND, NUM_MODULI>(
                    handle, CUBLAS_OP_N, Spq.op, fr, sq, sp,
                    alpha, Full + r + p * ldFull, ldFull, Spq.ptr, ldHerm,
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

template <typename THerm, typename TFull, typename TC, Backend BACKEND, unsigned NUM_MODULI>
std::vector<double> hemm_core(
    common::Handle_t handle,
    cublasSideMode_t side, cublasFillMode_t uplo,
    size_t rowsC, size_t colsC,
    const TC *alpha,
    const THerm *const Herm, size_t ldHerm,
    const TFull *const Full, size_t ldFull,
    const TC *beta,
    TC *const C, size_t ldc,
    bool fastmode,
    void *const work, void *const workHerm, void *const workFull,
    bool enable_skip_scalHerm, bool enable_skip_scalFull,
    bool skip_scalHerm, bool skip_scalFull,
    cudaStream_t stream //
) {
    static_assert(common::isComplex<THerm> && common::isComplex<TFull> && common::isComplex<TC>,
                  "THerm, TFull, and TC must be all complex");

    if (uplo == CUBLAS_FILL_MODE_FULL) {
        assert(false && "unsupported");
        return std::vector<double>(4, 0.0);
    }

    if (work == nullptr) {
        const size_t m = rowsC;
        const size_t n = colsC;
        const size_t k = (side == CUBLAS_SIDE_LEFT) ? rowsC : colsC;
        size_t wa = 0, wb = 0;
        const size_t total = workSize<common::isComplex<THerm>, BACKEND>(
            m, n, k, NUM_MODULI,
            side == CUBLAS_SIDE_LEFT ? enable_skip_scalHerm : enable_skip_scalFull,
            side == CUBLAS_SIDE_LEFT ? enable_skip_scalFull : enable_skip_scalHerm,
            &wa, &wb, fastmode);
        std::vector<double> timer(4, 0.0);
        timer[0] = static_cast<double>(total);
        timer[1] = static_cast<double>(side == CUBLAS_SIDE_LEFT ? wa : wb);
        timer[2] = static_cast<double>(side == CUBLAS_SIDE_LEFT ? wb : wa);
        return timer;
    }

    if (uplo == CUBLAS_FILL_MODE_UPPER) {
        return impl::run<CUBLAS_FILL_MODE_UPPER, THerm, TFull, TC, BACKEND, NUM_MODULI>(
            handle, side, rowsC, colsC, alpha, Herm, ldHerm, Full, ldFull, beta, C, ldc,
            fastmode, work, workHerm, workFull,
            enable_skip_scalHerm, enable_skip_scalFull, skip_scalHerm, skip_scalFull, stream);
    }
    return impl::run<CUBLAS_FILL_MODE_LOWER, THerm, TFull, TC, BACKEND, NUM_MODULI>(
        handle, side, rowsC, colsC, alpha, Herm, ldHerm, Full, ldFull, beta, C, ldc,
        fastmode, work, workHerm, workFull,
        enable_skip_scalHerm, enable_skip_scalFull, skip_scalHerm, skip_scalFull, stream);
}

} // namespace gemmul8::oz2::hemm
