#pragma once
#include "../common/common.hpp"
#include "../core/block_size.hpp"
#include "../core/blocking.hpp"
#include "../core/oz2_core.hpp"
#include "../gemm/worksize.hpp"
#include "worksize.hpp"

namespace gemmul8::oz2::trmm {

namespace impl {

template <cublasFillMode_t UPLO_T, cublasDiagType_t DIAG_T,
          typename TTri, typename TFull, typename TC,
          Backend BACKEND, unsigned NUM_MODULI>
std::vector<double> run(
    common::Handle_t handle,
    cublasSideMode_t side, cublasOperation_t trans,
    size_t rowsC, size_t colsC,
    const TC *alpha,
    const TTri *Tri, size_t ldTri,
    const TFull *Full, size_t ldFull,
    TC *C, size_t ldc,
    bool fastmode,
    void *work, void *workTri, void *workFull,
    bool enable_skip_scalTri, bool enable_skip_scalFull,
    bool skip_scalTri, bool skip_scalFull,
    cudaStream_t stream //
) {
    constexpr cublasFillMode_t UF  = CUBLAS_FILL_MODE_FULL;
    constexpr cublasDiagType_t DN  = CUBLAS_DIAG_NON_UNIT;
    constexpr common::MatStruct SF = common::MatStruct::Full;
    constexpr common::MatStruct ST = common::MatStruct::Triangular;

    const bool left               = side == CUBLAS_SIDE_LEFT;
    const bool memory_saving_mode = handle.config.memory_saving;
    const size_t limit            = handle.config.max_worksize;
    const bool memory_saving      = memory_saving_mode && limit > 0;
    const bool use_skipTri        = memory_saving_mode ? false : enable_skip_scalTri;
    const bool use_skipFull       = memory_saving_mode ? false : enable_skip_scalFull;
    const bool do_skipTri         = memory_saving_mode ? false : skip_scalTri;
    const bool do_skipFull        = memory_saving_mode ? false : skip_scalFull;

    const size_t s = left ? rowsC : colsC;
    const size_t f = left ? colsC : rowsC;
    const size_t m = rowsC, n = colsC, k = s;

    const size_t full_worksize = workSize<common::isComplex<TTri>, BACKEND>(
        m, n, k, NUM_MODULI,
        left ? use_skipTri : use_skipFull,
        left ? use_skipFull : use_skipTri,
        nullptr, nullptr, fastmode);

    if (!memory_saving || full_worksize <= limit) {
        if (left) {
            return core::oz2_core<Func::trmm, TTri, TFull, TC,
                                  BACKEND, NUM_MODULI, TC, TC,
                                  ST, SF, UPLO_T, UF, DIAG_T, DN>(
                handle, trans, CUBLAS_OP_N, m, n, k,
                alpha, Tri, ldTri, Full, ldFull, nullptr, C, ldc,
                fastmode, work, workTri, workFull,
                use_skipTri, use_skipFull, do_skipTri, do_skipFull, stream);
        }
        return core::oz2_core<Func::trmm, TFull, TTri, TC,
                              BACKEND, NUM_MODULI, TC, TC,
                              SF, ST, UF, UPLO_T, DN, DIAG_T>(
            handle, CUBLAS_OP_N, trans, m, n, k,
            alpha, Full, ldFull, Tri, ldTri, nullptr, C, ldc,
            fastmode, work, workFull, workTri,
            use_skipFull, use_skipTri, do_skipFull, do_skipTri, stream);
    }

    auto fits = [&](size_t sB, size_t fB, bool need_gemm) {
        const size_t sm = left ? sB : fB;
        const size_t sn = left ? fB : sB;

        if (workSize<common::isComplex<TTri>, BACKEND>(
                sm, sn, sB, NUM_MODULI, false, false, nullptr, nullptr, fastmode) > limit) {
            return false;
        }

        return !need_gemm ||
               gemm::workSize<common::isComplex<TTri>, BACKEND>(
                   sm, sn, sB, NUM_MODULI, false, false, nullptr, nullptr, fastmode) <= limit;
    };

    const core::BlockSize2D block = core::find_block_size_structured(s, f, true, fits);
    if (!block) {
        assert(false && "The workspace-size limit is too small.");
        return std::vector<double>(4, 0.0);
    }

    core::blocking::OneScalar<TC> one_storage;
    const TC *one = nullptr;
    if (s > block.xB) {
        one = one_storage.get(alpha, stream);
        if (!one) {
            assert(false && "Failed to create beta=1 scalar for blocked TRMM.");
            return std::vector<double>(4, 0.0);
        }
    }

    const cublasFillMode_t effective = core::blocking::effective_uplo(UPLO_T, trans);
    std::vector<double> timer(4, 0.0);

    core::blocking::for_each_structured_tile<true>(
        side, effective, s, f, block.xB, block.yB,
        [&](size_t q, size_t sq, size_t r, size_t fr) {
            const TTri *Tqq = Tri + q + q * ldTri;
            if (left) {
                const auto t = core::oz2_core<Func::trmm, TTri, TFull, TC,
                                              BACKEND, NUM_MODULI, TC, TC,
                                              ST, SF, UPLO_T, UF, DIAG_T, DN>(
                    handle, trans, CUBLAS_OP_N, sq, fr, sq,
                    alpha, Tqq, ldTri, Full + q + r * ldFull, ldFull,
                    nullptr, C + q + r * ldc, ldc,
                    fastmode, work, nullptr, nullptr,
                    false, false, false, false, stream);
                core::blocking::add_timer(timer, t);
            } else {
                const auto t = core::oz2_core<Func::trmm, TFull, TTri, TC,
                                              BACKEND, NUM_MODULI, TC, TC,
                                              SF, ST, UF, UPLO_T, DN, DIAG_T>(
                    handle, CUBLAS_OP_N, trans, fr, sq, sq,
                    alpha, Full + r + q * ldFull, ldFull, Tqq, ldTri,
                    nullptr, C + r + q * ldc, ldc,
                    fastmode, work, nullptr, nullptr,
                    false, false, false, false, stream);
                core::blocking::add_timer(timer, t);
            }
        },
        [&](size_t q, size_t sq, size_t p, size_t sp, size_t r, size_t fr) {
            if (left) {
                const TTri *Tqp = core::blocking::matrix_block_ptr(Tri, ldTri, trans, q, p);
                const auto t    = core::oz2_core<Func::gemm, TTri, TFull, TC, BACKEND, NUM_MODULI>(
                    handle, trans, CUBLAS_OP_N, sq, fr, sp,
                    alpha, Tqp, ldTri, Full + p + r * ldFull, ldFull,
                    one, C + q + r * ldc, ldc,
                    fastmode, work, nullptr, nullptr,
                    false, false, false, false, stream);
                core::blocking::add_timer(timer, t);
            } else {
                const TTri *Tpq = core::blocking::matrix_block_ptr(Tri, ldTri, trans, p, q);
                const auto t    = core::oz2_core<Func::gemm, TFull, TTri, TC, BACKEND, NUM_MODULI>(
                    handle, CUBLAS_OP_N, trans, fr, sq, sp,
                    alpha, Full + r + p * ldFull, ldFull, Tpq, ldTri,
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

template <typename TTri, typename TFull, typename TC, Backend BACKEND, unsigned NUM_MODULI>
std::vector<double> trmm_core(
    common::Handle_t handle,
    cublasSideMode_t side, cublasFillMode_t uplo,
    cublasOperation_t trans, cublasDiagType_t diag,
    size_t rowsC, size_t colsC,
    const TC *alpha,
    const TTri *const Tri, size_t ldTri,
    const TFull *const Full, size_t ldFull,
    TC *const C, size_t ldc,
    bool fastmode,
    void *const work, void *const workTri, void *const workFull,
    bool enable_skip_scalTri, bool enable_skip_scalFull,
    bool skip_scalTri, bool skip_scalFull,
    cudaStream_t stream //
) {
    static_assert(common::isComplex<TTri> == common::isComplex<TFull> &&
                      common::isComplex<TFull> == common::isComplex<TC>,
                  "TTri, TFull, and TC must be all real or all complex");

    if (uplo == CUBLAS_FILL_MODE_FULL) {
        assert(false && "unsupported");
        return std::vector<double>(4, 0.0);
    }

    if (work == nullptr) {
        const size_t m = rowsC, n = colsC, k = (side == CUBLAS_SIDE_LEFT) ? rowsC : colsC;
        size_t wa = 0, wb = 0;
        const size_t total = workSize<common::isComplex<TTri>, BACKEND>(
            m, n, k, NUM_MODULI,
            side == CUBLAS_SIDE_LEFT ? enable_skip_scalTri : enable_skip_scalFull,
            side == CUBLAS_SIDE_LEFT ? enable_skip_scalFull : enable_skip_scalTri,
            &wa, &wb, fastmode);
        std::vector<double> timer(4, 0.0);
        timer[0] = static_cast<double>(total);
        timer[1] = static_cast<double>(side == CUBLAS_SIDE_LEFT ? wa : wb);
        timer[2] = static_cast<double>(side == CUBLAS_SIDE_LEFT ? wb : wa);
        return timer;
    }

    if (uplo == CUBLAS_FILL_MODE_UPPER) {
        if (diag == CUBLAS_DIAG_NON_UNIT) {
            return impl::run<CUBLAS_FILL_MODE_UPPER, CUBLAS_DIAG_NON_UNIT, TTri, TFull, TC, BACKEND, NUM_MODULI>(
                handle, side, trans, rowsC, colsC, alpha, Tri, ldTri, Full, ldFull, C, ldc,
                fastmode, work, workTri, workFull, enable_skip_scalTri, enable_skip_scalFull,
                skip_scalTri, skip_scalFull, stream);
        }
        return impl::run<CUBLAS_FILL_MODE_UPPER, CUBLAS_DIAG_UNIT, TTri, TFull, TC, BACKEND, NUM_MODULI>(
            handle, side, trans, rowsC, colsC, alpha, Tri, ldTri, Full, ldFull, C, ldc,
            fastmode, work, workTri, workFull, enable_skip_scalTri, enable_skip_scalFull,
            skip_scalTri, skip_scalFull, stream);
    }
    if (diag == CUBLAS_DIAG_NON_UNIT) {
        return impl::run<CUBLAS_FILL_MODE_LOWER, CUBLAS_DIAG_NON_UNIT, TTri, TFull, TC, BACKEND, NUM_MODULI>(
            handle, side, trans, rowsC, colsC, alpha, Tri, ldTri, Full, ldFull, C, ldc,
            fastmode, work, workTri, workFull, enable_skip_scalTri, enable_skip_scalFull,
            skip_scalTri, skip_scalFull, stream);
    }
    return impl::run<CUBLAS_FILL_MODE_LOWER, CUBLAS_DIAG_UNIT, TTri, TFull, TC, BACKEND, NUM_MODULI>(
        handle, side, trans, rowsC, colsC, alpha, Tri, ldTri, Full, ldFull, C, ldc,
        fastmode, work, workTri, workFull, enable_skip_scalTri, enable_skip_scalFull,
        skip_scalTri, skip_scalFull, stream);
}

} // namespace gemmul8::oz2::trmm
