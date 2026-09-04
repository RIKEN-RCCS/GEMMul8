#pragma once
#include "../common/common.hpp"
#include "../common/self_hipify.hpp"
#include "../core/blocking.hpp"
#include "../gemm/gemm_core.hpp"
#include "blas.hpp"
#include "block_size.hpp"
#include "worksize.hpp"

namespace gemmul8::oz2::trsm {

template <typename TA, typename TB, Backend BACKEND, unsigned NUM_MODULI>
inline std::vector<double> trsm_right(
    common::Handle_t handle,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    cublasDiagType_t diag,
    size_t m, size_t n,
    const TB *alpha,
    const TA *const A, size_t lda,
    TB *const B, size_t ldb,
    bool fastmode,
    void *const work,
    cudaStream_t stream //
) {
    static_assert(std::is_same_v<TA, TB>, "trsm requires std::is_same_v<TA, TB>.");

    // Return workspace size.
    if (work == nullptr) {
        size_t workSize_total = workSize_right<TA, BACKEND>(m, n, NUM_MODULI, handle.config.block_size_trsm);
        std::vector<double> timer(4, 0.0);
        timer[0] = static_cast<double>(workSize_total);
        return timer;
    }

    cublasHandle_t trsm_handle{};
    bool created_trsm_handle = false;

    if (handle.kind == common::HandleKind::cuBLAS) {
        trsm_handle = handle.cublas;
    } else {
        cublasCreate(&trsm_handle);
        created_trsm_handle = true;

        cublasSetStream(trsm_handle, stream);

        const size_t lwork_blas = size_t(32) << 20;
        const size_t lwork_max  = handle.config.max_worksize;
        const size_t lwork_trsm = (handle.config.memory_saving && lwork_max > 0)
                                      ? std::min<size_t>(lwork_blas, lwork_max)
                                      : lwork_blas;
        cublasSetWorkspace(trsm_handle, work, lwork_trsm);

        const cublasPointerMode_t pointer_mode = undo_scaling::is_device_pointer(alpha)
                                                     ? CUBLAS_POINTER_MODE_DEVICE
                                                     : CUBLAS_POINTER_MODE_HOST;
        cublasSetPointerMode(trsm_handle, pointer_mode);
    }

    handle.arch     = 0;
    const size_t nB = size_t(block_size_trsm<TA, BACKEND>(n, handle.arch, handle.config.block_size_trsm));

    core::blocking::OneScalar<TB> one_storage;
    core::blocking::MinusOneScalar<TB> minus_one_storage;
    const TB *one       = nullptr;
    const TB *minus_one = nullptr;

    if (n > nB) {
        one       = one_storage.get(alpha, stream);
        minus_one = minus_one_storage.get(alpha, stream);

        if (!one || !minus_one) {
            assert(false && "Failed to create scalar constants for blocked TRSM.");
            one_storage.release(stream);
            minus_one_storage.release(stream);
            if (created_trsm_handle) cublasDestroy(trsm_handle);
            return std::vector<double>(4, 0.0);
        }
    }

    const unsigned num_events = unsigned((n + nB - 1) / nB) * 2u + 1u;
    std::vector<cudaEvent_t> events(num_events, nullptr);
    for (auto &ev : events) {
        cudaEventCreate(&ev);
    }

    const bool transN = (trans == CUBLAS_OP_N);
    const cublasFillMode_t eff_uplo =
        transN ? uplo : ((uplo == CUBLAS_FILL_MODE_LOWER) ? CUBLAS_FILL_MODE_UPPER : CUBLAS_FILL_MODE_LOWER);
    const bool forward = (eff_uplo == CUBLAS_FILL_MODE_UPPER);

    auto run_update =
        [&](const size_t col0,
            const size_t cols,
            const size_t j,
            const size_t jb,
            const TB *beta_gemm //
        ) {
            const TA *Ablk         = nullptr;
            cublasOperation_t op_A = CUBLAS_OP_N;

            if (transN) {
                Ablk = A + j + col0 * lda;
                op_A = CUBLAS_OP_N;
            } else {
                Ablk = A + col0 + j * lda;
                op_A = trans;
            }

            TB *const Bj   = B + j * ldb;
            TB *const Cblk = B + col0 * ldb;

            gemm::gemm_core<TB, TA, TB, BACKEND, NUM_MODULI>(
                handle, CUBLAS_OP_N, op_A, m, cols, jb,
                minus_one, Bj, ldb, Ablk, lda, beta_gemm, Cblk, ldb,
                fastmode, work, nullptr, nullptr, false, false, false, false, stream);
        };

    unsigned event_idx = 0;
    cudaEventRecord(events[event_idx], stream);

    if (forward) {

        for (size_t j = 0; j < n; j += nB) {
            const size_t jb = std::min<size_t>(nB, n - j);

            const TA *const Ajj = A + j + j * lda;
            TB *const Bj        = B + j * ldb;

            const bool first_block     = (j == 0);
            const TB *const alpha_trsm = first_block ? alpha : one;
            const TB *const beta_gemm  = first_block ? alpha : one;

            small_trsm<TB>(
                trsm_handle, CUBLAS_SIDE_RIGHT,
                uplo, trans, diag,
                static_cast<int>(m), static_cast<int>(jb),
                alpha_trsm,
                Ajj, static_cast<int>(lda),
                Bj, static_cast<int>(ldb));

            ++event_idx;
            cudaEventRecord(events[event_idx], stream);

            const size_t right = j + jb;
            if (right < n) {
                run_update(right, n - right, j, jb, beta_gemm);

                ++event_idx;
                cudaEventRecord(events[event_idx], stream);
            }
        }

    } else {

        for (size_t j_end = n; j_end > 0;) {
            const size_t jb = std::min<size_t>(nB, j_end);
            const size_t j  = j_end - jb;

            const TA *const Ajj = A + j + j * lda;
            TB *const Bj        = B + j * ldb;

            const bool first_block     = (j_end == n);
            const TB *const alpha_trsm = first_block ? alpha : one;
            const TB *const beta_gemm  = first_block ? alpha : one;

            small_trsm<TB>(
                trsm_handle,
                CUBLAS_SIDE_RIGHT, uplo, trans, diag,
                static_cast<int>(m), static_cast<int>(jb),
                alpha_trsm,
                Ajj, static_cast<int>(lda),
                Bj, static_cast<int>(ldb));

            ++event_idx;
            cudaEventRecord(events[event_idx], stream);

            if (j > 0) {
                run_update(0, j, j, jb, beta_gemm);

                ++event_idx;
                cudaEventRecord(events[event_idx], stream);
            }

            j_end = j;
        }
    }

    cudaEventSynchronize(events[event_idx]);

    std::vector<double> timer(4, 0.0);
    float ms      = 0.0f;
    int timer_idx = 0;

    for (unsigned i = 1; i <= event_idx; ++i) {
        cudaEventElapsedTime(&ms, events[i - 1], events[i]);
        timer[timer_idx] += double(ms) * 1.0e-3;
        timer_idx = 1 - timer_idx;
    }

    for (auto ev : events) {
        if (ev) cudaEventDestroy(ev);
    }
    one_storage.release(stream);
    minus_one_storage.release(stream);

    if (created_trsm_handle) {
        cublasDestroy(trsm_handle);
    }

    return timer;
}

} // namespace gemmul8::oz2::trsm
