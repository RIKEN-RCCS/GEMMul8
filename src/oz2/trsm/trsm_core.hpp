#pragma once
#include "../common/common.hpp"
#include "worksize.hpp"
#include "trsm_left.hpp"
#include "trsm_right.hpp"

namespace gemmul8::oz2::trsm {

template <typename TA, typename TB, Backend BACKEND, unsigned NUM_MODULI>
inline std::vector<double> trsm_core(
    common::Handle_t handle,
    cublasSideMode_t side, cublasFillMode_t uplo,
    cublasOperation_t trans, cublasDiagType_t diag,
    size_t m, size_t n,
    const TB *alpha,
    const TA *const A, size_t lda,
    TB *const B, size_t ldb,
    bool fastmode,
    void *const work,
    cudaStream_t stream //
) {
    static_assert(std::is_same_v<TA, TB>, "trsm requires std::is_same_v<TA, TB>.");

    if (work != nullptr && (m == 0 || n == 0)) return std::vector<double>(4, 0.0);

    if (uplo != CUBLAS_FILL_MODE_LOWER && uplo != CUBLAS_FILL_MODE_UPPER) {
        assert(false && "TRSM requires uplo = LOWER or UPPER.");
        return std::vector<double>(4, 0.0);
    }

    if (side == CUBLAS_SIDE_LEFT) {
        return trsm_left<TA, TB, BACKEND, NUM_MODULI>(
            handle, uplo, trans, diag, m, n, alpha, A, lda, B, ldb, fastmode, work, stream);
    }

    if (side == CUBLAS_SIDE_RIGHT) {
        return trsm_right<TA, TB, BACKEND, NUM_MODULI>(
            handle, uplo, trans, diag, m, n, alpha, A, lda, B, ldb, fastmode, work, stream);
    }

    assert(false && "TRSM requires side = LEFT or RIGHT.");
    return std::vector<double>(4, 0.0);
}

} // namespace gemmul8::oz2::trsm
