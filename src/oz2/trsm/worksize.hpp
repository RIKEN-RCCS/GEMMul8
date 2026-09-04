#pragma once
#include "../common/common.hpp"
#include "../core/worksize.hpp"
#include "block_size.hpp"

namespace gemmul8::oz2::trsm {

template <bool is_Complex, Backend BACKEND>
inline size_t workSize_update_gemm(size_t rows, size_t cols, size_t kk, unsigned NUM_MODULI) {
    constexpr size_t lwork_blas = size_t(32) << 20; // 32 MiB

    size_t lwork_gemm = 0;
    if (rows > 0 && cols > 0 && kk > 0) {
        constexpr common::MatMulKind KIND = common::MatMulKind::Gemm;

        lwork_gemm = core::workSize<is_Complex, BACKEND, KIND>(
            rows, cols, kk, NUM_MODULI, false, false, nullptr, nullptr, false);
    }

    return std::max<size_t>(lwork_gemm, lwork_blas);
}

template <typename T, Backend BACKEND>
inline size_t workSize_left(size_t m, size_t n, unsigned NUM_MODULI, int nB_override = 0) {
    int arch        = 0;
    const size_t nB = size_t(block_size_trsm<T, BACKEND>(m, arch, nB_override));
    const size_t jb = std::min<size_t>(nB, m);

    const size_t rows = (jb < m) ? (m - jb) : 0;
    return workSize_update_gemm<common::isComplex<T>, BACKEND>(rows, n, jb, NUM_MODULI);
}

template <typename T, Backend BACKEND>
inline size_t workSize_right(size_t m, size_t n, unsigned NUM_MODULI, int nB_override = 0) {
    int arch        = 0;
    const size_t nB = size_t(block_size_trsm<T, BACKEND>(n, arch, nB_override));
    const size_t jb = std::min<size_t>(nB, n);

    const size_t cols = (jb < n) ? (n - jb) : 0;
    return workSize_update_gemm<common::isComplex<T>, BACKEND>(m, cols, jb, NUM_MODULI);
}

} // namespace gemmul8::oz2::trsm
