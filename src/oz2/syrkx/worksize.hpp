#pragma once
#include "../common/common.hpp"
#include "../core/worksize.hpp"

namespace gemmul8::oz2::syrkx {

template <bool is_Complex, Backend BACKEND>
inline size_t workSize(
    size_t m, size_t n, size_t k, unsigned NUM_MODULI,
    bool enable_skip_scalA, bool enable_skip_scalB,
    size_t *workSizeA, size_t *workSizeB, bool fastmode //
) {
    constexpr common::MatMulKind KIND = common::MatMulKind::ATxB;
    return core::workSize<is_Complex, BACKEND, KIND>(
        m, n, k, NUM_MODULI, enable_skip_scalA, enable_skip_scalB, workSizeA, workSizeB, fastmode);
}

} // namespace gemmul8::oz2::syrkx
