#pragma once
#include "../common/common.hpp"
#include "../core/worksize.hpp"

namespace gemmul8::oz2::syrk {

template <bool is_Complex, Backend BACKEND>
inline size_t workSize(
    size_t m, size_t n, size_t k, unsigned NUM_MODULI,
    bool enable_skip_scalA, bool enable_skip_scalB,
    size_t *workSizeA, size_t *workSizeB, bool fastmode //
) {
    if (workSizeB != nullptr) *workSizeB = 0;
    return core::workSize_rk<is_Complex, BACKEND>(m, k, NUM_MODULI, workSizeA);
}

} // namespace gemmul8::oz2::syrk
