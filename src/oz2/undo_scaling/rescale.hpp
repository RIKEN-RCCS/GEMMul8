#pragma once
#include "../common/common.hpp"

namespace gemmul8::undo_scaling {

template <typename T> __device__ __forceinline__ T rescale_crt(T x, const int32_t sft);

template <>
__device__ __forceinline__ double rescale_crt<double>(double x, const int32_t sft) {
    const uint32_t hi = uint32_t(__double2hiint(x));
    const int32_t e   = int32_t((hi >> 20) & 0x7ffU);
    if (e == 0) return x;
    const int32_t new_e      = e + sft;
    const uint32_t sign_frac = hi & 0x800fffffU;
    if (uint32_t(new_e - 1) < 2046U) {
        return __hiloint2double(int(sign_frac | (uint32_t(new_e) << 20)), __double2loint(x));
    }
    if (new_e >= 2047) return __hiloint2double(int((hi & 0x80000000U) | 0x7ff00000U), 0);
    if (new_e < -52) return __hiloint2double(int(hi & 0x80000000U), 0);
    const double normal = __hiloint2double(int(sign_frac | (uint32_t(new_e + 54) << 20)), __double2loint(x));
    return normal * 0x1p-54;
}

template <>
__device__ __forceinline__ cuDoubleComplex rescale_crt<cuDoubleComplex>(cuDoubleComplex x, const int32_t sft) {
    return {rescale_crt(x.x, sft), rescale_crt(x.y, sft)};
}

} // namespace gemmul8::undo_scaling
