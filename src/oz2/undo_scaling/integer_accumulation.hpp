#pragma once
#include "accumulation.hpp"

namespace gemmul8::undo_scaling {

constexpr int crt_lattice_exponent(double q) {
    const uint64_t bits = std::bit_cast<uint64_t>(q);
    const uint64_t sig  = (bits & ((uint64_t(1) << 52) - 1)) | (uint64_t(1) << 52);
    return int((bits >> 52) & 0x7ff) - 1023 - 52 + std::countr_zero(sig);
}

template <Backend B, unsigned N, bool C>
struct crt_integer_traits {
    static constexpr int exponent =
        []<size_t... I>(std::index_sequence<I...>) {
            return std::min({crt_lattice_exponent(common::table::qPi_double2<B, N, I, C>().x)...});
        }(std::make_index_sequence<N>{});

    static constexpr double scale = std::bit_cast<double>(uint64_t(exponent + 1023) << 52);

    template <unsigned I>
    static constexpr int64_t coefficient = int64_t(common::table::qPi_double2<B, N, I, C>().x / scale);
};

#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ == 800 || __CUDA_ARCH__ == 900 || __CUDA_ARCH__ == 1000)
inline constexpr bool use_integer_crt_high = false; // FP64 is not slow
#elif defined(__HIP_DEVICE_COMPILE__) && (defined(__gfx908__) || defined(__gfx90a__) || defined(__gfx942__) || defined(__gfx950__))
inline constexpr bool use_integer_crt_high = false; // FP64 is not slow
#else
inline constexpr bool use_integer_crt_high = true; // FP64 is slow
#endif

template <Backend B, unsigned N, unsigned I = 0>
__device__ __forceinline__ void accumulate_real_integer_high(
    int64_t &hi,
    double &lo,
    const common::mid_t<B, false> *__restrict__ mid,
    const size_t inc //
) {
    using Q              = crt_integer_traits<B, N, false>;
    const int64_t c      = mid[I * inc];
    constexpr double qlo = common::table::qPi_double2<B, N, I>().y;
    if constexpr (I == 0) {
        hi = Q::template coefficient<I> * c;
        lo = qlo * double(c);
    } else {
        hi += Q::template coefficient<I> * c;
        lo = fma(qlo, double(c), lo);
    }
    if constexpr (I + 1U < N) accumulate_real_integer_high<B, N, I + 1U>(hi, lo, mid, inc);
}

template <Backend B, unsigned N>
__device__ __forceinline__ double2 accumulate_real_integer_high(
    const common::mid_t<B, false> *__restrict__ mid,
    const size_t inc //
) {
    using Q = crt_integer_traits<B, N, false>;
    int64_t hi;
    double lo;
    accumulate_real_integer_high<B, N>(hi, lo, mid, inc);
    return {double(hi) * Q::scale, lo};
}

template <Backend B, unsigned N, unsigned I = 0>
__device__ __forceinline__ void accumulate_complex_integer_high(
    int64_t &re,
    int64_t &im,
    double2 &lo,
    const common::mid_t<B, true> *__restrict__ mid,
    const size_t inc //
) {
    using Q              = crt_integer_traits<B, N, true>;
    const auto c         = mid[I * inc];
    constexpr double qlo = common::table::qPi_double2<B, N, I, true>().y;
    if constexpr (I == 0) {
        re = Q::template coefficient<I> * int64_t(c.x);
        im = Q::template coefficient<I> * int64_t(c.y);
        lo = {qlo * double(c.x), qlo * double(c.y)};
    } else {
        re += Q::template coefficient<I> * int64_t(c.x);
        im += Q::template coefficient<I> * int64_t(c.y);
        lo.x = fma(qlo, double(c.x), lo.x);
        lo.y = fma(qlo, double(c.y), lo.y);
    }
    if constexpr (I + 1U < N) accumulate_complex_integer_high<B, N, I + 1U>(re, im, lo, mid, inc);
}

template <Backend B, unsigned N>
__device__ __forceinline__ common::double2x2_t accumulate_complex_integer_high(
    const common::mid_t<B, true> *__restrict__ mid,
    const size_t inc //
) {
    using Q = crt_integer_traits<B, N, true>;
    int64_t re, im;
    double2 lo;
    accumulate_complex_integer_high<B, N>(re, im, lo, mid, inc);
    return {
        {double(re) * Q::scale, double(im) * Q::scale},
        lo
    };
}

} // namespace gemmul8::undo_scaling
