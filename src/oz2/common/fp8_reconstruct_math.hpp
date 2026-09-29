#pragma once
#include "fp8_plan.hpp"
#include <cstdint>

namespace gemmul8::common::fp8_plan {

template <int32_t Q>
__device__ __forceinline__ int32_t reduce_int_nowrap(int32_t x) {
    static_assert(Q >= 3);
    constexpr int32_t inv = int32_t((uint64_t(1) << 32) / Q);
    const int32_t k       = __mulhi(x, inv);
    return int32_t(uint32_t(x) - uint32_t(k) * uint32_t(Q));
}

__device__ __forceinline__ constexpr int64_t coeff_abs(int32_t c) {
    return c < 0 ? -int64_t(c) : int64_t(c);
}

__device__ __forceinline__ constexpr int64_t reconstruction_bound(info s, unsigned mask) {
    int64_t bound = 0;
    for (unsigned j = 0; j < s.products; ++j)
        bound += coeff_abs(s.coefficient[j]) * ((mask & (1U << j)) ? (1LL << 24) : 2LL * s.q[j]);
    return bound;
}

__device__ __forceinline__ constexpr unsigned integer_raw_mask(info s) {
    unsigned best = 0, count = 0;
    int64_t best_bound = INT64_MAX;
    for (unsigned mask = 0; mask < (1U << s.products); ++mask) {
        unsigned n = 0;
        for (unsigned j = 0; j < s.products; ++j) n += bool(mask & (1U << j));
        const int64_t bound = reconstruction_bound(s, mask);
        if (bound < int64_t(INT32_MAX) - s.p && (n > count || (n == count && bound < best_bound))) {
            best       = mask;
            count      = n;
            best_bound = bound;
        }
    }
    return best;
}

template <int32_t P, unsigned J>
__device__ __forceinline__ int32_t reconstruction_input(float f) {
    constexpr auto s = scheme<P>;
    const int32_t x  = static_cast<int32_t>(f);
    if constexpr (integer_raw_mask(s) & (1U << J)) return x;
    else return reduce_int_nowrap<s.q[J]>(x);
}

__device__ __forceinline__ constexpr int32_t small_inverse(int32_t x, int32_t q) {
    for (int32_t y = 1; y < q; ++y) {
        if ((x * y) % q == 1) return y > q / 2 ? y - q : y;
    }
    return 0;
}

struct signed_range { int64_t lo, hi; };

__device__ __forceinline__ constexpr int64_t range_abs(signed_range r) {
    return -r.lo > r.hi ? -r.lo : r.hi;
}

__device__ __forceinline__ constexpr signed_range remainder_range(int64_t bound, int32_t q) {
    const int64_t scale = 1LL << 32;
    const int64_t error = scale % q;
    return {-(bound * error / scale), q - 1 + (bound * error + scale - 1) / scale};
}

template <int32_t P> struct garner_plan {
    static constexpr auto s            = scheme<P>;
    static constexpr int32_t p01       = s.q[0] * s.q[1];
    static constexpr int32_t inv1      = small_inverse(s.q[0], s.q[1]);
    static constexpr int32_t inv2      = small_inverse(p01, s.q[2]);
    static constexpr signed_range r0   = remainder_range(1LL << 24, s.q[0]);
    static constexpr int64_t b1        = ((1LL << 24) + range_abs(r0)) * coeff_abs(inv1);
    static constexpr signed_range t1   = remainder_range(b1, s.q[1]);
    static constexpr signed_range r01  = {r0.lo + s.q[0] * t1.lo, r0.hi + s.q[0] * t1.hi};
    static constexpr int64_t b2        = ((1LL << 24) + range_abs(r01)) * coeff_abs(inv2);
    static constexpr signed_range t2   = remainder_range(b2, s.q[2]);
    static constexpr signed_range r012 = {r01.lo + p01 * t2.lo, r01.hi + p01 * t2.hi};
    static constexpr auto result       = s.products == 2 ? r01 : r012;
};

template <int32_t P, bool Center = true>
__device__ __forceinline__ int32_t reconstruct_garner(float c0, float c1, float c2) {
    using G = garner_plan<P>;
    static_assert(Center || (G::result.lo >= -2LL * P && G::result.hi <= 2LL * P));
    constexpr auto s = G::s;
    int32_t x        = reduce_int_nowrap<s.q[0]>(int32_t(c0));
    const int32_t t1 = reduce_int_nowrap<s.q[1]>((int32_t(c1) - x) * G::inv1);
    x += s.q[0] * t1;
    if constexpr (s.products == 3) {
        const int32_t t2 = reduce_int_nowrap<s.q[2]>((int32_t(c2) - x) * G::inv2);
        x += G::p01 * t2;
    }
    return (Center && x > P / 2) ? x - P : x;
}

template <int32_t P, bool Center = true>
__device__ __forceinline__ int32_t reconstruct(float c0, float c1, float c2) {
    constexpr auto s = scheme<P>;
    if constexpr (s.w == 0 && s.products >= 2 && !crt_scaled(P)) {
        using G = garner_plan<P>;

        constexpr bool first_power_of_two = (s.q[0] & (s.q[0] - 1)) == 0;
        if constexpr ((s.products == 2 && !first_power_of_two) ||
                      (s.products == 3 && coeff_abs(G::inv2) == 1)) {
            return reconstruct_garner<P, Center>(c0, c1, c2);
        }
    }
    int32_t x = s.coefficient[0] * reconstruction_input<P, 0>(c0);
    if constexpr (s.products > 1) {
        x += s.coefficient[1] * reconstruction_input<P, 1>(c1);
    }
    if constexpr (s.products > 2) {
        x += s.coefficient[2] * reconstruction_input<P, 2>(c2);
    }
    const int32_t r = reduce_int_nowrap<P>(x);
    return (Center && r > P / 2) ? r - P : r;
}

} // namespace gemmul8::common::fp8_plan
