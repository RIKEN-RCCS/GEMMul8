#pragma once
#include "mod_core.hpp"

namespace gemmul8::mod {

// A+ = Ar+s*Ai, A- = Ar-s*Ai, where s*s == -1 (mod p)
template <Backend BACKEND, unsigned IDX> struct complex_2m_traits {
    static constexpr int32_t p = common::table::moduli<BACKEND, IDX, true>;
    static_assert(IDX < 20U && p > 2 && (p & 1) != 0);
    static constexpr int32_t root = [] {
        for (int32_t s = 1; s <= p / 2; ++s) {
            if ((s * s + 1) % p == 0) return s;
        }
        return int32_t(0);
    }();
    static_assert(root != 0, "Complex modulus must admit sqrt(-1)");
    static constexpr int32_t half_root_mod = (root * ((p + 1) / 2)) % p;
    static constexpr int32_t half_root     = half_root_mod > p / 2 ? half_root_mod - p : half_root_mod;
};

template <unsigned IDX> struct complex_2m_dot_traits {
    using R                               = residue_traits<Backend::INT8, IDX, true>;
    static constexpr uint32_t p           = R::p;
    static constexpr uint32_t s           = uint32_t(complex_2m_traits<Backend::INT8, IDX>::root);
    static constexpr uint32_t w00         = (R::r00 * s) % p;
    static constexpr uint32_t w08         = (R::r08 * s) % p;
    static constexpr uint32_t w16         = (R::r16 * s) % p;
    static constexpr uint32_t w24         = (R::r24 * s) % p;
    static constexpr uint32_t w32         = (R::r32 * s) % p;
    static constexpr uint32_t w40         = (R::r40 * s) % p;
    static constexpr uint32_t w48         = (R::r48 * s) % p;
    static constexpr uint32_t w56         = (R::r56 * s) % p;
    static constexpr uint32_t w64         = (R::r64 * s) % p;
    static constexpr uint32_t weight_lo   = w00 | (w08 << 8) | (w16 << 16) | (w24 << 24);
    static constexpr uint32_t weight_hi   = w32 | (w40 << 8) | (w48 << 16) | (w56 << 24);
    static constexpr uint32_t neg_corr_32 = (p - w32) % p;
    static constexpr uint32_t neg_corr_64 = (p - w64) % p;
    static_assert(2ULL * (8ULL * 255ULL * (p - 1U) + p) < uint64_t(INT32_MAX));
};

template <Backend BACKEND, unsigned IDX, typename V>
__device__ __forceinline__ int2 project_complex_2m(const V v) {
    if constexpr (BACKEND == Backend::INT8 && has_u8_dot4 &&
                  (std::is_same_v<V, int2> || std::is_same_v<V, common::mant2_t>)) {
        using D = complex_2m_dot_traits<IDX>;
        int32_t ar, si;
        if constexpr (std::is_same_v<V, int2>) {
            ar = int32_t(reduce_i32_dot4<IDX, true>(v.x));
            si = int32_t(dot4_u8(uint32_t(v.y), D::weight_lo, v.y < 0 ? D::neg_corr_32 : 0U));
        } else {
            ar                = int32_t(reduce_i64_dot4<IDX, true>(v.x));
            const uint32_t hi = dot4_u8(uint32_t(v.y.hi), D::weight_hi, v.y.hi < 0 ? D::neg_corr_64 : 0U);
            si                = int32_t(dot4_u8(v.y.lo, D::weight_lo, hi));
        }
        return int2{mod_small<BACKEND, IDX, true>(ar + si),
                    mod_small<BACKEND, IDX, true>(ar - si)};
    } else if constexpr (BACKEND == Backend::FP8) {
        using R                  = complex_2m_traits<BACKEND, IDX>;
        const int32_t ar         = calc_mod_fp8_nowrap<IDX, true>(v.x);
        const int32_t ai         = calc_mod_fp8_nowrap<IDX, true>(v.y);
        const int32_t si         = R::root * ai;
        constexpr uint64_t bound = 2ULL * R::p * (1ULL + R::root);
        return int2{mod_bounded<BACKEND, IDX, true, bound>(ar + si),
                    mod_bounded<BACKEND, IDX, true, bound>(ar - si)};
    }
    using R          = complex_2m_traits<BACKEND, IDX>;
    const int32_t ar = calc_mod<BACKEND, IDX, true>(v.x);
    const int32_t ai = calc_mod<BACKEND, IDX, true>(v.y);
    const int32_t si = mod_small<BACKEND, IDX, true>(R::root * ai);
    return int2{wrapping<BACKEND, IDX, true>(ar + si),
                wrapping<BACKEND, IDX, true>(ar - si)};
}

template <Backend BACKEND, unsigned IDX, bool FLIP_IMAG = false>
__device__ __forceinline__ int2 reconstruct_complex_2m(const int32_t plus, const int32_t minus) {
    using R                          = complex_2m_traits<BACKEND, IDX>;
    const int32_t sum                = plus + minus;
    const int32_t re                 = (sum + (sum & 1) * (sum < 0 ? R::p : -R::p)) >> 1;
    const int32_t diff               = FLIP_IMAG ? plus - minus : minus - plus;
    constexpr uint64_t abs_half_root = R::half_root < 0 ? -int64_t(R::half_root) : R::half_root;
    constexpr uint64_t bound         = uint64_t(R::p - 1) * abs_half_root;
    const int32_t im                 = mod_bounded<BACKEND, IDX, true, bound>(R::half_root * diff);
    return int2{re, im};
}

} // namespace gemmul8::mod
