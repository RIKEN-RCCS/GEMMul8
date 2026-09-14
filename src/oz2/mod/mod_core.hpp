#pragma once
#include "../common/common.hpp"
#include "../common/table.hpp"

namespace gemmul8::mod {

template <Backend BACKEND, unsigned IDX, bool COMPLEX>
struct residue_traits {
    static constexpr uint32_t p = uint32_t(common::table::moduli<BACKEND, IDX, COMPLEX>);
    static_assert((BACKEND == Backend::INT8 && p <= 256U) ||
                  (BACKEND == Backend::FP8 && p <= 2401U));

    // 2^(8*i) mod p, i = 0,...,8
    static constexpr uint32_t r00 = 1U % p;
    static constexpr uint32_t r08 = (r00 << 8) % p;
    static constexpr uint32_t r16 = (r08 << 8) % p;
    static constexpr uint32_t r24 = (r16 << 8) % p;
    static constexpr uint32_t r32 = (r24 << 8) % p;
    static constexpr uint32_t r40 = (r32 << 8) % p;
    static constexpr uint32_t r48 = (r40 << 8) % p;
    static constexpr uint32_t r56 = (r48 << 8) % p;
    static constexpr uint32_t r64 = (r56 << 8) % p;

    static constexpr uint32_t neg_corr_32 = (r32 == 0U) ? 0U : (p - r32);
    static constexpr uint32_t neg_corr_64 = (r64 == 0U) ? 0U : (p - r64);

    static constexpr uint32_t p_inv = uint32_t(common::table::p_inv_32_v<BACKEND, IDX, COMPLEX>);

    // floor((2^32 mod p) * 2^32 / p)
    static constexpr uint32_t r32_inv = uint32_t((uint64_t(r32) << 32) / uint64_t(p));
};

#if defined(__CUDACC__) && !defined(__HIPCC__)

inline constexpr bool has_u8_dot4 = true;
__device__ __forceinline__ uint32_t dot4_u8(uint32_t a, uint32_t b, uint32_t c) { return __dp4a(a, b, c); }

#elif defined(__HIP_DEVICE_COMPILE__) &&                                     \
    (defined(__gfx906__) || defined(__gfx908__) || defined(__gfx90a__) ||    \
     defined(__gfx940__) || defined(__gfx941__) || defined(__gfx942__) ||    \
     defined(__gfx1011__) || defined(__gfx1012__) ||                         \
     defined(__gfx1030__) || defined(__gfx1031__) || defined(__gfx1032__) || \
     defined(__gfx1034__) || defined(__gfx1035__) || defined(__gfx1036__) || \
     defined(__gfx1100__) || defined(__gfx1101__) || defined(__gfx1102__) || \
     defined(__gfx1103__) || defined(__gfx1150__) || defined(__gfx1151__) || \
     defined(__gfx1152__) || defined(__gfx1200__) || defined(__gfx1201__))

inline constexpr bool has_u8_dot4 = true;
__device__ __forceinline__ uint32_t dot4_u8(uint32_t a, uint32_t b, uint32_t c) { return __builtin_amdgcn_udot4(a, b, c, false); }

#else

inline constexpr bool has_u8_dot4 = false;

__device__ __forceinline__ uint32_t dot4_u8(uint32_t a, uint32_t b, uint32_t c) {
    c += (a & 255U) * (b & 255U);
    c += ((a >> 8) & 255U) * ((b >> 8) & 255U);
    c += ((a >> 16) & 255U) * ((b >> 16) & 255U);
    c += ((a >> 24) & 255U) * ((b >> 24) & 255U);
    return c;
}

#endif

#if defined(__HIP_DEVICE_COMPILE__) &&                                       \
    (defined(__gfx906__) || defined(__gfx908__) || defined(__gfx90a__) ||    \
     defined(__gfx940__) || defined(__gfx941__) || defined(__gfx942__) ||    \
     defined(__gfx1011__) || defined(__gfx1012__) ||                         \
     defined(__gfx1030__) || defined(__gfx1031__) || defined(__gfx1032__) || \
     defined(__gfx1034__) || defined(__gfx1035__) || defined(__gfx1036__))

inline constexpr bool has_u16_dot2 = true;
using u16x2_t                      = unsigned short __attribute__((ext_vector_type(2)));
__device__ __forceinline__ uint32_t dot2_u16(uint32_t a, uint32_t b, uint32_t c) {
    const u16x2_t va = __builtin_bit_cast(u16x2_t, a);
    const u16x2_t vb = __builtin_bit_cast(u16x2_t, b);
    return __builtin_amdgcn_udot2(va, vb, c, false);
}

#else

inline constexpr bool has_u16_dot2 = false;
__device__ __forceinline__ uint32_t dot2_u16(uint32_t a, uint32_t b, uint32_t c) {
    c += (a & 65535U) * (b & 65535U);
    c += (a >> 16) * (b >> 16);
    return c;
}

#endif

template <unsigned IDX, bool COMPLEX>
struct int8_dot_traits : residue_traits<Backend::INT8, IDX, COMPLEX> {
    using R = residue_traits<Backend::INT8, IDX, COMPLEX>;
    static_assert(R::p <= 256U);

    static constexpr uint32_t weight_lo = R::r00 | (R::r08 << 8) | (R::r16 << 16) | (R::r24 << 24);
    static constexpr uint32_t weight_hi = R::r32 | (R::r40 << 8) | (R::r48 << 16) | (R::r56 << 24);

    static constexpr uint32_t dot4_constexpr(uint32_t a, uint32_t b, uint32_t c = 0U) {
        c += (a & 255U) * (b & 255U);
        c += ((a >> 8) & 255U) * ((b >> 8) & 255U);
        c += ((a >> 16) & 255U) * ((b >> 16) & 255U);
        c += ((a >> 24) & 255U) * ((b >> 24) & 255U);
        return c;
    }

    static constexpr uint64_t mant32_max = [] {
        const uint64_t positive = 255ULL * (R::r00 + R::r08 + R::r16) + 127ULL * R::r24;
        const uint64_t negative = 255ULL * (R::r00 + R::r08 + R::r16 + R::r24) + R::neg_corr_32;
        return positive > negative ? positive : negative;
    }();

    static constexpr uint64_t mant64_max = [] {
        const uint64_t positive = 255ULL * (R::r00 + R::r08 + R::r16 + R::r24 + R::r32 + R::r40 + R::r48) + 127ULL * R::r56;
        const uint64_t negative = 255ULL * (R::r00 + R::r08 + R::r16 + R::r24 + R::r32 + R::r40 + R::r48 + R::r56) + R::neg_corr_64;
        return positive > negative ? positive : negative;
    }();

    static constexpr uint32_t exp_first_max = [] {
        uint32_t vmax = 0U;
        for (unsigned bit = 0; bit < 64U; ++bit) {
            const bool hi      = bit >= 32U;
            const uint32_t val = uint32_t(1U) << (bit & 31U);
            const uint32_t raw = dot4_constexpr(val, hi ? weight_hi : weight_lo, 0U);
            vmax               = raw > vmax ? raw : vmax;
        }
        return vmax;
    }();

    static constexpr uint32_t exp_second_max = [] {
        uint32_t vmax = 0U;
        for (unsigned bit = 0; bit < 64U; ++bit) {
            const bool hi      = bit >= 32U;
            const uint32_t val = uint32_t(1U) << (bit & 31U);
            const uint32_t raw = dot4_constexpr(val, hi ? weight_hi : weight_lo, 0U);
            const uint32_t rem = dot4_constexpr(raw, weight_lo, 0U);
            vmax               = rem > vmax ? rem : vmax;
        }
        return vmax;
    }();

    static constexpr uint64_t fp32_scaled_mant_max = [] {
        const uint64_t positive = 128ULL * R::r00 + 255ULL * R::r08 + 255ULL * R::r16 + 127ULL * R::r24;
        const uint64_t negative = 128ULL * R::r00 + 255ULL * R::r08 + 255ULL * R::r16 + 192ULL * R::r24 + R::neg_corr_32;
        return positive > negative ? positive : negative;
    }();

    static constexpr uint64_t fp64_scaled_mant_max = [] {
        const uint64_t positive = 252ULL * R::r08 + 255ULL * (R::r16 + R::r24 + R::r32 + R::r40 + R::r48) + 127ULL * R::r56;
        const uint64_t negative = 252ULL * R::r08 + 255ULL * (R::r16 + R::r24 + R::r32 + R::r40 + R::r48) + 192ULL * R::r56 + R::neg_corr_64;
        return positive > negative ? positive : negative;
    }();

    static constexpr uint64_t i32_max = uint64_t(INT32_MAX);

    static constexpr bool fp32_one_exp_pass_safe =
        mant32_max <= i32_max && fp32_scaled_mant_max * uint64_t(exp_first_max) <= i32_max;

    static constexpr bool fp64_one_exp_pass_safe =
        mant64_max <= i32_max && fp64_scaled_mant_max * uint64_t(exp_first_max) <= i32_max;

    static_assert(mant32_max * uint64_t(exp_second_max) <= i32_max);
    static_assert(mant64_max * uint64_t(exp_second_max) <= i32_max);
};

template <unsigned IDX, bool COMPLEX>
struct fp8_large_traits : residue_traits<Backend::FP8, IDX, COMPLEX> {
    using R = residue_traits<Backend::FP8, IDX, COMPLEX>;

    static constexpr int32_t r32_centered = R::r32 > R::p / 2 ? int32_t(R::r32) - int32_t(R::p) : int32_t(R::r32);
    static constexpr uint64_t exp_max     = [] {
        uint64_t bound = 0;
        for (unsigned bit = 0; bit < 64U; ++bit) {
            const uint64_t val = uint64_t(1) << (bit & 31U);
            const uint64_t c   = bit >= 32U ? R::r32 : 1U;
            const uint64_t mu  = bit >= 32U ? R::r32_inv : R::p_inv;
            const uint64_t rem = val * c - ((val * mu) >> 32) * R::p;
            bound              = rem > bound ? rem : bound;
        }
        return bound;
    }();

    static constexpr int64_t hi_min           = -int64_t(R::r32 / 2U);
    static constexpr int64_t hi_max           = int64_t(R::p) + (R::r32 + 1U) / 2U - 1;
    static constexpr int64_t lo_max           = int64_t(R::p) + R::r32 - 1;
    static constexpr int64_t raw_min          = (r32_centered >= 0 ? hi_min : hi_max) * r32_centered;
    static constexpr int64_t raw_max          = (r32_centered >= 0 ? hi_max : hi_min) * r32_centered + lo_max;
    static constexpr uint64_t mant_max        = uint64_t(-raw_min > raw_max ? -raw_min : raw_max);
    static constexpr bool skip_mant_reduction = mant_max * exp_max <= uint64_t(INT32_MAX);
};

template <Backend BACKEND, unsigned IDX, bool COMPLEX>
__device__ __forceinline__ int32_t center_reduced(int32_t a) {
    constexpr int32_t p      = common::table::moduli<BACKEND, IDX, COMPLEX>;
    constexpr int32_t p_half = p / 2;
    return (a > p_half) ? (a - p) : a;
}

template <Backend BACKEND, unsigned IDX, bool COMPLEX>
__device__ __forceinline__ uint32_t reduce_exp_scaled(common::exp_t a) {
    using R           = residue_traits<BACKEND, IDX, COMPLEX>;
    const uint32_t c  = a.is_hi ? R::r32 : 1U;
    const uint32_t mu = a.is_hi ? R::r32_inv : R::p_inv;
    const uint32_t q  = __umulhi(a.val, mu);
    return a.val * c - q * R::p;
}

template <unsigned IDX, bool COMPLEX>
__device__ __forceinline__ uint32_t reduce_u32_dot4(uint32_t a) {
    using D = int8_dot_traits<IDX, COMPLEX>;
    return dot4_u8(a, D::weight_lo, 0U);
}

template <unsigned IDX, bool COMPLEX>
__device__ __forceinline__ uint32_t reduce_i32_dot4(int32_t a) {
    using D            = int8_dot_traits<IDX, COMPLEX>;
    const uint32_t acc = (a < 0) ? D::neg_corr_32 : 0U;
    return dot4_u8(uint32_t(a), D::weight_lo, acc);
}

template <unsigned IDX, bool COMPLEX>
__device__ __forceinline__ uint32_t reduce_i64_dot4(common::mant_t a) {
    using D               = int8_dot_traits<IDX, COMPLEX>;
    const uint32_t acc0   = (a.hi < 0) ? D::neg_corr_64 : 0U;
    const uint32_t rem_hi = dot4_u8(uint32_t(a.hi), D::weight_hi, acc0);
    return dot4_u8(a.lo, D::weight_lo, rem_hi);
}

template <unsigned IDX, bool COMPLEX, bool SECOND_PASS>
__device__ __forceinline__ uint32_t reduce_exp_dot4(common::exp_t a) {
    using D      = int8_dot_traits<IDX, COMPLEX>;
    uint32_t raw = dot4_u8(a.val, a.is_hi ? D::weight_hi : D::weight_lo, 0U);
    if constexpr (SECOND_PASS) {
        raw = dot4_u8(raw, D::weight_lo, 0U);
    }
    return raw;
}

template <unsigned IDX, bool COMPLEX>
__device__ __forceinline__ uint32_t reduce_i64_fp8_dot2(common::mant_t a) {
    using R                      = residue_traits<Backend::FP8, IDX, COMPLEX>;
    constexpr uint32_t weight_lo = R::r00 | (R::r16 << 16);
    constexpr uint32_t weight_hi = R::r32 | (R::r48 << 16);
    const uint32_t acc0          = (a.hi < 0) ? R::neg_corr_64 : 0U;
    const uint32_t lo            = dot2_u16(a.lo, weight_lo, acc0);
    return dot2_u16(uint32_t(a.hi), weight_hi, lo);
}

//------------------------------
// Calculate mod: a - round(a/p(j))*p(j)
//------------------------------

// return value in [-p/2, p/2]
template <Backend BACKEND, unsigned IDX, bool COMPLEX = false>
__device__ __forceinline__ int32_t wrapping(int32_t a) {
    constexpr int32_t p      = common::table::moduli<BACKEND, IDX, COMPLEX>;
    constexpr int32_t p_half = p / 2;
    return (a > p_half) ? (a - p) : ((a < -p_half) ? (a + p) : a);
}

template <Backend BACKEND, unsigned IDX, bool COMPLEX = false>
__device__ __forceinline__ uint32_t mod_small_nowrap_u32(uint32_t a) {
    constexpr uint32_t up = uint32_t(common::table::moduli<BACKEND, IDX, COMPLEX>);
    if constexpr (BACKEND == Backend::INT8 && up == 256U) {
        return a & 255U;
    }
    if constexpr (BACKEND == Backend::FP8 && up == 1024U) {
        return a & 1023U;
    }
    if constexpr (BACKEND == Backend::INT8 && has_u8_dot4) {
        return reduce_u32_dot4<IDX, COMPLEX>(a); // < 2^17
    }
    constexpr uint32_t p_inv_u32 = uint32_t(common::table::p_inv_32_v<BACKEND, IDX, COMPLEX>); // 2^32/p
    return a - up * __umulhi(a, p_inv_u32);
}

// |a| < 2^31 is guaranteed (#moduli <= common::threshold::S)
template <Backend BACKEND, unsigned IDX, bool COMPLEX = false>
__device__ __forceinline__ int32_t mod_small_nowrap(int32_t a) {
    constexpr uint32_t up = uint32_t(common::table::moduli<BACKEND, IDX, COMPLEX>);
    if constexpr (BACKEND == Backend::INT8 && up == 256U) {
        return a & 255;
    }
    if constexpr (BACKEND == Backend::FP8 && up == 1024U) {
        return a & 1023;
    }
    if constexpr (BACKEND == Backend::INT8 && has_u8_dot4) {
        return int32_t(reduce_i32_dot4<IDX, COMPLEX>(a)); // < 2^17
    }
    constexpr int32_t p         = common::table::moduli<BACKEND, IDX, COMPLEX>;
    constexpr int32_t p_inv_i32 = common::table::p_inv_32_v<BACKEND, IDX, COMPLEX>; // 2^32/p
    return int32_t(uint32_t(a) - uint32_t(p) * uint32_t(__mulhi(a, p_inv_i32)));
}

template <Backend BACKEND, unsigned IDX, bool COMPLEX = false>
__device__ __forceinline__ int32_t mod_small(int32_t a) {
    constexpr int32_t p = common::table::moduli<BACKEND, IDX, COMPLEX>;
    if constexpr (BACKEND == Backend::INT8 && p == 256) {
        return center_reduced<Backend::INT8, IDX, COMPLEX>(a & 255);
    }
    if constexpr (BACKEND == Backend::FP8 && p == 1024) {
        return center_reduced<Backend::FP8, IDX, COMPLEX>(a & 1023);
    }
    if constexpr (BACKEND == Backend::INT8 && p == 255 && has_u8_dot4) {
        const uint32_t s = dot4_u8(uint32_t(a), 0x01010101U, (a < 0) ? 254U : 0U);
        const uint32_t t = (s & 255U) + (s >> 8);
        return center_reduced<Backend::INT8, IDX, COMPLEX>(int32_t(t));
    }
    constexpr int32_t p_inv_i32 = common::table::p_inv_32_v<BACKEND, IDX, COMPLEX>; // 2^32/p
    const int32_t rem           = int32_t(uint32_t(a) - uint32_t(p) * uint32_t(__mulhi(a, p_inv_i32)));
    return center_reduced<BACKEND, IDX, COMPLEX>(rem);
}

template <Backend BACKEND, unsigned IDX, bool COMPLEX, uint64_t ABS_MAX>
__device__ __forceinline__ int32_t mod_bounded(int32_t a) {
    static_assert(ABS_MAX <= uint64_t(INT32_MAX));

    constexpr int32_t p           = common::table::moduli<BACKEND, IDX, COMPLEX>;
    constexpr uint64_t scale      = uint64_t(1) << 32;
    constexpr uint32_t m          = uint32_t((scale + p - 1) / p);
    constexpr uint64_t e          = uint64_t(m) * p - scale;
    constexpr uint64_t reciprocal = ((uint64_t(1) << 33) + p / 2) / p;
    constexpr int64_t error       = int64_t(reciprocal * p) - (int64_t(1) << 33);
    constexpr uint64_t abs_error  = error < 0 ? uint64_t(-error) : uint64_t(error);
    if constexpr ((p & 1) && COMPLEX && 2ULL * ABS_MAX * e + p < scale) {
        constexpr uint32_t bias = uint32_t((scale * (p / 2) + ABS_MAX * e) / p + 1);
        const uint32_t frac     = uint32_t(a) * m + bias;
        return int32_t(__umulhi(frac, uint32_t(p))) - p / 2;
    } else if constexpr ((p & 1) && reciprocal <= uint64_t(INT32_MAX) && ABS_MAX * abs_error < scale) {
        const int32_t twice_q = __mulhi(a, int32_t(reciprocal)) + 1;
        return int32_t(uint32_t(a) - uint32_t(twice_q >> 1) * uint32_t(p));
    } else {
        return mod_small<BACKEND, IDX, COMPLEX>(a);
    }
}

template <Backend BACKEND, unsigned IDX, bool COMPLEX = false>
__device__ __forceinline__ int32_t reduce_mant(common::mant_t a) {
    constexpr uint32_t up = uint32_t(common::table::moduli<BACKEND, IDX, COMPLEX>);
    if constexpr (BACKEND == Backend::INT8 && up == 256U) {
        return int32_t(a.lo & 255U);
    }
    if constexpr (BACKEND == Backend::FP8 && up == 1024U) {
        return int32_t(a.lo & 1023U);
    }
    if constexpr (BACKEND == Backend::INT8 && has_u8_dot4) {
        return int32_t(reduce_i64_dot4<IDX, COMPLEX>(a)); // < 2^18
    }
    if constexpr (BACKEND == Backend::FP8 && has_u16_dot2) {
        const uint32_t raw = reduce_i64_fp8_dot2<IDX, COMPLEX>(a);
        return int32_t(mod_small_nowrap_u32<BACKEND, IDX, COMPLEX>(raw));
    }
    using R               = residue_traits<BACKEND, IDX, COMPLEX>;
    const int32_t rem_hi  = mod_small_nowrap<BACKEND, IDX, COMPLEX>(a.hi);
    const uint32_t rem_lo = mod_small_nowrap_u32<BACKEND, IDX, COMPLEX>(a.lo);
    const int32_t raw     = rem_hi * int32_t(R::r32) + int32_t(rem_lo);
    if constexpr (BACKEND == Backend::INT8) {
        return raw;
    } else {
        return mod_small_nowrap<BACKEND, IDX, COMPLEX>(raw);
    }
}

template <Backend BACKEND, unsigned IDX, bool COMPLEX = false>
__device__ __forceinline__ int32_t reduce_mant_large(common::mant_t a) {
    constexpr uint32_t up = uint32_t(common::table::moduli<BACKEND, IDX, COMPLEX>);
    if constexpr (BACKEND == Backend::INT8 && up == 256U) {
        return int32_t(a.lo & 255U);
    }
    if constexpr (BACKEND == Backend::FP8 && up == 1024U) {
        return int32_t(a.lo & 1023U);
    }
    if constexpr (BACKEND == Backend::INT8 && has_u8_dot4) {
        return int32_t(reduce_i64_dot4<IDX, COMPLEX>(a));
    }
    if constexpr (BACKEND == Backend::FP8 && has_u16_dot2) {
        const uint32_t raw = reduce_i64_fp8_dot2<IDX, COMPLEX>(a);
        return int32_t(mod_small_nowrap_u32<BACKEND, IDX, COMPLEX>(raw));
    }
    using R               = residue_traits<BACKEND, IDX, COMPLEX>;
    const int32_t rem_hi  = mod_small_nowrap<BACKEND, IDX, COMPLEX>(a.hi);
    const uint32_t rem_lo = mod_small_nowrap_u32<BACKEND, IDX, COMPLEX>(a.lo);
    const int32_t raw     = rem_hi * int32_t(R::r32) + int32_t(rem_lo);
    return mod_small_nowrap<BACKEND, IDX, COMPLEX>(raw);
}

template <Backend BACKEND, unsigned IDX, bool COMPLEX = false, bool SECOND_PASS = true>
__device__ __forceinline__ int32_t reduce_exp(common::exp_t a) {
    constexpr uint32_t up = uint32_t(common::table::moduli<BACKEND, IDX, COMPLEX>);
    if constexpr (BACKEND == Backend::INT8 && up == 256U) {
        return int32_t((a.is_hi ? 0U : a.val) & 255U);
    }
    if constexpr (BACKEND == Backend::FP8 && up == 1024U) {
        return int32_t((a.is_hi ? 0U : a.val) & 1023U);
    }
    if constexpr (BACKEND == Backend::INT8 && has_u8_dot4) {
        return int32_t(reduce_exp_dot4<IDX, COMPLEX, SECOND_PASS>(a));
    }
    return int32_t(reduce_exp_scaled<BACKEND, IDX, COMPLEX>(a));
}

// |a| < 2^63 is guaranteed (common::threshold::S < #moduli <= common::threshold::M)
template <Backend BACKEND, unsigned IDX, bool COMPLEX = false>
__device__ __forceinline__ int32_t mod_middle(common::mant_t a) {
    constexpr int32_t p = common::table::moduli<BACKEND, IDX, COMPLEX>;
    if constexpr (BACKEND == Backend::INT8 && p == 256) {
        return center_reduced<Backend::INT8, IDX, COMPLEX>(int32_t(a.lo & 255U));
    }
    if constexpr (BACKEND == Backend::FP8 && p == 1024) {
        return center_reduced<Backend::FP8, IDX, COMPLEX>(int32_t(a.lo & 1023U));
    }
    const int32_t rem = reduce_mant<BACKEND, IDX, COMPLEX>(a);
    if constexpr (BACKEND == Backend::INT8) {
        if constexpr (has_u8_dot4) {
            using D = int8_dot_traits<IDX, COMPLEX>;
            return mod_bounded<BACKEND, IDX, COMPLEX, D::mant64_max>(rem);
        } else {
            return mod_small<BACKEND, IDX, COMPLEX>(rem);
        }
    } else {
        return center_reduced<BACKEND, IDX, COMPLEX>(rem);
    }
}

// |a| can be >= 2^63 (common::threshold::M < #moduli)
template <Backend BACKEND, unsigned IDX, bool COMPLEX = false, bool CENTER = true>
__device__ __forceinline__ int32_t mod_large(common::fp32_mant_exp a) {
    static_assert(CENTER || BACKEND == Backend::FP8);
    constexpr int32_t p = common::table::moduli<BACKEND, IDX, COMPLEX>;
    if constexpr (BACKEND == Backend::INT8 && p == 256) {
        const uint32_t exp_lo = a.exp.is_hi ? 0U : a.exp.val;
        const uint32_t prod   = uint32_t(a.mant) * exp_lo;
        return center_reduced<Backend::INT8, IDX, COMPLEX>(int32_t(prod & 255U));
    }
    if constexpr (BACKEND == Backend::FP8 && p == 1024) {
        const uint32_t exp_lo = a.exp.is_hi ? 0U : a.exp.val;
        const uint32_t prod   = uint32_t(a.mant) * exp_lo;
        if constexpr (CENTER) return center_reduced<Backend::FP8, IDX, COMPLEX>(int32_t(prod & 1023U));
        else return int32_t(prod & 1023U);
    }
    const int32_t rem1 = mod_small_nowrap<BACKEND, IDX, COMPLEX>(a.mant);
    if constexpr (BACKEND == Backend::INT8 && has_u8_dot4) {
        using D                    = int8_dot_traits<IDX, COMPLEX>;
        constexpr bool second_pass = !D::fp32_one_exp_pass_safe;
        const int32_t rem2         = reduce_exp<BACKEND, IDX, COMPLEX, second_pass>(a.exp);
        static_assert(second_pass || D::fp32_scaled_mant_max * uint64_t(D::exp_first_max) <= uint64_t(INT32_MAX));
        return mod_small<BACKEND, IDX, COMPLEX>(rem1 * rem2);
    } else {
        const int32_t rem2 = reduce_exp<BACKEND, IDX, COMPLEX>(a.exp);
        if constexpr (CENTER) return mod_small<BACKEND, IDX, COMPLEX>(rem1 * rem2);
        else return mod_small_nowrap<BACKEND, IDX, COMPLEX>(rem1 * rem2);
    }
}

template <Backend BACKEND, unsigned IDX, bool COMPLEX = false, bool CENTER = true>
__device__ __forceinline__ int32_t mod_large(common::fp64_mant_exp a) {
    static_assert(CENTER || BACKEND == Backend::FP8);
    constexpr int32_t p = common::table::moduli<BACKEND, IDX, COMPLEX>;
    if constexpr (BACKEND == Backend::INT8 && p == 256) {
        const uint32_t exp_lo = a.exp.is_hi ? 0U : a.exp.val;
        const uint32_t prod   = a.mant.lo * exp_lo;
        return center_reduced<Backend::INT8, IDX, COMPLEX>(int32_t(prod & 255U));
    }
    if constexpr (BACKEND == Backend::FP8 && p == 1024) {
        const uint32_t exp_lo = a.exp.is_hi ? 0U : a.exp.val;
        const uint32_t prod   = a.mant.lo * exp_lo;
        if constexpr (CENTER) return center_reduced<Backend::FP8, IDX, COMPLEX>(int32_t(prod & 1023U));
        else return int32_t(prod & 1023U);
    }

    const int32_t rem1 = [&] {
        if constexpr (BACKEND == Backend::FP8 && !has_u16_dot2) {
            using F = fp8_large_traits<IDX, COMPLEX>;
            if constexpr (F::skip_mant_reduction) {
                const int32_t hi  = mod_small_nowrap<BACKEND, IDX, COMPLEX>(a.mant.hi);
                const uint32_t lo = mod_small_nowrap_u32<BACKEND, IDX, COMPLEX>(a.mant.lo);
                return hi * F::r32_centered + int32_t(lo);
            }
        }
        return reduce_mant_large<BACKEND, IDX, COMPLEX>(a.mant);
    }();

    if constexpr (BACKEND == Backend::INT8 && has_u8_dot4) {
        using D                    = int8_dot_traits<IDX, COMPLEX>;
        constexpr bool second_pass = !D::fp64_one_exp_pass_safe;
        const int32_t rem2         = reduce_exp<BACKEND, IDX, COMPLEX, second_pass>(a.exp);
        static_assert(second_pass || D::fp64_scaled_mant_max * uint64_t(D::exp_first_max) <= uint64_t(INT32_MAX));
        return mod_small<BACKEND, IDX, COMPLEX>(rem1 * rem2);
    } else {
        const int32_t rem2 = reduce_exp<BACKEND, IDX, COMPLEX>(a.exp);
        if constexpr (CENTER) return mod_small<BACKEND, IDX, COMPLEX>(rem1 * rem2);
        else return mod_small_nowrap<BACKEND, IDX, COMPLEX>(rem1 * rem2);
    }
}

template <unsigned IDX, bool COMPLEX = false, typename V>
__device__ __forceinline__ int32_t calc_mod_fp8_nowrap(V a) {
    if constexpr (std::is_same_v<V, int32_t>) {
        return mod_small_nowrap<Backend::FP8, IDX, COMPLEX>(a);
    } else if constexpr (std::is_same_v<V, common::mant_t>) {
        return reduce_mant<Backend::FP8, IDX, COMPLEX>(a);
    } else {
        return mod_large<Backend::FP8, IDX, COMPLEX, false>(a);
    }
}

template <Backend BACKEND, unsigned IDX, bool COMPLEX = false> __device__ __forceinline__ int32_t calc_mod(int32_t a) { return mod_small<BACKEND, IDX, COMPLEX>(a); }
template <Backend BACKEND, unsigned IDX, bool COMPLEX = false> __device__ __forceinline__ int32_t calc_mod(common::mant_t a) { return mod_middle<BACKEND, IDX, COMPLEX>(a); }
template <Backend BACKEND, unsigned IDX, bool COMPLEX = false> __device__ __forceinline__ int32_t calc_mod(common::fp32_mant_exp a) { return mod_large<BACKEND, IDX, COMPLEX>(a); }
template <Backend BACKEND, unsigned IDX, bool COMPLEX = false> __device__ __forceinline__ int32_t calc_mod(common::fp64_mant_exp a) { return mod_large<BACKEND, IDX, COMPLEX>(a); }

} // namespace gemmul8::mod
