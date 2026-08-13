#pragma once
#include "../common/common.hpp"
#include "../common/table.hpp"

namespace gemmul8::mod {

//------------------------------
// Calculate mod: a - round(a/p(j))*p(j)
//------------------------------

// return value in [-p/2, p/2]
template <Backend BACKEND, unsigned IDX>
__device__ __forceinline__ int32_t wrapping(int32_t a) {
    constexpr int32_t p      = common::table::moduli<BACKEND, IDX>;
    constexpr int32_t p_half = p / 2;
    return (a > p_half) ? (a - p) : ((a < -p_half) ? (a + p) : a);
}

template <Backend BACKEND, unsigned IDX>
__device__ __forceinline__ uint32_t mod_small_nowrap_u32(uint32_t a) {
    constexpr uint32_t up = uint32_t(common::table::moduli<BACKEND, IDX>);
    if constexpr (BACKEND == Backend::INT8 && up == 256U) {
        return a & 255U;
    }
    if constexpr (BACKEND == Backend::FP8 && up == 1024U) {
        return a & 1023U;
    }
#if defined(__CUDACC__)
    if constexpr (BACKEND == Backend::INT8) {
        constexpr uint32_t pow2_00_mod = 1U % up;
        constexpr uint32_t pow2_08_mod = (pow2_00_mod << 8) % up;
        constexpr uint32_t pow2_16_mod = (pow2_08_mod << 8) % up;
        constexpr uint32_t pow2_24_mod = (pow2_16_mod << 8) % up;
        constexpr uint32_t weight      = pow2_00_mod | (pow2_08_mod << 8) | (pow2_16_mod << 16) | (pow2_24_mod << 24);
        return __dp4a(a, weight, 0U); // < 2^17
    }
#endif
    constexpr uint32_t p_inv_u32 = uint32_t(common::table::p_inv_32_v<BACKEND, IDX>); // 2^32/p
    const uint32_t rem           = a - up * __umulhi(a, p_inv_u32);
    return rem;
}

// |a| < 2^31 is guaranteed (#moduli <= common::threshold::S)
template <Backend BACKEND, unsigned IDX>
__device__ __forceinline__ int32_t mod_small_nowrap(int32_t a) {
    constexpr uint32_t up = uint32_t(common::table::moduli<BACKEND, IDX>);
    if constexpr (BACKEND == Backend::INT8 && up == 256U) {
        return a & 255;
    }
    if constexpr (BACKEND == Backend::FP8 && up == 1024U) {
        return a & 1023;
    }
#if defined(__CUDACC__)
    if constexpr (BACKEND == Backend::INT8) {
        constexpr uint32_t pow2_00_mod = 1U % up;
        constexpr uint32_t pow2_08_mod = (pow2_00_mod << 8) % up;
        constexpr uint32_t pow2_16_mod = (pow2_08_mod << 8) % up;
        constexpr uint32_t pow2_24_mod = (pow2_16_mod << 8) % up;
        constexpr uint32_t pow2_32_mod = (pow2_24_mod << 8) % up;
        constexpr uint32_t weight      = pow2_00_mod | (pow2_08_mod << 8) | (pow2_16_mod << 16) | (pow2_24_mod << 24);
        constexpr uint32_t neg_corr    = (pow2_32_mod == 0U) ? 0U : (up - pow2_32_mod);
        return int32_t(__dp4a(uint32_t(a), weight, (a < 0) ? neg_corr : 0U)); // < 2^17
    }
#endif
    constexpr int32_t p         = common::table::moduli<BACKEND, IDX>;
    constexpr int32_t p_inv_i32 = common::table::p_inv_32_v<BACKEND, IDX>; // 2^32/p
    const int32_t rem           = a - p * __mulhi(a, p_inv_i32);
    return rem;
}

template <Backend BACKEND, unsigned IDX>
__device__ __forceinline__ int32_t mod_small(int32_t a) {
    constexpr int32_t p = common::table::moduli<BACKEND, IDX>;
    if constexpr (BACKEND == Backend::INT8 && p == 256) {
        return wrapping<Backend::INT8, IDX>(a & 255);
    }
    if constexpr (BACKEND == Backend::FP8 && p == 1024) {
        return wrapping<Backend::FP8, IDX>(a & 1023);
    }
#if defined(__CUDACC__)
    if constexpr (BACKEND == Backend::INT8 && p == 255) {
        const uint32_t s = __dp4a(uint32_t(a), 0x01010101U, (a < 0) ? 254U : 0U);
        const uint32_t t = (s & 255U) + ((s >> 8));
        return wrapping<Backend::INT8, IDX>(int32_t(t));
    }
#endif
    constexpr int32_t p_inv_i32 = common::table::p_inv_32_v<BACKEND, IDX>; // 2^32/p
    const int32_t rem           = a - p * __mulhi(a, p_inv_i32);
    return wrapping<BACKEND, IDX>(rem);
}

template <Backend BACKEND, unsigned IDX>
__device__ __forceinline__ int32_t reduce_mant(common::mant_t a) {
    constexpr uint32_t up = uint32_t(common::table::moduli<BACKEND, IDX>);
    if constexpr (BACKEND == Backend::INT8 && up == 256U) {
        return int32_t(a.lo & 255U);
    }
    if constexpr (BACKEND == Backend::FP8 && up == 1024U) {
        return int32_t(a.lo & 1023U);
    }
#if defined(__CUDACC__)
    if constexpr (BACKEND == Backend::INT8) {
        constexpr uint32_t pow2_00_mod = 1U % up;
        constexpr uint32_t pow2_08_mod = (pow2_00_mod << 8) % up;
        constexpr uint32_t pow2_16_mod = (pow2_08_mod << 8) % up;
        constexpr uint32_t pow2_24_mod = (pow2_16_mod << 8) % up;
        constexpr uint32_t pow2_32_mod = (pow2_24_mod << 8) % up;
        constexpr uint32_t pow2_40_mod = (pow2_32_mod << 8) % up;
        constexpr uint32_t pow2_48_mod = (pow2_40_mod << 8) % up;
        constexpr uint32_t pow2_56_mod = (pow2_48_mod << 8) % up;
        constexpr uint32_t pow2_64_mod = (pow2_56_mod << 8) % up;
        constexpr uint32_t neg_corr    = (pow2_64_mod == 0U) ? 0U : (up - pow2_64_mod);
        constexpr uint32_t weight_lo   = pow2_00_mod | (pow2_08_mod << 8) | (pow2_16_mod << 16) | (pow2_24_mod << 24);
        constexpr uint32_t weight_hi   = pow2_32_mod | (pow2_40_mod << 8) | (pow2_48_mod << 16) | (pow2_56_mod << 24);

        const uint32_t acc0   = (a.hi < 0) ? neg_corr : 0U;
        const uint32_t rem_hi = __dp4a(uint32_t(a.hi), weight_hi, acc0);
        const uint32_t raw    = __dp4a(a.lo, weight_lo, rem_hi);
        return int32_t(raw); // < 2^18
    }
#endif
    constexpr uint64_t p          = uint64_t(common::table::moduli<BACKEND, IDX>);
    constexpr int32_t pow2_32_mod = int32_t((uint64_t(1) << 32) % p);
    const int32_t rem_hi          = mod_small_nowrap<BACKEND, IDX>(a.hi);
    const uint32_t rem_lo         = mod_small_nowrap_u32<BACKEND, IDX>(a.lo);
    const int32_t raw             = rem_hi * pow2_32_mod + int32_t(rem_lo);
    if constexpr (BACKEND == Backend::INT8) {
        return raw;
    } else {
        return mod_small_nowrap<BACKEND, IDX>(raw);
    }
}

template <Backend BACKEND, unsigned IDX>
__device__ __forceinline__ int32_t reduce_mant_large(common::mant_t a) {
    constexpr uint32_t up = uint32_t(common::table::moduli<BACKEND, IDX>);
    if constexpr (BACKEND == Backend::INT8 && up == 256U) {
        return int32_t(a.lo & 255U);
    }
    if constexpr (BACKEND == Backend::FP8 && up == 1024U) {
        return int32_t(a.lo & 1023U);
    }
#if defined(__CUDACC__)
    if constexpr (BACKEND == Backend::INT8) {
        constexpr uint32_t pow2_00_mod = 1U % up;
        constexpr uint32_t pow2_08_mod = (pow2_00_mod << 8) % up;
        constexpr uint32_t pow2_16_mod = (pow2_08_mod << 8) % up;
        constexpr uint32_t pow2_24_mod = (pow2_16_mod << 8) % up;
        constexpr uint32_t pow2_32_mod = (pow2_24_mod << 8) % up;
        constexpr uint32_t pow2_40_mod = (pow2_32_mod << 8) % up;
        constexpr uint32_t pow2_48_mod = (pow2_40_mod << 8) % up;
        constexpr uint32_t pow2_56_mod = (pow2_48_mod << 8) % up;
        constexpr uint32_t pow2_64_mod = (pow2_56_mod << 8) % up;
        constexpr uint32_t neg_corr    = (pow2_64_mod == 0U) ? 0U : (up - pow2_64_mod);
        constexpr uint32_t weight_lo   = pow2_00_mod | (pow2_08_mod << 8) | (pow2_16_mod << 16) | (pow2_24_mod << 24);
        constexpr uint32_t weight_hi   = pow2_32_mod | (pow2_40_mod << 8) | (pow2_48_mod << 16) | (pow2_56_mod << 24);

        const uint32_t acc0   = (a.hi < 0) ? neg_corr : 0U;
        const uint32_t rem_hi = __dp4a(uint32_t(a.hi), weight_hi, acc0);
        const uint32_t raw    = __dp4a(a.lo, weight_lo, rem_hi);
        return int32_t(__dp4a(raw, weight_lo, 0U)); // < 2^15
    }
#endif
    constexpr uint64_t p          = uint64_t(common::table::moduli<BACKEND, IDX>);
    constexpr int32_t pow2_32_mod = int32_t((uint64_t(1) << 32) % p);
    const int32_t rem_hi          = mod_small_nowrap<BACKEND, IDX>(a.hi);
    const uint32_t rem_lo         = mod_small_nowrap_u32<BACKEND, IDX>(a.lo);
    const int32_t raw             = rem_hi * pow2_32_mod + int32_t(rem_lo);
    if constexpr (BACKEND == Backend::INT8) {
        return raw;
    } else {
        return mod_small_nowrap<BACKEND, IDX>(raw);
    }
}

template <Backend BACKEND, unsigned IDX>
__device__ __forceinline__ int32_t reduce_exp(common::exp_t a) {
    constexpr uint32_t up = uint32_t(common::table::moduli<BACKEND, IDX>);
    if constexpr (BACKEND == Backend::INT8 && up == 256U) {
        return int32_t((a.is_hi ? 0 : a.val) & 255U);
    }
    if constexpr (BACKEND == Backend::FP8 && up == 1024U) {
        return int32_t((a.is_hi ? 0 : a.val) & 1023U);
    }
#if defined(__CUDACC__)
    if constexpr (BACKEND == Backend::INT8) {
        constexpr uint32_t pow2_00_mod = 1U % up;
        constexpr uint32_t pow2_08_mod = (pow2_00_mod << 8) % up;
        constexpr uint32_t pow2_16_mod = (pow2_08_mod << 8) % up;
        constexpr uint32_t pow2_24_mod = (pow2_16_mod << 8) % up;
        constexpr uint32_t pow2_32_mod = (pow2_24_mod << 8) % up;
        constexpr uint32_t pow2_40_mod = (pow2_32_mod << 8) % up;
        constexpr uint32_t pow2_48_mod = (pow2_40_mod << 8) % up;
        constexpr uint32_t pow2_56_mod = (pow2_48_mod << 8) % up;
        constexpr uint32_t weight_lo   = pow2_00_mod | (pow2_08_mod << 8) | (pow2_16_mod << 16) | (pow2_24_mod << 24);
        constexpr uint32_t weight_hi   = pow2_32_mod | (pow2_40_mod << 8) | (pow2_48_mod << 16) | (pow2_56_mod << 24);
        const uint32_t raw             = __dp4a(a.val, a.is_hi ? weight_hi : weight_lo, 0U);
        return int32_t(__dp4a(raw, weight_lo, 0U)); // < 2^15
    }
#endif
    constexpr uint64_t p          = uint64_t(common::table::moduli<BACKEND, IDX>);
    constexpr int32_t pow2_32_mod = int32_t((uint64_t(1) << 32) % p);
    const int32_t sft             = a.is_hi ? pow2_32_mod : 1;
    const int32_t raw             = int32_t(mod_small_nowrap_u32<BACKEND, IDX>(a.val)) * sft;
    if constexpr (BACKEND == Backend::INT8) {
        return raw;
    } else {
        return mod_small_nowrap<BACKEND, IDX>(raw);
    }
}

// |a| < 2^63 is guaranteed (common::threshold::S < #moduli <= common::threshold::M)
template <Backend BACKEND, unsigned IDX>
__device__ __forceinline__ int32_t mod_middle(common::mant_t a) {
    constexpr int32_t p = common::table::moduli<BACKEND, IDX>;
    if constexpr (BACKEND == Backend::INT8 && p == 256) {
        return wrapping<Backend::INT8, IDX>(int32_t(a.lo & 255U));
    }
    if constexpr (BACKEND == Backend::FP8 && p == 1024) {
        return wrapping<Backend::FP8, IDX>(int32_t(a.lo & 1023U));
    }
    const int32_t rem = reduce_mant<BACKEND, IDX>(a);
    if constexpr (BACKEND == Backend::INT8) {
        return mod_small<BACKEND, IDX>(rem);
    } else {
        return wrapping<BACKEND, IDX>(rem);
    }
}

// |a| can be >= 2^63 (common::threshold::M < #moduli)
template <Backend BACKEND, unsigned IDX>
__device__ __forceinline__ int32_t mod_large(common::fp32_mant_exp a) {
    constexpr int32_t p = common::table::moduli<BACKEND, IDX>;
    if constexpr (BACKEND == Backend::INT8 && p == 256) {
        const uint32_t exp_lo = a.exp.is_hi ? 0U : a.exp.val;
        const uint32_t prod   = a.mant * exp_lo;
        return wrapping<Backend::INT8, IDX>(int32_t(prod & 255U));
    }
    if constexpr (BACKEND == Backend::FP8 && p == 1024) {
        const uint32_t exp_lo = a.exp.is_hi ? 0U : a.exp.val;
        const uint32_t prod   = a.mant * exp_lo;
        return wrapping<Backend::FP8, IDX>(int32_t(prod & 1023U));
    }
    const int32_t rem1 = mod_small_nowrap<BACKEND, IDX>(a.mant);
    const int32_t rem2 = reduce_exp<BACKEND, IDX>(a.exp);
    return mod_small<BACKEND, IDX>(rem1 * rem2);
}

template <Backend BACKEND, unsigned IDX>
__device__ __forceinline__ int32_t mod_large(common::fp64_mant_exp a) {
    constexpr int32_t p = common::table::moduli<BACKEND, IDX>;
    if constexpr (BACKEND == Backend::INT8 && p == 256) {
        const uint32_t exp_lo = a.exp.is_hi ? 0U : a.exp.val;
        const uint32_t prod   = a.mant.lo * exp_lo;
        return wrapping<Backend::INT8, IDX>(int32_t(prod & 255U));
    }
    if constexpr (BACKEND == Backend::FP8 && p == 1024) {
        const uint32_t exp_lo = a.exp.is_hi ? 0U : a.exp.val;
        const uint32_t prod   = a.mant.lo * exp_lo;
        return wrapping<Backend::FP8, IDX>(int32_t(prod & 1023U));
    }
    const int32_t rem1 = reduce_mant_large<BACKEND, IDX>(a.mant);
    const int32_t rem2 = reduce_exp<BACKEND, IDX>(a.exp);
    return mod_small<BACKEND, IDX>(rem1 * rem2);
}

template <Backend BACKEND, unsigned IDX> __device__ __forceinline__ int32_t calc_mod(int32_t a) { return mod_small<BACKEND, IDX>(a); }
template <Backend BACKEND, unsigned IDX> __device__ __forceinline__ int32_t calc_mod(common::mant_t a) { return mod_middle<BACKEND, IDX>(a); }
template <Backend BACKEND, unsigned IDX> __device__ __forceinline__ int32_t calc_mod(common::fp32_mant_exp a) { return mod_large<BACKEND, IDX>(a); }
template <Backend BACKEND, unsigned IDX> __device__ __forceinline__ int32_t calc_mod(common::fp64_mant_exp a) { return mod_large<BACKEND, IDX>(a); }

} // namespace gemmul8::mod
