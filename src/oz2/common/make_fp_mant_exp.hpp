#pragma once
#include "common.hpp"

namespace gemmul8::common {

//------------------------------
// Transform a floating-point integer into mant*exp where exp = 2^x >= 1
//------------------------------

__device__ __forceinline__ fp32_mant_exp make_fp32_mant_exp(float a) {

    fp32_mant_exp d;

    // bit pattern
    const uint32_t bits = __float_as_uint(a);
    const uint32_t sign = bits >> 31;
    const uint32_t e    = (bits >> 23) & 0xFFu; // exponent bits
    const uint32_t frac = bits & 0x7FFFFFu;     // fraction bits

    // a == 0
    if (e == 0u) {
        d.mant = 0;
        d.exp  = common::exp_t{};
        return d;
    }

    const int32_t unbiased = int32_t(e) - 127;  // floor(log2(|a|)) for normalized integer a
    const uint32_t sig     = (1u << 23) | frac; // 24-bit significand
    const int32_t out_exp  = max(unbiased - 30, 0);

    int32_t mant;
    if (unbiased > 30) {
        mant = (int32_t)(sig << 7);
    } else {
        const int32_t sh   = unbiased - 23;
        const uint32_t mag = (sh >= 0) ? (sig << sh) : (sig >> (-sh));
        mant               = (int32_t)mag;
    }

    const bool is_hi = (out_exp >= 32);
    d.mant           = sign ? -mant : mant;
    d.exp            = common::exp_t(uint32_t(1u << (out_exp - (is_hi ? 32 : 0))), is_hi);

    return d;
}

template <typename T>
__device__ __forceinline__ T make_fp64_mant_exp(double a);

template <>
__device__ __forceinline__ fp64_mant_exp make_fp64_mant_exp<fp64_mant_exp>(double a) {
    fp64_mant_exp d;

    const uint64_t bits = (uint64_t)__double_as_longlong(a);
    const uint64_t sign = bits >> 63;
    const uint64_t e    = (bits >> 52) & 0x7FFull;      // exponent bits
    const uint64_t frac = bits & ((1ull << 52) - 1ull); // fraction bits

    if (e == 0ull) {
        d.mant = common::mant_t{};
        d.exp  = common::exp_t{};
        return d;
    }

    const int32_t unbiased = int32_t(e) - 1023;   // floor(log2(|a|)) for normalized integer a
    const uint64_t sig     = (1ull << 52) | frac; // 53-bit significand
    const int32_t out_exp  = max(unbiased - 62, 0);

    int64_t mant;
    if (unbiased > 62) {
        mant = (int64_t)(sig << 10);
    } else {
        const int32_t sh   = unbiased - 52;
        const uint64_t mag = (sh >= 0) ? (sig << sh) : (sig >> (-sh));
        mant               = (int64_t)mag;
    }

    const bool is_hi = (out_exp >= 32);
    d.mant           = common::mant_t(sign ? -mant : mant);
    d.exp            = common::exp_t(uint32_t(1u << (out_exp - (is_hi ? 32 : 0))), is_hi);

    return d;
}

__device__ __forceinline__ fp32_mant_exp2 make_fp32_mant_exp2(cuFloatComplex v) {
    fp32_mant_exp2 d;
    d.x = make_fp32_mant_exp(v.x);
    d.y = make_fp32_mant_exp(v.y);
    return d;
}

template <typename T>
__device__ __forceinline__ T make_fp64_mant_exp2(cuDoubleComplex v);

template <>
__device__ __forceinline__ fp64_mant_exp2 make_fp64_mant_exp2<fp64_mant_exp2>(cuDoubleComplex v) {
    fp64_mant_exp2 d;
    d.x = make_fp64_mant_exp<fp64_mant_exp>(v.x);
    d.y = make_fp64_mant_exp<fp64_mant_exp>(v.y);
    return d;
}

template <typename V>
__device__ __forceinline__ V make_fp32_mant_exp_as(float v) {
    if constexpr (std::is_same_v<V, fp64_mant_exp>) {
        return make_fp64_mant_exp<fp64_mant_exp>(static_cast<double>(v));
    } else {
        static_assert(std::is_same_v<V, fp32_mant_exp>);
        return make_fp32_mant_exp(v);
    }
}

template <typename V>
__device__ __forceinline__ V make_fp32_mant_exp2_as(cuFloatComplex v) {
    V d;
    using Part = decltype(d.x);
    d.x        = make_fp32_mant_exp_as<Part>(v.x);
    d.y        = make_fp32_mant_exp_as<Part>(v.y);
    return d;
}

} // namespace gemmul8::common
