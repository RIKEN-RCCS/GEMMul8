#pragma once
#include "../../common/common.hpp"

namespace gemmul8::scaling::general {

__device__ __forceinline__ int32_t trunc_scalbn_to_i32_bit(double x, const int32_t sft) {
    const uint64_t bits = static_cast<uint64_t>(__double_as_longlong(x));

    const uint64_t sign = bits >> 63;
    const uint64_t e    = (bits >> 52) & 0x7ffULL;
    const uint64_t frac = bits & ((1ULL << 52) - 1ULL);

    if (e == 0ULL) {
        return __double2int_rz(scalbn(x, sft));
    }

    const int32_t exp = int32_t(e) - 1023 + sft;

    if (exp < 0) {
        return 0;
    }

    const uint64_t sig = (1ULL << 52) | frac;

    uint32_t mag;
    if (exp >= 52) {
        mag = uint32_t(sig << (exp - 52));
    } else {
        mag = uint32_t(sig >> (52 - exp));
    }

    const int32_t y = sign ? -static_cast<int32_t>(mag) : static_cast<int32_t>(mag);
    return y;
}

__device__ __forceinline__ int32_t trunc_scalbn_to_i32_bit(float x, const int32_t sft) {
    const uint32_t bits = __float_as_uint(x);

    const uint32_t sign = bits >> 31;
    const uint32_t e    = (bits >> 23) & 0xffU;
    const uint32_t frac = bits & ((1U << 23) - 1U);

    if (e == 0U) {
        return __float2int_rz(scalbnf(x, sft));
    }

    const int32_t exp = int32_t(e) - 127 + sft;

    if (exp < 0) {
        return 0;
    }

    const uint32_t sig = (1U << 23) | frac;

    uint32_t mag;
    if (exp >= 23) {
        mag = sig << (exp - 23);
    } else {
        mag = sig >> (23 - exp);
    }

    const int32_t y = sign ? -static_cast<int32_t>(mag) : static_cast<int32_t>(mag);
    return y;
}

__device__ __forceinline__ common::mant_t trunc_scalbn_to_i64_bit(double x, const int32_t sft) {
    const uint64_t bits = static_cast<uint64_t>(__double_as_longlong(x));

    const uint64_t sign = bits >> 63;
    const uint64_t e    = (bits >> 52) & 0x7ffULL;
    const uint64_t frac = bits & ((1ULL << 52) - 1ULL);

    if (e == 0ULL) {
        return common::mant_t(__double2ll_rz(scalbn(x, sft)));
    }

    const int32_t exp = int32_t(e) - 1023 + sft;

    if (exp < 0) {
        return common::mant_t{};
    }

    const uint64_t sig = (1ULL << 52) | frac;

    uint64_t mag;
    if (exp >= 52) {
        mag = sig << (exp - 52);
    } else {
        mag = sig >> (52 - exp);
    }

    const int64_t y = sign ? -static_cast<int64_t>(mag) : static_cast<int64_t>(mag);
    return common::mant_t(y);
}

__device__ __forceinline__ common::mant_t trunc_scalbn_to_i64_bit(float x, const int32_t sft) {
    const uint32_t bits = __float_as_uint(x);

    const uint32_t sign = bits >> 31;
    const uint32_t e    = (bits >> 23) & 0xffU;
    const uint32_t frac = bits & ((1U << 23) - 1U);

    if (e == 0U) {
        return common::mant_t(__float2ll_rz(scalbnf(x, sft)));
    }

    const int32_t exp = int32_t(e) - 127 + sft;

    if (exp < 0) {
        return common::mant_t{};
    }

    const uint32_t sig = (1U << 23) | frac;

    uint64_t mag;
    if (exp >= 23) {
        mag = uint64_t(sig) << (exp - 23);
    } else {
        mag = uint64_t(sig >> (23 - exp));
    }

    const int64_t y = sign ? -static_cast<int64_t>(mag) : static_cast<int64_t>(mag);
    return common::mant_t(y);
}

// trunc(scalbn(in, sft))
template <typename T>
__device__ __forceinline__ T trunc_scalbn_to_fp(T in, const int32_t sft) {
    return common::Ttrunc<T>(common::Tscalbn<T>(in, sft));
}

template <typename T>
__device__ __forceinline__ int32_t trunc_scalbn_to_i32(T in, const int32_t sft) {
    // return common::__fp2int_rz<T>(common::Tscalbn<T>(in, sft));
    return trunc_scalbn_to_i32_bit(in, sft);
}

template <typename T>
__device__ __forceinline__ common::mant_t trunc_scalbn_to_i64(T in, const int32_t sft) {
    // return common::mant_t(common::__fp2ll_rz<T>(common::Tscalbn<T>(in, sft)));
    return trunc_scalbn_to_i64_bit(in, sft);
}

template <bool cast_flag, typename T, Backend BACKEND, unsigned NUM_MODULI, typename Enable = void> struct trunc_scalbn;

// NUM_MODULI <= common::threshold<BACKEND>::S, real
template <typename T, Backend BACKEND, unsigned NUM_MODULI>
struct trunc_scalbn<true, T, BACKEND, NUM_MODULI,
                    std::enable_if_t<(NUM_MODULI <= common::threshold<BACKEND>::S &&
                                      !common::isComplex<T>)>> {
    __device__ __forceinline__ static int32_t run(T in, int32_t sft) { return trunc_scalbn_to_i32<T>(in, sft); }
    __device__ __forceinline__ static int32_t cast(int32_t in) { return in; }
};
template <typename T, Backend BACKEND, unsigned NUM_MODULI>
struct trunc_scalbn<false, T, BACKEND, NUM_MODULI,
                    std::enable_if_t<(NUM_MODULI <= common::threshold<BACKEND>::S &&
                                      !common::isComplex<T>)>> {
    __device__ __forceinline__ static T run(T in, int32_t sft) { return common::Tscalbn<T>(in, sft); }
    __device__ __forceinline__ static int32_t cast(T in) { return common::__fp2int_rz<T>(in); }
};

// common::threshold<BACKEND>::S < NUM_MODULI <= common::threshold<BACKEND>::M, real
template <typename T, Backend BACKEND, unsigned NUM_MODULI>
struct trunc_scalbn<true, T, BACKEND, NUM_MODULI,
                    std::enable_if_t<(common::threshold<BACKEND>::S < NUM_MODULI &&
                                      NUM_MODULI <= common::threshold<BACKEND>::M &&
                                      !common::isComplex<T>)>> {
    __device__ __forceinline__ static common::mant_t run(T in, int32_t sft) { return trunc_scalbn_to_i64<T>(in, sft); }
    __device__ __forceinline__ static common::mant_t cast(common::mant_t in) { return in; }
};
template <typename T, Backend BACKEND, unsigned NUM_MODULI>
struct trunc_scalbn<false, T, BACKEND, NUM_MODULI,
                    std::enable_if_t<(common::threshold<BACKEND>::S < NUM_MODULI &&
                                      NUM_MODULI <= common::threshold<BACKEND>::M &&
                                      !common::isComplex<T>)>> {
    __device__ __forceinline__ static T run(T in, int32_t sft) { return common::Tscalbn<T>(in, sft); }
    __device__ __forceinline__ static common::mant_t cast(T in) { return common::mant_t(common::__fp2ll_rz<T>(in)); }
};

// common::threshold<BACKEND>::M < NUM_MODULI, real
template <typename T, Backend BACKEND, unsigned NUM_MODULI>
struct trunc_scalbn<true, T, BACKEND, NUM_MODULI,
                    std::enable_if_t<(common::threshold<BACKEND>::M < NUM_MODULI &&
                                      !common::isComplex<T>)>> {
    __device__ __forceinline__ static T run(T in, int32_t sft) { return trunc_scalbn_to_fp<T>(in, sft); }
    __device__ __forceinline__ static T cast(T in) { return in; }
};
template <typename T, Backend BACKEND, unsigned NUM_MODULI>
struct trunc_scalbn<false, T, BACKEND, NUM_MODULI,
                    std::enable_if_t<(common::threshold<BACKEND>::M < NUM_MODULI &&
                                      !common::isComplex<T>)>> {
    __device__ __forceinline__ static T run(T in, int32_t sft) { return common::Tscalbn<T>(in, sft); }
    __device__ __forceinline__ static T cast(T in) { return common::Ttrunc<T>(in); }
};

// NUM_MODULI <= common::threshold<BACKEND, true>::S, complex
template <typename T, Backend BACKEND, unsigned NUM_MODULI>
struct trunc_scalbn<true, T, BACKEND, NUM_MODULI,
                    std::enable_if_t<(NUM_MODULI <= common::threshold<BACKEND, true>::S &&
                                      common::isComplex<T>)>> {
    __device__ __forceinline__ static int2 run(T in, int32_t sft) {
        int2 r;
        r.x = trunc_scalbn_to_i32<common::underlying_t<T>>(in.x, sft);
        r.y = trunc_scalbn_to_i32<common::underlying_t<T>>(in.y, sft);
        return r;
    }
    __device__ __forceinline__ static int2 cast(int2 in) { return in; }
};
template <typename T, Backend BACKEND, unsigned NUM_MODULI>
struct trunc_scalbn<false, T, BACKEND, NUM_MODULI,
                    std::enable_if_t<(NUM_MODULI <= common::threshold<BACKEND, true>::S &&
                                      common::isComplex<T>)>> {
    __device__ __forceinline__ static T run(T in, int32_t sft) { return common::Tscalbn<T>(in, sft); }
    __device__ __forceinline__ static int2 cast(T in) {
        int2 r;
        r.x = common::__fp2int_rz<common::underlying_t<T>>(in.x);
        r.y = common::__fp2int_rz<common::underlying_t<T>>(in.y);
        return r;
    }
};

// common::threshold<BACKEND, true>::S < NUM_MODULI <= common::threshold<BACKEND, true>::M, complex
template <typename T, Backend BACKEND, unsigned NUM_MODULI>
struct trunc_scalbn<true, T, BACKEND, NUM_MODULI,
                    std::enable_if_t<(common::threshold<BACKEND, true>::S < NUM_MODULI &&
                                      NUM_MODULI <= common::threshold<BACKEND, true>::M &&
                                      common::isComplex<T>)>> {
    __device__ __forceinline__ static common::mant2_t run(T in, int32_t sft) {
        common::mant2_t r;
        r.x = trunc_scalbn_to_i64<common::underlying_t<T>>(in.x, sft);
        r.y = trunc_scalbn_to_i64<common::underlying_t<T>>(in.y, sft);
        return r;
    }
    __device__ __forceinline__ static common::mant2_t cast(common::mant2_t in) { return in; }
};
template <typename T, Backend BACKEND, unsigned NUM_MODULI>
struct trunc_scalbn<false, T, BACKEND, NUM_MODULI,
                    std::enable_if_t<(common::threshold<BACKEND, true>::S < NUM_MODULI &&
                                      NUM_MODULI <= common::threshold<BACKEND, true>::M &&
                                      common::isComplex<T>)>> {
    __device__ __forceinline__ static T run(T in, int32_t sft) { return common::Tscalbn<T>(in, sft); }
    __device__ __forceinline__ static common::mant2_t cast(T in) {
        common::mant2_t r;
        r.x = common::mant_t(common::__fp2ll_rz<common::underlying_t<T>>(in.x));
        r.y = common::mant_t(common::__fp2ll_rz<common::underlying_t<T>>(in.y));
        return r;
    }
};

// common::threshold<BACKEND, true>::M < NUM_MODULI, complex
template <typename T, Backend BACKEND, unsigned NUM_MODULI>
struct trunc_scalbn<true, T, BACKEND, NUM_MODULI,
                    std::enable_if_t<(common::threshold<BACKEND, true>::M < NUM_MODULI &&
                                      common::isComplex<T>)>> {
    __device__ __forceinline__ static T run(T in, int32_t sft) {
        T r;
        r.x = trunc_scalbn_to_fp<common::underlying_t<T>>(in.x, sft);
        r.y = trunc_scalbn_to_fp<common::underlying_t<T>>(in.y, sft);
        return r;
    }
    __device__ __forceinline__ static T cast(T in) { return in; }
};
template <typename T, Backend BACKEND, unsigned NUM_MODULI>
struct trunc_scalbn<false, T, BACKEND, NUM_MODULI,
                    std::enable_if_t<(common::threshold<BACKEND, true>::M < NUM_MODULI &&
                                      common::isComplex<T>)>> {
    __device__ __forceinline__ static T run(T in, int32_t sft) { return common::Tscalbn<T>(in, sft); }
    __device__ __forceinline__ static T cast(T in) {
        T r;
        r.x = common::Ttrunc<common::underlying_t<T>>(in.x);
        r.y = common::Ttrunc<common::underlying_t<T>>(in.y);
        return r;
    }
};

} // namespace gemmul8::scaling::general
