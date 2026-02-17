#pragma once

template <typename T> __forceinline__ __device__ int8_t trunc_scalbn_8i(T in, const int sft) {
    using I = int_t<T>;

    // fabs (sign clear)
    I bits = __fp_as_int<T>(in);
    bits &= ~extract_sign<T>(bits);
    if (bits == 0) return int8_t(0);

    int exp = (int)extract_exp<T>(bits);
    I frac  = extract_significand<T>(bits);

    I mant;
    int e;

    if (exp) {
        // normal
        mant = frac | (I(1) << fp<T>::prec);
        e    = exp - fp<T>::bias;
    } else {
        // subnormal (only place using clz)
        int k = countlz<T>(frac) - (fp<T>::bits - fp<T>::prec);
        mant  = (frac << k) | (I(1) << fp<T>::prec);
        e     = (1 - fp<T>::bias) - k;
    }

    // apply scalbn
    e += sft;

    // integer exponent alignment
    int shift = fp<T>::prec - e;

    // integer result (no rounding needed)
    if (shift <= 0) return int8_t(mant << (-shift));

    // small values
    if (shift >= fp<T>::prec + 1) return int8_t(1);

    // compute ceil
    I mask     = (I(1) << shift) - 1;
    I floor    = mant >> shift;
    I has_frac = (mant & mask) != 0;

    return int8_t(floor + has_frac);
};

template <typename T> __device__ __forceinline__ __nv_fp8_e4m3 fp8_e4m3_ru(T a) {
    __nv_fp8_e4m3 r               = __nv_fp8_e4m3(a);
    const T y                     = T(r);
    const __nv_fp8_storage_t bits = (__nv_fp8_storage_t)(unsigned(r.__x) + unsigned(y < a));
    r.__x                         = bits;
    return r;
}

template <gemmul8::Backend backend, typename T> __forceinline__ __device__ upperBound_t<backend, T> upperBound_lo(T in, const int sft);

template <> __forceinline__ __device__ int8_t upperBound_lo<gemmul8::Backend::INT8, float>(float in, const int sft) {
    return trunc_scalbn_8i<float>(in, sft);
}
template <> __forceinline__ __device__ int8_t upperBound_lo<gemmul8::Backend::INT8, double>(double in, const int sft) {
    return trunc_scalbn_8i<double>(in, sft);
}
template <> __forceinline__ __device__ char2 upperBound_lo<gemmul8::Backend::INT8, cuFloatComplex>(cuFloatComplex in, const int sft) {
    char2 out;
    out.x = trunc_scalbn_8i<float>(in.x, sft);
    out.y = trunc_scalbn_8i<float>(in.y, sft);
    return out;
}
template <> __forceinline__ __device__ char2 upperBound_lo<gemmul8::Backend::INT8, cuDoubleComplex>(cuDoubleComplex in, const int sft) {
    char2 out;
    out.x = trunc_scalbn_8i<double>(in.x, sft);
    out.y = trunc_scalbn_8i<double>(in.y, sft);
    return out;
}

template <> __forceinline__ __device__ __nv_fp8_e4m3 upperBound_lo<gemmul8::Backend::FP8, float>(float in, const int sft) {
    return fp8_e4m3_ru<float>(scalbnf(fabsf(in), sft));
}
template <> __forceinline__ __device__ __nv_fp8_e4m3 upperBound_lo<gemmul8::Backend::FP8, double>(double in, const int sft) {
    return fp8_e4m3_ru<double>(scalbn(fabs(in), sft));
}
template <> __forceinline__ __device__ fp8x2_e4m3 upperBound_lo<gemmul8::Backend::FP8, cuFloatComplex>(cuFloatComplex in, const int sft) {
    fp8x2_e4m3 out;
    out.x = fp8_e4m3_ru<float>(scalbnf(fabsf(in.x), sft));
    out.y = fp8_e4m3_ru<float>(scalbnf(fabsf(in.y), sft));
    return out;
}
template <> __forceinline__ __device__ fp8x2_e4m3 upperBound_lo<gemmul8::Backend::FP8, cuDoubleComplex>(cuDoubleComplex in, const int sft) {
    fp8x2_e4m3 out;
    out.x = fp8_e4m3_ru<double>(scalbn(fabs(in.x), sft));
    out.y = fp8_e4m3_ru<double>(scalbn(fabs(in.y), sft));
    return out;
}

//------------------------------
// Return trunc(scalbn(in, sft))
//------------------------------
template <typename T> __forceinline__ __device__ T trunc_scalbn_to_fp(T in, const int sft) {
    using I        = int_t<T>;
    I bits         = __fp_as_int<T>(in);
    const I sign   = extract_sign<T>(bits);
    int exp_biased = (int)extract_exp<T>(bits);
    I significand  = extract_significand<T>(bits);

    if (exp_biased != 0) {
        exp_biased += sft;
        if (exp_biased < fp<T>::bias) {
            return __int_as_fp<T>(sign);
        }
        if (exp_biased >= (fp<T>::bias + fp<T>::prec)) {
            bits = sign | (I(exp_biased) << fp<T>::prec) | significand;
            return __int_as_fp<T>(bits);
        }

        significand |= (I(1) << fp<T>::prec);
        int chop_bits = (fp<T>::bias + fp<T>::prec) - exp_biased;
        I mask        = I(-1) << chop_bits;
        significand   = extract_significand<T>(significand & mask);
        bits          = sign | (I)exp_biased << fp<T>::prec | significand;
        return __int_as_fp<T>(bits);
    }

    if (significand == 0) {
        return in;
    }

    int lz = (fp<T>::bits - fp<T>::prec) - countlz<T>(significand);
    int e  = lz + sft;

    if (e < fp<T>::bias) {
        return __int_as_fp<T>(sign);
    }

    I frac_full = (significand << (2 - lz)) ^ (I(1) << fp<T>::prec);
    I mask      = I(-1) << max(I(fp<T>::bias + fp<T>::prec - e), I(0));
    bits        = sign | I(e) << fp<T>::prec | (frac_full & mask);
    return __int_as_fp<T>(bits);
}

template <typename T> __forceinline__ __device__ int32_t trunc_scalbn_to_i32(T in, const int sft) {
    using I = int_t<T>;
    using U = uint_t<T>;

    const I bits         = __fp_as_int<T>(in);
    const bool neg       = (extract_sign<T>(bits) != I(0));
    const int exp_biased = (int)extract_exp<T>(bits);
    const U frac         = (U)extract_significand<T>(bits);

    // 0
    if (exp_biased == 0 && frac == 0) return 0;

    // make normalized significand S & unbiased exponent E
    U S;
    int E;
    if (exp_biased != 0) {
        // normal: value = (1.frac)*2^(exp-bias) = S * 2^(E - prec)
        S = (U(1) << fp<T>::prec) | frac; // S in [2^prec, 2^(prec+1))
        E = exp_biased - fp<T>::bias;
    } else {
        // subnormal: value = frac * 2^(1-bias-prec), frac != 0
        const int lz  = countlz<T>((I)frac);  // frac != 0
        const int msb = fp<T>::bits - 1 - lz; // 0..prec-1
        const int sh  = fp<T>::prec - msb;    // 1..prec
        S             = (U)(frac << sh);      // MSB -> bit prec
        E             = (1 - fp<T>::bias) - sh;
    }

    // scalbn
    const int E_scaled = E + sft;

    // |in|*2^sft < 1
    if (E_scaled < 0) return 0;

    // m = floor(|in|*2^sft)
    uint64_t mag = (uint64_t)S;
    if (E_scaled > fp<T>::prec) {
        const int sh = E_scaled - fp<T>::prec;
        mag <<= sh;
    } else {
        const int sh = fp<T>::prec - E_scaled;
        mag >>= sh;
    }

    const int32_t x = (int32_t)mag;
    return neg ? -x : x;
}

template <typename T> __forceinline__ __device__ int64_t trunc_scalbn_to_i64(T in, const int sft) {
    using I = int_t<T>;
    using U = uint_t<T>;

    const I bits         = __fp_as_int<T>(in);
    const bool neg       = (extract_sign<T>(bits) != I(0));
    const int exp_biased = (int)extract_exp<T>(bits);
    const U frac         = (U)extract_significand<T>(bits);

    // 0
    if (exp_biased == 0 && frac == 0) return 0;

    // make normalized significand S & unbiased exponent E
    U S;
    int E;
    if (exp_biased != 0) {
        // normal: value = (1.frac)*2^(exp-bias) = S * 2^(E - prec)
        S = (U(1) << fp<T>::prec) | frac; // S in [2^prec, 2^(prec+1))
        E = exp_biased - fp<T>::bias;
    } else {
        // subnormal: value = frac * 2^(1-bias-prec), frac != 0
        const int lz  = countlz<T>((I)frac);  // frac != 0
        const int msb = fp<T>::bits - 1 - lz; // 0..prec-1
        const int sh  = fp<T>::prec - msb;    // 1..prec
        S             = (U)(frac << sh);      // MSB -> bit prec
        E             = (1 - fp<T>::bias) - sh;
    }

    // scalbn
    const int E_scaled = E + sft;

    // |in|*2^sft < 1
    if (E_scaled < 0) return 0;

    // m = floor(|in|*2^sft)
    uint64_t mag = (uint64_t)S;
    if (E_scaled > fp<T>::prec) {
        const int sh = E_scaled - fp<T>::prec;
        mag <<= sh;
    } else {
        const int sh = fp<T>::prec - E_scaled;
        mag >>= sh;
    }

    const int64_t x = (int64_t)mag;
    return neg ? -x : x;
}

template <gemmul8::Backend backend, typename T, int N, typename Enable = void> struct trunc_scalbn;

// real
template <gemmul8::Backend backend, typename T, int N> struct trunc_scalbn<backend, T, N, std::enable_if_t<(N <= threshold<backend>::S && !isComplex<T>)>> {
    __device__ __forceinline__ static int32_t run(T in, int sft) {
        return trunc_scalbn_to_i32<T>(in, sft);
    }
};
template <gemmul8::Backend backend, typename T, int N> struct trunc_scalbn<backend, T, N, std::enable_if_t<(threshold<backend>::S < N && N <= threshold<backend>::M && !isComplex<T>)>> {
    __device__ __forceinline__ static int64_t run(T in, int sft) {
        return trunc_scalbn_to_i64<T>(in, sft);
    }
};
template <gemmul8::Backend backend, typename T, int N> struct trunc_scalbn<backend, T, N, std::enable_if_t<(threshold<backend>::M < N && !isComplex<T>)>> {
    __device__ __forceinline__ static T run(T in, int sft) {
        return trunc_scalbn_to_fp<T>(in, sft);
    }
};

// complex
template <gemmul8::Backend backend, typename T, int N> struct trunc_scalbn<backend, T, N, std::enable_if_t<(N <= threshold<backend>::S && isComplex<T>)>> {
    __device__ __forceinline__ static int2 run(T in, int sft) {
        int2 r;
        r.x = trunc_scalbn_to_i32<underlying_t<T>>(in.x, sft);
        r.y = trunc_scalbn_to_i32<underlying_t<T>>(in.y, sft);
        return r;
    }
};
template <gemmul8::Backend backend, typename T, int N> struct trunc_scalbn<backend, T, N, std::enable_if_t<(threshold<backend>::S < N && N <= threshold<backend>::M && isComplex<T>)>> {
    __device__ __forceinline__ static int64x2_t run(T in, int sft) {
        int64x2_t r;
        r.x = trunc_scalbn_to_i64<underlying_t<T>>(in.x, sft);
        r.y = trunc_scalbn_to_i64<underlying_t<T>>(in.y, sft);
        return r;
    }
};
template <gemmul8::Backend backend, typename T, int N> struct trunc_scalbn<backend, T, N, std::enable_if_t<(threshold<backend>::M < N && isComplex<T>)>> {
    __device__ __forceinline__ static T run(T in, int sft) {
        T r;
        r.x = trunc_scalbn_to_fp<underlying_t<T>>(in.x, sft);
        r.y = trunc_scalbn_to_fp<underlying_t<T>>(in.y, sft);
        return r;
    }
};
