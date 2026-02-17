#pragma once

//------------------------------
// Calculate mod: a - round(a/p(j))*p(j)
//------------------------------

// return value in [-p/2, p/2]
template <gemmul8::Backend backend, int IDX> __forceinline__ __device__ int32_t wrapping(int32_t a) {
    constexpr int32_t p      = table::moduli<backend, IDX>;
    constexpr int32_t p_half = p / 2;
    return (a > p_half) ? (a - p) : ((a < -p_half) ? (a + p) : a);
}

// |a| < 2^31 is guaranteed (#moduli <= threshold::S)
template <gemmul8::Backend backend, int IDX> __forceinline__ __device__ int32_t mod_small(int32_t a) {
    constexpr int32_t p        = table::moduli<backend, IDX>;
    constexpr int32_t p_inv_32 = int32_t(4294967296ULL / uint64_t(p)); // 2^32/p

    const int32_t rem = a - p * __mulhi(a, p_inv_32);
    return wrapping<backend, IDX>(rem);
}
template <gemmul8::Backend backend, int IDX> __forceinline__ __device__ int32_t mod_small_nowrap(int32_t a) {
    constexpr int32_t p        = table::moduli<backend, IDX>;
    constexpr int32_t p_inv_32 = int32_t(4294967296ULL / uint64_t(p)); // 2^32/p

    const int32_t rem = a - p * __mulhi(a, p_inv_32);
    return rem;
}

// |a| < 2^63 is guaranteed (threshold::S < #moduli <= threshold::M)
template <gemmul8::Backend backend, int IDX> __forceinline__ __device__ int32_t mod_middle(int64_t a) {
    constexpr int32_t p        = table::moduli<backend, IDX>;
    constexpr int64_t p_inv_64 = int64_t(18446744073709551615ULL / uint64_t(p)) + (18446744073709551615u % uint64_t(p) == uint64_t(p - 1)); // 2^64/p

    const int32_t rem = int32_t(a - p * __mul64hi(a, p_inv_64));
    return wrapping<backend, IDX>(rem);
}

// |a| can be >= 2^63 (threshold::M < #moduli <= threshold::L)
template <gemmul8::Backend backend, int IDX> __forceinline__ __device__ int32_t mod_large(fp_mant_exp<float> a) {
    constexpr int32_t p        = table::moduli<backend, IDX>;
    constexpr int32_t p_inv_32 = (int32_t)(4294967296ULL / uint64_t(p)); // 2^32/p

    const int32_t rem1 = a.mant - p * __mulhi(a.mant, p_inv_32);
    const int32_t rem2 = table::get_mod_pow2<backend, IDX>(a.exp);
    return mod_small<backend, IDX>(rem1 * rem2);
}
template <gemmul8::Backend backend, int IDX> __forceinline__ __device__ int32_t mod_large(fp_mant_exp<double> a) {
    constexpr int32_t p        = table::moduli<backend, IDX>;
    constexpr int64_t p_inv_64 = int64_t(18446744073709551615ULL / uint64_t(p)) + (18446744073709551615u % uint64_t(p) == uint64_t(p - 1)); // 2^64/p

    const int32_t rem1 = int32_t(a.mant - p * __mul64hi(a.mant, p_inv_64));
    const int32_t rem2 = table::get_mod_pow2<backend, IDX>(a.exp);
    return mod_small<backend, IDX>(rem1 * rem2);
}

// mod 256 for gemmul8::Backend::INT8 (static_cast<int8_t> performs wrapping)
template <> __forceinline__ __device__ int32_t mod_small<gemmul8::Backend::INT8, 0>(int32_t a) {
    return wrapping<gemmul8::Backend::INT8, 0>(a & 255);
}
template <> __forceinline__ __device__ int32_t mod_small_nowrap<gemmul8::Backend::INT8, 0>(int32_t a) {
    return a & 255;
}
template <> __forceinline__ __device__ int32_t mod_middle<gemmul8::Backend::INT8, 0>(int64_t a) {
    return wrapping<gemmul8::Backend::INT8, 0>(int32_t(a & 255LL));
}
template <> __forceinline__ __device__ int32_t mod_large<gemmul8::Backend::INT8, 0>(fp_mant_exp<float> a) {
    if (a.exp >= 8) return 0;
    return wrapping<gemmul8::Backend::INT8, 0>((a.mant << a.exp) & 255);
}
template <> __forceinline__ __device__ int32_t mod_large<gemmul8::Backend::INT8, 0>(fp_mant_exp<double> a) {
    if (a.exp >= 8) return 0;
    return wrapping<gemmul8::Backend::INT8, 0>(static_cast<int32_t>((a.mant << a.exp) & 255LL));
}

// mod 1024 for gemmul8::Backend::FP8
template <> __forceinline__ __device__ int32_t mod_small<gemmul8::Backend::FP8, 1>(int32_t a) {
    return wrapping<gemmul8::Backend::FP8, 1>(a & 1023);
}
template <> __forceinline__ __device__ int32_t mod_small_nowrap<gemmul8::Backend::FP8, 1>(int32_t a) {
    return a & 1023;
}
template <> __forceinline__ __device__ int32_t mod_middle<gemmul8::Backend::FP8, 1>(int64_t a) {
    return wrapping<gemmul8::Backend::FP8, 1>(int32_t(a & 1023LL));
}
template <> __forceinline__ __device__ int32_t mod_large<gemmul8::Backend::FP8, 1>(fp_mant_exp<float> a) {
    if (a.exp >= 10) return 0;
    return wrapping<gemmul8::Backend::FP8, 1>((a.mant << a.exp) & 1023);
}
template <> __forceinline__ __device__ int32_t mod_large<gemmul8::Backend::FP8, 1>(fp_mant_exp<double> a) {
    if (a.exp >= 10) return 0;
    return wrapping<gemmul8::Backend::FP8, 1>((a.mant << a.exp) & 1023LL);
}

template <gemmul8::Backend backend, int IDX> __forceinline__ __device__ int32_t mod_select(int32_t a) { return mod_small<backend, IDX>(a); }
template <gemmul8::Backend backend, int IDX> __forceinline__ __device__ int32_t mod_select(int64_t a) { return mod_middle<backend, IDX>(a); }
template <gemmul8::Backend backend, int IDX> __forceinline__ __device__ int32_t mod_select(fp_mant_exp<float> a) { return mod_large<backend, IDX>(a); }
template <gemmul8::Backend backend, int IDX> __forceinline__ __device__ int32_t mod_select(fp_mant_exp<double> a) { return mod_large<backend, IDX>(a); }

//------------------------------
// Transform 3 float integers
//     isSqr<p> => sqrt(p_i)*(C0 + C1) + C2     with int16, where C0=Ahi*Blo, C1=Alo*Bhi, C2=Alo*Blo
//    !isSqr<p> => 256*C0 + 16*(C2-C0-C1) + C1 with C0=Ahi*Bhi, C1=Alo*Blo, C2=(Ahi+Alo)*(Bhi+Blo)
// into 1 INT32
//------------------------------
template <int IDX> __forceinline__ __device__ int32_t mod_f32x3_2_i32(
    const float C0, const float C1, const float C2 //
) {
    const int32_t c0 = __float2int_rn(C0);
    const int32_t c1 = __float2int_rn(C1);
    const int32_t c2 = __float2int_rn(C2);

    const int32_t r0 = mod_small_nowrap<gemmul8::Backend::FP8, IDX>(c0);
    const int32_t r1 = mod_small_nowrap<gemmul8::Backend::FP8, IDX>(c1);
    const int32_t r2 = mod_small_nowrap<gemmul8::Backend::FP8, IDX>(c2);

    if constexpr (IDX < table::not_Karatsuba) {

        constexpr int32_t sqrt_p = table::sqrt_moduli<IDX>;
        const int32_t t          = sqrt_p * (r0 + r1) + r2;
        const int32_t rem        = mod_small<gemmul8::Backend::FP8, IDX>(t);
        return rem;

    } else {

        const int32_t t   = (r0 * 256) + ((r2 - r0 - r1) * 16) + r1;
        const int32_t rem = mod_middle<gemmul8::Backend::FP8, IDX>(t);
        return rem;
    }
}
template <int IDX> __forceinline__ __device__ int32_t mod_f32x3_2_i32_nowrap(
    const float C0, const float C1, const float C2 //
) {
    const int32_t c0 = __float2int_rn(C0);
    const int32_t c1 = __float2int_rn(C1);
    const int32_t c2 = __float2int_rn(C2);

    const int32_t r0 = mod_small_nowrap<gemmul8::Backend::FP8, IDX>(c0);
    const int32_t r1 = mod_small_nowrap<gemmul8::Backend::FP8, IDX>(c1);
    const int32_t r2 = mod_small_nowrap<gemmul8::Backend::FP8, IDX>(c2);

    if constexpr (IDX < table::not_Karatsuba) {

        constexpr int32_t sqrt_p = table::sqrt_moduli<IDX>;
        const int32_t rem        = sqrt_p * (r0 + r1) + r2;
        return rem;

    } else {

        const int32_t rem = (r0 * 256) + ((r2 - r0 - r1) * 16) + r1;
        return rem;
    }
}

//------------------------------
// Transform a 9-bit integer into 2*FP8
//------------------------------
// eror-free conversion a to {b.x, b.y, b.z} s.t. a = 16*b.x + b.y & b.z = b.x + b.y
__forceinline__ __device__ fp8x3_e4m3 make_fp8x3(int32_t a) {
    const uint32_t ua   = (uint32_t)a;
    const uint32_t sign = ua >> 31;                          // 0: a>=0, 1: a<0
    const uint32_t absu = (ua ^ (0u - sign)) + sign;         // |a| as uint32
    const uint32_t q    = (absu + 15u) >> 4;                 // ceil(|a|/16)
    const int32_t bx    = (sign ? -(int32_t)q : (int32_t)q); // sign(a)*ceil(|a|/16)
    const int32_t by    = a - 16 * bx;

    fp8x3_e4m3 b;
    b.x = __nv_fp8_e4m3(bx);      // without error
    b.y = __nv_fp8_e4m3(by);      // without error
    b.z = __nv_fp8_e4m3(bx + by); // without error

    return b;
}

// eror-free conversion a to {b.x, b.y} s.t. a = sqrt(p)*b.x + b.y for p = int^2
template <int IDX> __forceinline__ __device__ fp8x2_e4m3 make_fp8x2(int32_t a) {
    constexpr float inv_sqrtp = 1.0f / static_cast<float>(table::sqrt_moduli<IDX>); // 1/sqrt(p)
    const float a_f           = __int2float_rz(a);                                  // without error
    const float q             = a_f * inv_sqrtp;                                    // a/sqrt(p)
    const float bx_f          = rintf(q);                                           // round(a/sqrt(p))
    constexpr float m_sqrtp   = -static_cast<float>(table::sqrt_moduli<IDX>);       // -sqrt(p)
    const float by_f          = fmaf(m_sqrtp, bx_f, a_f);                           // a - sqrt(p) * round(a/sqrt(p)

    fp8x2_e4m3 b;
    b.x = __nv_fp8_e4m3(bx_f); // without error
    b.y = __nv_fp8_e4m3(by_f); // without error

    return b;
}

//------------------------------
// Transform a floating-point integer into mant*2^exp with exp >= 0
//------------------------------
__forceinline__ __device__ fp_mant_exp<float> make_fp_mant_exp(float a) {
    // fp_mant_exp<float> d;
    // d.exp  = max(Tilogb<float>(a) - 30, 0);
    // d.mant = __fp2int_rz<float>(Tscalbn<float>(a, -d.exp));
    // return d;

    using T = float;
    using I = int_t<T>;
    fp_mant_exp<T> d;

    // bit pattern
    const uint32_t bits = __float_as_uint(a);
    const uint32_t sign = bits >> 31;
    const uint32_t e    = (bits >> 23) & 0xFFu; // exponent bits
    const uint32_t frac = bits & 0x7FFFFFu;     // fraction bits

    // a == 0
    if (e == 0u) {
        d.exp  = 0;
        d.mant = 0;
        return d;
    }

    const int unbiased = int(e) - 127;      // floor(log2(|a|)) for normalized integer a
    const uint32_t sig = (1u << 23) | frac; // 24-bit significand
    const int out_exp  = max(unbiased - 30, 0);

    int32_t mant;
    if (unbiased > 30) {
        mant = (int32_t)(sig << 7);
    } else {
        const int sh       = unbiased - 23;
        const uint32_t mag = (sh >= 0) ? (sig << sh) : (sig >> (-sh));
        mant               = (int32_t)mag;
    }

    d.exp  = out_exp;
    d.mant = sign ? -mant : mant;
    return d;
}

__forceinline__ __device__ fp_mant_exp<double> make_fp_mant_exp(double a) {
    // fp_mant_exp<double> d;
    // d.exp  = max(Tilogb<double>(a) - 62, 0);
    // d.mant = __fp2int_rz<double>(Tscalbn<double>(a, -d.exp));
    // return d;

    fp_mant_exp<double> d;

    const uint64_t bits = (uint64_t)__double_as_longlong(a);
    const uint64_t sign = bits >> 63;
    const uint64_t e    = (bits >> 52) & 0x7FFull;      // exponent bits
    const uint64_t frac = bits & ((1ull << 52) - 1ull); // fraction bits

    if (e == 0ull) {
        d.exp  = 0;
        d.mant = 0;
        return d;
    }

    const int unbiased = int(e) - 1023;       // floor(log2(|a|)) for normalized integer a
    const uint64_t sig = (1ull << 52) | frac; // 53-bit significand
    const int out_exp  = max(unbiased - 62, 0);

    int64_t mant;
    if (unbiased > 62) {
        mant = (int64_t)(sig << 10);
    } else {
        const int sh       = unbiased - 52;
        const uint64_t mag = (sh >= 0) ? (sig << sh) : (sig >> (-sh));
        mant               = (int64_t)mag;
    }

    d.exp  = out_exp;
    d.mant = sign ? -mant : mant;
    return d;
}

__forceinline__ __device__ fp_mant_exp2<float> make_fp_mant_exp2(cuFloatComplex v) {
    fp_mant_exp2<float> d;
    d.x = make_fp_mant_exp(v.x);
    d.y = make_fp_mant_exp(v.y);
    return d;
}

__forceinline__ __device__ fp_mant_exp2<double> make_fp_mant_exp2(cuDoubleComplex v) {
    fp_mant_exp2<double> d;
    d.x = make_fp_mant_exp(v.x);
    d.y = make_fp_mant_exp(v.y);
    return d;
}

//------------------------------
// Launcher of mod functions
//------------------------------

//=====
// INT8
//=====
// launcher for V in {int32_t, int64_t, float, double}
template <int IDX, typename V> __forceinline__ __device__ void mod_launcher(
    int8_t *__restrict__ out,
    V v //
) {
    out[0] = static_cast<int8_t>(mod_select<gemmul8::Backend::INT8, IDX>(v));
}

// launcher for V in {int32_t, int64_t, float, double}
template <int IDX, typename V> __forceinline__ __device__ void mod_launcher(
    char4 *__restrict__ out,
    V v0, V v1, V v2, V v3 //
) {
    char4 rem;
    rem.x  = static_cast<int8_t>(mod_select<gemmul8::Backend::INT8, IDX>(v0));
    rem.y  = static_cast<int8_t>(mod_select<gemmul8::Backend::INT8, IDX>(v1));
    rem.z  = static_cast<int8_t>(mod_select<gemmul8::Backend::INT8, IDX>(v2));
    rem.w  = static_cast<int8_t>(mod_select<gemmul8::Backend::INT8, IDX>(v3));
    out[0] = rem;
}

// launcher for V in {int2, int64x2_t, float2, double2}
template <int IDX, typename V> __forceinline__ __device__ void mod_launcher(
    int8_t *__restrict__ out_r,
    int8_t *__restrict__ out_i,
    int8_t *__restrict__ out_ri,
    V v //
) {
    const int8_t vr = static_cast<int8_t>(mod_select<gemmul8::Backend::INT8, IDX>(v.x));
    const int8_t vi = static_cast<int8_t>(mod_select<gemmul8::Backend::INT8, IDX>(v.y));
    out_r[0]        = vr;
    out_i[0]        = vi;
    out_ri[0]       = static_cast<int8_t>(wrapping<gemmul8::Backend::INT8, IDX>(int(vr) + int(vi)));
}

// launcher for V in {int2, int64x2_t, float2, double2}
template <int IDX, typename V> __forceinline__ __device__ void mod_launcher(
    char4 *__restrict__ out_r,
    char4 *__restrict__ out_i,
    char4 *__restrict__ out_ri,
    V v0, V v1, V v2, V v3 //
) {
    char4 rem_r;
    rem_r.x  = static_cast<int8_t>(mod_select<gemmul8::Backend::INT8, IDX>(v0.x));
    rem_r.y  = static_cast<int8_t>(mod_select<gemmul8::Backend::INT8, IDX>(v1.x));
    rem_r.z  = static_cast<int8_t>(mod_select<gemmul8::Backend::INT8, IDX>(v2.x));
    rem_r.w  = static_cast<int8_t>(mod_select<gemmul8::Backend::INT8, IDX>(v3.x));
    out_r[0] = rem_r;

    char4 rem_i;
    rem_i.x  = static_cast<int8_t>(mod_select<gemmul8::Backend::INT8, IDX>(v0.y));
    rem_i.y  = static_cast<int8_t>(mod_select<gemmul8::Backend::INT8, IDX>(v1.y));
    rem_i.z  = static_cast<int8_t>(mod_select<gemmul8::Backend::INT8, IDX>(v2.y));
    rem_i.w  = static_cast<int8_t>(mod_select<gemmul8::Backend::INT8, IDX>(v3.y));
    out_i[0] = rem_i;

    char4 rem_ri;
    rem_ri.x  = static_cast<int8_t>(wrapping<gemmul8::Backend::INT8, IDX>(int(rem_r.x) + int(rem_i.x)));
    rem_ri.y  = static_cast<int8_t>(wrapping<gemmul8::Backend::INT8, IDX>(int(rem_r.y) + int(rem_i.y)));
    rem_ri.z  = static_cast<int8_t>(wrapping<gemmul8::Backend::INT8, IDX>(int(rem_r.z) + int(rem_i.z)));
    rem_ri.w  = static_cast<int8_t>(wrapping<gemmul8::Backend::INT8, IDX>(int(rem_r.w) + int(rem_i.w)));
    out_ri[0] = rem_ri;
}

//=====
// FP8
//=====
// launcher for V in {int32_t, int64_t, float, double}
template <int IDX, typename V> __forceinline__ __device__ void mod_launcher(
    __nv_fp8_e4m3 *__restrict__ out,
    size_t next,
    V v //
) {
    if constexpr (IDX < table::not_Karatsuba) {

        const fp8x2_e4m3 r = make_fp8x2<IDX>(mod_select<gemmul8::Backend::FP8, IDX>(v));

        out[0]    = r.x;
        out[next] = r.y;

    } else {

        const fp8x3_e4m3 rem = make_fp8x3(mod_select<gemmul8::Backend::FP8, IDX>(v));

        out[0]        = rem.x;
        out[next]     = rem.y;
        out[next * 2] = rem.z;
    }
}

// launcher for V in {int32_t, int64_t, float, double}
template <int IDX, typename V> __forceinline__ __device__ void mod_launcher(
    __nv_fp8x4_e4m3 *__restrict__ out,
    size_t next,
    V v0, V v1, V v2, V v3 //
) {
    if constexpr (IDX < table::not_Karatsuba) {

        const fp8x2_e4m3 rem0 = make_fp8x2<IDX>(mod_select<gemmul8::Backend::FP8, IDX>(v0));
        const fp8x2_e4m3 rem1 = make_fp8x2<IDX>(mod_select<gemmul8::Backend::FP8, IDX>(v1));
        const fp8x2_e4m3 rem2 = make_fp8x2<IDX>(mod_select<gemmul8::Backend::FP8, IDX>(v2));
        const fp8x2_e4m3 rem3 = make_fp8x2<IDX>(mod_select<gemmul8::Backend::FP8, IDX>(v3));

        out[0]    = concat(rem0.x, rem1.x, rem2.x, rem3.x);
        out[next] = concat(rem0.y, rem1.y, rem2.y, rem3.y);

    } else {

        const fp8x3_e4m3 r0 = make_fp8x3(mod_select<gemmul8::Backend::FP8, IDX>(v0));
        const fp8x3_e4m3 r1 = make_fp8x3(mod_select<gemmul8::Backend::FP8, IDX>(v1));
        const fp8x3_e4m3 r2 = make_fp8x3(mod_select<gemmul8::Backend::FP8, IDX>(v2));
        const fp8x3_e4m3 r3 = make_fp8x3(mod_select<gemmul8::Backend::FP8, IDX>(v3));

        out[0]        = concat(r0.x, r1.x, r2.x, r3.x);
        out[next]     = concat(r0.y, r1.y, r2.y, r3.y);
        out[next * 2] = concat(r0.z, r1.z, r2.z, r3.z);
    }
}

// V in {int2, int64x2_t, float2, double2}
template <int IDX, typename V> __forceinline__ __device__ void mod_launcher(
    __nv_fp8_e4m3 *__restrict__ out_r,
    __nv_fp8_e4m3 *__restrict__ out_i,
    __nv_fp8_e4m3 *__restrict__ out_ri,
    size_t next,
    V v //
) {
    if constexpr (IDX < table::not_Karatsuba) {

        const int32_t rem_r      = mod_select<gemmul8::Backend::FP8, IDX>(v.x);
        const fp8x2_e4m3 rem_rx2 = make_fp8x2<IDX>(rem_r);
        out_r[0]                 = rem_rx2.x;
        out_r[next]              = rem_rx2.y;

        const int32_t rem_i      = mod_select<gemmul8::Backend::FP8, IDX>(v.y);
        const fp8x2_e4m3 rem_ix2 = make_fp8x2<IDX>(rem_i);
        out_i[0]                 = rem_ix2.x;
        out_i[next]              = rem_ix2.y;

        const int32_t rem_ri      = wrapping<gemmul8::Backend::FP8, IDX>(rem_r + rem_i);
        const fp8x2_e4m3 rem_rix2 = make_fp8x2<IDX>(rem_ri);
        out_ri[0]                 = rem_rix2.x;
        out_ri[next]              = rem_rix2.y;

    } else {

        const int32_t rem_r      = mod_select<gemmul8::Backend::FP8, IDX>(v.x);
        const fp8x3_e4m3 rem_rx3 = make_fp8x3(rem_r);
        out_r[0]                 = rem_rx3.x;
        out_r[next]              = rem_rx3.y;
        out_r[next * 2]          = rem_rx3.z;

        const int32_t rem_i      = mod_select<gemmul8::Backend::FP8, IDX>(v.y);
        const fp8x3_e4m3 rem_ix3 = make_fp8x3(rem_i);
        out_i[0]                 = rem_ix3.x;
        out_i[next]              = rem_ix3.y;
        out_i[next * 2]          = rem_ix3.z;

        const int32_t rem_ri      = wrapping<gemmul8::Backend::FP8, IDX>(rem_r + rem_i);
        const fp8x3_e4m3 rem_rix3 = make_fp8x3(rem_ri);
        out_ri[0]                 = rem_rix3.x;
        out_ri[next]              = rem_rix3.y;
        out_ri[next * 2]          = rem_rix3.z;
    }
}

// V in {int2, int64x2_t, float2, double2}
template <int IDX, typename V> __forceinline__ __device__ void mod_launcher(
    __nv_fp8x4_e4m3 *__restrict__ out_r,
    __nv_fp8x4_e4m3 *__restrict__ out_i,
    __nv_fp8x4_e4m3 *__restrict__ out_ri,
    size_t next,
    V v0, V v1, V v2, V v3 //
) {
    if constexpr (IDX < table::not_Karatsuba) {

        const int32_t rem_r0 = mod_select<gemmul8::Backend::FP8, IDX>(v0.x);
        const int32_t rem_r1 = mod_select<gemmul8::Backend::FP8, IDX>(v1.x);
        const int32_t rem_r2 = mod_select<gemmul8::Backend::FP8, IDX>(v2.x);
        const int32_t rem_r3 = mod_select<gemmul8::Backend::FP8, IDX>(v3.x);

        const fp8x2_e4m3 rem_r0x2 = make_fp8x2<IDX>(rem_r0);
        const fp8x2_e4m3 rem_r1x2 = make_fp8x2<IDX>(rem_r1);
        const fp8x2_e4m3 rem_r2x2 = make_fp8x2<IDX>(rem_r2);
        const fp8x2_e4m3 rem_r3x2 = make_fp8x2<IDX>(rem_r3);

        out_r[0]    = concat(rem_r0x2.x, rem_r1x2.x, rem_r2x2.x, rem_r3x2.x);
        out_r[next] = concat(rem_r0x2.y, rem_r1x2.y, rem_r2x2.y, rem_r3x2.y);

        const int32_t rem_i0 = mod_select<gemmul8::Backend::FP8, IDX>(v0.y);
        const int32_t rem_i1 = mod_select<gemmul8::Backend::FP8, IDX>(v1.y);
        const int32_t rem_i2 = mod_select<gemmul8::Backend::FP8, IDX>(v2.y);
        const int32_t rem_i3 = mod_select<gemmul8::Backend::FP8, IDX>(v3.y);

        const fp8x2_e4m3 rem_i0x2 = make_fp8x2<IDX>(rem_i0);
        const fp8x2_e4m3 rem_i1x2 = make_fp8x2<IDX>(rem_i1);
        const fp8x2_e4m3 rem_i2x2 = make_fp8x2<IDX>(rem_i2);
        const fp8x2_e4m3 rem_i3x2 = make_fp8x2<IDX>(rem_i3);

        out_i[0]    = concat(rem_i0x2.x, rem_i1x2.x, rem_i2x2.x, rem_i3x2.x);
        out_i[next] = concat(rem_i0x2.y, rem_i1x2.y, rem_i2x2.y, rem_i3x2.y);

        const int32_t rem_ri0 = wrapping<gemmul8::Backend::FP8, IDX>(rem_r0 + rem_i0);
        const int32_t rem_ri1 = wrapping<gemmul8::Backend::FP8, IDX>(rem_r1 + rem_i1);
        const int32_t rem_ri2 = wrapping<gemmul8::Backend::FP8, IDX>(rem_r2 + rem_i2);
        const int32_t rem_ri3 = wrapping<gemmul8::Backend::FP8, IDX>(rem_r3 + rem_i3);

        const fp8x2_e4m3 rem_ri0x2 = make_fp8x2<IDX>(rem_ri0);
        const fp8x2_e4m3 rem_ri1x2 = make_fp8x2<IDX>(rem_ri1);
        const fp8x2_e4m3 rem_ri2x2 = make_fp8x2<IDX>(rem_ri2);
        const fp8x2_e4m3 rem_ri3x2 = make_fp8x2<IDX>(rem_ri3);

        out_ri[0]    = concat(rem_ri0x2.x, rem_ri1x2.x, rem_ri2x2.x, rem_ri3x2.x);
        out_ri[next] = concat(rem_ri0x2.y, rem_ri1x2.y, rem_ri2x2.y, rem_ri3x2.y);

    } else {

        const int32_t rem_r0 = mod_select<gemmul8::Backend::FP8, IDX>(v0.x);
        const int32_t rem_r1 = mod_select<gemmul8::Backend::FP8, IDX>(v1.x);
        const int32_t rem_r2 = mod_select<gemmul8::Backend::FP8, IDX>(v2.x);
        const int32_t rem_r3 = mod_select<gemmul8::Backend::FP8, IDX>(v3.x);

        const fp8x3_e4m3 rem_r0x3 = make_fp8x3(rem_r0);
        const fp8x3_e4m3 rem_r1x3 = make_fp8x3(rem_r1);
        const fp8x3_e4m3 rem_r2x3 = make_fp8x3(rem_r2);
        const fp8x3_e4m3 rem_r3x3 = make_fp8x3(rem_r3);

        out_r[0]        = concat(rem_r0x3.x, rem_r1x3.x, rem_r2x3.x, rem_r3x3.x);
        out_r[next]     = concat(rem_r0x3.y, rem_r1x3.y, rem_r2x3.y, rem_r3x3.y);
        out_r[next * 2] = concat(rem_r0x3.z, rem_r1x3.z, rem_r2x3.z, rem_r3x3.z);

        const int32_t rem_i0 = mod_select<gemmul8::Backend::FP8, IDX>(v0.y);
        const int32_t rem_i1 = mod_select<gemmul8::Backend::FP8, IDX>(v1.y);
        const int32_t rem_i2 = mod_select<gemmul8::Backend::FP8, IDX>(v2.y);
        const int32_t rem_i3 = mod_select<gemmul8::Backend::FP8, IDX>(v3.y);

        const fp8x3_e4m3 rem_i0x3 = make_fp8x3(rem_i0);
        const fp8x3_e4m3 rem_i1x3 = make_fp8x3(rem_i1);
        const fp8x3_e4m3 rem_i2x3 = make_fp8x3(rem_i2);
        const fp8x3_e4m3 rem_i3x3 = make_fp8x3(rem_i3);

        out_i[0]        = concat(rem_i0x3.x, rem_i1x3.x, rem_i2x3.x, rem_i3x3.x);
        out_i[next]     = concat(rem_i0x3.y, rem_i1x3.y, rem_i2x3.y, rem_i3x3.y);
        out_i[next * 2] = concat(rem_i0x3.z, rem_i1x3.z, rem_i2x3.z, rem_i3x3.z);

        const int32_t rem_ri0 = wrapping<gemmul8::Backend::FP8, IDX>(rem_r0 + rem_i0);
        const int32_t rem_ri1 = wrapping<gemmul8::Backend::FP8, IDX>(rem_r1 + rem_i1);
        const int32_t rem_ri2 = wrapping<gemmul8::Backend::FP8, IDX>(rem_r2 + rem_i2);
        const int32_t rem_ri3 = wrapping<gemmul8::Backend::FP8, IDX>(rem_r3 + rem_i3);

        const fp8x3_e4m3 rem_ri0x3 = make_fp8x3(rem_ri0);
        const fp8x3_e4m3 rem_ri1x3 = make_fp8x3(rem_ri1);
        const fp8x3_e4m3 rem_ri2x3 = make_fp8x3(rem_ri2);
        const fp8x3_e4m3 rem_ri3x3 = make_fp8x3(rem_ri3);

        out_ri[0]        = concat(rem_ri0x3.x, rem_ri1x3.x, rem_ri2x3.x, rem_ri3x3.x);
        out_ri[next]     = concat(rem_ri0x3.y, rem_ri1x3.y, rem_ri2x3.y, rem_ri3x3.y);
        out_ri[next * 2] = concat(rem_ri0x3.z, rem_ri1x3.z, rem_ri2x3.z, rem_ri3x3.z);
    }
}

//------------------------------
// Interface
//------------------------------
#define GEMMUL8_FOR_EACH_IDX_20(M) \
    M(0)                           \
    M(1)                           \
    M(2)                           \
    M(3)                           \
    M(4)                           \
    M(5)                           \
    M(6)                           \
    M(7)                           \
    M(8)                           \
    M(9)                           \
    M(10)                          \
    M(11)                          \
    M(12)                          \
    M(13)                          \
    M(14)                          \
    M(15)                          \
    M(16)                          \
    M(17)                          \
    M(18)                          \
    M(19)

#define GEMMUL8_INT8_RUN_SCALAR_STEP(I) \
    if constexpr (num_moduli > (I)) {   \
        mod_launcher<(I), V>(out, v);   \
        out += inc;                     \
    }

#define GEMMUL8_INT8_RUN_VEC4_STEP(I)              \
    if constexpr (num_moduli > (I)) {              \
        mod_launcher<(I), V>(out, v0, v1, v2, v3); \
        out += inc;                                \
    }

#define GEMMUL8_INT8_RUN_CPLX_SCALAR_STEP(I)       \
    if constexpr (num_moduli > (I)) {              \
        mod_launcher<(I), V>(out0, out1, out2, v); \
        out0 += inc;                               \
        out1 += inc;                               \
        out2 += inc;                               \
    }

#define GEMMUL8_INT8_RUN_CPLX_VEC4_STEP(I)                      \
    if constexpr (num_moduli > (I)) {                           \
        mod_launcher<(I), V>(out0, out1, out2, v0, v1, v2, v3); \
        out0 += inc;                                            \
        out1 += inc;                                            \
        out2 += inc;                                            \
    }

#define GEMMUL8_FP8_RUN_SCALAR_STEP(I)                       \
    if constexpr (num_moduli > (I)) {                        \
        mod_launcher<(I), V>(out, inc, v);                   \
        out += (((I) < table::not_Karatsuba) ? 2 : 3) * inc; \
    }

#define GEMMUL8_FP8_RUN_VEC4_STEP(I)                         \
    if constexpr (num_moduli > (I)) {                        \
        mod_launcher<(I), V>(out, inc, v0, v1, v2, v3);      \
        out += (((I) < table::not_Karatsuba) ? 2 : 3) * inc; \
    }

#define GEMMUL8_FP8_RUN_CPLX_SCALAR_STEP(I)                               \
    if constexpr (num_moduli > (I)) {                                     \
        mod_launcher<(I), V>(out0, out1, out2, inc, v);                   \
        const size_t step = (((I) < table::not_Karatsuba) ? 2 : 3) * inc; \
        out0 += step;                                                     \
        out1 += step;                                                     \
        out2 += step;                                                     \
    }

#define GEMMUL8_FP8_RUN_CPLX_VEC4_STEP(I)                                 \
    if constexpr (num_moduli > (I)) {                                     \
        mod_launcher<(I), V>(out0, out1, out2, inc, v0, v1, v2, v3);      \
        const size_t step = (((I) < table::not_Karatsuba) ? 2 : 3) * inc; \
        out0 += step;                                                     \
        out1 += step;                                                     \
        out2 += step;                                                     \
    }

// interface for general num_moduli
template <int num_moduli, typename V> struct ModUnroll {

    //=====
    // INT8
    //=====
    template <int IDX = 0> __inline__ __device__ static void run(
        int8_t *__restrict__ out,
        size_t inc,
        V v //
    ) {
        if constexpr (isComplex<V>) return;
        GEMMUL8_FOR_EACH_IDX_20(GEMMUL8_INT8_RUN_SCALAR_STEP);
    }

    template <int IDX = 0> __inline__ __device__ static void run(
        char4 *__restrict__ out,
        size_t inc,
        V v0, V v1, V v2, V v3 //
    ) {
        if constexpr (isComplex<V>) return;
        GEMMUL8_FOR_EACH_IDX_20(GEMMUL8_INT8_RUN_VEC4_STEP);
    }

    template <int IDX = 0> __inline__ __device__ static void run(
        int8_t *__restrict__ out0, int8_t *__restrict__ out1, int8_t *__restrict__ out2,
        size_t inc,
        V v //
    ) {
        if constexpr (!isComplex<V>) return;
        GEMMUL8_FOR_EACH_IDX_20(GEMMUL8_INT8_RUN_CPLX_SCALAR_STEP);
    }

    template <int IDX = 0> __inline__ __device__ static void run(
        char4 *__restrict__ out0, char4 *__restrict__ out1, char4 *__restrict__ out2,
        size_t inc,
        V v0, V v1, V v2, V v3 //
    ) {
        if constexpr (!isComplex<V>) return;
        GEMMUL8_FOR_EACH_IDX_20(GEMMUL8_INT8_RUN_CPLX_VEC4_STEP);
    }

    //=====
    // FP8
    //=====
    template <int IDX = 0> __inline__ __device__ static void run(
        __nv_fp8_e4m3 *__restrict__ out,
        size_t inc,
        V v //
    ) {
        if constexpr (isComplex<V>) return;
        GEMMUL8_FOR_EACH_IDX_20(GEMMUL8_FP8_RUN_SCALAR_STEP);
    }

    template <int IDX = 0> __inline__ __device__ static void run(
        __nv_fp8x4_e4m3 *__restrict__ out,
        size_t inc,
        V v0, V v1, V v2, V v3 //
    ) {
        if constexpr (isComplex<V>) return;
        GEMMUL8_FOR_EACH_IDX_20(GEMMUL8_FP8_RUN_VEC4_STEP);
    }

    template <int IDX = 0> __inline__ __device__ static void run(
        __nv_fp8_e4m3 *__restrict__ out0, __nv_fp8_e4m3 *__restrict__ out1, __nv_fp8_e4m3 *__restrict__ out2,
        size_t inc,
        V v //
    ) {
        if constexpr (!isComplex<V>) return;
        GEMMUL8_FOR_EACH_IDX_20(GEMMUL8_FP8_RUN_CPLX_SCALAR_STEP);
    }

    template <int IDX = 0> __inline__ __device__ static void run(
        __nv_fp8x4_e4m3 *__restrict__ out0, __nv_fp8x4_e4m3 *__restrict__ out1, __nv_fp8x4_e4m3 *__restrict__ out2,
        size_t inc,
        V v0, V v1, V v2, V v3 //
    ) {
        if constexpr (!isComplex<V>) return;
        GEMMUL8_FOR_EACH_IDX_20(GEMMUL8_FP8_RUN_CPLX_VEC4_STEP);
    }
};

#undef GEMMUL8_FOR_EACH_IDX_20
#undef GEMMUL8_INT8_RUN_SCALAR_STEP
#undef GEMMUL8_INT8_RUN_VEC4_STEP
#undef GEMMUL8_INT8_RUN_CPLX_SCALAR_STEP
#undef GEMMUL8_INT8_RUN_CPLX_VEC4_STEP
#undef GEMMUL8_FP8_RUN_SCALAR_STEP
#undef GEMMUL8_FP8_RUN_VEC4_STEP
#undef GEMMUL8_FP8_RUN_CPLX_SCALAR_STEP
#undef GEMMUL8_FP8_RUN_CPLX_VEC4_STEP

// interface for num_moduli > threshold::M & V = float
template <int num_moduli> struct ModUnroll<num_moduli, float> {

    //=====
    // INT8
    //=====
    __forceinline__ __device__ static void run(int8_t *__restrict__ out, size_t inc, float v) {
        fp_mant_exp<float> d = make_fp_mant_exp(v);
        ModUnroll<num_moduli, fp_mant_exp<float>>::run(out, inc, d);
    }

    __forceinline__ __device__ static void run(char4 *__restrict__ out, size_t inc, float v0, float v1, float v2, float v3) {
        fp_mant_exp<float> d0 = make_fp_mant_exp(v0);
        fp_mant_exp<float> d1 = make_fp_mant_exp(v1);
        fp_mant_exp<float> d2 = make_fp_mant_exp(v2);
        fp_mant_exp<float> d3 = make_fp_mant_exp(v3);
        ModUnroll<num_moduli, fp_mant_exp<float>>::run(out, inc, d0, d1, d2, d3);
    }

    //=====
    // FP8
    //=====
    __forceinline__ __device__ static void run(__nv_fp8_e4m3 *__restrict__ out, size_t inc, float v) {
        fp_mant_exp<float> d = make_fp_mant_exp(v);
        ModUnroll<num_moduli, fp_mant_exp<float>>::run(out, inc, d);
    }

    __forceinline__ __device__ static void run(__nv_fp8x4_e4m3 *__restrict__ out, size_t inc, float v0, float v1, float v2, float v3) {
        fp_mant_exp<float> d0 = make_fp_mant_exp(v0);
        fp_mant_exp<float> d1 = make_fp_mant_exp(v1);
        fp_mant_exp<float> d2 = make_fp_mant_exp(v2);
        fp_mant_exp<float> d3 = make_fp_mant_exp(v3);
        ModUnroll<num_moduli, fp_mant_exp<float>>::run(out, inc, d0, d1, d2, d3);
    }
};

// interface for num_moduli > threshold::M & V = double
template <int num_moduli> struct ModUnroll<num_moduli, double> {

    //=====
    // INT8
    //=====
    __forceinline__ __device__ static void run(int8_t *__restrict__ out, size_t inc, double v) {
        fp_mant_exp<double> d = make_fp_mant_exp(v);
        ModUnroll<num_moduli, fp_mant_exp<double>>::run(out, inc, d);
    }

    __forceinline__ __device__ static void run(char4 *__restrict__ out, size_t inc, double v0, double v1, double v2, double v3) {
        fp_mant_exp<double> d0 = make_fp_mant_exp(v0);
        fp_mant_exp<double> d1 = make_fp_mant_exp(v1);
        fp_mant_exp<double> d2 = make_fp_mant_exp(v2);
        fp_mant_exp<double> d3 = make_fp_mant_exp(v3);
        ModUnroll<num_moduli, fp_mant_exp<double>>::run(out, inc, d0, d1, d2, d3);
    }

    //=====
    // FP8
    //=====
    __forceinline__ __device__ static void run(__nv_fp8_e4m3 *__restrict__ out, size_t inc, double v) {
        fp_mant_exp<double> d = make_fp_mant_exp(v);
        ModUnroll<num_moduli, fp_mant_exp<double>>::run(out, inc, d);
    }

    __forceinline__ __device__ static void run(__nv_fp8x4_e4m3 *__restrict__ out, size_t inc, double v0, double v1, double v2, double v3) {
        fp_mant_exp<double> d0 = make_fp_mant_exp(v0);
        fp_mant_exp<double> d1 = make_fp_mant_exp(v1);
        fp_mant_exp<double> d2 = make_fp_mant_exp(v2);
        fp_mant_exp<double> d3 = make_fp_mant_exp(v3);
        ModUnroll<num_moduli, fp_mant_exp<double>>::run(out, inc, d0, d1, d2, d3);
    }
};

// interface for num_moduli > threshold::M & V = cuFloatComplex
template <int num_moduli> struct ModUnroll<num_moduli, cuFloatComplex> {

    //=====
    // INT8
    //=====
    __forceinline__ __device__ static void run(int8_t *__restrict__ out0, int8_t *__restrict__ out1, int8_t *__restrict__ out2, size_t inc, cuFloatComplex v) {
        fp_mant_exp2<float> d = make_fp_mant_exp2(v);
        ModUnroll<num_moduli, fp_mant_exp2<float>>::run(out0, out1, out2, inc, d);
    }

    __forceinline__ __device__ static void run(char4 *__restrict__ out0, char4 *__restrict__ out1, char4 *__restrict__ out2,
                                               size_t inc, cuFloatComplex v0, cuFloatComplex v1, cuFloatComplex v2, cuFloatComplex v3) {
        fp_mant_exp2<float> d0 = make_fp_mant_exp2(v0);
        fp_mant_exp2<float> d1 = make_fp_mant_exp2(v1);
        fp_mant_exp2<float> d2 = make_fp_mant_exp2(v2);
        fp_mant_exp2<float> d3 = make_fp_mant_exp2(v3);
        ModUnroll<num_moduli, fp_mant_exp2<float>>::run(out0, out1, out2, inc, d0, d1, d2, d3);
    }

    //=====
    // FP8
    //=====
    __forceinline__ __device__ static void run(__nv_fp8_e4m3 *__restrict__ out0, __nv_fp8_e4m3 *__restrict__ out1, __nv_fp8_e4m3 *__restrict__ out2, size_t inc, cuFloatComplex v) {
        fp_mant_exp2<float> d = make_fp_mant_exp2(v);
        ModUnroll<num_moduli, fp_mant_exp2<float>>::run(out0, out1, out2, inc, d);
    }

    __forceinline__ __device__ static void run(__nv_fp8x4_e4m3 *__restrict__ out0, __nv_fp8x4_e4m3 *__restrict__ out1, __nv_fp8x4_e4m3 *__restrict__ out2,
                                               size_t inc, cuFloatComplex v0, cuFloatComplex v1, cuFloatComplex v2, cuFloatComplex v3) {
        fp_mant_exp2<float> d0 = make_fp_mant_exp2(v0);
        fp_mant_exp2<float> d1 = make_fp_mant_exp2(v1);
        fp_mant_exp2<float> d2 = make_fp_mant_exp2(v2);
        fp_mant_exp2<float> d3 = make_fp_mant_exp2(v3);
        ModUnroll<num_moduli, fp_mant_exp2<float>>::run(out0, out1, out2, inc, d0, d1, d2, d3);
    }
};

// interface for num_moduli > threshold::M & V = cuDoubleComplex
template <int num_moduli> struct ModUnroll<num_moduli, cuDoubleComplex> {

    //=====
    // INT8
    //=====
    __forceinline__ __device__ static void run(int8_t *__restrict__ out0, int8_t *__restrict__ out1, int8_t *__restrict__ out2,
                                               size_t inc, cuDoubleComplex v) {
        fp_mant_exp2<double> d = make_fp_mant_exp2(v);
        ModUnroll<num_moduli, fp_mant_exp2<double>>::run(out0, out1, out2, inc, d);
    }

    __forceinline__ __device__ static void run(char4 *__restrict__ out0, char4 *__restrict__ out1, char4 *__restrict__ out2,
                                               size_t inc, cuDoubleComplex v0, cuDoubleComplex v1, cuDoubleComplex v2, cuDoubleComplex v3) {
        fp_mant_exp2<double> d0 = make_fp_mant_exp2(v0);
        fp_mant_exp2<double> d1 = make_fp_mant_exp2(v1);
        fp_mant_exp2<double> d2 = make_fp_mant_exp2(v2);
        fp_mant_exp2<double> d3 = make_fp_mant_exp2(v3);
        ModUnroll<num_moduli, fp_mant_exp2<double>>::run(out0, out1, out2, inc, d0, d1, d2, d3);
    }

    //=====
    // FP8
    //=====
    __forceinline__ __device__ static void run(__nv_fp8_e4m3 *__restrict__ out0, __nv_fp8_e4m3 *__restrict__ out1, __nv_fp8_e4m3 *__restrict__ out2,
                                               size_t inc, cuDoubleComplex v) {
        fp_mant_exp2<double> d = make_fp_mant_exp2(v);
        ModUnroll<num_moduli, fp_mant_exp2<double>>::run(out0, out1, out2, inc, d);
    }

    __forceinline__ __device__ static void run(__nv_fp8x4_e4m3 *__restrict__ out0, __nv_fp8x4_e4m3 *__restrict__ out1, __nv_fp8x4_e4m3 *__restrict__ out2,
                                               size_t inc, cuDoubleComplex v0, cuDoubleComplex v1, cuDoubleComplex v2, cuDoubleComplex v3) {
        fp_mant_exp2<double> d0 = make_fp_mant_exp2(v0);
        fp_mant_exp2<double> d1 = make_fp_mant_exp2(v1);
        fp_mant_exp2<double> d2 = make_fp_mant_exp2(v2);
        fp_mant_exp2<double> d3 = make_fp_mant_exp2(v3);
        ModUnroll<num_moduli, fp_mant_exp2<double>>::run(out0, out1, out2, inc, d0, d1, d2, d3);
    }
};
