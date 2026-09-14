#pragma once
#include "common.hpp"
#include "table.hpp"
#include "fp8_limb_selection.hpp"
#include "fp8_residue_lut.hpp"

namespace gemmul8::common {

namespace make_f8 {

//==========
// Constant unsigned division by a small compile-time constant
//==========
template <uint32_t D>
inline constexpr uint32_t ceildiv_u32 = uint32_t(((uint64_t(1) << 32) + D - 1) / D);

template <uint32_t D>
__device__ __forceinline__ uint32_t div_small_u32(const uint32_t x) {
#if defined(__CUDACC__) && !defined(__HIPCC__)
    uint32_t q;
    asm("mul.hi.u32 %0, %1, %2;"
        : "=r"(q)
        : "r"(x), "n"(ceildiv_u32<D>));
    return q;
#else
    return x / D;
#endif
}

//==========
// Small integer helpers
//==========
__device__ __forceinline__ uint32_t is_e4m3_hole_mag(const uint32_t x) {
    return uint32_t((x & 0x11u) == 0x11u);
}

template <uint32_t MAX_ABS>
__device__ __forceinline__ uint32_t valid_mag_pair(const uint32_t x, const uint32_t y) {
    static_assert(MAX_ABS <= 32u);
    const uint32_t m    = (x > y) ? x : y;
    const uint32_t hole = is_e4m3_hole_mag(x) | is_e4m3_hole_mag(y);
    return uint32_t(m <= MAX_ABS) & (hole ^ 1u);
}

template <uint32_t MAX_ABS>
__device__ __forceinline__ uint32_t valid_i32_triple(const int32_t x, const int32_t y, const int32_t z) {
    static_assert(MAX_ABS <= 32u);
    const uint32_t ax   = uint32_t(abs(x));
    const uint32_t ay   = uint32_t(abs(y));
    const uint32_t az   = uint32_t(abs(z));
    const uint32_t mxy  = (ax > ay) ? ax : ay;
    const uint32_t m    = (mxy > az) ? mxy : az;
    const uint32_t hole = is_e4m3_hole_mag(ax) | is_e4m3_hole_mag(ay) | is_e4m3_hole_mag(az);
    return uint32_t(m <= MAX_ABS) & (hole ^ 1u);
}

__device__ __forceinline__ int32_t select_i32(const uint32_t pred, const int32_t x, const int32_t y) {
    return pred ? x : y;
}

//==========
// FP8 conversion
//==========
#if defined(__CUDACC__) && !defined(__HIPCC__)

__device__ __forceinline__ uint16_t cvt_e4m3x2(const float lo, const float hi) {
    uint16_t r;
    asm("cvt.rn.satfinite.e4m3x2.f32 %0, %2, %1;"
        : "=h"(r)
        : "f"(lo), "f"(hi));
    return r;
}

__device__ __forceinline__ uint16_t cvt_e4m3x1_lo(const float x) {
    uint16_t r;
    asm("cvt.rn.satfinite.e4m3x2.f32 "
        "%0, 0f00000000, %1;"
        : "=h"(r)
        : "f"(x));
    return r;
}

__device__ __forceinline__ fp8x2_e4m3 make_fp8x2_raw(const uint16_t raw) {
    fp8x2_e4m3 out;
    out.x.__x = static_cast<__nv_fp8_storage_t>(raw);
    out.y.__x = static_cast<__nv_fp8_storage_t>(raw >> 8);
    return out;
}

__device__ __forceinline__ fp8x2_e4m3 convert_i32x2(const int32_t x0, const int32_t x1, const int32_t sign_mask) {
    const float f0    = __int2float_rn(x0);
    const float f1    = __int2float_rn(x1);
    uint16_t raw      = cvt_e4m3x2(f0, f1);
    const uint16_t sm = uint16_t(uint32_t(sign_mask) & 0x8080u);
    raw ^= sm;
    return make_fp8x2_raw(raw);
}

__device__ __forceinline__ fp8x3_e4m3 convert_i32x3(const int32_t x0, const int32_t x1, const int32_t x2, const int32_t sign_mask) {
    const float f0    = __int2float_rn(x0);
    const float f1    = __int2float_rn(x1);
    const float f2    = __int2float_rn(x2);
    uint16_t raw01    = cvt_e4m3x2(f0, f1);
    uint16_t raw2     = cvt_e4m3x1_lo(f2);
    const uint16_t sm = uint16_t(uint32_t(sign_mask) & 0x8080u);
    raw01 ^= sm;
    raw2 ^= sm;

    fp8x3_e4m3 out;
    out.x.__x = static_cast<__nv_fp8_storage_t>(raw01);
    out.y.__x = static_cast<__nv_fp8_storage_t>(raw01 >> 8);
    out.z.__x = static_cast<__nv_fp8_storage_t>(raw2);
    return out;
}

__device__ __forceinline__ fp8x2_e4m3 convert_f32x2(const float x0, const float x1) {
    return make_fp8x2_raw(cvt_e4m3x2(x0, x1));
}

#else

__device__ __forceinline__ fp8x2_e4m3 convert_i32x2(const int32_t x0, const int32_t x1, const int32_t sign_mask) {
    const int32_t s = sign_mask | 1;
    fp8x2_e4m3 out;
    out.x = __nv_fp8_e4m3(s * x0);
    out.y = __nv_fp8_e4m3(s * x1);
    return out;
}

__device__ __forceinline__ fp8x3_e4m3 convert_i32x3(const int32_t x0, const int32_t x1, const int32_t x2, const int32_t sign_mask) {
    const int32_t s = sign_mask | 1;
    fp8x3_e4m3 out;
    out.x = __nv_fp8_e4m3(s * x0);
    out.y = __nv_fp8_e4m3(s * x1);
    out.z = __nv_fp8_e4m3(s * x2);
    return out;
}

__device__ __forceinline__ fp8x2_e4m3 convert_f32x2(const float x0, const float x1) {
    fp8x2_e4m3 out;
    out.x = __nv_fp8_e4m3(x0);
    out.y = __nv_fp8_e4m3(x1);
    return out;
}

#endif

//==========
// Maximum limb magnitude of the selected Karatsuba scheme
//==========
template <int32_t P, table::KaratsubaType TYPE>
inline constexpr uint32_t kara_max_abs = [] {
    if constexpr (TYPE == table::KaratsubaType::BASE49) {
        if constexpr (P == 997) return 24u;
        else if constexpr (P == 937) return 28u;
        else return 32u;
    } else if constexpr (TYPE == table::KaratsubaType::BASE33_SUM) {
        return 30u;
    } else if constexpr (TYPE == table::KaratsubaType::BASE32_DIFF) {
        return 32u;
    } else if constexpr (TYPE == table::KaratsubaType::BASE32_SUM) {
        if constexpr (P == 709) return 20u;
        else if constexpr (P == 769 || P == 761) return 22u;
        else if constexpr (P == 911 || P == 907 || P == 905 || P == 901 || P == 757) return 26u;
        else if constexpr (P == 853) return 30u;
        else return 32u;
    } else {
        return 32u;
    }
}();

//==========
// Square modulus
// a == BASE*q + r (mod p), p = BASE^2, BASE > 32
//==========
template <int32_t BASE>
__device__ __forceinline__ void decompose_square_large(const uint32_t u, int32_t &q, int32_t &r) {
    static_assert(32 < BASE && BASE <= 49);

    q = int32_t(div_small_u32<uint32_t(BASE)>(u + 16u));
    r = int32_t(u) - BASE * q;

    const int32_t fix_r = int32_t(uint32_t(r > 16) & uint32_t(r & 1));
    q += fix_r;
    r -= BASE * fix_r;

    const int32_t fix_q = int32_t(is_e4m3_hole_mag(uint32_t(q)));
    q -= BASE * fix_q;
}

//==========
// BASE32_SUM / BASE32_DIFF
// a == 32*a0 + a1 (mod p)
// SUM : a2 = a0 + a1
// DIFF: a2 = a1 - a0
//==========
template <int32_t P, table::KaratsubaType TYPE>
__device__ __forceinline__ void decompose_base32(const uint32_t u, int32_t &a0, int32_t &a1, int32_t &a2) {
    static_assert(TYPE == table::KaratsubaType::BASE32_SUM || TYPE == table::KaratsubaType::BASE32_DIFF);
    constexpr bool SUM         = TYPE == table::KaratsubaType::BASE32_SUM;
    constexpr uint32_t MAX_ABS = kara_max_abs<P, TYPE>;
    constexpr int32_t H        = (P + 16) >> 5;
    constexpr int32_t T        = P - (H << 5);
    const int32_t q            = int32_t(u >> 5);
    const int32_t r            = int32_t(u & 31u);

    // u = 32*q + r
    const int32_t b0 = q;
    const int32_t b1 = r;
    const int32_t b2 = SUM ? (r + q) : (r - q);
    uint32_t b2_mag;
    if constexpr (SUM) {
        b2_mag = uint32_t(b2);
    } else {
        b2_mag = uint32_t(abs(b2));
    }
    const uint32_t ok_B = valid_mag_pair<MAX_ABS>(uint32_t(r), b2_mag);

    // u = 32*(q+1) + (r-32)
    const int32_t c0      = q + 1;
    const int32_t c1      = r - 32;
    const int32_t c2      = SUM ? (r + q - 31) : (r - q - 33);
    const uint32_t c1_mag = uint32_t(32 - r);
    uint32_t c2_mag;
    if constexpr (SUM) {
        c2_mag = uint32_t(abs(c2));
    } else {
        c2_mag = uint32_t(33 + q - r);
    }
    const uint32_t ok_C  = valid_mag_pair<MAX_ABS>(c1_mag, c2_mag);
    const uint32_t use_C = (ok_B ^ 1u) & ok_C;

    a0 = select_i32(use_C, c0, b0);
    a1 = select_i32(use_C, c1, b1);

    const uint32_t need_A = (ok_B | ok_C) ^ 1u;

    // P = 32*H + T.
    constexpr int32_t A_CARRY = (P == 853 || P == 793 || P == 797) ? 1 : 0;
    const int32_t aa0         = q - H + A_CARRY;
    const int32_t aa1         = r - T - 32 * A_CARRY;

    if constexpr (P != 797) {
        a0 = select_i32(need_A, aa0, a0);
        a1 = select_i32(need_A, aa1, a1);
    } else {
        const int32_t aa2    = SUM ? (aa0 + aa1) : (aa1 - aa0);
        const uint32_t ok_A  = valid_i32_triple<MAX_ABS>(aa0, aa1, aa2);
        const uint32_t use_A = need_A & ok_A;
        const uint32_t use_D = need_A & (ok_A ^ 1u);
        const int32_t d0     = q + H;
        const int32_t d1     = r + T;

        a0 = select_i32(use_A, aa0, a0);
        a1 = select_i32(use_A, aa1, a1);
        a0 = select_i32(use_D, d0, a0);
        a1 = select_i32(use_D, d1, a1);
    }

    a2 = SUM ? (a0 + a1) : (a1 - a0);
}

//==========
// BASE33_SUM
// a == 33*a0 + a1 (mod p) == 32*a0 + a2 (mod p)
// a2 = a0 + a1
//==========
template <int32_t P>
__device__ __forceinline__ void decompose_base33_sum(const uint32_t u, int32_t &a0, int32_t &a1, int32_t &a2) {
    static_assert(P == 1033);

    constexpr uint32_t MAX_ABS = 30u;
    const int32_t q            = int32_t(u >> 5);
    const int32_t r            = int32_t(u & 31u);

    const int32_t b0    = q;
    const int32_t b2    = r;
    const int32_t b1    = r - q;
    const uint32_t ok_B = valid_mag_pair<MAX_ABS>(uint32_t(abs(b1)), uint32_t(r));

    const int32_t c0     = q + 1;
    const int32_t c2     = r - 32;
    const uint32_t ok_C  = valid_mag_pair<MAX_ABS>(uint32_t(33 + q - r), uint32_t(32 - r));
    const uint32_t use_C = (ok_B ^ 1u) & ok_C;

    a0 = select_i32(use_C, c0, b0);
    a2 = select_i32(use_C, c2, b2);

    const uint32_t need_A = (ok_B | ok_C) ^ 1u;

    const int32_t aa0 = q - 31;
    const int32_t aa2 = r - 41;

    a0 = select_i32(need_A, aa0, a0);
    a2 = select_i32(need_A, aa2, a2);

    a1 = a2 - a0;
}

//==========
// BASE31_DIFF
// a == 31*a0 + a1 (mod p) == 32*a0 + a2 (mod p)
// a2 = a1-a0
//==========
template <int32_t P>
__device__ __forceinline__ void decompose_base31_diff(const uint32_t u, int32_t &a0, int32_t &a1, int32_t &a2) {

    constexpr int32_t H = (P + 16) / 32;
    constexpr int32_t T = P - 32 * H;
    const int32_t q     = int32_t(u >> 5);
    const int32_t r     = int32_t(u & 31u);

    // B
    const int32_t b0    = q;
    const int32_t b2    = r;
    const int32_t b1    = b0 + b2;
    const uint32_t ok_B = valid_i32_triple<32u>(b0, b1, b2);

    // C
    const int32_t c0     = q + 1;
    const int32_t c2     = r - 32;
    const int32_t c1     = c0 + c2;
    const uint32_t ok_C  = valid_i32_triple<32u>(c0, c1, c2);
    const uint32_t use_C = (ok_B ^ 1u) & ok_C;

    a0                    = select_i32(use_C, c0, b0);
    a2                    = select_i32(use_C, c2, b2);
    const uint32_t need_A = (ok_B | ok_C) ^ 1u;

    // A: u-p
    const int32_t aa0    = q - H;
    const int32_t aa2    = r - T;
    const int32_t aa1    = aa0 + aa2;
    const uint32_t ok_A  = valid_i32_triple<32u>(aa0, aa1, aa2);
    const uint32_t use_A = need_A & ok_A;

    // D: u+p
    const int32_t d0     = q + H;
    const int32_t d2     = r + T;
    const uint32_t use_D = need_A & (ok_A ^ 1u);

    a0 = select_i32(use_A, aa0, a0);
    a2 = select_i32(use_A, aa2, a2);
    a0 = select_i32(use_D, d0, a0);
    a2 = select_i32(use_D, d2, a2);

    a1 = a0 + a2;
}

//==========
// BASE49
// a == 49*a0 + a1 (mod p) == 48*a0 + a2 (mod p)
// a2 = a0+a1
//==========
template <int32_t P>
__device__ __forceinline__ void decompose_base49(const uint32_t u, int32_t &a0, int32_t &a1, int32_t &a2) {

    constexpr uint32_t MAX_ABS = kara_max_abs<P, table::KaratsubaType::BASE49>;
    constexpr int32_t H        = (P + 24) / 48;
    constexpr int32_t T        = P - 48 * H;

    const int32_t q = int32_t(div_small_u32<48u>(u));
    const int32_t r = int32_t(u) - 48 * q;

    const int32_t b0    = q;
    const int32_t b2    = r;
    const int32_t b1    = r - q;
    const uint32_t ok_B = valid_mag_pair<MAX_ABS>(uint32_t(abs(b1)), uint32_t(r));

    const int32_t c0     = q + 1;
    const int32_t c2     = r - 48;
    const uint32_t ok_C  = valid_mag_pair<MAX_ABS>(uint32_t(49 + q - r), uint32_t(48 - r));
    const uint32_t use_C = (ok_B ^ 1u) & ok_C;

    a0 = select_i32(use_C, c0, b0);
    a2 = select_i32(use_C, c2, b2);

    const uint32_t need_A = (ok_B | ok_C) ^ 1u;

    const int32_t aa0    = q + 1 - H;
    const int32_t aa2    = r - T - 48;
    const int32_t aa1    = aa2 - aa0;
    const uint32_t ok_A  = valid_i32_triple<MAX_ABS>(aa0, aa1, aa2);
    const uint32_t use_A = need_A & ok_A;

    const uint32_t use_D = need_A & (ok_A ^ 1u);
    const int32_t d0     = q + H;
    const int32_t d2     = r + T;

    a0 = select_i32(use_A, aa0, a0);
    a2 = select_i32(use_A, aa2, a2);
    a0 = select_i32(use_D, d0, a0);
    a2 = select_i32(use_D, d2, a2);

    a1 = a2 - a0;
}

} // namespace make_f8

//==========
// Square modulus
// a == BASE*a0 + a1 (mod BASE^2)
//==========
template <unsigned IDX, bool COMPLEX = false>
__device__ __forceinline__ fp8x2_e4m3 make_fp8x2(const int32_t a) {
    constexpr int32_t BASE = table::sqrt_moduli<IDX, COMPLEX>;
    static_assert(IDX < 20U);
    static_assert(table::kara_type<IDX, COMPLEX> == table::KaratsubaType::NONE);
    static_assert(BASE > 0 && BASE <= 49);
    static_assert(BASE * BASE == table::moduli<Backend::FP8, IDX, COMPLEX>);

#if GEMMUL8_FP8_RESIDUE_LUT && defined(__CUDACC__) && !defined(__HIPCC__)
    const uint32_t raw = make_f8::lookup_residue<IDX, COMPLEX, 32U>(a);
    fp8x2_e4m3 out;
    out.x.__x = uint8_t(raw);
    out.y.__x = uint8_t(raw >> 8);
    return out;
#endif

    if constexpr (BASE <= 32) {
        constexpr float BASE_F     = float(BASE);
        constexpr float INV_BASE_F = 1.0f / float(BASE);
        const float af             = __int2float_rn(a);
        const float q              = rintf(af * INV_BASE_F);
        const float r              = __fmaf_rn(-BASE_F, q, af);
        return make_f8::convert_f32x2(q, r);
    } else {
        const int32_t sign_mask = a >> 31;
        const uint32_t u        = uint32_t(abs(a));
        int32_t q;
        int32_t r;
        make_f8::decompose_square_large<BASE>(u, q, r);
        return make_f8::convert_i32x2(q, r, sign_mask);
    }
}

//==========
// Karatsuba with |limb| <= 32.
//==========
template <unsigned IDX, bool COMPLEX = false>
__device__ __forceinline__ void decompose_fp8x3(const uint32_t u, int32_t &a0, int32_t &a1, int32_t &a2) {
    constexpr int32_t P = table::moduli<Backend::FP8, IDX, COMPLEX>;
    constexpr auto TYPE = table::kara_type<IDX, COMPLEX>;
    static_assert(IDX < 20U);
    static_assert(TYPE != table::KaratsubaType::NONE);
    static_assert(make_f8::kara_max_abs<P, TYPE> <= 32U);

    using Plan = make_f8::limb_selection<P, TYPE, make_f8::kara_max_abs<P, TYPE>>;
    if constexpr (TYPE != table::KaratsubaType::BASE31_DIFF && Plan::masks.valid) {
        const int32_t q = Plan::base49 ? int32_t(make_f8::div_small_u32<48U>(u)) : int32_t(u >> 5);
        const int32_t r = int32_t(u) - Plan::base * q;
        make_f8::select_limbs<Plan>(q, r, a0, a1, a2);
    } else if constexpr (TYPE == table::KaratsubaType::BASE32_SUM || TYPE == table::KaratsubaType::BASE32_DIFF) {
        make_f8::decompose_base32<P, TYPE>(u, a0, a1, a2);
    } else if constexpr (TYPE == table::KaratsubaType::BASE33_SUM) {
        make_f8::decompose_base33_sum<P>(u, a0, a1, a2);
    } else if constexpr (TYPE == table::KaratsubaType::BASE31_DIFF) {
        make_f8::decompose_base31_diff<P>(u, a0, a1, a2);
    } else {
        static_assert(TYPE == table::KaratsubaType::BASE49);
        make_f8::decompose_base49<P>(u, a0, a1, a2);
    }
}

template <unsigned IDX, bool COMPLEX = false>
__device__ __forceinline__ fp8x3_e4m3 make_fp8x3(const int32_t a) {
#if GEMMUL8_FP8_RESIDUE_LUT && defined(__CUDACC__) && !defined(__HIPCC__)
    constexpr int32_t p      = table::moduli<Backend::FP8, IDX, COMPLEX>;
    constexpr auto type      = table::kara_type<IDX, COMPLEX>;
    constexpr uint32_t bound = make_f8::kara_max_abs<p, type>;
    if constexpr (make_f8::residue_bytes_host<IDX, COMPLEX, bound>.valid) {
        const uint32_t raw = make_f8::lookup_residue<IDX, COMPLEX, bound>(a);
        fp8x3_e4m3 out;
        out.x.__x = uint8_t(raw);
        out.y.__x = uint8_t(raw >> 8);
        out.z.__x = uint8_t(raw >> 16);
        return out;
    }
#endif
    int32_t a0, a1, a2;
    decompose_fp8x3<IDX, COMPLEX>(uint32_t(abs(a)), a0, a1, a2);
    return make_f8::convert_i32x3(a0, a1, a2, a >> 31);
}

namespace make_f8 {

__device__ __forceinline__ uint32_t signs4(int32_t v0, int32_t v1, int32_t v2, int32_t v3) {
    return ((uint32_t(v0) >> 24) & 0x00000080U) |
           ((uint32_t(v1) >> 16) & 0x00008000U) |
           ((uint32_t(v2) >> 8) & 0x00800000U) | (uint32_t(v3) & 0x80000000U);
}

__device__ __forceinline__ __nv_fp8x4_e4m3 convert_f32x4(float x0, float x1, float x2, float x3, uint32_t signs = 0U) {
    __nv_fp8x4_e4m3 out;
#if defined(__CUDACC__) && !defined(__HIPCC__)
    out.__x = (uint32_t(cvt_e4m3x2(x0, x1)) | (uint32_t(cvt_e4m3x2(x2, x3)) << 16)) ^ signs;
#else
    out = common::concat(__nv_fp8_e4m3(x0), __nv_fp8_e4m3(x1), __nv_fp8_e4m3(x2), __nv_fp8_e4m3(x3));
    out.__x ^= signs;
#endif
    return out;
}

__device__ __forceinline__ __nv_fp8x4_e4m3 convert_i32x4(int32_t x0, int32_t x1, int32_t x2, int32_t x3, uint32_t signs) {
    return convert_f32x4(__int2float_rn(x0), __int2float_rn(x1), __int2float_rn(x2), __int2float_rn(x3), signs);
}

} // namespace make_f8

template <unsigned IDX, bool COMPLEX = false>
__device__ __forceinline__ void make_fp8x3_vec4(int32_t v0, int32_t v1, int32_t v2, int32_t v3,
                                                __nv_fp8x4_e4m3 &x, __nv_fp8x4_e4m3 &y, __nv_fp8x4_e4m3 &z) {
#if GEMMUL8_FP8_RESIDUE_LUT && defined(__CUDACC__) && !defined(__HIPCC__)
    constexpr int32_t p      = table::moduli<Backend::FP8, IDX, COMPLEX>;
    constexpr auto type      = table::kara_type<IDX, COMPLEX>;
    constexpr uint32_t bound = make_f8::kara_max_abs<p, type>;
    if constexpr (make_f8::residue_bytes_host<IDX, COMPLEX, bound>.valid) {
        const uint32_t a = make_f8::lookup_residue<IDX, COMPLEX, bound>(v0);
        const uint32_t b = make_f8::lookup_residue<IDX, COMPLEX, bound>(v1);
        const uint32_t c = make_f8::lookup_residue<IDX, COMPLEX, bound>(v2);
        const uint32_t d = make_f8::lookup_residue<IDX, COMPLEX, bound>(v3);
        make_f8::transpose_residue_bytes<true>(a, b, c, d, x.__x, y.__x, z.__x);
        return;
    }
#endif
    int32_t a0, b0, c0, a1, b1, c1, a2, b2, c2, a3, b3, c3;
    decompose_fp8x3<IDX, COMPLEX>(uint32_t(abs(v0)), a0, b0, c0);
    decompose_fp8x3<IDX, COMPLEX>(uint32_t(abs(v1)), a1, b1, c1);
    decompose_fp8x3<IDX, COMPLEX>(uint32_t(abs(v2)), a2, b2, c2);
    decompose_fp8x3<IDX, COMPLEX>(uint32_t(abs(v3)), a3, b3, c3);
    const uint32_t signs = make_f8::signs4(v0, v1, v2, v3);
    x                    = make_f8::convert_i32x4(a0, a1, a2, a3, signs);
    y                    = make_f8::convert_i32x4(b0, b1, b2, b3, signs);
    z                    = make_f8::convert_i32x4(c0, c1, c2, c3, signs);
}

template <unsigned IDX, bool COMPLEX = false>
__device__ __forceinline__ void make_fp8x2_vec4(int32_t v0, int32_t v1, int32_t v2, int32_t v3,
                                                __nv_fp8x4_e4m3 &x, __nv_fp8x4_e4m3 &y) {
    constexpr int32_t BASE = table::sqrt_moduli<IDX, COMPLEX>;
    static_assert(table::kara_type<IDX, COMPLEX> == table::KaratsubaType::NONE);
#if GEMMUL8_FP8_RESIDUE_LUT && defined(__CUDACC__) && !defined(__HIPCC__)
    const uint32_t a = make_f8::lookup_residue<IDX, COMPLEX, 32U>(v0);
    const uint32_t b = make_f8::lookup_residue<IDX, COMPLEX, 32U>(v1);
    const uint32_t c = make_f8::lookup_residue<IDX, COMPLEX, 32U>(v2);
    const uint32_t d = make_f8::lookup_residue<IDX, COMPLEX, 32U>(v3);
    uint32_t unused;
    make_f8::transpose_residue_bytes<false>(a, b, c, d, x.__x, y.__x, unused);
    return;
#endif
    if constexpr (BASE <= 32) {
        constexpr float inv = 1.0f / float(BASE);
        const float a0 = __int2float_rn(v0), a1 = __int2float_rn(v1);
        const float a2 = __int2float_rn(v2), a3 = __int2float_rn(v3);
        const float q0 = rintf(a0 * inv), q1 = rintf(a1 * inv);
        const float q2 = rintf(a2 * inv), q3 = rintf(a3 * inv);
        x = make_f8::convert_f32x4(q0, q1, q2, q3);
        y = make_f8::convert_f32x4(__fmaf_rn(-float(BASE), q0, a0), __fmaf_rn(-float(BASE), q1, a1),
                                   __fmaf_rn(-float(BASE), q2, a2), __fmaf_rn(-float(BASE), q3, a3));
    } else {
        int32_t q0, r0, q1, r1, q2, r2, q3, r3;
        make_f8::decompose_square_large<BASE>(uint32_t(abs(v0)), q0, r0);
        make_f8::decompose_square_large<BASE>(uint32_t(abs(v1)), q1, r1);
        make_f8::decompose_square_large<BASE>(uint32_t(abs(v2)), q2, r2);
        make_f8::decompose_square_large<BASE>(uint32_t(abs(v3)), q3, r3);
        const uint32_t signs = make_f8::signs4(v0, v1, v2, v3);
        x                    = make_f8::convert_i32x4(q0, q1, q2, q3, signs);
        y                    = make_f8::convert_i32x4(r0, r1, r2, r3, signs);
    }
}

} // namespace gemmul8::common
