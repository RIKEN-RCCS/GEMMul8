#pragma once
#include "common.hpp"
#include "table.hpp"

namespace gemmul8::common {

namespace make_f8 {

template <uint32_t D>
inline constexpr uint32_t div_magic_u32 = uint32_t(((uint64_t(1) << 32) + D - 1) / D);

template <uint32_t D>
__device__ __forceinline__ uint32_t div_small_u32(const uint32_t x) {
#if defined(__CUDACC__) && !defined(__HIPCC__)
    uint32_t q;
    asm("mul.hi.u32 %0, %1, %2;"
        : "=r"(q)
        : "r"(x), "n"(div_magic_u32<D>));
    return q;
#else
    return x / D;
#endif
}

template <uint32_t MASK>
__device__ __forceinline__ uint32_t bfe_bit(const uint32_t pos) {
#if defined(__CUDACC__) && !defined(__HIPCC__)
    uint32_t r;
    asm("bfe.u32 %0, %1, %2, 1;"
        : "=r"(r)
        : "n"(MASK), "r"(pos));
    return r;
#else
    return (uint32_t(pos < 32u) * ((MASK >> (pos & 31u)) & 1u));
#endif
}

template <uint32_t LUT>
__device__ __forceinline__ int32_t bfe_s8(const uint32_t pos) {
#if defined(__CUDACC__) && !defined(__HIPCC__)
    int32_t r;
    asm("bfe.s32 %0, %1, %2, 8;"
        : "=r"(r)
        : "n"(LUT), "r"(pos));
    return r;
#else
    const uint32_t x = (LUT >> (pos & 31u)) & 0xffu;
    return int32_t(int8_t(x));
#endif
}

template <int32_t X0, int32_t X1, int32_t X2, int32_t X3>
inline constexpr uint32_t pack_s8x4_v = (uint32_t(uint8_t(X0))) |
                                        (uint32_t(uint8_t(X1)) << 8) |
                                        (uint32_t(uint8_t(X2)) << 16) |
                                        (uint32_t(uint8_t(X3)) << 24);

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

__device__ __forceinline__ fp8x2_e4m3 make_fp8x2_raw(uint16_t raw) {
    fp8x2_e4m3 out;
    out.x.__x = static_cast<__nv_fp8_storage_t>(raw);
    out.y.__x = static_cast<__nv_fp8_storage_t>(raw >> 8);
    return out;
}
#endif

__device__ __forceinline__ fp8x2_e4m3 convert_i32x2(const int32_t x0, const int32_t x1, const int32_t sign_mask) {
#if defined(__CUDACC__) && !defined(__HIPCC__)
    const float f0    = __int2float_rn(x0);
    const float f1    = __int2float_rn(x1);
    uint16_t raw      = cvt_e4m3x2(f0, f1);
    const uint16_t sm = uint16_t(uint32_t(sign_mask) & 0x8080u);
    raw ^= sm;
    return make_fp8x2_raw(raw);
#else
    const int32_t s = sign_mask | 1;
    fp8x2_e4m3 out;
    out.x = __nv_fp8_e4m3(s * x0);
    out.y = __nv_fp8_e4m3(s * x1);
    return out;
#endif
}

__device__ __forceinline__ fp8x3_e4m3 convert_i32x3(const int32_t x0, const int32_t x1, const int32_t x2, const int32_t sign_mask) {
#if defined(__CUDACC__) && !defined(__HIPCC__)
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
#else
    const int32_t s = sign_mask | 1;
    fp8x3_e4m3 out;
    out.x = __nv_fp8_e4m3(s * x0);
    out.y = __nv_fp8_e4m3(s * x1);
    out.z = __nv_fp8_e4m3(s * x2);
    return out;
#endif
}

__device__ __forceinline__ fp8x2_e4m3 convert_f32x2(const float x0, const float x1) {
#if defined(__CUDACC__) && !defined(__HIPCC__)
    return make_fp8x2_raw(cvt_e4m3x2(x0, x1));
#else
    fp8x2_e4m3 out;
    out.x = __nv_fp8_e4m3(x0);
    out.y = __nv_fp8_e4m3(x1);
    return out;
#endif
}

inline constexpr uint32_t FP8_POSITIVE_HOLES = 0xAAAA0000u;

template <int32_t BASE>
__device__ __forceinline__ void decompose_square_large(const uint32_t u, int32_t &q, int32_t &r) {
    static_assert(BASE > 32);
    static_assert(BASE <= 49);

    q = int32_t(div_small_u32<uint32_t(BASE)>(u + 16u));
    r = int32_t(u) - BASE * q;

    const int32_t fix_r = int32_t(bfe_bit<FP8_POSITIVE_HOLES>(uint32_t(r)));
    q += fix_r;
    r -= BASE * fix_r;

    const int32_t fix_q = int32_t(bfe_bit<FP8_POSITIVE_HOLES>(uint32_t(q)));
    q -= BASE * fix_q;
}

inline constexpr uint32_t KAR_BAD_C = 0xAAAAFFFFu;
inline constexpr uint32_t KAR_BAD_B = 0xFFFEAAAAu;

template <unsigned IDX>
inline constexpr uint32_t KAR_BAD_A0 = IDX == 6U    ? 0x000000AAu
                                       : IDX == 7U  ? 0x0000002Au
                                       : IDX == 8U  ? 0x00000015u
                                       : IDX == 10U ? 0x0000000Au
                                       : IDX == 11U ? 0x0000000Au
                                       : IDX == 13U ? 0x00000005u
                                       : IDX == 14U ? 0x00000005u
                                       : IDX == 15U ? 0x00000002u
                                       : IDX == 16U ? 0x00000002u
                                       : IDX == 17U ? 0x00000001u
                                       : IDX == 18U ? 0x00000001u
                                                    : 0u;

template <unsigned IDX>
inline constexpr bool KAR_A0_ONLY = IDX == 6U ||
                                    IDX == 7U ||
                                    IDX == 10U ||
                                    IDX == 11U;

template <unsigned IDX>
inline constexpr uint32_t KAR_BAD_A1 = IDX == 8U    ? 0xAAA00000u
                                       : IDX == 13U ? 0xA8000000u
                                       : IDX == 14U ? 0xAA800000u
                                       : IDX == 15U ? 0x00001555u
                                       : IDX == 16U ? 0x00000155u
                                       : IDX == 17U ? 0x00002AAAu
                                       : IDX == 18U ? 0x000002AAu
                                                    : 0u;

} // namespace make_f8

// ============================================================
// a == BASE*a0 + a1  (mod BASE^2)
// ============================================================
template <unsigned IDX>
__device__ __forceinline__ fp8x2_e4m3 make_fp8x2(const int32_t a) {
    constexpr int32_t BASE = table::sqrt_moduli<IDX>;

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

// ============================================================
// a == 49*a0 + a1  (mod p) == 48*a0 + a2  (mod p)
// a2 = a0 + a1
// ============================================================
template <unsigned IDX>
__device__ __forceinline__ fp8x3_e4m3 make_fp8x3(const int32_t a) {
    constexpr int32_t P = table::moduli<Backend::FP8, IDX>;

    constexpr int32_t H = (P + 24) / 48;
    constexpr int32_t T = P - 48 * H;

    const int32_t sign_mask = a >> 31;
    const uint32_t u        = uint32_t(abs(a));

    const int32_t q = int32_t(make_f8::div_small_u32<48u>(u));
    const int32_t r = int32_t(u) - 48 * q;
    const int32_t d = r - q;

    const uint32_t bad_C = make_f8::bfe_bit<make_f8::KAR_BAD_C>(uint32_t(r)) |
                           make_f8::bfe_bit<make_f8::KAR_BAD_C>(uint32_t(d - 1));
    const uint32_t bad_B = make_f8::bfe_bit<make_f8::KAR_BAD_B>(uint32_t(r - 16));
    uint32_t bad_A       = make_f8::bfe_bit<make_f8::KAR_BAD_A0<IDX>>(uint32_t(q));
    if constexpr (!make_f8::KAR_A0_ONLY<IDX>) {
        bad_A |= make_f8::bfe_bit<make_f8::KAR_BAD_A1<IDX>>(uint32_t(d));
    }

    const uint32_t cb         = bad_C & bad_B;
    const uint32_t k          = bad_C + cb + (cb & bad_A);
    constexpr uint32_t LUT_A0 = make_f8::pack_s8x4_v<+1, 0, 1 - H, H>;
    constexpr uint32_t LUT_A2 = make_f8::pack_s8x4_v<-48, 0, -T - 48, T>;
    const uint32_t bitpos     = k << 3;
    const int32_t a0          = q + make_f8::bfe_s8<LUT_A0>(bitpos);
    const int32_t a2          = r + make_f8::bfe_s8<LUT_A2>(bitpos);
    const int32_t a1          = a2 - a0;

    return make_f8::convert_i32x3(a0, a1, a2, sign_mask);
}

} // namespace gemmul8::common