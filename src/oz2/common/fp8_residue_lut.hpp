#pragma once
#include "fp8_limb_selection.hpp"

#ifndef GEMMUL8_FP8_RESIDUE_LUT
    #define GEMMUL8_FP8_RESIDUE_LUT 1
#endif

namespace gemmul8::common::make_f8 {

constexpr uint32_t encode_limb(int32_t x) {
    const uint32_t a = uint32_t(x < 0 ? -x : x);
    if (a == 0) return 0;
    unsigned e = 0;
    for (unsigned t = a; t > 1; t >>= 1) ++e;
    return (x < 0 ? 0x80U : 0U) | ((e + 7U) << 3) | (((a << 3) >> e) & 7U);
}

template <unsigned SIZE> struct residue_bytes {
    uint32_t value[SIZE]{};
    bool valid = true;
};

template <unsigned IDX, bool COMPLEX>
inline constexpr bool use_signed_square_lut =
    table::kara_type<IDX, COMPLEX> == table::KaratsubaType::NONE &&
    table::sqrt_moduli<IDX, COMPLEX> != 0;

struct square_limbs {
    int32_t a;
    int32_t b;
};

template <unsigned IDX, bool COMPLEX>
constexpr square_limbs make_square_limbs(int32_t u) {
    constexpr int32_t base = table::sqrt_moduli<IDX, COMPLEX>;
    int32_t a, b;
    if constexpr (base <= 32) {
        a               = u / base;
        const int32_t r = u - base * a;
        if (2 * r > base || (2 * r == base && (a & 1))) ++a;
        b = u - base * a;
    } else {
        a = (u + 16) / base;
        b = u - base * a;
        if (b > 16 && (b & 1)) {
            ++a;
            b -= base;
        }
        if ((uint32_t(a) & 0x11U) == 0x11U) a -= base;
    }
    return {a, b};
}

template <unsigned SIZE> struct signed_square_bytes {
    uint16_t value[SIZE]{};
    bool valid = true;
};

template <unsigned IDX, bool COMPLEX>
constexpr auto build_signed_square_table() {
    constexpr int32_t p    = table::moduli<Backend::FP8, IDX, COMPLEX>;
    constexpr int32_t half = p / 2;
    constexpr int32_t base = table::sqrt_moduli<IDX, COMPLEX>;
    static_assert(use_signed_square_lut<IDX, COMPLEX>);

    signed_square_bytes<unsigned(2 * half + 1)> result;
    for (int32_t v = -half; v <= half; ++v) {
        const int32_t u   = v < 0 ? -v : v;
        const auto [a, b] = make_square_limbs<IDX, COMPLEX>(u);

        for (const int32_t limb : {a, b}) {
            const uint32_t mag = uint32_t(limb < 0 ? -limb : limb);
            if (mag > 32U || (mag & 0x11U) == 0x11U) result.valid = false;
        }

        const uint32_t raw = encode_limb(a) | (encode_limb(b) << 8);
        uint32_t sign      = v < 0 ? 0x00008080U : 0U;
        if constexpr (base <= 32) {
            if ((raw & 0xff00U) == 0) sign &= ~0x8000U;
        }
        result.value[unsigned(v + half)] = uint16_t(raw ^ sign);
    }
    return result;
}

template <unsigned IDX, bool COMPLEX>
inline constexpr auto signed_square_table_host = build_signed_square_table<IDX, COMPLEX>();

#if defined(__CUDACC__) && !defined(__HIPCC__)
template <unsigned IDX, bool COMPLEX>
static __device__ const auto signed_square_table_device = signed_square_table_host<IDX, COMPLEX>;
#endif

template <unsigned IDX, bool COMPLEX, uint32_t MAX_ABS>
constexpr auto build_residue_bytes() {
    constexpr int32_t p    = table::moduli<Backend::FP8, IDX, COMPLEX>;
    constexpr auto type    = table::kara_type<IDX, COMPLEX>;
    constexpr int32_t base = table::sqrt_moduli<IDX, COMPLEX>;
    residue_bytes<unsigned(p / 2 + 1)> result;
    for (int32_t u = 0; u <= p / 2; ++u) {
        int32_t a, b, c = 0;
        if constexpr (type == table::KaratsubaType::NONE) {
            const auto limbs = make_square_limbs<IDX, COMPLEX>(u);
            a                = limbs.a;
            b                = limbs.b;
        } else {
            using Plan   = limb_selection<p, type, MAX_ABS>;
            const auto s = Plan::classify(u);
            if (type == table::KaratsubaType::BASE31_DIFF || s.choice > 3U) {
                result.valid = false;
                continue;
            }
            a = s.f[0] + Plan::delta0[s.choice];
            b = s.f[1] + Plan::delta1[s.choice];
            if constexpr (Plan::base49 || Plan::base33) {
                c = b;
                b -= a;
            } else c = Plan::sum32 ? a + b : b - a;
        }
        for (const int32_t limb : {a, b, c}) {
            const uint32_t mag = uint32_t(limb < 0 ? -limb : limb);
            if (mag > 32U || (mag & 0x11U) == 0x11U) result.valid = false;
        }
        result.value[u] = encode_limb(a) | (encode_limb(b) << 8) | (encode_limb(c) << 16);
    }
    return result;
}

template <unsigned IDX, bool COMPLEX, uint32_t MAX_ABS>
inline constexpr auto residue_bytes_host = build_residue_bytes<IDX, COMPLEX, MAX_ABS>();

#if defined(__CUDACC__) && !defined(__HIPCC__)
template <unsigned IDX, bool COMPLEX, uint32_t MAX_ABS>
static __device__ const auto residue_bytes_device = residue_bytes_host<IDX, COMPLEX, MAX_ABS>;
#endif

template <unsigned IDX, bool COMPLEX, uint32_t MAX_ABS>
__device__ __forceinline__ uint32_t lookup_residue(int32_t v) {
    if constexpr (use_signed_square_lut<IDX, COMPLEX>) {
        static_assert(signed_square_table_host<IDX, COMPLEX>.valid);
        constexpr int32_t half = table::moduli<Backend::FP8, IDX, COMPLEX> / 2;
        const uint32_t i       = uint32_t(v + half);
#if defined(__CUDACC__) && !defined(__HIPCC__)
        return uint32_t(__ldg(&signed_square_table_device<IDX, COMPLEX>.value[i]));
#else
        return uint32_t(signed_square_table_host<IDX, COMPLEX>.value[i]);
#endif
    } else {
        static_assert(residue_bytes_host<IDX, COMPLEX, MAX_ABS>.valid);
        const uint32_t u = uint32_t(abs(v));
#if defined(__CUDACC__) && !defined(__HIPCC__)
        const uint32_t raw = __ldg(&residue_bytes_device<IDX, COMPLEX, MAX_ABS>.value[u]);
#else
        const uint32_t raw = residue_bytes_host<IDX, COMPLEX, MAX_ABS>.value[u];
#endif
        uint32_t sign = uint32_t(v >> 31) & 0x00808080U;
        if constexpr (table::kara_type<IDX, COMPLEX> == table::KaratsubaType::NONE &&
                      table::sqrt_moduli<IDX, COMPLEX> <= 32) {
            if ((raw & 0xff00U) == 0) sign &= ~0x8000U;
        }
        return raw ^ sign;
    }
}

__device__ __forceinline__ uint32_t permute_bytes(uint32_t a, uint32_t b, uint32_t selector) {
#if defined(__CUDACC__) && !defined(__HIPCC__)
    return __byte_perm(a, b, selector);
#else
    const uint64_t bits = uint64_t(a) | (uint64_t(b) << 32);
    uint32_t out        = 0;
    for (unsigned j = 0; j < 4; ++j) out |= uint32_t((bits >> (8 * ((selector >> (4 * j)) & 7))) & 255) << (8 * j);
    return out;
#endif
}

template <bool THREE>
__device__ __forceinline__ void transpose_residue_bytes(
    uint32_t a, uint32_t b, uint32_t c, uint32_t d,
    uint32_t &x, uint32_t &y, uint32_t &z //
) {
    const uint32_t ab = permute_bytes(a, b, 0x5140);
    const uint32_t cd = permute_bytes(c, d, 0x5140);
    x                 = permute_bytes(ab, cd, 0x5410);
    y                 = permute_bytes(ab, cd, 0x7632);
    if constexpr (THREE) {
        const uint32_t ab2 = permute_bytes(a, b, 0x0062);
        const uint32_t cd2 = permute_bytes(c, d, 0x0062);
        z                  = permute_bytes(ab2, cd2, 0x5410);
    }
}

} // namespace gemmul8::common::make_f8
