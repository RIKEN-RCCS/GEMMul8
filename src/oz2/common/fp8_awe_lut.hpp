#pragma once
#include "fp8_plan.hpp"
#include "fp8_residue_lut.hpp"

namespace gemmul8::common::fp8_plan {

constexpr uint32_t format_byte(uint32_t x) {
#if FP8_FNUZ
    return (x & 0x7fU) == 0 ? 0U : x + 8U;
#else
    return x;
#endif
}

template <unsigned N, typename T = uint32_t> struct packed_lut { T value[N]{}; };

constexpr int congruence_limb(int x, int q) {
    x %= q;
    if (x < 0) x += q;
    if (x > q / 2) x -= q;
    if ((x > 16 || x < -16) && (x & 1)) x += x < 0 ? q : -q;
    return x;
}

template <int32_t P> constexpr auto make_congruence_encoding() {
    using T = std::conditional_t<scheme<P>.products == 2, uint16_t, uint32_t>;
    packed_lut<P, T> out{};
    for (int i = 0; i < P; ++i) {
        const int a = i - P / 2;
        int u = 0, v = 0, w = 0;
        if constexpr (P == 343) {
            // a = 7*u - 6*v; u == v (mod 7), hence a*b = 7*U - 6*V.
            int best = INT32_MAX;
            for (int x = -32; x <= 32; ++x) {
                if ((x > 16 || x < -16) && (x & 1)) continue;
                const int y = ((7 * x - a) * 286 % P + P) % P;
                const int z = y > P / 2 ? y - P : y;
                if (z < -32 || z > 32 || ((z < -16 || z > 16) && (z & 1)) || (x - z) % 7) continue;
                const int ax = x < 0 ? -x : x, az = z < 0 ? -z : z;
                const int score = 129 * (ax > az ? ax : az) + ax + az;
                if (score < best) {
                    best = score;
                    u    = x;
                    v    = z;
                }
            }
            if (best == INT32_MAX) {
                out.value[i] = T(-1);
                continue;
            }
        } else {
            static_assert(P == 19305 || P == 11583);
            // a = u + r*digit, 3 | r
            constexpr int r = -scheme<P>.w;
            u               = congruence_limb(a, r);
            const int digit = (a - u) / r;
            v               = congruence_limb(22 * digit + 12 * a, 33);
            w               = congruence_limb(13 * (u + digit) + 27 * a, 39);
        }
        out.value[i] = T(format_byte(make_f8::encode_limb(u)) |
                         (format_byte(make_f8::encode_limb(v)) << 8) |
                         (format_byte(make_f8::encode_limb(w)) << 16));
    }
    return out;
}

// signed (a,b,a+b) table
template <int32_t P> constexpr auto make_encoding() {
    constexpr auto s = scheme<P>;
    struct pair { int a, b, score = INT32_MAX; } best[P / 2 + 1]{};
    const auto abs_i = [](int x) { return x < 0 ? -x : x; };
    const auto valid = [&](int x) { return abs_i(x) <= 32 && (abs_i(x) <= 16 || !(x & 1)); };
    for (int a = -32; a <= 32; ++a) {
        if (!valid(a)) continue;
        for (int b = -32; b <= 32; ++b) {
            if (!valid(b) || !valid(a + b)) continue;
            const int r = ((a + s.w * b) % P + P) % P;
            if (r > P / 2) continue;
            const int aa    = abs_i(a);
            const int bb    = abs_i(b);
            const int cc    = abs_i(a + b);
            const int ab    = aa > bb ? aa : bb;
            const int score = 129 * (ab > cc ? ab : cc) + aa + bb + cc;
            if (score < best[r].score) best[r] = {a, b, score};
        }
    }

    // x = h(a + wb): h is the least positive root of h^2 = c0 + c2 (mod P)
    int h = 1, inverse = 1;
    while ((h * h - s.coefficient[0] - s.coefficient[2]) % P != 0) ++h;
    while ((inverse * h) % P != 1) ++inverse;
    packed_lut<P> out{};
    for (int i = 0; i < P; ++i) {
        int y = ((i - P / 2) * inverse) % P;
        if (y > P / 2) y -= P;
        if (y < -P / 2) y += P;
        const auto v = best[abs_i(y)];
        if (v.score == INT32_MAX) {
            out.value[i] = UINT32_MAX;
            continue;
        }
        const int a = y < 0 ? -v.a : v.a, b = y < 0 ? -v.b : v.b;
        out.value[i] = format_byte(make_f8::encode_limb(a)) |
                       (format_byte(make_f8::encode_limb(b)) << 8) |
                       (format_byte(make_f8::encode_limb(a + b)) << 16);
    }
    return out;
}

template <int32_t P> struct encoding {
    static constexpr auto data = [] {
        if constexpr (scheme<P>.w < 0) return make_congruence_encoding<P>();
        else return make_encoding<P>();
    }();
    static_assert(
        [] {
            if constexpr (scheme<P>.w < 0) {
                for (const auto x : data.value) {
                    for (unsigned j = 0; j < scheme<P>.products; ++j)
                        if (((x >> (8 * j)) & 0x7fU) > 0x60U + 8U * FP8_FNUZ) return false;
                }
            }
            return true;
        }(),
        "The congruence encoding must cover every residue with exact FP8 limbs");
};

template <int32_t P> constexpr auto make_wide_encoding() {
    using T = std::conditional_t<scheme<P>.products == 2, uint16_t, uint32_t>;
    packed_lut<2 * P + 1, T> out{};
    for (int i = 0; i <= 2 * P; ++i) out.value[i] = encoding<P>::data.value[i % P];
    return out;
}

template <int32_t P> inline constexpr auto wide_host = make_wide_encoding<P>();

template <int32_t P> static __device__ const auto encoding_device = wide_host<P>;

} // namespace gemmul8::common::fp8_plan
