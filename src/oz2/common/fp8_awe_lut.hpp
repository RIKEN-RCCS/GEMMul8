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

template <unsigned N> struct packed_lut { uint32_t value[N]{}; };

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
    static constexpr auto data = make_encoding<P>();
};

template <int32_t P> constexpr auto make_wide_encoding() {
    packed_lut<2 * P + 1> out{};
    for (int i = 0; i <= 2 * P; ++i) out.value[i] = encoding<P>::data.value[i % P];
    return out;
}

template <int32_t P> inline constexpr auto wide_host = make_wide_encoding<P>();

template <int32_t P> static __device__ const auto encoding_device = wide_host<P>;

} // namespace gemmul8::common::fp8_plan
