#pragma once
#include "../common/common.hpp"
#include "../common/make_f8.hpp"
#include "../common/fp8_awe_lut.hpp"
#include "complex_2m.hpp"

namespace gemmul8::mod {

template <unsigned IDX> struct fp8_exp_lut {
    uint16_t value[66];
};

template <unsigned IDX> constexpr auto build_fp8_exp_lut() {
    fp8_exp_lut<IDX> out{};

    constexpr uint32_t p = common::table::moduli<Backend::FP8, IDX>;
    static_assert(p >= 3 && p <= 32767);

    uint32_t x = 1;
    for (unsigned i = 0; i < 64; ++i) {
        out.value[33U * (i / 32U) + 31U - i % 32U] = uint16_t(x);
        x                                          = x * 2 % p;
    }
    return out;
}

template <unsigned IDX> inline constexpr auto fp8_exp_host          = build_fp8_exp_lut<IDX>();
template <unsigned IDX> static __device__ const auto fp8_exp_device = fp8_exp_host<IDX>;

template <unsigned IDX>
__device__ __forceinline__ int32_t fp8_exp_lookup(common::exp_t exp) {
    const unsigned i = unsigned(__clz(exp.val)) + (exp.is_hi ? 33U : 0U);
    return __ldg(&fp8_exp_device<IDX>.value[i]);
}

template <unsigned IDX, typename V>
__device__ __forceinline__ int32_t fp8_encode_nowrap(V v) {
    if constexpr (std::is_same_v<V, common::fp32_mant_exp> || std::is_same_v<V, common::fp64_mant_exp>) {
        constexpr int32_t p = common::table::moduli<Backend::FP8, IDX>;
        const int32_t mant  = [&] {
            if constexpr (std::is_same_v<V, common::fp32_mant_exp>) {
                return mod_small_nowrap<Backend::FP8, IDX>(v.mant);
            } else {
                using F = fp8_large_traits<IDX, false>;
                if constexpr (F::mant_max * uint64_t(p - 1) <= uint64_t(INT32_MAX)) {
                    const int32_t hi  = mod_small_nowrap<Backend::FP8, IDX>(v.mant.hi);
                    const uint32_t lo = mod_small_nowrap_u32<Backend::FP8, IDX>(v.mant.lo);
                    return hi * F::r32_centered + int32_t(lo);
                } else {
                    return reduce_mant_large<Backend::FP8, IDX>(v.mant);
                }
            }
        }();
        return mod_small_nowrap<Backend::FP8, IDX>(mant * fp8_exp_lookup<IDX>(v.exp));
    } else {
        return calc_mod_fp8_nowrap<IDX>(v);
    }
}

template <unsigned IDX, typename V>
__device__ __forceinline__ int32_t fp8_encode_centered(V v) {
    return center_reduced<Backend::FP8, IDX, false>(fp8_encode_nowrap<IDX>(v));
}

#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ == 890 || (__CUDA_ARCH__ >= 1000 && __CUDA_ARCH__ < 1200))

template <typename V>
inline constexpr bool fp8_unsigned_index =
    std::is_same_v<V, common::fp32_mant_exp> || std::is_same_v<V, common::fp64_mant_exp>;

#else

template <typename V>
inline constexpr bool fp8_unsigned_index = false;

#endif

template <int32_t P, bool BiasPointer = false, bool UnsignedIndex = false>
__device__ __forceinline__ uint32_t fp8_awe_bytes(int32_t r) {
    if constexpr (UnsignedIndex) {
        return __ldg(common::fp8_plan::encoding_device<P>.value + uint32_t(r + P / 2));
    } else if constexpr (BiasPointer) {
        return __ldg(common::fp8_plan::encoding_device<P>.value + P / 2 + r);
    } else {
        return __ldg(&common::fp8_plan::encoding_device<P>.value[r + P / 2]);
    }
}

template <int32_t Q>
__device__ __forceinline__ int32_t fp8_small_limb(int32_t r) {
    static_assert(Q >= 2 && Q <= 49);
    if constexpr (Q > 33) {
        const int32_t odd = r & 1;
        return r + ((r < -16 && odd) ? Q : 0) - ((r > 16 && odd) ? Q : 0);
    }
    return r;
}

template <int32_t Q> struct fp8_small_weights : residue_traits<Backend::FP8, unsigned(Q) + 20U, false> {
    using R                            = residue_traits<Backend::FP8, unsigned(Q) + 20U, false>;
    static constexpr uint32_t lo       = R::r00 | (R::r08 << 8) | (R::r16 << 16) | (R::r24 << 24);
    static constexpr uint32_t hi       = R::r32 | (R::r40 << 8) | (R::r48 << 16) | (R::r56 << 24);
    static constexpr uint64_t max_mant = 8ULL * 255ULL * (Q - 1) + Q;
    static constexpr uint64_t max_exp  = Q - 1;
    static_assert(max_mant * max_exp < uint64_t(INT32_MAX));
};

template <int32_t Q>
__device__ __forceinline__ uint32_t fp8_small_mant(int32_t v) {
    using W = fp8_small_weights<Q>;
    return dot4_u8(uint32_t(v), W::lo, v < 0 ? W::neg_corr_32 : 0U);
}

template <int32_t Q>
__device__ __forceinline__ uint32_t fp8_small_mant(common::mant_t v) {
    using W           = fp8_small_weights<Q>;
    const uint32_t hi = dot4_u8(uint32_t(v.hi), W::hi, v.hi < 0 ? W::neg_corr_64 : 0U);
    return dot4_u8(v.lo, W::lo, hi);
}

template <int32_t Q, typename V>
__device__ __forceinline__ int32_t fp8_small_mod(V v) {
    constexpr unsigned id = unsigned(Q) + 20U;
    if constexpr ((Q & (Q - 1)) == 0 && (std::is_same_v<V, int32_t> || std::is_same_v<V, common::mant_t>)) {
        static_assert(Q >= 2);
        constexpr unsigned shift = [] {
            unsigned s = 32U;
            for (uint32_t q = uint32_t(Q); q > 1U; q >>= 1U) --s;
            return s;
        }();
        const uint32_t low = [&] {
            if constexpr (std::is_same_v<V, int32_t>) return uint32_t(v);
            else return v.lo;
        }();
        return int32_t(low << shift) >> shift;
    } else if constexpr (std::is_same_v<V, int32_t>) {
        return fp8_encode_centered<id>(v);
    } else {
        using W = fp8_small_weights<Q>;
        if constexpr (std::is_same_v<V, common::mant_t>) {
            return mod_bounded<Backend::FP8, id, false, W::max_mant>(int32_t(fp8_small_mant<Q>(v)));
        } else {
            const uint32_t mant = fp8_small_mant<Q>(v.mant);
            const uint32_t exp  = fp8_exp_lookup<id>(v.exp);
            return mod_bounded<Backend::FP8, id, false, W::max_mant *(Q - 1)>(int32_t(mant * exp));
        }
    }
}

template <int32_t Q0, int32_t Q1, int32_t U0 = 1, int32_t U1 = 1>
struct pair_lut { uint16_t value[2 * Q0 * Q1 + 1]; };

template <int32_t Q0, int32_t Q1, int32_t U0 = 1, int32_t U1 = 1>
constexpr auto build_pair_lut() {
    constexpr int32_t p = Q0 * Q1;
    static_assert(Q0 >= 2 && Q0 <= 49 && Q1 >= 2 && Q1 <= 49 && p <= 2401);
    pair_lut<Q0, Q1, U0, U1> out{};
    for (int32_t x = -p / 2; x <= 2 * p - p / 2; ++x) {
        uint32_t packed = 0;
        for (unsigned j = 0; j < 2; ++j) {
            const int32_t q = j == 0 ? Q0 : Q1;
            int32_t r       = ((x % q) * (j == 0 ? U0 : U1)) % q;
            if (r > q / 2) r -= q;
            else if (r < -q / 2) r += q;
            if ((r > 16 || r < -16) && (r & 1)) r += r < 0 ? q : -q;
            packed |= common::fp8_plan::format_byte(common::make_f8::encode_limb(r)) << (8 * j);
        }
        out.value[x + p / 2] = uint16_t(packed);
    }
    return out;
}

template <int32_t Q0, int32_t Q1, int32_t U0 = 1, int32_t U1 = 1>
inline constexpr auto pair_host = build_pair_lut<Q0, Q1, U0, U1>();

template <int32_t Q0, int32_t Q1, int32_t U0 = 1, int32_t U1 = 1>
static __device__ const auto pair_device = pair_host<Q0, Q1, U0, U1>;

template <int32_t Q0, int32_t Q1, int32_t U0 = 1, int32_t U1 = 1,
          bool BiasPointer = false, bool UnsignedIndex = false>
__device__ __forceinline__ uint32_t pair_bytes(int32_t r) {
    if constexpr (UnsignedIndex) {
        return __ldg(pair_device<Q0, Q1, U0, U1>.value + uint32_t(r + Q0 * Q1 / 2));
    } else if constexpr (BiasPointer) {
        return __ldg(pair_device<Q0, Q1, U0, U1>.value + Q0 * Q1 / 2 + r);
    } else {
        return __ldg(&pair_device<Q0, Q1, U0, U1>.value[r + Q0 * Q1 / 2]);
    }
}

template <int32_t P, typename V>
__device__ __forceinline__ uint32_t crt_pair(V v) {
    constexpr auto s      = common::fp8_plan::scheme<P>;
    constexpr unsigned id = unsigned(s.q[0] * s.q[1]) + 20U;
    return pair_bytes<s.q[0], s.q[1],
                      common::fp8_plan::crt_unscale(s.p, 0),
                      common::fp8_plan::crt_unscale(s.p, 1),
                      (!std::is_same_v<V, int32_t> || (P != 59163 && P != 2303 && P != 1376 && P != 899 && P != 575)),
                      fp8_unsigned_index<V>>(
        fp8_encode_nowrap<id>(v));
}

template <int32_t P, unsigned J, typename V>
__device__ __forceinline__ int32_t fp8_crt_limb(V v) {
    constexpr int32_t q = common::fp8_plan::scheme<P>.q[J];
    return fp8_small_limb<q>(fp8_small_mod<q>(v));
}

template <int32_t P, typename V>
__device__ __forceinline__ int2 fp8_projection(V v) {
    constexpr auto s       = common::fp8_plan::scheme<P>;
    constexpr int32_t p    = s.w != 0 ? P : s.q[0] * s.q[1];
    constexpr unsigned id  = unsigned(p) + 20U;
    constexpr int32_t root = s.root_minus_one % p;
    static_assert((root * root + 1) % p == 0);
    const int32_t ar         = fp8_encode_nowrap<id>(v.x);
    const int32_t ai         = fp8_encode_nowrap<id>(v.y);
    const int32_t si         = root * ai;
    constexpr uint64_t bound = 2ULL * p * (1ULL + root);
    return {mod_bounded<Backend::FP8, id, true, bound>(ar + si),
            mod_bounded<Backend::FP8, id, true, bound>(ar - si)};
}

template <int32_t Q, int32_t Root> constexpr auto build_fp8_projection_lut() {
    common::fp8_plan::packed_lut<Q * Q, uint16_t> out{};
    for (int32_t a = -Q / 2; a <= Q / 2; ++a) {
        for (int32_t b = -Q / 2; b <= Q / 2; ++b) {
            const int32_t x = common::fp8_plan::congruence_limb(a + Root * b, Q);
            const int32_t y = common::fp8_plan::congruence_limb(a - Root * b, Q);

            out.value[(a + Q / 2) * Q + b + Q / 2] = uint16_t(
                common::fp8_plan::format_byte(common::make_f8::encode_limb(x)) |
                (common::fp8_plan::format_byte(common::make_f8::encode_limb(y)) << 8));
        }
    }
    return out;
}

template <int32_t Q, int32_t Root>
static __device__ const auto fp8_projection_device = build_fp8_projection_lut<Q, Root>();

template <int32_t P, unsigned J, typename V>
__device__ __forceinline__ uint32_t fp8_crt_projection(V v) {
    constexpr auto s    = common::fp8_plan::scheme<P>;
    constexpr int32_t q = s.q[J], root = s.root_minus_one % q;
    static_assert(root != 0 && (root * root + 1) % q == 0);
    const int32_t ar = fp8_small_mod<q>(v.x), ai = fp8_small_mod<q>(v.y);
    return __ldg(&fp8_projection_device<q, root>.value[(ar + q / 2) * q + ai + q / 2]);
}

__device__ __forceinline__ __nv_fp8_e4m3 fp8_raw(uint32_t raw) {
    __nv_fp8_e4m3 out;
    out.__x = static_cast<decltype(out.__x)>(raw);
    return out;
}

template <unsigned Products = 3>
__device__ __forceinline__ void store_awe_bytes(
    __nv_fp8_e4m3 *out,
    size_t next,
    uint32_t raw //
) {
    out[0]    = fp8_raw(raw);
    out[next] = fp8_raw(raw >> 8);
    if constexpr (Products == 3) out[2 * next] = fp8_raw(raw >> 16);
}

template <unsigned Products = 3>
__device__ __forceinline__ void store_awe_bytes(
    __nv_fp8x4_e4m3 *out,
    size_t next,
    uint32_t a, uint32_t b, uint32_t c, uint32_t d //
) {
    __nv_fp8x4_e4m3 x, y, z;
    common::make_f8::transpose_residue_bytes<Products == 3>(a, b, c, d, x.__x, y.__x, z.__x);
    out[0]    = x;
    out[next] = y;
    if constexpr (Products == 3) out[2 * next] = z;
}

template <unsigned IDX, typename V>
__device__ __forceinline__ void mod_launch(
    __nv_fp8_e4m3 *out,
    size_t next,
    V v //
) {
    static_assert(IDX >= 20U);
    constexpr int32_t p = int32_t(IDX - 20U);
    constexpr auto s    = common::fp8_plan::scheme<p>;
    if constexpr (s.w != 0) {
        store_awe_bytes<s.products>(out, next, fp8_awe_bytes<p, true, fp8_unsigned_index<V>>(fp8_encode_nowrap<IDX>(v)));
    } else {
        if constexpr (s.products >= 2) {
            const uint32_t raw = crt_pair<p>(v);
            out[0]             = fp8_raw(raw);
            out[next]          = fp8_raw(raw >> 8);
        } else {
            out[0] = __nv_fp8_e4m3(fp8_crt_limb<p, 0>(v));
        }
        if constexpr (s.products == 3) out[2 * next] = __nv_fp8_e4m3(fp8_crt_limb<p, 2>(v));
    }
}

template <unsigned IDX, typename V>
__device__ __forceinline__ void mod_launch(
    __nv_fp8x4_e4m3 *out,
    size_t next,
    V v0, V v1, V v2, V v3 //
) {
    static_assert(IDX >= 20U);
    constexpr int32_t p = int32_t(IDX - 20U);
    constexpr auto s    = common::fp8_plan::scheme<p>;
    if constexpr (s.w != 0) {
        const uint32_t a = fp8_awe_bytes<p, true, fp8_unsigned_index<V>>(fp8_encode_nowrap<IDX>(v0));
        const uint32_t b = fp8_awe_bytes<p, true, fp8_unsigned_index<V>>(fp8_encode_nowrap<IDX>(v1));
        const uint32_t c = fp8_awe_bytes<p, true, fp8_unsigned_index<V>>(fp8_encode_nowrap<IDX>(v2));
        const uint32_t d = fp8_awe_bytes<p, true, fp8_unsigned_index<V>>(fp8_encode_nowrap<IDX>(v3));
        store_awe_bytes<s.products>(out, next, a, b, c, d);
    } else {
        if constexpr (s.products >= 2) {
            __nv_fp8x4_e4m3 x, y;
            uint32_t unused;
            common::make_f8::transpose_residue_bytes<false>(
                crt_pair<p>(v0),
                crt_pair<p>(v1),
                crt_pair<p>(v2),
                crt_pair<p>(v3),
                x.__x, y.__x, unused);
            out[0]    = x;
            out[next] = y;
        } else {
            out[0] = common::make_f8::convert_i32x4(
                fp8_crt_limb<p, 0>(v0),
                fp8_crt_limb<p, 0>(v1),
                fp8_crt_limb<p, 0>(v2),
                fp8_crt_limb<p, 0>(v3),
                0U);
        }
        if constexpr (s.products == 3) {
            out[2 * next] = common::make_f8::convert_i32x4(
                fp8_crt_limb<p, 2>(v0),
                fp8_crt_limb<p, 2>(v1),
                fp8_crt_limb<p, 2>(v2),
                fp8_crt_limb<p, 2>(v3),
                0U);
        }
    }
}

struct fp8_three_words { uint32_t x, y, z; };

template <int32_t P, typename V>
__device__ __forceinline__ fp8_three_words fp8_three_pair(V v) {
    constexpr auto s = common::fp8_plan::scheme<P>;
    static_assert(s.products == 2 && s.w == 0);
    constexpr unsigned id = unsigned(P) + 20U;
    const int32_t ar = fp8_encode_nowrap<id>(v.x), ai = fp8_encode_nowrap<id>(v.y);
    const int32_t sum = mod_small_nowrap<Backend::FP8, id>(ar + ai);
    const auto bytes  = [](int32_t r) {
        return pair_bytes<s.q[0], s.q[1],
                          common::fp8_plan::crt_unscale(P, 0),
                          common::fp8_plan::crt_unscale(P, 1)>(r);
    };
    return {bytes(ar), bytes(ai), bytes(sum)};
}

template <int32_t P> struct fp8_three_lut { uint32_t value[P * P]; };

template <int32_t P> constexpr auto build_fp8_three_lut() {
    fp8_three_lut<P> out{};
    for (int32_t a = -P / 2; a <= P / 2; ++a) {
        for (int32_t b = -P / 2; b <= P / 2; ++b) {
            uint32_t packed = 0;
            for (unsigned j = 0; j < 3; ++j) {
                int32_t x = j == 0   ? a
                            : j == 1 ? b
                                     : a + b;

                if (x > P / 2) {
                    x -= P;
                } else if (x < -P / 2) {
                    x += P;
                }

                if ((x > 16 || x < -16) && (x & 1)) {
                    x += x < 0 ? P : -P;
                }

                packed |= common::fp8_plan::format_byte(common::make_f8::encode_limb(x)) << (8 * j);
            }
            out.value[(a + P / 2) * P + b + P / 2] = packed;
        }
    }
    return out;
}

template <int32_t P>
static __device__ const auto fp8_three_device = build_fp8_three_lut<P>();

template <int32_t P, typename V>
__device__ __forceinline__ uint32_t fp8_three_direct(V v) {
    static_assert(P == 33 || P == 49);
    const int32_t ar = fp8_small_mod<P>(v.x), ai = fp8_small_mod<P>(v.y);
    return __ldg(&fp8_three_device<P>.value[(ar + P / 2) * P + ai + P / 2]);
}

template <unsigned IDX, typename V>
__device__ __forceinline__ void mod_launch(__nv_fp8_e4m3 *plus, __nv_fp8_e4m3 *minus, size_t next, V v) {
    static_assert(IDX >= 20U);
    constexpr int32_t p = int32_t(IDX - 20U);
    constexpr auto s    = common::fp8_plan::scheme<p>;
    if constexpr (s.root_minus_one == 0) {
        if constexpr (s.products == 2) {
            const auto r = fp8_three_pair<p>(v);
            store_awe_bytes<2>(plus, next, r.x);
            store_awe_bytes<2>(minus, next, r.y);
            store_awe_bytes<2>(minus + 2 * next, next, r.z);
        } else {
            const auto r = fp8_three_direct<p>(v);
            plus[0]      = fp8_raw(r);
            minus[0]     = fp8_raw(r >> 8);
            minus[next]  = fp8_raw(r >> 16);
        }
    } else if constexpr (s.w != 0) {
        const int2 r = fp8_projection<p>(v);
        store_awe_bytes(plus, next, fp8_awe_bytes<p>(r.x));
        store_awe_bytes(minus, next, fp8_awe_bytes<p>(r.y));
    } else {
        if constexpr (s.products >= 2) {
            const int2 r     = fp8_projection<p>(v);
            const uint32_t a = pair_bytes<s.q[0], s.q[1],
                                          common::fp8_plan::crt_unscale(s.p, 0),
                                          common::fp8_plan::crt_unscale(s.p, 1)>(r.x);
            const uint32_t b = pair_bytes<s.q[0], s.q[1],
                                          common::fp8_plan::crt_unscale(s.p, 0),
                                          common::fp8_plan::crt_unscale(s.p, 1)>(r.y);
            plus[0]          = fp8_raw(a);
            plus[next]       = fp8_raw(a >> 8);
            minus[0]         = fp8_raw(b);
            minus[next]      = fp8_raw(b >> 8);
        }
        if constexpr (s.products != 2) {
            constexpr unsigned J = s.products == 3 ? 2U : 0U;
            const uint32_t r     = fp8_crt_projection<p, J>(v);
            plus[J * next]       = fp8_raw(r);
            minus[J * next]      = fp8_raw(r >> 8);
        }
    }
}

template <unsigned IDX, typename V>
__device__ __forceinline__ void mod_launch(
    __nv_fp8x4_e4m3 *plus, __nv_fp8x4_e4m3 *minus,
    size_t next, V v0, V v1, V v2, V v3 //
) {
    static_assert(IDX >= 20U);
    constexpr int32_t p = int32_t(IDX - 20U);
    constexpr auto s    = common::fp8_plan::scheme<p>;
    if constexpr (s.root_minus_one == 0) {
        if constexpr (s.products == 2) {
            const auto a = fp8_three_pair<p>(v0), b = fp8_three_pair<p>(v1);
            const auto c = fp8_three_pair<p>(v2), d = fp8_three_pair<p>(v3);
            store_awe_bytes<2>(plus, next, a.x, b.x, c.x, d.x);
            store_awe_bytes<2>(minus, next, a.y, b.y, c.y, d.y);
            store_awe_bytes<2>(minus + 2 * next, next, a.z, b.z, c.z, d.z);
        } else {
            const auto a = fp8_three_direct<p>(v0), b = fp8_three_direct<p>(v1);
            const auto c = fp8_three_direct<p>(v2), d = fp8_three_direct<p>(v3);
            __nv_fp8x4_e4m3 x, y, z;
            common::make_f8::transpose_residue_bytes<true>(a, b, c, d, x.__x, y.__x, z.__x);
            plus[0]     = x;
            minus[0]    = y;
            minus[next] = z;
        }
    } else if constexpr (s.w != 0) {
        auto encode = [](V v) {
            const int2 r = fp8_projection<p>(v);
            return uint2{fp8_awe_bytes<p>(r.x), fp8_awe_bytes<p>(r.y)};
        };
        const uint2 a = encode(v0), b = encode(v1), c = encode(v2), d = encode(v3);
        store_awe_bytes(plus, next, a.x, b.x, c.x, d.x);
        store_awe_bytes(minus, next, a.y, b.y, c.y, d.y);
    } else {
        if constexpr (s.products >= 2) {
            const int2 a = fp8_projection<p>(v0), b = fp8_projection<p>(v1);
            const int2 c = fp8_projection<p>(v2), d = fp8_projection<p>(v3);
            __nv_fp8x4_e4m3 x, y;
            uint32_t unused;
            common::make_f8::transpose_residue_bytes<false>(
                pair_bytes<s.q[0], s.q[1],
                           common::fp8_plan::crt_unscale(s.p, 0),
                           common::fp8_plan::crt_unscale(s.p, 1)>(a.x),
                pair_bytes<s.q[0], s.q[1],
                           common::fp8_plan::crt_unscale(s.p, 0),
                           common::fp8_plan::crt_unscale(s.p, 1)>(b.x),
                pair_bytes<s.q[0], s.q[1],
                           common::fp8_plan::crt_unscale(s.p, 0),
                           common::fp8_plan::crt_unscale(s.p, 1)>(c.x),
                pair_bytes<s.q[0], s.q[1],
                           common::fp8_plan::crt_unscale(s.p, 0),
                           common::fp8_plan::crt_unscale(s.p, 1)>(d.x),
                x.__x, y.__x, unused);
            plus[0]    = x;
            plus[next] = y;
            common::make_f8::transpose_residue_bytes<false>(
                pair_bytes<s.q[0], s.q[1],
                           common::fp8_plan::crt_unscale(s.p, 0),
                           common::fp8_plan::crt_unscale(s.p, 1)>(a.y),
                pair_bytes<s.q[0], s.q[1],
                           common::fp8_plan::crt_unscale(s.p, 0),
                           common::fp8_plan::crt_unscale(s.p, 1)>(b.y),
                pair_bytes<s.q[0], s.q[1],
                           common::fp8_plan::crt_unscale(s.p, 0),
                           common::fp8_plan::crt_unscale(s.p, 1)>(c.y),
                pair_bytes<s.q[0], s.q[1],
                           common::fp8_plan::crt_unscale(s.p, 0),
                           common::fp8_plan::crt_unscale(s.p, 1)>(d.y),
                x.__x, y.__x, unused);
            minus[0]    = x;
            minus[next] = y;
        }
        if constexpr (s.products != 2) {
            constexpr unsigned J = s.products == 3 ? 2U : 0U;
            const uint32_t a = fp8_crt_projection<p, J>(v0), b = fp8_crt_projection<p, J>(v1);
            const uint32_t c = fp8_crt_projection<p, J>(v2), d = fp8_crt_projection<p, J>(v3);
            __nv_fp8x4_e4m3 x, y;
            uint32_t unused;
            common::make_f8::transpose_residue_bytes<false>(a, b, c, d, x.__x, y.__x, unused);
            plus[J * next]  = x;
            minus[J * next] = y;
        }
    }
}

} // namespace gemmul8::mod
