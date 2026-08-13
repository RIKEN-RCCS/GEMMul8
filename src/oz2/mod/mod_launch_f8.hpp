#pragma once
#include "../common/common.hpp"
#include "../common/make_f8.hpp"
#include "mod_core.hpp"

namespace gemmul8::mod {

// launcher for V in {int32_t, int64_t, float, double}
template <unsigned IDX, typename V> __device__ __forceinline__ void mod_launch(
    __nv_fp8_e4m3 *__restrict__ out,
    size_t next,
    V v //
) {
    if constexpr (common::table::sqrt_moduli<IDX> != 0) {

        const common::fp8x2_e4m3 r = common::make_fp8x2<IDX>(calc_mod<Backend::FP8, IDX>(v));

        out[0]    = r.x;
        out[next] = r.y;

    } else {

        const common::fp8x3_e4m3 rem = common::make_fp8x3<IDX>(calc_mod<Backend::FP8, IDX>(v));

        out[0]        = rem.x;
        out[next]     = rem.y;
        out[next * 2] = rem.z;
    }
}

// launcher for V in {int32_t, int64_t, float, double}
template <unsigned IDX, typename V> __device__ __forceinline__ void mod_launch(
    __nv_fp8x4_e4m3 *__restrict__ out,
    size_t next,
    V v0, V v1, V v2, V v3 //
) {
    if constexpr (common::table::sqrt_moduli<IDX> != 0) {

        const common::fp8x2_e4m3 rem0 = common::make_fp8x2<IDX>(calc_mod<Backend::FP8, IDX>(v0));
        const common::fp8x2_e4m3 rem1 = common::make_fp8x2<IDX>(calc_mod<Backend::FP8, IDX>(v1));
        const common::fp8x2_e4m3 rem2 = common::make_fp8x2<IDX>(calc_mod<Backend::FP8, IDX>(v2));
        const common::fp8x2_e4m3 rem3 = common::make_fp8x2<IDX>(calc_mod<Backend::FP8, IDX>(v3));

        out[0]    = common::concat(rem0.x, rem1.x, rem2.x, rem3.x);
        out[next] = common::concat(rem0.y, rem1.y, rem2.y, rem3.y);

    } else {

        const common::fp8x3_e4m3 r0 = common::make_fp8x3<IDX>(calc_mod<Backend::FP8, IDX>(v0));
        const common::fp8x3_e4m3 r1 = common::make_fp8x3<IDX>(calc_mod<Backend::FP8, IDX>(v1));
        const common::fp8x3_e4m3 r2 = common::make_fp8x3<IDX>(calc_mod<Backend::FP8, IDX>(v2));
        const common::fp8x3_e4m3 r3 = common::make_fp8x3<IDX>(calc_mod<Backend::FP8, IDX>(v3));

        out[0]        = common::concat(r0.x, r1.x, r2.x, r3.x);
        out[next]     = common::concat(r0.y, r1.y, r2.y, r3.y);
        out[next * 2] = common::concat(r0.z, r1.z, r2.z, r3.z);
    }
}

// V in {int2, int64x2_t, float2, double2}
template <unsigned IDX, typename V> __device__ __forceinline__ void mod_launch(
    __nv_fp8_e4m3 *__restrict__ out_r,
    __nv_fp8_e4m3 *__restrict__ out_i,
    __nv_fp8_e4m3 *__restrict__ out_ri,
    size_t next,
    V v //
) {
    if constexpr (common::table::sqrt_moduli<IDX> != 0) {

        const int32_t rem_r              = calc_mod<Backend::FP8, IDX>(v.x);
        const common::fp8x2_e4m3 rem_rx2 = common::make_fp8x2<IDX>(rem_r);
        out_r[0]                         = rem_rx2.x;
        out_r[next]                      = rem_rx2.y;

        const int32_t rem_i              = calc_mod<Backend::FP8, IDX>(v.y);
        const common::fp8x2_e4m3 rem_ix2 = common::make_fp8x2<IDX>(rem_i);
        out_i[0]                         = rem_ix2.x;
        out_i[next]                      = rem_ix2.y;

        const int32_t rem_ri              = wrapping<Backend::FP8, IDX>(rem_r + rem_i);
        const common::fp8x2_e4m3 rem_rix2 = common::make_fp8x2<IDX>(rem_ri);
        out_ri[0]                         = rem_rix2.x;
        out_ri[next]                      = rem_rix2.y;

    } else {

        const int32_t rem_r              = calc_mod<Backend::FP8, IDX>(v.x);
        const common::fp8x3_e4m3 rem_rx3 = common::make_fp8x3<IDX>(rem_r);
        out_r[0]                         = rem_rx3.x;
        out_r[next]                      = rem_rx3.y;
        out_r[next * 2]                  = rem_rx3.z;

        const int32_t rem_i              = calc_mod<Backend::FP8, IDX>(v.y);
        const common::fp8x3_e4m3 rem_ix3 = common::make_fp8x3<IDX>(rem_i);
        out_i[0]                         = rem_ix3.x;
        out_i[next]                      = rem_ix3.y;
        out_i[next * 2]                  = rem_ix3.z;

        const int32_t rem_ri              = wrapping<Backend::FP8, IDX>(rem_r + rem_i);
        const common::fp8x3_e4m3 rem_rix3 = common::make_fp8x3<IDX>(rem_ri);
        out_ri[0]                         = rem_rix3.x;
        out_ri[next]                      = rem_rix3.y;
        out_ri[next * 2]                  = rem_rix3.z;
    }
}

// V in {int2, int64x2_t, float2, double2}
template <unsigned IDX, typename V> __device__ __forceinline__ void mod_launch(
    __nv_fp8x4_e4m3 *__restrict__ out_r,
    __nv_fp8x4_e4m3 *__restrict__ out_i,
    __nv_fp8x4_e4m3 *__restrict__ out_ri,
    size_t next,
    V v0, V v1, V v2, V v3 //
) {
    if constexpr (common::table::sqrt_moduli<IDX> != 0) {

        const int32_t rem_r0 = calc_mod<Backend::FP8, IDX>(v0.x);
        const int32_t rem_r1 = calc_mod<Backend::FP8, IDX>(v1.x);
        const int32_t rem_r2 = calc_mod<Backend::FP8, IDX>(v2.x);
        const int32_t rem_r3 = calc_mod<Backend::FP8, IDX>(v3.x);

        const common::fp8x2_e4m3 rem_r0x2 = common::make_fp8x2<IDX>(rem_r0);
        const common::fp8x2_e4m3 rem_r1x2 = common::make_fp8x2<IDX>(rem_r1);
        const common::fp8x2_e4m3 rem_r2x2 = common::make_fp8x2<IDX>(rem_r2);
        const common::fp8x2_e4m3 rem_r3x2 = common::make_fp8x2<IDX>(rem_r3);

        out_r[0]    = common::concat(rem_r0x2.x, rem_r1x2.x, rem_r2x2.x, rem_r3x2.x);
        out_r[next] = common::concat(rem_r0x2.y, rem_r1x2.y, rem_r2x2.y, rem_r3x2.y);

        const int32_t rem_i0 = calc_mod<Backend::FP8, IDX>(v0.y);
        const int32_t rem_i1 = calc_mod<Backend::FP8, IDX>(v1.y);
        const int32_t rem_i2 = calc_mod<Backend::FP8, IDX>(v2.y);
        const int32_t rem_i3 = calc_mod<Backend::FP8, IDX>(v3.y);

        const common::fp8x2_e4m3 rem_i0x2 = common::make_fp8x2<IDX>(rem_i0);
        const common::fp8x2_e4m3 rem_i1x2 = common::make_fp8x2<IDX>(rem_i1);
        const common::fp8x2_e4m3 rem_i2x2 = common::make_fp8x2<IDX>(rem_i2);
        const common::fp8x2_e4m3 rem_i3x2 = common::make_fp8x2<IDX>(rem_i3);

        out_i[0]    = common::concat(rem_i0x2.x, rem_i1x2.x, rem_i2x2.x, rem_i3x2.x);
        out_i[next] = common::concat(rem_i0x2.y, rem_i1x2.y, rem_i2x2.y, rem_i3x2.y);

        const int32_t rem_ri0 = wrapping<Backend::FP8, IDX>(rem_r0 + rem_i0);
        const int32_t rem_ri1 = wrapping<Backend::FP8, IDX>(rem_r1 + rem_i1);
        const int32_t rem_ri2 = wrapping<Backend::FP8, IDX>(rem_r2 + rem_i2);
        const int32_t rem_ri3 = wrapping<Backend::FP8, IDX>(rem_r3 + rem_i3);

        const common::fp8x2_e4m3 rem_ri0x2 = common::make_fp8x2<IDX>(rem_ri0);
        const common::fp8x2_e4m3 rem_ri1x2 = common::make_fp8x2<IDX>(rem_ri1);
        const common::fp8x2_e4m3 rem_ri2x2 = common::make_fp8x2<IDX>(rem_ri2);
        const common::fp8x2_e4m3 rem_ri3x2 = common::make_fp8x2<IDX>(rem_ri3);

        out_ri[0]    = common::concat(rem_ri0x2.x, rem_ri1x2.x, rem_ri2x2.x, rem_ri3x2.x);
        out_ri[next] = common::concat(rem_ri0x2.y, rem_ri1x2.y, rem_ri2x2.y, rem_ri3x2.y);

    } else {

        const int32_t rem_r0 = calc_mod<Backend::FP8, IDX>(v0.x);
        const int32_t rem_r1 = calc_mod<Backend::FP8, IDX>(v1.x);
        const int32_t rem_r2 = calc_mod<Backend::FP8, IDX>(v2.x);
        const int32_t rem_r3 = calc_mod<Backend::FP8, IDX>(v3.x);

        const common::fp8x3_e4m3 rem_r0x3 = common::make_fp8x3<IDX>(rem_r0);
        const common::fp8x3_e4m3 rem_r1x3 = common::make_fp8x3<IDX>(rem_r1);
        const common::fp8x3_e4m3 rem_r2x3 = common::make_fp8x3<IDX>(rem_r2);
        const common::fp8x3_e4m3 rem_r3x3 = common::make_fp8x3<IDX>(rem_r3);

        out_r[0]        = common::concat(rem_r0x3.x, rem_r1x3.x, rem_r2x3.x, rem_r3x3.x);
        out_r[next]     = common::concat(rem_r0x3.y, rem_r1x3.y, rem_r2x3.y, rem_r3x3.y);
        out_r[next * 2] = common::concat(rem_r0x3.z, rem_r1x3.z, rem_r2x3.z, rem_r3x3.z);

        const int32_t rem_i0 = calc_mod<Backend::FP8, IDX>(v0.y);
        const int32_t rem_i1 = calc_mod<Backend::FP8, IDX>(v1.y);
        const int32_t rem_i2 = calc_mod<Backend::FP8, IDX>(v2.y);
        const int32_t rem_i3 = calc_mod<Backend::FP8, IDX>(v3.y);

        const common::fp8x3_e4m3 rem_i0x3 = common::make_fp8x3<IDX>(rem_i0);
        const common::fp8x3_e4m3 rem_i1x3 = common::make_fp8x3<IDX>(rem_i1);
        const common::fp8x3_e4m3 rem_i2x3 = common::make_fp8x3<IDX>(rem_i2);
        const common::fp8x3_e4m3 rem_i3x3 = common::make_fp8x3<IDX>(rem_i3);

        out_i[0]        = common::concat(rem_i0x3.x, rem_i1x3.x, rem_i2x3.x, rem_i3x3.x);
        out_i[next]     = common::concat(rem_i0x3.y, rem_i1x3.y, rem_i2x3.y, rem_i3x3.y);
        out_i[next * 2] = common::concat(rem_i0x3.z, rem_i1x3.z, rem_i2x3.z, rem_i3x3.z);

        const int32_t rem_ri0 = wrapping<Backend::FP8, IDX>(rem_r0 + rem_i0);
        const int32_t rem_ri1 = wrapping<Backend::FP8, IDX>(rem_r1 + rem_i1);
        const int32_t rem_ri2 = wrapping<Backend::FP8, IDX>(rem_r2 + rem_i2);
        const int32_t rem_ri3 = wrapping<Backend::FP8, IDX>(rem_r3 + rem_i3);

        const common::fp8x3_e4m3 rem_ri0x3 = common::make_fp8x3<IDX>(rem_ri0);
        const common::fp8x3_e4m3 rem_ri1x3 = common::make_fp8x3<IDX>(rem_ri1);
        const common::fp8x3_e4m3 rem_ri2x3 = common::make_fp8x3<IDX>(rem_ri2);
        const common::fp8x3_e4m3 rem_ri3x3 = common::make_fp8x3<IDX>(rem_ri3);

        out_ri[0]        = common::concat(rem_ri0x3.x, rem_ri1x3.x, rem_ri2x3.x, rem_ri3x3.x);
        out_ri[next]     = common::concat(rem_ri0x3.y, rem_ri1x3.y, rem_ri2x3.y, rem_ri3x3.y);
        out_ri[next * 2] = common::concat(rem_ri0x3.z, rem_ri1x3.z, rem_ri2x3.z, rem_ri3x3.z);
    }
}

} // namespace gemmul8::mod
