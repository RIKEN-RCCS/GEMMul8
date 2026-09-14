#pragma once
#include "../common/common.hpp"
#include "../common/make_f8.hpp"
#include "complex_2m.hpp"

namespace gemmul8::mod {

template <unsigned IDX, bool COMPLEX = false> __device__ __forceinline__ void store_fp8_residue(
    __nv_fp8_e4m3 *__restrict__ out,
    size_t next,
    int32_t v //
) {
    if constexpr (common::table::kara_type<IDX, COMPLEX> == common::table::KaratsubaType::NONE) {

        const common::fp8x2_e4m3 r = common::make_fp8x2<IDX, COMPLEX>(v);

        out[0]    = r.x;
        out[next] = r.y;

    } else {

        const common::fp8x3_e4m3 rem = common::make_fp8x3<IDX, COMPLEX>(v);

        out[0]        = rem.x;
        out[next]     = rem.y;
        out[next * 2] = rem.z;
    }
}

// launcher for V in {int32_t, int64_t, float, double}
template <unsigned IDX, bool COMPLEX = false> __device__ __forceinline__ void store_fp8_residue(
    __nv_fp8x4_e4m3 *__restrict__ out,
    size_t next,
    int32_t v0, int32_t v1, int32_t v2, int32_t v3 //
) {
#if defined(__CUDACC__) && !defined(__HIPCC__)
    __nv_fp8x4_e4m3 x, y, z;
    if constexpr (common::table::kara_type<IDX, COMPLEX> == common::table::KaratsubaType::NONE) {
        common::make_fp8x2_vec4<IDX, COMPLEX>(v0, v1, v2, v3, x, y);
    } else {
        common::make_fp8x3_vec4<IDX, COMPLEX>(v0, v1, v2, v3, x, y, z);
        out[next * 2] = z;
    }
    out[0]    = x;
    out[next] = y;
#else
    if constexpr (common::table::kara_type<IDX, COMPLEX> == common::table::KaratsubaType::NONE) {

        const common::fp8x2_e4m3 rem0 = common::make_fp8x2<IDX, COMPLEX>(v0);
        const common::fp8x2_e4m3 rem1 = common::make_fp8x2<IDX, COMPLEX>(v1);
        const common::fp8x2_e4m3 rem2 = common::make_fp8x2<IDX, COMPLEX>(v2);
        const common::fp8x2_e4m3 rem3 = common::make_fp8x2<IDX, COMPLEX>(v3);

        out[0]    = common::concat(rem0.x, rem1.x, rem2.x, rem3.x);
        out[next] = common::concat(rem0.y, rem1.y, rem2.y, rem3.y);

    } else {

        const common::fp8x3_e4m3 r0 = common::make_fp8x3<IDX, COMPLEX>(v0);
        const common::fp8x3_e4m3 r1 = common::make_fp8x3<IDX, COMPLEX>(v1);
        const common::fp8x3_e4m3 r2 = common::make_fp8x3<IDX, COMPLEX>(v2);
        const common::fp8x3_e4m3 r3 = common::make_fp8x3<IDX, COMPLEX>(v3);

        out[0]        = common::concat(r0.x, r1.x, r2.x, r3.x);
        out[next]     = common::concat(r0.y, r1.y, r2.y, r3.y);
        out[next * 2] = common::concat(r0.z, r1.z, r2.z, r3.z);
    }
#endif
}

template <unsigned IDX, typename V> __device__ __forceinline__ void mod_launch(
    __nv_fp8_e4m3 *__restrict__ out, size_t next, V v) {
    store_fp8_residue<IDX>(out, next, calc_mod<Backend::FP8, IDX>(v));
}

template <unsigned IDX, typename V> __device__ __forceinline__ void mod_launch(
    __nv_fp8x4_e4m3 *__restrict__ out, size_t next, V v0, V v1, V v2, V v3) {
    store_fp8_residue<IDX>(out, next,
                           calc_mod<Backend::FP8, IDX>(v0), calc_mod<Backend::FP8, IDX>(v1),
                           calc_mod<Backend::FP8, IDX>(v2), calc_mod<Backend::FP8, IDX>(v3));
}

template <unsigned IDX, typename V> __device__ __forceinline__ void mod_launch(
    __nv_fp8_e4m3 *__restrict__ out_plus, __nv_fp8_e4m3 *__restrict__ out_minus,
    size_t next, V v) {
    const int2 r = project_complex_2m<Backend::FP8, IDX>(v);
    store_fp8_residue<IDX, true>(out_plus, next, r.x);
    store_fp8_residue<IDX, true>(out_minus, next, r.y);
}

template <unsigned IDX, typename V> __device__ __forceinline__ void mod_launch(
    __nv_fp8x4_e4m3 *__restrict__ out_plus, __nv_fp8x4_e4m3 *__restrict__ out_minus,
    size_t next, V v0, V v1, V v2, V v3) {
    const int2 r0 = project_complex_2m<Backend::FP8, IDX>(v0);
    const int2 r1 = project_complex_2m<Backend::FP8, IDX>(v1);
    const int2 r2 = project_complex_2m<Backend::FP8, IDX>(v2);
    const int2 r3 = project_complex_2m<Backend::FP8, IDX>(v3);
    store_fp8_residue<IDX, true>(out_plus, next, r0.x, r1.x, r2.x, r3.x);
    store_fp8_residue<IDX, true>(out_minus, next, r0.y, r1.y, r2.y, r3.y);
}

} // namespace gemmul8::mod
