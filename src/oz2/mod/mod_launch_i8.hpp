#pragma once
#include "../common/common.hpp"
#include "complex_2m.hpp"

namespace gemmul8::mod {

// launcher for V in {int32_t, int64_t, float, double}
template <unsigned IDX, typename V> __device__ __forceinline__ void mod_launch(
    int8_t *__restrict__ out,
    V v //
) {
    out[0] = static_cast<int8_t>(calc_mod<Backend::INT8, IDX>(v));
}

// launcher for V in {int32_t, int64_t, float, double}
template <unsigned IDX, typename V> __device__ __forceinline__ void mod_launch(
    char4 *__restrict__ out,
    V v0, V v1, V v2, V v3 //
) {
    char4 rem;
    rem.x  = static_cast<int8_t>(calc_mod<Backend::INT8, IDX>(v0));
    rem.y  = static_cast<int8_t>(calc_mod<Backend::INT8, IDX>(v1));
    rem.z  = static_cast<int8_t>(calc_mod<Backend::INT8, IDX>(v2));
    rem.w  = static_cast<int8_t>(calc_mod<Backend::INT8, IDX>(v3));
    out[0] = rem;
}

// The two streams are the centered 2M residues, not real/imaginary parts.
template <unsigned IDX, typename V> __device__ __forceinline__ void mod_launch(
    int8_t *__restrict__ out_plus, int8_t *__restrict__ out_minus, V v) {
    const int2 r = project_complex_2m<Backend::INT8, IDX>(v);
    out_plus[0]  = static_cast<int8_t>(r.x);
    out_minus[0] = static_cast<int8_t>(r.y);
}

template <unsigned IDX, typename V> __device__ __forceinline__ void mod_launch(
    char4 *__restrict__ out_plus, char4 *__restrict__ out_minus,
    V v0, V v1, V v2, V v3) {
    const int2 r0 = project_complex_2m<Backend::INT8, IDX>(v0);
    const int2 r1 = project_complex_2m<Backend::INT8, IDX>(v1);
    const int2 r2 = project_complex_2m<Backend::INT8, IDX>(v2);
    const int2 r3 = project_complex_2m<Backend::INT8, IDX>(v3);
    out_plus[0]   = char4{static_cast<int8_t>(r0.x), static_cast<int8_t>(r1.x),
                          static_cast<int8_t>(r2.x), static_cast<int8_t>(r3.x)};
    out_minus[0]  = char4{static_cast<int8_t>(r0.y), static_cast<int8_t>(r1.y),
                          static_cast<int8_t>(r2.y), static_cast<int8_t>(r3.y)};
}

} // namespace gemmul8::mod
