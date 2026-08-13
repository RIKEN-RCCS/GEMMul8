#pragma once
#include "../common/common.hpp"
#include "mod_launch_i8.hpp"
#include "mod_launch_f8.hpp"
#include "../common/make_fp_mant_exp.hpp"

namespace gemmul8::mod {

//------------------------------
// Interface
//------------------------------
#define GEMMUL8_FOR_EACH_IDX_20(M) \
    M(0U)                          \
    M(1U)                          \
    M(2U)                          \
    M(3U)                          \
    M(4U)                          \
    M(5U)                          \
    M(6U)                          \
    M(7U)                          \
    M(8U)                          \
    M(9U)                          \
    M(10U)                         \
    M(11U)                         \
    M(12U)                         \
    M(13U)                         \
    M(14U)                         \
    M(15U)                         \
    M(16U)                         \
    M(17U)                         \
    M(18U)                         \
    M(19U)

#define GEMMUL8_INT8_RUN_SCALAR_STEP(I) \
    if constexpr (NUM_MODULI > (I)) {   \
        mod_launch<(I), V>(out, v);     \
        out += inc;                     \
    }

#define GEMMUL8_INT8_RUN_VEC4_STEP(I)            \
    if constexpr (NUM_MODULI > (I)) {            \
        mod_launch<(I), V>(out, v0, v1, v2, v3); \
        out += inc;                              \
    }

#define GEMMUL8_INT8_RUN_CPLX_SCALAR_STEP(I)     \
    if constexpr (NUM_MODULI > (I)) {            \
        mod_launch<(I), V>(out0, out1, out2, v); \
        out0 += inc;                             \
        out1 += inc;                             \
        out2 += inc;                             \
    }

#define GEMMUL8_INT8_RUN_CPLX_VEC4_STEP(I)                    \
    if constexpr (NUM_MODULI > (I)) {                         \
        mod_launch<(I), V>(out0, out1, out2, v0, v1, v2, v3); \
        out0 += inc;                                          \
        out1 += inc;                                          \
        out2 += inc;                                          \
    }

#define GEMMUL8_FP8_RUN_SCALAR_STEP(I)                              \
    if constexpr (NUM_MODULI > (I)) {                               \
        mod_launch<(I), V>(out, inc, v);                            \
        out += ((common::table::sqrt_moduli<I> == 0) ? 3 : 2) * inc; \
    }

#define GEMMUL8_FP8_RUN_VEC4_STEP(I)                                \
    if constexpr (NUM_MODULI > (I)) {                               \
        mod_launch<(I), V>(out, inc, v0, v1, v2, v3);               \
        out += ((common::table::sqrt_moduli<I> == 0) ? 3 : 2) * inc; \
    }

#define GEMMUL8_FP8_RUN_CPLX_SCALAR_STEP(I)                                      \
    if constexpr (NUM_MODULI > (I)) {                                            \
        mod_launch<(I), V>(out0, out1, out2, inc, v);                            \
        const size_t step = ((common::table::sqrt_moduli<I> == 0) ? 3 : 2) * inc; \
        out0 += step;                                                            \
        out1 += step;                                                            \
        out2 += step;                                                            \
    }

#define GEMMUL8_FP8_RUN_CPLX_VEC4_STEP(I)                                        \
    if constexpr (NUM_MODULI > (I)) {                                            \
        mod_launch<(I), V>(out0, out1, out2, inc, v0, v1, v2, v3);               \
        const size_t step = ((common::table::sqrt_moduli<I> == 0) ? 3 : 2) * inc; \
        out0 += step;                                                            \
        out1 += step;                                                            \
        out2 += step;                                                            \
    }

// interface for general NUM_MODULI
template <unsigned NUM_MODULI, typename V> struct ModUnroll {

    //=====
    // INT8
    //=====
    __device__ __forceinline__ static void run(
        int8_t *__restrict__ out,
        size_t inc,
        V v //
    ) {
        if constexpr (common::isComplex<V>) return;
        GEMMUL8_FOR_EACH_IDX_20(GEMMUL8_INT8_RUN_SCALAR_STEP);
    }

    __device__ __forceinline__ static void run(
        char4 *__restrict__ out,
        size_t inc,
        V v0, V v1, V v2, V v3 //
    ) {
        if constexpr (common::isComplex<V>) return;
        GEMMUL8_FOR_EACH_IDX_20(GEMMUL8_INT8_RUN_VEC4_STEP);
    }

    __device__ __forceinline__ static void run(
        int8_t *__restrict__ out0,
        int8_t *__restrict__ out1,
        int8_t *__restrict__ out2,
        size_t inc,
        V v //
    ) {
        if constexpr (!common::isComplex<V>) return;
        GEMMUL8_FOR_EACH_IDX_20(GEMMUL8_INT8_RUN_CPLX_SCALAR_STEP);
    }

    __device__ __forceinline__ static void run(
        char4 *__restrict__ out0,
        char4 *__restrict__ out1,
        char4 *__restrict__ out2,
        size_t inc,
        V v0, V v1, V v2, V v3 //
    ) {
        if constexpr (!common::isComplex<V>) return;
        GEMMUL8_FOR_EACH_IDX_20(GEMMUL8_INT8_RUN_CPLX_VEC4_STEP);
    }

    //=====
    // FP8
    //=====
    __device__ __forceinline__ static void run(
        __nv_fp8_e4m3 *__restrict__ out,
        size_t inc,
        V v //
    ) {
        if constexpr (common::isComplex<V>) return;
        GEMMUL8_FOR_EACH_IDX_20(GEMMUL8_FP8_RUN_SCALAR_STEP);
    }

    __device__ __forceinline__ static void run(
        __nv_fp8x4_e4m3 *__restrict__ out,
        size_t inc,
        V v0, V v1, V v2, V v3 //
    ) {
        if constexpr (common::isComplex<V>) return;
        GEMMUL8_FOR_EACH_IDX_20(GEMMUL8_FP8_RUN_VEC4_STEP);
    }

    __device__ __forceinline__ static void run(
        __nv_fp8_e4m3 *__restrict__ out0,
        __nv_fp8_e4m3 *__restrict__ out1,
        __nv_fp8_e4m3 *__restrict__ out2,
        size_t inc,
        V v //
    ) {
        if constexpr (!common::isComplex<V>) return;
        GEMMUL8_FOR_EACH_IDX_20(GEMMUL8_FP8_RUN_CPLX_SCALAR_STEP);
    }

    __device__ __forceinline__ static void run(
        __nv_fp8x4_e4m3 *__restrict__ out0,
        __nv_fp8x4_e4m3 *__restrict__ out1,
        __nv_fp8x4_e4m3 *__restrict__ out2,
        size_t inc,
        V v0, V v1, V v2, V v3 //
    ) {
        if constexpr (!common::isComplex<V>) return;
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

// interface for NUM_MODULI > common::threshold::M & V = float
template <unsigned NUM_MODULI> struct ModUnroll<NUM_MODULI, float> {

    //=====
    // INT8
    //=====
    __device__ __forceinline__ static void run(
        int8_t *__restrict__ out,
        size_t inc,
        float v //
    ) {
        common::fp32_mant_exp d = common::make_fp32_mant_exp(v);
        ModUnroll<NUM_MODULI, common::fp32_mant_exp>::run(out, inc, d);
    }

    __device__ __forceinline__ static void run(
        char4 *__restrict__ out,
        size_t inc,
        float v0, float v1, float v2, float v3 //
    ) {
        common::fp32_mant_exp d0 = common::make_fp32_mant_exp(v0);
        common::fp32_mant_exp d1 = common::make_fp32_mant_exp(v1);
        common::fp32_mant_exp d2 = common::make_fp32_mant_exp(v2);
        common::fp32_mant_exp d3 = common::make_fp32_mant_exp(v3);
        ModUnroll<NUM_MODULI, common::fp32_mant_exp>::run(out, inc, d0, d1, d2, d3);
    }

    //=====
    // FP8
    //=====
    __device__ __forceinline__ static void run(
        __nv_fp8_e4m3 *__restrict__ out,
        size_t inc,
        float v //
    ) {
        common::fp32_mant_exp d = common::make_fp32_mant_exp(v);
        ModUnroll<NUM_MODULI, common::fp32_mant_exp>::run(out, inc, d);
    }

    __device__ __forceinline__ static void run(
        __nv_fp8x4_e4m3 *__restrict__ out,
        size_t inc,
        float v0, float v1, float v2, float v3 //
    ) {
        common::fp32_mant_exp d0 = common::make_fp32_mant_exp(v0);
        common::fp32_mant_exp d1 = common::make_fp32_mant_exp(v1);
        common::fp32_mant_exp d2 = common::make_fp32_mant_exp(v2);
        common::fp32_mant_exp d3 = common::make_fp32_mant_exp(v3);
        ModUnroll<NUM_MODULI, common::fp32_mant_exp>::run(out, inc, d0, d1, d2, d3);
    }
};

// interface for NUM_MODULI > common::threshold::M & V = double
template <unsigned NUM_MODULI> struct ModUnroll<NUM_MODULI, double> {

    //=====
    // INT8
    //=====
    __device__ __forceinline__ static void run(
        int8_t *__restrict__ out,
        size_t inc,
        double v //
    ) {
        using VT = common::fp64_mant_exp;

        VT d = common::make_fp64_mant_exp<VT>(v);
        ModUnroll<NUM_MODULI, VT>::run(out, inc, d);
    }

    __device__ __forceinline__ static void run(
        char4 *__restrict__ out,
        size_t inc,
        double v0, double v1, double v2, double v3 //
    ) {
        using VT = common::fp64_mant_exp;

        VT d0 = common::make_fp64_mant_exp<VT>(v0);
        VT d1 = common::make_fp64_mant_exp<VT>(v1);
        VT d2 = common::make_fp64_mant_exp<VT>(v2);
        VT d3 = common::make_fp64_mant_exp<VT>(v3);
        ModUnroll<NUM_MODULI, VT>::run(out, inc, d0, d1, d2, d3);
    }

    //=====
    // FP8
    //=====
    __device__ __forceinline__ static void run(
        __nv_fp8_e4m3 *__restrict__ out,
        size_t inc,
        double v //
    ) {
        using VT = common::fp64_mant_exp;

        VT d = common::make_fp64_mant_exp<VT>(v);
        ModUnroll<NUM_MODULI, VT>::run(out, inc, d);
    }

    __device__ __forceinline__ static void run(
        __nv_fp8x4_e4m3 *__restrict__ out,
        size_t inc,
        double v0, double v1, double v2, double v3 //
    ) {
        using VT = common::fp64_mant_exp;

        VT d0 = common::make_fp64_mant_exp<VT>(v0);
        VT d1 = common::make_fp64_mant_exp<VT>(v1);
        VT d2 = common::make_fp64_mant_exp<VT>(v2);
        VT d3 = common::make_fp64_mant_exp<VT>(v3);
        ModUnroll<NUM_MODULI, VT>::run(out, inc, d0, d1, d2, d3);
    }
};

// interface for NUM_MODULI > common::threshold::M & V = cuFloatComplex
template <unsigned NUM_MODULI> struct ModUnroll<NUM_MODULI, cuFloatComplex> {

    //=====
    // INT8
    //=====
    __device__ __forceinline__ static void run(
        int8_t *__restrict__ out0,
        int8_t *__restrict__ out1,
        int8_t *__restrict__ out2,
        size_t inc,
        cuFloatComplex v //
    ) {
        common::fp32_mant_exp2 d = common::make_fp32_mant_exp2(v);
        ModUnroll<NUM_MODULI, common::fp32_mant_exp2>::run(out0, out1, out2, inc, d);
    }

    __device__ __forceinline__ static void run(
        char4 *__restrict__ out0,
        char4 *__restrict__ out1,
        char4 *__restrict__ out2,
        size_t inc,
        cuFloatComplex v0,
        cuFloatComplex v1,
        cuFloatComplex v2,
        cuFloatComplex v3 //
    ) {
        common::fp32_mant_exp2 d0 = common::make_fp32_mant_exp2(v0);
        common::fp32_mant_exp2 d1 = common::make_fp32_mant_exp2(v1);
        common::fp32_mant_exp2 d2 = common::make_fp32_mant_exp2(v2);
        common::fp32_mant_exp2 d3 = common::make_fp32_mant_exp2(v3);
        ModUnroll<NUM_MODULI, common::fp32_mant_exp2>::run(out0, out1, out2, inc, d0, d1, d2, d3);
    }

    //=====
    // FP8
    //=====
    __device__ __forceinline__ static void run(
        __nv_fp8_e4m3 *__restrict__ out0,
        __nv_fp8_e4m3 *__restrict__ out1,
        __nv_fp8_e4m3 *__restrict__ out2,
        size_t inc,
        cuFloatComplex v //
    ) {
        common::fp32_mant_exp2 d = common::make_fp32_mant_exp2(v);
        ModUnroll<NUM_MODULI, common::fp32_mant_exp2>::run(out0, out1, out2, inc, d);
    }

    __device__ __forceinline__ static void run(
        __nv_fp8x4_e4m3 *__restrict__ out0,
        __nv_fp8x4_e4m3 *__restrict__ out1,
        __nv_fp8x4_e4m3 *__restrict__ out2,
        size_t inc,
        cuFloatComplex v0,
        cuFloatComplex v1,
        cuFloatComplex v2,
        cuFloatComplex v3 //
    ) {
        common::fp32_mant_exp2 d0 = common::make_fp32_mant_exp2(v0);
        common::fp32_mant_exp2 d1 = common::make_fp32_mant_exp2(v1);
        common::fp32_mant_exp2 d2 = common::make_fp32_mant_exp2(v2);
        common::fp32_mant_exp2 d3 = common::make_fp32_mant_exp2(v3);
        ModUnroll<NUM_MODULI, common::fp32_mant_exp2>::run(out0, out1, out2, inc, d0, d1, d2, d3);
    }
};

// interface for NUM_MODULI > common::threshold::M & V = cuDoubleComplex
template <unsigned NUM_MODULI> struct ModUnroll<NUM_MODULI, cuDoubleComplex> {

    //=====
    // INT8
    //=====
    __device__ __forceinline__ static void run(
        int8_t *__restrict__ out0,
        int8_t *__restrict__ out1,
        int8_t *__restrict__ out2,
        size_t inc,
        cuDoubleComplex v //
    ) {
        using VT = common::fp64_mant_exp2;

        VT d = common::make_fp64_mant_exp2<VT>(v);
        ModUnroll<NUM_MODULI, VT>::run(out0, out1, out2, inc, d);
    }

    __device__ __forceinline__ static void run(
        char4 *__restrict__ out0,
        char4 *__restrict__ out1,
        char4 *__restrict__ out2,
        size_t inc,
        cuDoubleComplex v0,
        cuDoubleComplex v1,
        cuDoubleComplex v2,
        cuDoubleComplex v3 //
    ) {
        using VT = common::fp64_mant_exp2;

        VT d0 = common::make_fp64_mant_exp2<VT>(v0);
        VT d1 = common::make_fp64_mant_exp2<VT>(v1);
        VT d2 = common::make_fp64_mant_exp2<VT>(v2);
        VT d3 = common::make_fp64_mant_exp2<VT>(v3);
        ModUnroll<NUM_MODULI, VT>::run(out0, out1, out2, inc, d0, d1, d2, d3);
    }

    //=====
    // FP8
    //=====
    __device__ __forceinline__ static void run(
        __nv_fp8_e4m3 *__restrict__ out0,
        __nv_fp8_e4m3 *__restrict__ out1,
        __nv_fp8_e4m3 *__restrict__ out2,
        size_t inc, cuDoubleComplex v //
    ) {
        using VT = common::fp64_mant_exp2;

        VT d = common::make_fp64_mant_exp2<VT>(v);
        ModUnroll<NUM_MODULI, VT>::run(out0, out1, out2, inc, d);
    }

    __device__ __forceinline__ static void run(
        __nv_fp8x4_e4m3 *__restrict__ out0,
        __nv_fp8x4_e4m3 *__restrict__ out1,
        __nv_fp8x4_e4m3 *__restrict__ out2,
        size_t inc,
        cuDoubleComplex v0,
        cuDoubleComplex v1,
        cuDoubleComplex v2,
        cuDoubleComplex v3 //
    ) {
        using VT = common::fp64_mant_exp2;

        VT d0 = common::make_fp64_mant_exp2<VT>(v0);
        VT d1 = common::make_fp64_mant_exp2<VT>(v1);
        VT d2 = common::make_fp64_mant_exp2<VT>(v2);
        VT d3 = common::make_fp64_mant_exp2<VT>(v3);
        ModUnroll<NUM_MODULI, VT>::run(out0, out1, out2, inc, d0, d1, d2, d3);
    }
};

template <Backend BACKEND, unsigned NUM_MODULI, typename Out>
struct ModUnrollFillZero {

    __device__ __forceinline__ static void run(
        Out *__restrict__ out,
        size_t inc //
    ) {
        constexpr unsigned num_mat = common::table::num_mat_v<BACKEND, NUM_MODULI>;
        const Out zero             = common::Tconst<Out>::zero();

#pragma unroll
        for (unsigned i = 0; i < num_mat; ++i) {
            out[i * inc] = zero;
        }
    }

    __device__ __forceinline__ static void run(
        Out *__restrict__ out0,
        Out *__restrict__ out1,
        Out *__restrict__ out2,
        size_t inc //
    ) {
        constexpr unsigned num_mat = common::table::num_mat_v<BACKEND, NUM_MODULI>;
        const Out zero             = common::Tconst<Out>::zero();

#pragma unroll
        for (unsigned i = 0; i < num_mat; ++i) {
            out0[i * inc] = zero;
            out1[i * inc] = zero;
            out2[i * inc] = zero;
        }
    }
};

} // namespace gemmul8::mod
