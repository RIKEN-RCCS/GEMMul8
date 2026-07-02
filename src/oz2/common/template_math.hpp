#pragma once

namespace gemmul8::common {

//------------------------------
// abs(x)
//------------------------------
template <typename T> __device__ __forceinline__ T Tabs(T in);
template <> __device__ __forceinline__ double Tabs<double>(double in) { return fabs(in); }
template <> __device__ __forceinline__ float Tabs<float>(float in) { return fabsf(in); }
template <> __device__ __forceinline__ int32_t Tabs<int32_t>(int32_t in) { return abs(in); }
template <> __device__ __forceinline__ cuDoubleComplex Tabs<cuDoubleComplex>(cuDoubleComplex in) { return cuDoubleComplex{fabs(in.x), fabs(in.y)}; }
template <> __device__ __forceinline__ cuFloatComplex Tabs<cuFloatComplex>(cuFloatComplex in) { return cuFloatComplex{fabsf(in.x), fabsf(in.y)}; }

//------------------------------
// x+y
//------------------------------
template <typename T> __device__ __forceinline__ T Tadd(T a, T b) { return a + b; }
template <> __device__ __forceinline__ cuDoubleComplex Tadd<cuDoubleComplex>(cuDoubleComplex a, cuDoubleComplex b) { return cuDoubleComplex{a.x + b.x, a.y + b.y}; }
template <> __device__ __forceinline__ cuFloatComplex Tadd<cuFloatComplex>(cuFloatComplex a, cuFloatComplex b) { return cuFloatComplex{a.x + b.x, a.y + b.y}; }

//------------------------------
// x-y
//------------------------------
template <typename T> __device__ __forceinline__ T Tsub(T a, T b) { return a - b; }
template <> __device__ __forceinline__ cuDoubleComplex Tsub<cuDoubleComplex>(cuDoubleComplex a, cuDoubleComplex b) { return cuDoubleComplex{a.x - b.x, a.y - b.y}; }
template <> __device__ __forceinline__ cuFloatComplex Tsub<cuFloatComplex>(cuFloatComplex a, cuFloatComplex b) { return cuFloatComplex{a.x - b.x, a.y - b.y}; }

//------------------------------
// x*y
//------------------------------
template <typename T1, typename T2 = T1> __device__ __forceinline__ T2 Tmul(T1 a, T2 b) { return a * b; }
template <> __device__ __forceinline__ cuDoubleComplex Tmul<cuDoubleComplex, cuDoubleComplex>(cuDoubleComplex a, cuDoubleComplex b) {
    cuDoubleComplex prod;
    prod.x = fma(-a.y, b.y, a.x * b.x);
    prod.y = fma(a.y, b.x, a.x * b.y);
    return prod;
}
template <> __device__ __forceinline__ cuDoubleComplex Tmul<double, cuDoubleComplex>(double a, cuDoubleComplex b) { return cuDoubleComplex{a * b.x, a * b.y}; }
template <> __device__ __forceinline__ cuDoubleComplex Tmul<float, cuDoubleComplex>(float a, cuDoubleComplex b) { return Tmul<double, cuDoubleComplex>(double(a), b); }
template <> __device__ __forceinline__ cuFloatComplex Tmul<cuFloatComplex, cuFloatComplex>(cuFloatComplex a, cuFloatComplex b) {
    cuFloatComplex prod;
    prod.x = fmaf(-a.y, b.y, a.x * b.x);
    prod.y = fmaf(a.y, b.x, a.x * b.y);
    return prod;
}
template <> __device__ __forceinline__ cuFloatComplex Tmul<float, cuFloatComplex>(float a, cuFloatComplex b) { return cuFloatComplex{a * b.x, a * b.y}; }
template <> __device__ __forceinline__ cuFloatComplex Tmul<double, cuFloatComplex>(double a, cuFloatComplex b) { return Tmul<float, cuFloatComplex>(double(a), b); }

//------------------------------
// -x
//------------------------------
template <typename T> __device__ __forceinline__ T Tneg(T a) { return -a; }
template <> __device__ __forceinline__ cuDoubleComplex Tneg<cuDoubleComplex>(cuDoubleComplex a) { return cuDoubleComplex{-a.x, -a.y}; }
template <> __device__ __forceinline__ cuFloatComplex Tneg<cuFloatComplex>(cuFloatComplex a) { return cuFloatComplex{-a.x, -a.y}; }

//------------------------------
// x^2 + y
//------------------------------
template <typename T> __device__ __forceinline__ underlying_t<T> Tsqr_add_ru(T in1, underlying_t<T> in2);
template <> __device__ __forceinline__ double Tsqr_add_ru<double>(double in1, double in2) { return __fma_ru(in1, in1, in2); }
template <> __device__ __forceinline__ float Tsqr_add_ru<float>(float in1, float in2) { return __fmaf_ru(in1, in1, in2); }
template <> __device__ __forceinline__ double Tsqr_add_ru<cuDoubleComplex>(cuDoubleComplex in1, double in2) { return __fma_ru(in1.y, in1.y, __fma_ru(in1.x, in1.x, in2)); }
template <> __device__ __forceinline__ float Tsqr_add_ru<cuFloatComplex>(cuFloatComplex in1, float in2) { return __fmaf_ru(in1.y, in1.y, __fmaf_ru(in1.x, in1.x, in2)); }

//------------------------------
// x+y in round-up mode
//------------------------------
template <typename T> __device__ __forceinline__ T __Tadd_ru(T in1, T in2);
template <> __device__ __forceinline__ double __Tadd_ru<double>(double in1, double in2) { return __dadd_ru(in1, in2); }
template <> __device__ __forceinline__ float __Tadd_ru<float>(float in1, float in2) { return __fadd_ru(in1, in2); }

//------------------------------
// a*x + b*y
//------------------------------
template <typename T, typename Ta, typename Tb> __device__ __forceinline__ T Taxpby(Ta a, T x, Tb b, T y);
template <> __device__ __forceinline__ double Taxpby<double, double, double>(double a, double x, double b, double y) {
    return fma(b, y, a * x);
}
template <> __device__ __forceinline__ float Taxpby<float, float, float>(float a, float x, float b, float y) {
    return fmaf(b, y, a * x);
}
template <> __device__ __forceinline__ cuDoubleComplex Taxpby<cuDoubleComplex, cuDoubleComplex, cuDoubleComplex>(
    cuDoubleComplex a, cuDoubleComplex x, cuDoubleComplex b, cuDoubleComplex y) {
    cuDoubleComplex out;
    out.x = fma(-b.y, y.y, fma(b.x, y.x, fma(-a.y, x.y, a.x * x.x)));
    out.y = fma(b.y, y.x, fma(b.x, y.y, fma(a.y, x.x, a.x * x.y)));
    return out;
}
template <> __device__ __forceinline__ cuFloatComplex Taxpby<cuFloatComplex, cuFloatComplex, cuFloatComplex>(
    cuFloatComplex a, cuFloatComplex x, cuFloatComplex b, cuFloatComplex y) {
    cuFloatComplex out;
    out.x = fmaf(-b.y, y.y, fmaf(b.x, y.x, fmaf(-a.y, x.y, a.x * x.x)));
    out.y = fmaf(b.y, y.x, fmaf(b.x, y.y, fmaf(a.y, x.x, a.x * x.y)));
    return out;
}
template <> __device__ __forceinline__ cuDoubleComplex Taxpby<cuDoubleComplex, double, double>(
    double a, cuDoubleComplex x, double b, cuDoubleComplex y) {
    cuDoubleComplex out;
    out.x = fma(b, y.x, a * x.x);
    out.y = fma(b, y.y, a * x.y);
    return out;
}
template <> __device__ __forceinline__ cuFloatComplex Taxpby<cuFloatComplex, float, float>(
    float a, cuFloatComplex x, float b, cuFloatComplex y) {
    cuFloatComplex out;
    out.x = fmaf(b, y.x, a * x.x);
    out.y = fmaf(b, y.y, a * x.y);
    return out;
}
template <> __device__ __forceinline__ cuDoubleComplex Taxpby<cuDoubleComplex, cuDoubleComplex, double>(
    cuDoubleComplex a, cuDoubleComplex x, double b, cuDoubleComplex y) {
    cuDoubleComplex out;
    out.x = fma(b, y.x, fma(-a.y, x.y, a.x * x.x));
    out.y = fma(b, y.y, fma(a.y, x.x, a.x * x.y));
    return out;
}
template <> __device__ __forceinline__ cuFloatComplex Taxpby<cuFloatComplex, cuFloatComplex, float>(
    cuFloatComplex a, cuFloatComplex x, float b, cuFloatComplex y) {
    cuFloatComplex out;
    out.x = fmaf(b, y.x, fmaf(-a.y, x.y, a.x * x.x));
    out.y = fmaf(b, y.y, fmaf(a.y, x.x, a.x * x.y));
    return out;
}

//------------------------------
// a*x + y
//------------------------------
template <typename T, typename Ta> __device__ __forceinline__ T Taxpy(Ta a, T x, T y);
template <> __device__ __forceinline__ double Taxpy<double, double>(double a, double x, double y) {
    return fma(a, x, y);
}
template <> __device__ __forceinline__ float Taxpy<float, float>(float a, float x, float y) {
    return fmaf(a, x, y);
}
template <> __device__ __forceinline__ cuDoubleComplex Taxpy<cuDoubleComplex, cuDoubleComplex>(
    cuDoubleComplex a, cuDoubleComplex x, cuDoubleComplex y) {
    cuDoubleComplex out;
    out.x = fma(-a.y, x.y, fma(a.x, x.x, y.x));
    out.y = fma(a.y, x.x, fma(a.x, x.y, y.y));
    return out;
}
template <> __device__ __forceinline__ cuFloatComplex Taxpy<cuFloatComplex, cuFloatComplex>(
    cuFloatComplex a, cuFloatComplex x, cuFloatComplex y) {
    cuFloatComplex out;
    out.x = fmaf(-a.y, x.y, fmaf(a.x, x.x, y.x));
    out.y = fmaf(a.y, x.x, fmaf(a.x, x.y, y.y));
    return out;
}
template <> __device__ __forceinline__ cuDoubleComplex Taxpy<cuDoubleComplex, double>(
    double a, cuDoubleComplex x, cuDoubleComplex y) {
    cuDoubleComplex out;
    out.x = fma(a, x.x, y.x);
    out.y = fma(a, x.y, y.y);
    return out;
}
template <> __device__ __forceinline__ cuFloatComplex Taxpy<cuFloatComplex, float>(
    float a, cuFloatComplex x, cuFloatComplex y) {
    cuFloatComplex out;
    out.x = fmaf(a, x.x, y.x);
    out.y = fmaf(a, x.y, y.y);
    return out;
}

//------------------------------
// x*2^y
//------------------------------
template <typename T> __device__ __forceinline__ T Tscalbn(T in, int32_t sft);
template <> __device__ __forceinline__ double Tscalbn<double>(double in, int32_t sft) { return scalbn(in, sft); }
template <> __device__ __forceinline__ float Tscalbn<float>(float in, int32_t sft) { return scalbnf(in, sft); }
template <> __device__ __forceinline__ cuDoubleComplex Tscalbn<cuDoubleComplex>(cuDoubleComplex in, int32_t sft) { return cuDoubleComplex{scalbn(in.x, sft), scalbn(in.y, sft)}; }
template <> __device__ __forceinline__ cuFloatComplex Tscalbn<cuFloatComplex>(cuFloatComplex in, int32_t sft) { return cuFloatComplex{scalbnf(in.x, sft), scalbnf(in.y, sft)}; }

//------------------------------
// rint(x) = round(x)
//------------------------------
template <typename T> __device__ __forceinline__ T Trint(T in);
template <> __device__ __forceinline__ double Trint<double>(double in) { return rint(in); }
template <> __device__ __forceinline__ cuDoubleComplex Trint<cuDoubleComplex>(cuDoubleComplex in) { return cuDoubleComplex{rint(in.x), rint(in.y)}; }

//------------------------------
// trunc(x)
//------------------------------
template <typename T> __device__ __forceinline__ T Ttrunc(T in);
template <> __device__ __forceinline__ double Ttrunc<double>(double in) { return trunc(in); }
template <> __device__ __forceinline__ float Ttrunc<float>(float in) { return truncf(in); }
template <> __device__ __forceinline__ cuDoubleComplex Ttrunc<cuDoubleComplex>(cuDoubleComplex in) { return cuDoubleComplex{trunc(in.x), trunc(in.y)}; }
template <> __device__ __forceinline__ cuFloatComplex Ttrunc<cuFloatComplex>(cuFloatComplex in) { return cuFloatComplex{truncf(in.x), truncf(in.y)}; }

//------------------------------
// ilogb(x) = floor(log2(x))
//------------------------------
template <typename T> __device__ __forceinline__ int32_t Tilogb(T in);
template <> __device__ __forceinline__ int32_t Tilogb<double>(double in) { return (in == 0.0) ? 0 : ilogb(in); }
template <> __device__ __forceinline__ int32_t Tilogb<float>(float in) { return (in == 0.0F) ? 0 : ilogbf(in); }

//------------------------------
// max(x,y)
//------------------------------
template <typename T> __device__ __forceinline__ underlying_t<T> Tmax(T in1, underlying_t<T> in2);
template <> __device__ __forceinline__ double Tmax<double>(double in1, double in2) { return max(in1, in2); }
template <> __device__ __forceinline__ float Tmax<float>(float in1, float in2) { return max(in1, in2); }
template <> __device__ __forceinline__ int32_t Tmax<int32_t>(int32_t in1, int32_t in2) { return max(in1, in2); }
template <> __device__ __forceinline__ double Tmax<cuDoubleComplex>(cuDoubleComplex in1, double in2) { return max(max(in1.x, in1.y), in2); }
template <> __device__ __forceinline__ float Tmax<cuFloatComplex>(cuFloatComplex in1, float in2) { return max(max(in1.x, in1.y), in2); }

//------------------------------
// Cast Tin to Tout
//------------------------------
template <typename Tin, typename Tout> __device__ __forceinline__ Tout Tcast(Tin in) { return Tout(in); };
template <> __device__ __forceinline__ double Tcast<double, double>(double in) { return in; }
template <> __device__ __forceinline__ float Tcast<double, float>(double in) { return __double2float_rn(in); }
template <> __device__ __forceinline__ cuDoubleComplex Tcast<cuDoubleComplex, cuDoubleComplex>(cuDoubleComplex in) { return in; }
template <> __device__ __forceinline__ cuFloatComplex Tcast<cuDoubleComplex, cuFloatComplex>(cuDoubleComplex in) { return cuFloatComplex{__double2float_rn(in.x), __double2float_rn(in.y)}; }
template <> __device__ __forceinline__ cuDoubleComplex Tcast<char2, cuDoubleComplex>(char2 in) { return cuDoubleComplex{double(in.x), double(in.y)}; }
template <> __device__ __forceinline__ cuDoubleComplex Tcast<short2, cuDoubleComplex>(short2 in) { return cuDoubleComplex{double(in.x), double(in.y)}; }
template <> __device__ __forceinline__ double Tcast<uint8_t, double>(uint8_t in) { return double(in); }

template <typename Tin> __device__ __forceinline__ float float_ru(Tin in);
template <> __device__ __forceinline__ float float_ru<double>(double in) { return __double2float_ru(in); };
template <> __device__ __forceinline__ float float_ru<float>(float in) { return in; };

//------------------------------
// static_cast (fp -> int)
//------------------------------
template <typename T> __device__ __forceinline__ int32_type<T> __fp2int_rz(T in);
template <> __device__ __forceinline__ int32_t __fp2int_rz<double>(double in) { return __double2int_rz(in); }
template <> __device__ __forceinline__ int32_t __fp2int_rz<float>(float in) { return __float2int_rz(in); }
template <> __device__ __forceinline__ int2 __fp2int_rz<cuDoubleComplex>(cuDoubleComplex in) { return int2{__double2int_rz(in.x), __double2int_rz(in.y)}; }
template <> __device__ __forceinline__ int2 __fp2int_rz<cuFloatComplex>(cuFloatComplex in) { return int2{__float2int_rz(in.x), __float2int_rz(in.y)}; }

template <typename T> __device__ __forceinline__ int64_type<T> __fp2ll_rz(T in);
template <> __device__ __forceinline__ int64_t __fp2ll_rz<double>(double in) { return __double2ll_rz(in); }
template <> __device__ __forceinline__ int64_t __fp2ll_rz<float>(float in) { return __float2ll_rz(in); }
template <> __device__ __forceinline__ int64x2_t __fp2ll_rz<cuDoubleComplex>(cuDoubleComplex in) { return int64x2_t{__double2ll_rz(in.x), __double2ll_rz(in.y)}; }
template <> __device__ __forceinline__ int64x2_t __fp2ll_rz<cuFloatComplex>(cuFloatComplex in) { return int64x2_t{__float2ll_rz(in.x), __float2ll_rz(in.y)}; }

template <typename T> __device__ __forceinline__ int32_type<T> __fp2int_ru(T in);
template <> __device__ __forceinline__ int32_t __fp2int_ru<double>(double in) { return __double2int_ru(in); }
template <> __device__ __forceinline__ int32_t __fp2int_ru<float>(float in) { return __float2int_ru(in); }
template <> __device__ __forceinline__ int2 __fp2int_ru<cuDoubleComplex>(cuDoubleComplex in) { return int2{__double2int_ru(in.x), __double2int_ru(in.y)}; }
template <> __device__ __forceinline__ int2 __fp2int_ru<cuFloatComplex>(cuFloatComplex in) { return int2{__float2int_ru(in.x), __float2int_ru(in.y)}; }

//------------------------------
// Retrieves the complex conjugate of a complex number.
//------------------------------
template <typename T, bool CONJ> __device__ __forceinline__ T conj(T in) { return in; };
template <> __device__ __forceinline__ cuDoubleComplex conj<cuDoubleComplex, true>(cuDoubleComplex in) { return cuDoubleComplex{in.x, -in.y}; };
template <> __device__ __forceinline__ cuFloatComplex conj<cuFloatComplex, true>(cuFloatComplex in) { return cuFloatComplex{in.x, -in.y}; };

//------------------------------
// Warp reduction (max)
//------------------------------
template <typename T, int32_t width = 32> __device__ __forceinline__ T inner_warp_max(T amax) {
#pragma unroll
    for (unsigned i = width / 2; i > 0; i >>= 1) {
        amax = max(amax, __shfl_down_sync(0xFFFFFFFFu, amax, i, width));
    }
    return amax;
}

//------------------------------
// Warp reduction (sum in round-up mode)
//------------------------------
template <typename T, int32_t width = 32> __device__ __forceinline__ T inner_warp_sum(T sum) {
#pragma unroll
    for (unsigned i = width / 2; i > 0; i >>= 1) {
        sum = __Tadd_ru<T>(sum, __shfl_down_sync(0xFFFFFFFFu, sum, i, width));
    }
    return sum;
}
template <> __device__ __forceinline__ double inner_warp_sum<double, 32>(double sum) {
    sum = __dadd_ru(sum, __shfl_down_sync(0xFFFFFFFFu, sum, 16u, 32));
    sum = __dadd_ru(sum, __shfl_down_sync(0xFFFFFFFFu, sum, 8u, 32));
    sum = __dadd_ru(sum, __shfl_down_sync(0xFFFFFFFFu, sum, 4u, 32));
    sum = __dadd_ru(sum, __shfl_down_sync(0xFFFFFFFFu, sum, 2u, 32));
    sum = __dadd_ru(sum, __shfl_down_sync(0xFFFFFFFFu, sum, 1u, 32));
    return sum;
}
template <> __device__ __forceinline__ float inner_warp_sum<float, 32>(float sum) {
    sum = __fadd_ru(sum, __shfl_down_sync(0xFFFFFFFFu, sum, 16u, 32));
    sum = __fadd_ru(sum, __shfl_down_sync(0xFFFFFFFFu, sum, 8u, 32));
    sum = __fadd_ru(sum, __shfl_down_sync(0xFFFFFFFFu, sum, 4u, 32));
    sum = __fadd_ru(sum, __shfl_down_sync(0xFFFFFFFFu, sum, 2u, 32));
    sum = __fadd_ru(sum, __shfl_down_sync(0xFFFFFFFFu, sum, 1u, 32));
    return sum;
}

} // namespace gemmul8::common
