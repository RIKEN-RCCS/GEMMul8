#pragma once

//------------------------------
// abs(x)
//------------------------------
template <typename T> __forceinline__ __device__ T Tabs(T in);
template <> __forceinline__ __device__ double Tabs<double>(double in) { return fabs(in); }
template <> __forceinline__ __device__ float Tabs<float>(float in) { return fabsf(in); }
template <> __forceinline__ __device__ int32_t Tabs<int32_t>(int32_t in) { return abs(in); }
template <> __forceinline__ __device__ cuDoubleComplex Tabs<cuDoubleComplex>(cuDoubleComplex in) { return make_cuDoubleComplex(fabs(in.x), fabs(in.y)); }
template <> __forceinline__ __device__ cuFloatComplex Tabs<cuFloatComplex>(cuFloatComplex in) { return make_cuFloatComplex(fabsf(in.x), fabsf(in.y)); }

//------------------------------
// x+y
//------------------------------
template <typename T> __forceinline__ __device__ T Tadd(T a, T b) { return a + b; }
template <> __forceinline__ __device__ cuDoubleComplex Tadd<cuDoubleComplex>(cuDoubleComplex a, cuDoubleComplex b) { return cuCadd(a, b); }
template <> __forceinline__ __device__ cuFloatComplex Tadd<cuFloatComplex>(cuFloatComplex a, cuFloatComplex b) { return cuCaddf(a, b); }

//------------------------------
// x-y
//------------------------------
template <typename T> __forceinline__ __device__ T Tsub(T a, T b) { return a - b; }
template <> __forceinline__ __device__ cuDoubleComplex Tsub<cuDoubleComplex>(cuDoubleComplex a, cuDoubleComplex b) { return cuCsub(a, b); }
template <> __forceinline__ __device__ cuFloatComplex Tsub<cuFloatComplex>(cuFloatComplex a, cuFloatComplex b) { return cuCsubf(a, b); }

//------------------------------
// x*y
//------------------------------
template <typename T1, typename T2 = T1> __forceinline__ __device__ T2 Tmul(T1 a, T2 b) { return a * b; }
template <> __forceinline__ __device__ cuDoubleComplex Tmul<cuDoubleComplex, cuDoubleComplex>(cuDoubleComplex a, cuDoubleComplex b) { return cuCmul(a, b); }
template <> __forceinline__ __device__ cuFloatComplex Tmul<cuFloatComplex, cuFloatComplex>(cuFloatComplex a, cuFloatComplex b) { return cuCmulf(a, b); }
template <> __forceinline__ __device__ cuDoubleComplex Tmul<double, cuDoubleComplex>(double a, cuDoubleComplex b) { return make_cuDoubleComplex(a * b.x, a * b.y); }

//------------------------------
// -x
//------------------------------
template <typename T> __forceinline__ __device__ T Tneg(T a) { return -a; }
template <> __forceinline__ __device__ cuDoubleComplex Tneg<cuDoubleComplex>(cuDoubleComplex a) { return make_cuDoubleComplex(-a.x, -a.y); }
template <> __forceinline__ __device__ cuFloatComplex Tneg<cuFloatComplex>(cuFloatComplex a) { return make_cuFloatComplex(-a.x, -a.y); }

//------------------------------
// x^2 + y
//------------------------------
template <typename T> __forceinline__ __device__ underlying_t<T> Tsqr_add_ru(T in1, underlying_t<T> in2);
template <> __forceinline__ __device__ double Tsqr_add_ru<double>(double in1, double in2) { return __fma_ru(in1, in1, in2); }
template <> __forceinline__ __device__ float Tsqr_add_ru<float>(float in1, float in2) { return __fmaf_ru(in1, in1, in2); }
template <> __forceinline__ __device__ double Tsqr_add_ru<cuDoubleComplex>(cuDoubleComplex in1, double in2) { return __fma_ru(in1.y, in1.y, __fma_ru(in1.x, in1.x, in2)); }
template <> __forceinline__ __device__ float Tsqr_add_ru<cuFloatComplex>(cuFloatComplex in1, float in2) { return __fmaf_ru(in1.y, in1.y, __fmaf_ru(in1.x, in1.x, in2)); }

//------------------------------
// x+y in round-up mode
//------------------------------
template <typename T> __forceinline__ __device__ T __Tadd_ru(T in1, T in2);
template <> __forceinline__ __device__ double __Tadd_ru<double>(double in1, double in2) { return __dadd_ru(in1, in2); }
template <> __forceinline__ __device__ float __Tadd_ru<float>(float in1, float in2) { return __fadd_ru(in1, in2); }

//------------------------------
// a*x + b*y
//------------------------------
template <typename T> __forceinline__ __device__ T Taxpby_scal(T a, T x, T b, T y);
template <> __forceinline__ __device__ double Taxpby_scal<double>(double a, double x, double b, double y) { return fma(b, y, a * x); }
template <> __forceinline__ __device__ float Taxpby_scal<float>(float a, float x, float b, float y) { return fmaf(b, y, a * x); }
template <> __device__ __forceinline__ cuDoubleComplex Taxpby_scal<cuDoubleComplex>(cuDoubleComplex a, cuDoubleComplex x, cuDoubleComplex b, cuDoubleComplex y) {
    double2 out;
    out.x = fma(-b.y, y.y, fma(b.x, y.x, fma(-a.y, x.y, a.x * x.x))); // a.x*x.x - a.y*x.y + b.x*y.x - b.y*y.y
    out.y = fma(b.y, y.x, fma(b.x, y.y, fma(a.y, x.x, a.x * x.y)));   // a.x*x.y + a.y*x.x + b.x*y.y + b.y*y.x
    return out;
}
template <> __device__ __forceinline__ cuFloatComplex Taxpby_scal<cuFloatComplex>(cuFloatComplex a, cuFloatComplex x, cuFloatComplex b, cuFloatComplex y) {
    float2 out;
    out.x = fmaf(-b.y, y.y, fmaf(b.x, y.x, fmaf(-a.y, x.y, a.x * x.x))); // a.x*x.x - a.y*x.y + b.x*y.x - b.y*y.y
    out.y = fmaf(b.y, y.x, fmaf(b.x, y.y, fmaf(a.y, x.x, a.x * x.y)));   // a.x*x.y + a.y*x.x + b.x*y.y + b.y*y.x
    return out;
}

//------------------------------
// x*2^y
//------------------------------
template <typename T> __forceinline__ __device__ T Tscalbn(T in, int sft);
template <> __forceinline__ __device__ double Tscalbn<double>(double in, int sft) { return scalbn(in, sft); }
template <> __forceinline__ __device__ float Tscalbn<float>(float in, int sft) { return scalbnf(in, sft); }
template <> __forceinline__ __device__ cuDoubleComplex Tscalbn<cuDoubleComplex>(cuDoubleComplex in, int sft) { return make_cuDoubleComplex(scalbn(in.x, sft), scalbn(in.y, sft)); }
template <> __forceinline__ __device__ cuFloatComplex Tscalbn<cuFloatComplex>(cuFloatComplex in, int sft) { return make_cuFloatComplex(scalbnf(in.x, sft), scalbnf(in.y, sft)); }

//------------------------------
// rint(x) = round(x)
//------------------------------
template <typename T> __forceinline__ __device__ T Trint(T in);
template <> __forceinline__ __device__ double Trint<double>(double in) { return rint(in); }
template <> __forceinline__ __device__ cuDoubleComplex Trint<cuDoubleComplex>(cuDoubleComplex in) { return make_cuDoubleComplex(rint(in.x), rint(in.y)); }

//------------------------------
// ilogb(x) = floor(log2(x))
//------------------------------
template <typename T> __forceinline__ __device__ int Tilogb(T in);
template <> __forceinline__ __device__ int Tilogb<double>(double in) { return (in == 0.0) ? 0 : ilogb(in); }
template <> __forceinline__ __device__ int Tilogb<float>(float in) { return (in == 0.0F) ? 0 : ilogbf(in); }
//------------------------------
// max(x,y)
//------------------------------
template <typename T> __forceinline__ __device__ underlying_t<T> Tmax(T in1, underlying_t<T> in2);
template <> __forceinline__ __device__ double Tmax<double>(double in1, double in2) { return max(in1, in2); }
template <> __forceinline__ __device__ float Tmax<float>(float in1, float in2) { return max(in1, in2); }
template <> __forceinline__ __device__ int32_t Tmax<int32_t>(int32_t in1, int32_t in2) { return max(in1, in2); }
template <> __forceinline__ __device__ double Tmax<cuDoubleComplex>(cuDoubleComplex in1, double in2) { return max(max(in1.x, in1.y), in2); }
template <> __forceinline__ __device__ float Tmax<cuFloatComplex>(cuFloatComplex in1, float in2) { return max(max(in1.x, in1.y), in2); }

//------------------------------
// Cast Tin to Tout
//------------------------------
template <typename Tin, typename Tout> __forceinline__ __device__ Tout Tcast(Tin in);
template <> __forceinline__ __device__ double Tcast<double, double>(double in) { return in; }
template <> __forceinline__ __device__ float Tcast<double, float>(double in) { return __double2float_rn(in); }
template <> __forceinline__ __device__ cuDoubleComplex Tcast<cuDoubleComplex, cuDoubleComplex>(cuDoubleComplex in) { return in; }
template <> __forceinline__ __device__ cuFloatComplex Tcast<cuDoubleComplex, cuFloatComplex>(cuDoubleComplex in) { return make_cuFloatComplex(__double2float_rn(in.x), __double2float_rn(in.y)); }
template <> __forceinline__ __device__ cuDoubleComplex Tcast<uchar2, cuDoubleComplex>(uchar2 in) { return make_cuDoubleComplex(static_cast<double>(in.x), static_cast<double>(in.y)); }
template <> __forceinline__ __device__ double Tcast<uint8_t, double>(uint8_t in) { return static_cast<double>(in); }

//------------------------------
// static_cast (fp -> int)
//------------------------------
template <typename T> __forceinline__ __device__ int_t<T> __fp2int_rz(T in);
template <> __forceinline__ __device__ int_t<double> __fp2int_rz<double>(double in) { return __double2ll_rz(in); }
template <> __forceinline__ __device__ int_t<float> __fp2int_rz<float>(float in) { return __float2int_rz(in); }

//------------------------------
// reinterpret_cast (fp -> int)
//------------------------------
template <typename T> __forceinline__ __device__ int_t<T> __fp_as_int(T in);
template <> __forceinline__ __device__ int_t<double> __fp_as_int<double>(double in) { return __double_as_longlong(in); }
template <> __forceinline__ __device__ int_t<float> __fp_as_int<float>(float in) { return __float_as_int(in); }

//------------------------------
// reinterpret_cast (fp <- int)
//------------------------------
template <typename T> __forceinline__ __device__ fp_t<T> __int_as_fp(int_t<T> in);
template <> __forceinline__ __device__ fp_t<double> __int_as_fp<double>(int_t<double> in) { return __longlong_as_double(in); }
template <> __forceinline__ __device__ fp_t<float> __int_as_fp<float>(int_t<float> in) { return __int_as_float(in); }

//------------------------------
// Extract sign part
//------------------------------
template <typename T> __forceinline__ __device__ int_t<T> extract_sign(int_t<T> in);
template <> __forceinline__ __device__ int_t<double> extract_sign<double>(int_t<double> in) { return in & 0x8000000000000000LL; }
template <> __forceinline__ __device__ int_t<float> extract_sign<float>(int_t<float> in) { return in & 0x80000000; }

//------------------------------
// Extract exponent part
//------------------------------
template <typename T> __forceinline__ __device__ int_t<T> extract_exp(int_t<T> in);
template <> __forceinline__ __device__ int_t<double> extract_exp<double>(int_t<double> in) { return (in >> 52) & 0x7FF; }
template <> __forceinline__ __device__ int_t<float> extract_exp<float>(int_t<float> in) { return (in >> 23) & 0xFF; }

//------------------------------
// Extract significand part
//------------------------------
template <typename T> __forceinline__ __device__ int_t<T> extract_significand(int_t<T> in);
template <> __forceinline__ __device__ int_t<double> extract_significand<double>(int_t<double> in) { return in & 0x000FFFFFFFFFFFFFULL; }
template <> __forceinline__ __device__ int_t<float> extract_significand<float>(int_t<float> in) { return in & 0x007FFFFF; }

//------------------------------
// Count the number of consecutive high-order zero bits
//------------------------------
template <typename T> __forceinline__ __device__ int32_t countlz(int_t<T> in);
template <> __forceinline__ __device__ int32_t countlz<double>(int_t<double> in) { return __clzll(in); }
template <> __forceinline__ __device__ int32_t countlz<float>(int_t<float> in) { return __clz(in); }

//------------------------------
// Retrieves the complex conjugate of a complex number.
//------------------------------
template <typename T, bool CONJ> __forceinline__ __device__ T conj(T in) { return in; };
template <> __forceinline__ __device__ cuDoubleComplex conj<cuDoubleComplex, true>(cuDoubleComplex in) { return make_cuDoubleComplex(in.x, -in.y); };
template <> __forceinline__ __device__ cuFloatComplex conj<cuFloatComplex, true>(cuFloatComplex in) { return make_cuFloatComplex(in.x, -in.y); };

//------------------------------
// Warp reduction (max)
//------------------------------
template <typename T, int width = 32> __forceinline__ __device__ T inner_warp_max(T amax) {
#pragma unroll
    for (unsigned i = width / 2; i > 0; i >>= 1) {
        amax = max(amax, __shfl_down_sync(0xFFFFFFFFu, amax, i, width));
    }
    return amax;
}

//------------------------------
// Warp reduction (sum in round-up mode)
//------------------------------
template <typename T, int width = 32> __forceinline__ __device__ T inner_warp_sum(T sum) {
#pragma unroll
    for (unsigned i = width / 2; i > 0; i >>= 1) {
        sum = __Tadd_ru<T>(sum, __shfl_down_sync(0xFFFFFFFFu, sum, i, width));
    }
    return sum;
}
template <> __forceinline__ __device__ double inner_warp_sum<double, 32>(double sum) {
    sum = __dadd_ru(sum, __shfl_down_sync(0xFFFFFFFFu, sum, 16u, 32));
    sum = __dadd_ru(sum, __shfl_down_sync(0xFFFFFFFFu, sum, 8u, 32));
    sum = __dadd_ru(sum, __shfl_down_sync(0xFFFFFFFFu, sum, 4u, 32));
    sum = __dadd_ru(sum, __shfl_down_sync(0xFFFFFFFFu, sum, 2u, 32));
    sum = __dadd_ru(sum, __shfl_down_sync(0xFFFFFFFFu, sum, 1u, 32));
    return sum;
}
template <> __forceinline__ __device__ float inner_warp_sum<float, 32>(float sum) {
    sum = __fadd_ru(sum, __shfl_down_sync(0xFFFFFFFFu, sum, 16u, 32));
    sum = __fadd_ru(sum, __shfl_down_sync(0xFFFFFFFFu, sum, 8u, 32));
    sum = __fadd_ru(sum, __shfl_down_sync(0xFFFFFFFFu, sum, 4u, 32));
    sum = __fadd_ru(sum, __shfl_down_sync(0xFFFFFFFFu, sum, 2u, 32));
    sum = __fadd_ru(sum, __shfl_down_sync(0xFFFFFFFFu, sum, 1u, 32));
    return sum;
}
