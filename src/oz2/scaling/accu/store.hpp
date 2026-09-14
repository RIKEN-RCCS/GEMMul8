#pragma once
#include "../general/roundup.hpp"

namespace gemmul8::scaling::accu {

template <Backend b, bool COMPLEX> struct upperBound_impl;
template <> struct upperBound_impl<Backend::INT8, false> { using type = int8_t; };
template <> struct upperBound_impl<Backend::INT8, true> { using type = char2; };
template <> struct upperBound_impl<Backend::FP8, false> { using type = __nv_fp8_e4m3; };
template <> struct upperBound_impl<Backend::FP8, true> { using type = common::fp8x2_e4m3; };
template <typename T, Backend b> using upperBound_t = typename upperBound_impl<b, common::isComplex<T>>::type;

// ceil(in * 2^sft) for in >= 0
template <typename T>
__device__ __forceinline__ int8_t ceil_scalbn_8i(T in, const int32_t sft) {
    return static_cast<int8_t>(common::__fp2int_ru<T>(common::Tscalbn<T>(common::Tabs<T>(in), sft)));
}

// round_up(in) in FP8_E4M3
template <typename T>
__device__ __forceinline__ __nv_fp8_e4m3 fp8_e4m3_ru(T in) {
    __nv_fp8_e4m3 r               = __nv_fp8_e4m3(in);
    const T y                     = T(r);
    const __nv_fp8_storage_t bits = (__nv_fp8_storage_t)(unsigned(r.__x) + unsigned(y < in));
    r.__x                         = bits;
    return r;
}

template <typename T, Backend BACKEND>
__device__ __forceinline__ upperBound_t<T, BACKEND> upperBound_lo(T in, const int32_t sft);

template <> __device__ __forceinline__ int8_t upperBound_lo<float, Backend::INT8>(float in, const int32_t sft) {
    return ceil_scalbn_8i<float>(in, sft);
}

template <> __device__ __forceinline__ int8_t upperBound_lo<double, Backend::INT8>(double in, const int32_t sft) {
    return ceil_scalbn_8i<double>(in, sft);
}

template <> __device__ __forceinline__ char2 upperBound_lo<cuFloatComplex, Backend::INT8>(cuFloatComplex in, const int32_t sft) {
    char2 out;
    out.x = ceil_scalbn_8i<float>(in.x, sft);
    out.y = ceil_scalbn_8i<float>(in.y, sft);
    return out;
}

template <> __device__ __forceinline__ char2 upperBound_lo<cuDoubleComplex, Backend::INT8>(cuDoubleComplex in, const int32_t sft) {
    char2 out;
    out.x = ceil_scalbn_8i<double>(in.x, sft);
    out.y = ceil_scalbn_8i<double>(in.y, sft);
    return out;
}

template <> __device__ __forceinline__ __nv_fp8_e4m3 upperBound_lo<float, Backend::FP8>(float in, const int32_t sft) {
    return fp8_e4m3_ru<float>(scalbnf(fabsf(in), sft));
}

template <> __device__ __forceinline__ __nv_fp8_e4m3 upperBound_lo<double, Backend::FP8>(double in, const int32_t sft) {
    return fp8_e4m3_ru<double>(scalbn(fabs(in), sft));
}

template <> __device__ __forceinline__ common::fp8x2_e4m3 upperBound_lo<cuFloatComplex, Backend::FP8>(cuFloatComplex in, const int32_t sft) {
    common::fp8x2_e4m3 out;
    out.x = fp8_e4m3_ru<float>(scalbnf(fabsf(in.x), sft));
    out.y = fp8_e4m3_ru<float>(scalbnf(fabsf(in.y), sft));
    return out;
}

template <> __device__ __forceinline__ common::fp8x2_e4m3 upperBound_lo<cuDoubleComplex, Backend::FP8>(cuDoubleComplex in, const int32_t sft) {
    common::fp8x2_e4m3 out;
    out.x = fp8_e4m3_ru<double>(scalbn(fabs(in.x), sft));
    out.y = fp8_e4m3_ru<double>(scalbn(fabs(in.y), sft));
    return out;
}

__device__ __forceinline__ int8_t sub_ru_8bit(int8_t a, int8_t b) {
    return int8_t(a - b);
}
__device__ __forceinline__ __nv_fp8_e4m3 sub_ru_8bit(__nv_fp8_e4m3 a, __nv_fp8_e4m3 b) {
    return fp8_e4m3_ru(fabsf(float(a) - float(b)));
}

template <typename T, Backend BACKEND>
__device__ __forceinline__ void extract_store_one_real(
    common::matptr_t<common::low_t<BACKEND>, false> A_lo,
    const size_t idx,
    const T a,
    const int32_t sft //
) {
    const auto v   = upperBound_lo<T, BACKEND>(a, sft);
    A_lo.ptr0[idx] = v;
}

template <typename T, Backend BACKEND>
__device__ __forceinline__ void extract_store_one_complex(
    common::matptr_t<common::low_t<BACKEND>, true> A_lo,
    const size_t idx,
    const T a,
    const int32_t sft //
) {
    const auto v = upperBound_lo<T, BACKEND>(a, sft);

    A_lo.ptr0[idx] = v.x;                   // Re
    A_lo.ptr1[idx] = v.y;                   // Im
    A_lo.ptr2[idx] = sub_ru_8bit(v.x, v.y); // Re - Im
}

template <typename T, Backend BACKEND>
__device__ __forceinline__ void extract_store_one(
    common::matptr_t<common::low_t<BACKEND>, common::isComplex<T>> A_lo,
    const size_t idx,
    const T a,
    const int32_t sft //
) {
    if constexpr (common::isComplex<T>) {
        extract_store_one_complex<T, BACKEND>(A_lo, idx, a, sft);
    } else {
        extract_store_one_real<T, BACKEND>(A_lo, idx, a, sft);
    }
}

template <typename T, Backend BACKEND>
__device__ __forceinline__ void extract_colwise_store4_real(
    common::lowx4_t<BACKEND> *__restrict__ out,
    const unsigned i,
    const T a0, const T a1, const T a2, const T a3,
    const int32_t sft //
) {
    using ValT = upperBound_t<T, BACKEND>;

    const ValT v0 = upperBound_lo<T, BACKEND>(a0, sft);
    const ValT v1 = upperBound_lo<T, BACKEND>(a1, sft);
    const ValT v2 = upperBound_lo<T, BACKEND>(a2, sft);
    const ValT v3 = upperBound_lo<T, BACKEND>(a3, sft);

    out[i] = common::concat(v0, v1, v2, v3);
}

template <typename T, Backend BACKEND>
__device__ __forceinline__ void extract_colwise_store4_complex(
    common::lowx4_t<BACKEND> *__restrict__ out_1,
    common::lowx4_t<BACKEND> *__restrict__ out_2,
    common::lowx4_t<BACKEND> *__restrict__ out_3,
    const unsigned i,
    const T a0, const T a1, const T a2, const T a3,
    const int32_t sft //
) {
    using ValT = upperBound_t<T, BACKEND>;

    const ValT v0 = upperBound_lo<T, BACKEND>(a0, sft);
    const ValT v1 = upperBound_lo<T, BACKEND>(a1, sft);
    const ValT v2 = upperBound_lo<T, BACKEND>(a2, sft);
    const ValT v3 = upperBound_lo<T, BACKEND>(a3, sft);

    out_1[i] = common::concat(v0.x, v1.x, v2.x, v3.x);                          // Re
    out_2[i] = common::concat(v0.y, v1.y, v2.y, v3.y);                          // Im
    out_3[i] = common::concat(sub_ru_8bit(v0.x, v0.y), sub_ru_8bit(v1.x, v1.y), // INT8: Re - Im; FP8: abs(Re - Im), rounded up
                              sub_ru_8bit(v2.x, v2.y), sub_ru_8bit(v3.x, v3.y));
}

} // namespace gemmul8::scaling::accu
