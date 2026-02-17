#pragma once

//==========
// Helpers for FP8
//==========
struct fp8x3_e4m3 {
    __nv_fp8_e4m3 x, y, z;
};
struct fp8x2_e4m3 {
    __nv_fp8_e4m3 x, y;
};
static __device__ __forceinline__ __nv_fp8x4_e4m3 concat(__nv_fp8_e4m3 a0, __nv_fp8_e4m3 a1, __nv_fp8_e4m3 a2, __nv_fp8_e4m3 a3) {
    uchar4 b;
    b.x = *reinterpret_cast<uint8_t *>(&a0.__x);
    b.y = *reinterpret_cast<uint8_t *>(&a1.__x);
    b.z = *reinterpret_cast<uint8_t *>(&a2.__x);
    b.w = *reinterpret_cast<uint8_t *>(&a3.__x);
    __nv_fp8x4_e4m3 v;
    v.__x = *reinterpret_cast<__nv_fp8x4_storage_t *>(&b);
    return v;
}
static __device__ __forceinline__ char4 concat(int8_t a0, int8_t a1, int8_t a2, int8_t a3) {
    return char4(a0, a1, a2, a3);
}

//==========
// Backend traits (type mapping)
//==========
template <gemmul8::Backend> struct backend_traits;
template <> struct backend_traits<gemmul8::Backend::INT8> {
    using low   = int8_t;
    using lowx2 = char2;
    using lowx4 = char4;

    using mid   = int8_t;
    using midx2 = char2;
    using midx4 = char4;
    using midx8 = short4; // char2*4

    using hi   = int32_t;
    using hix4 = int4;
};
template <> struct backend_traits<gemmul8::Backend::FP8> {
    using low   = __nv_fp8_e4m3;
    using lowx2 = fp8x2_e4m3;
    using lowx4 = __nv_fp8x4_e4m3;

    using mid   = int16_t;
    using midx2 = short2;
    using midx4 = short4;
    using midx8 = int4; // short2*4

    using hi   = float;
    using hix4 = float4;
};
template <gemmul8::Backend b> using low_t   = typename backend_traits<b>::low;
template <gemmul8::Backend b> using lowx2_t = typename backend_traits<b>::lowx2;
template <gemmul8::Backend b> using lowx4_t = typename backend_traits<b>::lowx4;
template <gemmul8::Backend b> using mid_t   = typename backend_traits<b>::mid;
template <gemmul8::Backend b> using midx2_t = typename backend_traits<b>::midx2;
template <gemmul8::Backend b> using midx4_t = typename backend_traits<b>::midx4;
template <gemmul8::Backend b> using midx8_t = typename backend_traits<b>::midx8;
template <gemmul8::Backend b> using hi_t    = typename backend_traits<b>::hi;
template <gemmul8::Backend b> using hix4_t  = typename backend_traits<b>::hix4;

//==========
// types of scaled values
//==========
template <typename T> struct fp_mant_exp;
template <> struct fp_mant_exp<float> {
    int32_t exp;
    int32_t mant;
};
template <> struct fp_mant_exp<double> {
    int32_t exp;
    int64_t mant;
};
template <typename T> struct fp_mant_exp2 {
    fp_mant_exp<T> x; // real
    fp_mant_exp<T> y; // imag
};

struct __align__(16) int64x2_t {
    int64_t x;
    int64_t y;
};

//==========
// Check if type is complex
//==========
template <typename T> inline constexpr bool isComplex             = false;
template <> inline constexpr bool isComplex<cuDoubleComplex>      = true;
template <> inline constexpr bool isComplex<cuFloatComplex>       = true;
template <> inline constexpr bool isComplex<int2>                 = true;
template <> inline constexpr bool isComplex<int64x2_t>            = true;
template <> inline constexpr bool isComplex<fp_mant_exp2<float>>  = true;
template <> inline constexpr bool isComplex<fp_mant_exp2<double>> = true;

//==========
// same-size maps
//==========
template <typename T> using int_t  = std::conditional_t<(sizeof(T) == 8), int64_t, int32_t>;
template <typename T> using uint_t = std::conditional_t<(sizeof(T) == 8), uint64_t, uint32_t>;
template <typename T> using fp_t   = std::conditional_t<(sizeof(T) == 8), double, float>;

//==========
// Map type to underlying scalar type
//==========
template <typename T> struct underlying_type {
    using type = T;
};
template <> struct underlying_type<int2> {
    using type = int32_t;
};
template <> struct underlying_type<int64x2_t> {
    using type = int64_t;
};
template <> struct underlying_type<cuFloatComplex> {
    using type = float;
};
template <> struct underlying_type<cuDoubleComplex> {
    using type = double;
};
template <typename T> using underlying_t = typename underlying_type<T>::type;

//==========
// Output type of upperBound_lo<backend, T>
//==========
template <gemmul8::Backend b, bool complex> struct upperBound_impl;
template <> struct upperBound_impl<gemmul8::Backend::INT8, false> {
    using type = int8_t;
};
template <> struct upperBound_impl<gemmul8::Backend::INT8, true> {
    using type = char2;
};
template <> struct upperBound_impl<gemmul8::Backend::FP8, false> {
    using type = __nv_fp8_e4m3;
};
template <> struct upperBound_impl<gemmul8::Backend::FP8, true> {
    using type = fp8x2_e4m3;
};
template <gemmul8::Backend b, typename T> using upperBound_t = typename upperBound_impl<b, isComplex<T>>::type;

//==========
// Number of bits for extraction
//==========
template <gemmul8::Backend b> inline constexpr int maxUFP = (b == gemmul8::Backend::FP8) ? 6 : 5;

//==========
// Floating-point traits
//==========
template <typename T> struct fp;
template <> struct fp<double> {
    static constexpr int32_t bias = 1023;
    static constexpr int32_t prec = 52;
    static constexpr int32_t bits = 64;
};
template <> struct fp<float> {
    static constexpr int32_t bias = 127;
    static constexpr int32_t prec = 23;
    static constexpr int32_t bits = 32;
};

//==========
// Type-specific constants
//==========
template <typename T> struct Tconst {
    __device__ __host__ __forceinline__ static constexpr T zero() { return static_cast<T>(0); }
    __device__ __host__ __forceinline__ static constexpr T one() { return static_cast<T>(1.0); }
    __device__ __host__ __forceinline__ static constexpr T mone() { return static_cast<T>(-1.0); }
};

template <> struct Tconst<double2> {
    __device__ __host__ __forceinline__ static constexpr double2 zero() { return {0.0, 0.0}; }
    __device__ __host__ __forceinline__ static constexpr double2 one() { return {1.0, 0.0}; }
    __device__ __host__ __forceinline__ static constexpr double2 mone() { return {-1.0, 0.0}; }
};

template <> struct Tconst<float2> {
    __device__ __host__ __forceinline__ static constexpr float2 zero() { return {0.0f, 0.0f}; }
    __device__ __host__ __forceinline__ static constexpr float2 one() { return {1.0f, 0.0f}; }
    __device__ __host__ __forceinline__ static constexpr float2 mone() { return {-1.0f, 0.0f}; }
};

template <> struct Tconst<char2> {
    __device__ __host__ __forceinline__ static constexpr char2 zero() {
        return {static_cast<int8_t>(0), static_cast<int8_t>(0)};
    }
};
