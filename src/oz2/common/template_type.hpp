#pragma once
#include "include.hpp"
#include "threshold.hpp"

namespace gemmul8::common {

//==========
// matrix structure
//==========
enum class MatStruct {
    Full,
    Triangular,
    Symmetric,
    Hermitian
};

template <MatStruct S> inline constexpr bool is_full       = (S == MatStruct::Full);
template <MatStruct S> inline constexpr bool is_triangular = (S == MatStruct::Triangular);
template <MatStruct S> inline constexpr bool is_symmetric  = (S == MatStruct::Symmetric);
template <MatStruct S> inline constexpr bool is_hermitian  = (S == MatStruct::Hermitian);
template <MatStruct S> inline constexpr bool needs_uplo    = S == MatStruct::Triangular ||
                                                             S == MatStruct::Symmetric ||
                                                             S == MatStruct::Hermitian;

enum class MatMulKind : unsigned char {
    Gemm,
    ATxB,
    ATxA,
    AHxA,
    TrmmLeft,
    TrmmRight,
    Trtrmm
};

//==========
// workspace pointer set
//==========
template <typename T, bool COMPLEX> struct matptr_type;
template <typename T> struct matptr_type<T, true> {
    T *ptr0, *ptr1, *ptr2;
    __forceinline__ void shift(const size_t offset) {
        ptr0 += offset;
        ptr1 += offset;
        if (ptr2 != nullptr) ptr2 += offset;
    }
};
template <typename T> struct matptr_type<T, false> {
    T *ptr0;
    __forceinline__ void shift(const size_t offset) {
        ptr0 += offset;
    }
};
template <typename T, bool COMPLEX> using matptr_t = matptr_type<T, COMPLEX>;

template <typename T, bool COMPLEX, unsigned PARTS = COMPLEX ? 3U : 1U>
__forceinline__ matptr_t<T, COMPLEX> make_matptr(T *const base, const size_t offset) {
    if constexpr (COMPLEX) {
        static_assert(PARTS == 2U || PARTS == 3U);
        return matptr_t<T, true>{base, base + offset, PARTS == 3U ? base + offset * 2 : nullptr};
    } else {
        return matptr_t<T, false>{base};
    }
}

//==========
// Helpers for FP8
//==========
struct fp8x2_e4m3 {
    __nv_fp8_e4m3 x, y;
};

struct fp8x3_e4m3 {
    __nv_fp8_e4m3 x, y, z;
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
    return char4{a0, a1, a2, a3};
}

//==========
// Backend traits (type mapping)
//==========
template <Backend> struct backend_traits;
template <> struct backend_traits<Backend::INT8> {
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
template <> struct backend_traits<Backend::FP8> {
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
template <Backend b> using low_t           = typename backend_traits<b>::low;
template <Backend b> using lowx2_t         = typename backend_traits<b>::lowx2;
template <Backend b> using lowx4_t         = typename backend_traits<b>::lowx4;
template <Backend b> using mid_t_real      = typename backend_traits<b>::mid;
template <Backend b> using mid_t_complex   = typename backend_traits<b>::midx2;
template <Backend b> using midx4_t_real    = typename backend_traits<b>::midx4;
template <Backend b> using midx4_t_complex = typename backend_traits<b>::midx8;
template <Backend b> using hi_t            = typename backend_traits<b>::hi;
template <Backend b> using hix4_t          = typename backend_traits<b>::hix4;

template <Backend b, bool COMPLEX = false> using mid_t   = std::conditional_t<COMPLEX, mid_t_complex<b>, mid_t_real<b>>;
template <Backend b, bool COMPLEX = false> using midx4_t = std::conditional_t<COMPLEX, midx4_t_complex<b>, midx4_t_real<b>>;

//==========
// types of scaled values
//==========
struct __align__(16) int64x2_t {
    int64_t x;
    int64_t y;
};

struct __align__(16) double2x2_t {
    double2 hi;
    double2 lo;
};

struct mant_t {
    int32_t hi;
    uint32_t lo;
    __device__ __forceinline__ constexpr mant_t() : hi(0), lo(0u) {}
    __device__ __forceinline__ constexpr mant_t(int64_t a)
        : hi(int32_t(uint64_t(a) >> 32)), lo(uint32_t(uint64_t(a))) {}
};

struct mant2_t {
    mant_t x; // real
    mant_t y; // imag
    __device__ __forceinline__ constexpr mant2_t() : x(mant_t{}), y(mant_t{}) {}
    __device__ __forceinline__ constexpr mant2_t(int64x2_t a) : x(a.x), y(a.y) {}
};

struct exp_t {
    uint32_t val;
    bool is_hi;
    __device__ __forceinline__ constexpr exp_t() : val(0u), is_hi(false) {}
    __device__ __forceinline__ constexpr exp_t(uint32_t val_, bool is_hi_) : val(val_), is_hi(is_hi_) {}
};

struct fp32_mant_exp {
    int32_t mant;
    exp_t exp;
};

struct fp64_mant_exp {
    mant_t mant;
    exp_t exp;
};

struct fp32_mant_exp2 {
    fp32_mant_exp x; // real
    fp32_mant_exp y; // imag
};

struct fp64_mant_exp2 {
    fp64_mant_exp x; // real
    fp64_mant_exp y; // imag
};

//==========
// Check if type is complex
//==========
template <typename T> inline constexpr bool isComplex             = false;
template <> inline constexpr bool isComplex<cuDoubleComplex>      = true;
template <> inline constexpr bool isComplex<cuFloatComplex>       = true;
template <> inline constexpr bool isComplex<int2>                 = true;
template <> inline constexpr bool isComplex<int64x2_t>            = true;
template <> inline constexpr bool isComplex<mant2_t>              = true;
template <> inline constexpr bool isComplex<fp32_mant_exp2>       = true;
template <> inline constexpr bool isComplex<fp64_mant_exp2>       = true;

//==========
// same-size maps
//==========
template <typename T> using int32_type = std::conditional_t<(isComplex<T>), int2, int32_t>;
template <typename T> using int64_type = std::conditional_t<(isComplex<T>), int64x2_t, int64_t>;

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
    __device__ __host__ __forceinline__ static constexpr T zero() { return T(0); }
    __device__ __host__ __forceinline__ static constexpr T one() { return T(1.0); }
    __device__ __host__ __forceinline__ static constexpr T mone() { return T(-1.0); }
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

template <> struct Tconst<char4> {
    __device__ __host__ __forceinline__ static constexpr char4 zero() { return char4{}; }
};

template <> struct Tconst<char2> {
    __device__ __host__ __forceinline__ static constexpr char2 zero() { return char2{}; }
};

template <> struct Tconst<int8_t> {
    __device__ __host__ __forceinline__ static constexpr int8_t zero() { return int8_t{}; }
};

template <> struct Tconst<__nv_fp8x4_e4m3> {
    __device__ __forceinline__ static constexpr __nv_fp8x4_e4m3 zero() { return __nv_fp8x4_e4m3{}; }
};

template <> struct Tconst<__nv_fp8_e4m3> {
    __device__ __forceinline__ static constexpr __nv_fp8_e4m3 zero() { return __nv_fp8_e4m3{}; }
};

//==========
// 4-vector type for aligned vector loads
//==========
template <typename T> struct vec4_type;
template <>
struct vec4_type<float> {
    using type                        = float4;
    static constexpr size_t alignment = 16ULL;
};

template <>
struct vec4_type<double> {
#if defined(__CUDACC_VER_MAJOR__) && (__CUDACC_VER_MAJOR__ >= 13) && !defined(__HIPCC__)
    using type                        = double4_32a;
    static constexpr size_t alignment = 32ULL;
#else
    using type                        = double4;
    static constexpr size_t alignment = 32ULL;
#endif
};

template <typename T> using vec4_t = typename vec4_type<T>::type;

template <typename T> inline constexpr size_t vec4_alignment_v = vec4_type<T>::alignment;

template <typename T>
__host__ __device__ __forceinline__ bool is_aligned_as(const void *ptr) {
    if constexpr (isComplex<T>) {
        return false;
    } else {
        return (reinterpret_cast<uintptr_t>(ptr) & (vec4_alignment_v<T> - 1ULL)) == 0ULL;
    }
}

} // namespace gemmul8::common
