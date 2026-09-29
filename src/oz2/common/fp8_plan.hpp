#pragma once
#include <cstdint>

namespace gemmul8::common::fp8_plan {

struct info {
    int32_t p;
    unsigned products;
    int32_t w;
    int32_t q[3];
    int32_t coefficient[3];
    unsigned max_abs;
    int32_t root_minus_one;
};

constexpr unsigned modulus_index(unsigned p) { return 20U + p; }

constexpr info get(int32_t p) {
    switch (p) {
    case 59163:
        return {
            59163,
            3U,
            0,
            {   41,     39,    37},
            {-7215, -15170, 22386},
            24U,
            0
        };
    case 2303:
        return {
            2303,
            2U,
            0,
            { 49, 47, 1},
            {-47, 49, 0},
            32U,
            0
        };
    case 1376:
        return {
            1376,
            2U,
            0,
            {  43,  32, 1},
            {-128, 129, 0},
            26U,
            0
        };
    case 899:
        return {
            899,
            2U,
            0,
            { 31, 29, 1},
            {-29, 62, 0},
            15U,
            0
        };
    case 575:
        return {
            575,
            2U,
            0,
            { 25, 23, 1},
            {-46, 25, 0},
            12U,
            0
        };
    case 1283:
        return {
            1283,
            3U,
            201,
            {1283, 1283, 1283},
            {  64,  -34,  -13},
            32U,
            0
        };
    case 1279:
        return {
            1279,
            3U,
            599,
            {1279, 1279, 1279},
            {  32,   17,   15},
            32U,
            0
        };
    case 1249:
        return {
            1249,
            3U,
            293,
            {1249, 1249, 1249},
            {  17,   15,  -47},
            30U,
            585
        };
    case 323:
        return {
            323,
            2U,
            0,
            {19, 17, 1},
            {17, 19, 0},
            9U,
            0
        };
    case 1223:
        return {
            1223,
            3U,
            286,
            {1223, 1223, 1223},
            {  47,   11,  -30},
            32U,
            0
        };
    case 1201:
        return {
            1201,
            3U,
            128,
            {1201, 1201, 1201},
            {  47,  -11,  -19},
            32U,
            49
        };
    case 1193:
        return {
            1193,
            3U,
            186,
            {1193, 1193, 1193},
            {  45,  -19,  -13},
            30U,
            186
        };
    case 1181:
        return {
            1181,
            3U,
            185,
            {1181, 1181, 1181},
            {  51,   13,   45},
            30U,
            243
        };
    case 1177:
        return {
            1177,
            3U,
            313,
            {1177, 1177, 1177},
            { -15,  -13,   49},
            32U,
            0
        };
    case 1171:
        return {
            1171,
            3U,
            137,
            {1171, 1171, 1171},
            { -60,   23,   26},
            32U,
            0
        };
    case 1163:
        return {
            1163,
            3U,
            104,
            {1163, 1163, 1163},
            {  34,  -47,   56},
            30U,
            0
        };
    case 1153:
        return {
            1153,
            3U,
            140,
            {1153, 1153, 1153},
            {  25,  -41,    8},
            32U,
            140
        };
    case 1151:
        return {
            1151,
            3U,
            205,
            {1151, 1151, 1151},
            {  28,   15,   17},
            28U,
            0
        };
    case 1129:
        return {
            1129,
            3U,
            176,
            {1129, 1129, 1129},
            { -45,   17,   13},
            28U,
            168
        };
    case 1123:
        return {
            1123,
            3U,
            55,
            {1123, 1123, 1123},
            {  21,  -32,   41},
            30U,
            0
        };
    case 65231:
        return {
            65231,
            3U,
            0,
            {   43,   41,    37},
            {27306, 7955, 29971},
            26U,
            0
        };
    case 992:
        return {
            992,
            2U,
            0,
            { 32, 31, 1},
            {-31, 32, 0},
            16U,
            0
        };
    case 783:
        return {
            783,
            2U,
            0,
            { 29,  27, 1},
            {-54, -29, 0},
            14U,
            0
        };
    case 4199:
        return {
            4199,
            3U,
            0,
            {  19,  17,   13},
            {1768, 494, 1938},
            9U,
            0
        };
    case 11:
        return {
            11,
            1U,
            0,
            {11, 1, 1},
            { 1, 0, 0},
            5U,
            0
        };
    case 18241:
        return {
            18241,
            3U,
            0,
            {   37,    29,    17},
            {-1479, -8177, -8584},
            20U,
            191
        };
    case 1025:
        return {
            1025,
            2U,
            0,
            { 41,  25, 1},
            {-25, -41, 0},
            24U,
            32
        };
    case 1313:
        return {
            1313,
            3U,
            349,
            {1313, 1313, 1313},
            {  49,  -32,   15},
            32U,
            515
        };
    case 1517:
        return {
            1517,
            2U,
            0,
            { 41, 37, 1},
            {-37, 41, 0},
            24U,
            401
        };
    case 725:
        return {
            725,
            2U,
            0,
            {25,  29, 1},
            {29, -25, 0},
            14U,
            157
        };
    case 1117:
        return {
            1117,
            3U,
            472,
            {1117, 1117, 1117},
            { -26,  -15,   64},
            28U,
            214
        };
    case 1109:
        return {
            1109,
            3U,
            406,
            {1109, 1109, 1109},
            { -30,  -19,  -11},
            32U,
            354
        };
    case 1097:
        return {
            1097,
            3U,
            189,
            {1097, 1097, 1097},
            {  35,  -33,   29},
            32U,
            341
        };
    case 1093:
        return {
            1093,
            3U,
            282,
            {1093, 1093, 1093},
            { -31,   -2,   35},
            32U,
            530
        };
    case 1069:
        return {
            1069,
            3U,
            87,
            {1069, 1069, 1069},
            { -12,  -25,   37},
            28U,
            249
        };
    case 6409:
        return {
            6409,
            3U,
            0,
            {   29,   17,   13},
            {-1768, 2262, -493},
            14U,
            684
        };
    case 25:
        return {
            25,
            1U,
            0,
            {25, 1, 1},
            { 1, 0, 0},
            12U,
            7
        };
    case 1061:
        return {
            1061,
            3U,
            50,
            {1061, 1061, 1061},
            { -64,   17,   22},
            26U,
            103
        };
    case 1049:
        return {
            1049,
            3U,
            142,
            {1049, 1049, 1049},
            {  22,   23,   82},
            28U,
            426
        };
    case 1033:
        return {
            1033,
            3U,
            178,
            {1033, 1033, 1033},
            {  29,    3,  -35},
            30U,
            355
        };
    case 1021:
        return {
            1021,
            3U,
            50,
            {1021, 1021, 1021},
            {  41,   -8,  -21},
            28U,
            374
        };
    case 1013:
        return {
            1013,
            3U,
            71,
            {1013, 1013, 1013},
            {  29,  -33,   14},
            28U,
            45
        };
    case 1009:
        return {
            1009,
            3U,
            82,
            {1009, 1009, 1009},
            {  37,   -7,  -25},
            28U,
            469
        };
    default: return {};
    }
}

template <int32_t P> inline constexpr info scheme = get(P);

__host__ __device__ constexpr int32_t crt_unscale(int32_t p, unsigned j) {
    switch (p) {
    case 2303: return j == 0 ? 5 : 27;
    case 899: return j == 0 ? 4 : 15;
    case 575: return j == 0 ? 13 : 14;
    case 323: return j == 0 ? 16 : 3;
    case 783: return j == 0 ? 15 : 11;
    case 1025: return j == 0 ? 31 : 17;
    case 1517: return j == 0 ? 21 : 19;
    case 725: return j == 0 ? 13 : 15;
    default: return 1;
    }
}

__host__ __device__ constexpr bool crt_scaled(int32_t p) {
    return crt_unscale(p, 0) != 1 || crt_unscale(p, 1) != 1;
}

inline constexpr int32_t real_prefix[]    = {59163, 2303, 1376, 899, 575, 1283, 1279, 1249, 323, 1223, 1201, 1193, 1181, 1177, 1171, 1163, 1153, 1151, 1129, 1123};
inline constexpr int32_t real_20[]        = {65231, 2303, 992, 783, 575, 4199, 1283, 1279, 1249, 1223, 1201, 1193, 1181, 1171, 1163, 1153, 1151, 1129, 1123, 11};
inline constexpr int32_t complex_small[]  = {18241, 1025, 1313, 1249, 1201, 1193};
inline constexpr int32_t complex_medium[] = {1517, 725, 1313, 1249, 1201, 1193, 1181, 1153, 1129, 1117, 1109, 1097, 1093, 1069};
inline constexpr int32_t complex_large[]  = {6409, 1517, 25, 1249, 1201, 1193, 1181, 1153, 1129, 1117, 1109, 1097, 1093, 1069, 1061, 1049, 1033, 1021, 1013, 1009};

constexpr int32_t modulus(unsigned n, unsigned i, bool complex) {
    if (!complex) return n == 20U ? real_20[i] : real_prefix[i];
    return n <= 6U    ? complex_small[i]
           : n <= 14U ? complex_medium[i]
                      : complex_large[i];
}

constexpr unsigned index(unsigned n, unsigned i, bool complex) {
    return modulus_index(unsigned(modulus(n, i, complex)));
}

constexpr unsigned planes(unsigned n, unsigned first, unsigned count, bool complex) {
    unsigned sum = 0;
    for (unsigned i = first; i < first + count; ++i) sum += get(modulus(n, i, complex)).products;
    return sum;
}

template <unsigned N, bool Complex>
inline constexpr unsigned plane_count = planes(N, 0U, N, Complex);

constexpr unsigned k_block(int32_t p, bool first, unsigned alignment = 256U) {
    const auto s          = get(p);
    const unsigned budget = (1U << 24) - (first ? 0U : unsigned(p / 2));
    return (budget / (s.max_abs * s.max_abs) / alignment) * alignment;
}

} // namespace gemmul8::common::fp8_plan

#define GEMMUL8_FP8_FOR_EACH_MODULUS(M) \
    M(59163)                            \
    M(2303)                             \
    M(1376)                             \
    M(899)                              \
    M(575)                              \
    M(1283)                             \
    M(1279)                             \
    M(1249)                             \
    M(323)                              \
    M(1223)                             \
    M(1201)                             \
    M(1193)                             \
    M(1181)                             \
    M(1177)                             \
    M(1171)                             \
    M(1163)                             \
    M(1153)                             \
    M(1151)                             \
    M(1129)                             \
    M(1123)                             \
    M(65231)                            \
    M(992)                              \
    M(783)                              \
    M(4199)                             \
    M(11)                               \
    M(18241)                            \
    M(1025)                             \
    M(1313)                             \
    M(1517)                             \
    M(725)                              \
    M(1117)                             \
    M(1109)                             \
    M(1097)                             \
    M(1093)                             \
    M(1069)                             \
    M(6409)                             \
    M(25)                               \
    M(1061)                             \
    M(1049)                             \
    M(1033)                             \
    M(1021)                             \
    M(1013)                             \
    M(1009)

#define GEMMUL8_FP8_FOR_EACH_COMPLEX_MODULUS(M) \
    M(1249)                                     \
    M(1201)                                     \
    M(1193)                                     \
    M(1181)                                     \
    M(1153)                                     \
    M(1129)                                     \
    M(18241)                                    \
    M(1025)                                     \
    M(1313)                                     \
    M(1517)                                     \
    M(725)                                      \
    M(1117)                                     \
    M(1109)                                     \
    M(1097)                                     \
    M(1093)                                     \
    M(1069)                                     \
    M(6409)                                     \
    M(25)                                       \
    M(1061)                                     \
    M(1049)                                     \
    M(1033)                                     \
    M(1021)                                     \
    M(1013)                                     \
    M(1009)
