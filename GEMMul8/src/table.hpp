#pragma once
// Constants for num_moduli <= 20

namespace table {

//==========
// moduli
//==========
template <gemmul8::Backend backend, int IDX> inline constexpr int32_t moduli = 0;

// INT8: moduli
template <> inline constexpr int32_t moduli<gemmul8::Backend::INT8, 0>  = 256;
template <> inline constexpr int32_t moduli<gemmul8::Backend::INT8, 1>  = 255;
template <> inline constexpr int32_t moduli<gemmul8::Backend::INT8, 2>  = 253;
template <> inline constexpr int32_t moduli<gemmul8::Backend::INT8, 3>  = 251;
template <> inline constexpr int32_t moduli<gemmul8::Backend::INT8, 4>  = 247;
template <> inline constexpr int32_t moduli<gemmul8::Backend::INT8, 5>  = 241;
template <> inline constexpr int32_t moduli<gemmul8::Backend::INT8, 6>  = 239;
template <> inline constexpr int32_t moduli<gemmul8::Backend::INT8, 7>  = 233;
template <> inline constexpr int32_t moduli<gemmul8::Backend::INT8, 8>  = 229;
template <> inline constexpr int32_t moduli<gemmul8::Backend::INT8, 9>  = 227;
template <> inline constexpr int32_t moduli<gemmul8::Backend::INT8, 10> = 223;
template <> inline constexpr int32_t moduli<gemmul8::Backend::INT8, 11> = 217;
template <> inline constexpr int32_t moduli<gemmul8::Backend::INT8, 12> = 211;
template <> inline constexpr int32_t moduli<gemmul8::Backend::INT8, 13> = 199;
template <> inline constexpr int32_t moduli<gemmul8::Backend::INT8, 14> = 197;
template <> inline constexpr int32_t moduli<gemmul8::Backend::INT8, 15> = 193;
template <> inline constexpr int32_t moduli<gemmul8::Backend::INT8, 16> = 191;
template <> inline constexpr int32_t moduli<gemmul8::Backend::INT8, 17> = 181;
template <> inline constexpr int32_t moduli<gemmul8::Backend::INT8, 18> = 179;
template <> inline constexpr int32_t moduli<gemmul8::Backend::INT8, 19> = 173;

// FP8: moduli
template <> inline constexpr int32_t moduli<gemmul8::Backend::FP8, 0>  = 1089; // NOT Karatsuba
template <> inline constexpr int32_t moduli<gemmul8::Backend::FP8, 1>  = 1024; // NOT Karatsuba
template <> inline constexpr int32_t moduli<gemmul8::Backend::FP8, 2>  = 961;  // NOT Karatsuba
template <> inline constexpr int32_t moduli<gemmul8::Backend::FP8, 3>  = 841;  // NOT Karatsuba
template <> inline constexpr int32_t moduli<gemmul8::Backend::FP8, 4>  = 625;  // NOT Karatsuba
template <> inline constexpr int32_t moduli<gemmul8::Backend::FP8, 5>  = 529;  // NOT Karatsuba
template <> inline constexpr int32_t moduli<gemmul8::Backend::FP8, 6>  = 481;  // Karatsuba
template <> inline constexpr int32_t moduli<gemmul8::Backend::FP8, 7>  = 479;  // Karatsuba
template <> inline constexpr int32_t moduli<gemmul8::Backend::FP8, 8>  = 469;  // Karatsuba
template <> inline constexpr int32_t moduli<gemmul8::Backend::FP8, 9>  = 467;  // Karatsuba
template <> inline constexpr int32_t moduli<gemmul8::Backend::FP8, 10> = 463;  // Karatsuba
template <> inline constexpr int32_t moduli<gemmul8::Backend::FP8, 11> = 461;  // Karatsuba
template <> inline constexpr int32_t moduli<gemmul8::Backend::FP8, 12> = 457;  // Karatsuba
template <> inline constexpr int32_t moduli<gemmul8::Backend::FP8, 13> = 449;  // Karatsuba
template <> inline constexpr int32_t moduli<gemmul8::Backend::FP8, 14> = 443;  // Karatsuba
template <> inline constexpr int32_t moduli<gemmul8::Backend::FP8, 15> = 439;  // Karatsuba
template <> inline constexpr int32_t moduli<gemmul8::Backend::FP8, 16> = 433;  // Karatsuba
template <> inline constexpr int32_t moduli<gemmul8::Backend::FP8, 17> = 431;  // Karatsuba
template <> inline constexpr int32_t moduli<gemmul8::Backend::FP8, 18> = 421;  // Karatsuba
template <> inline constexpr int32_t moduli<gemmul8::Backend::FP8, 19> = 419;  // Karatsuba

// FP8: sqrt(moduli)
template <int IDX> inline constexpr int sqrt_moduli = 33;
template <> inline constexpr int sqrt_moduli<0>     = 33;
template <> inline constexpr int sqrt_moduli<1>     = 32;
template <> inline constexpr int sqrt_moduli<2>     = 31;
template <> inline constexpr int sqrt_moduli<3>     = 29;
template <> inline constexpr int sqrt_moduli<4>     = 25;
template <> inline constexpr int sqrt_moduli<5>     = 23;

inline constexpr int not_Karatsuba = 6;

//==========
// number of matrices for workspace of A/B
//==========
template <gemmul8::Backend backend> constexpr unsigned num_mat(unsigned num_moduli) {
    if constexpr (backend == gemmul8::Backend::INT8) return num_moduli;
    else {
        if (num_moduli <= not_Karatsuba) return 2 * num_moduli;
        else return 2 * not_Karatsuba + 3 * (num_moduli - 5);
    }
};

//==========
// P[i] = -1 * prod(p[0],...,p[i+1]) in double-double
//==========
namespace INT8 {
constexpr double2 P[19] = {
    { -0x1.fe00000000000p+15,   0x0.0000000000000p+0}, // -p[0]*p[1]
    { -0x1.f806000000000p+23,   0x0.0000000000000p+0}, // -p[0]*p[1]*p[2]
    { -0x1.ee2de20000000p+31,   0x0.0000000000000p+0}, // -p[0]*p[1]*p[2]*p[3]
    { -0x1.dcce450e00000p+39,   0x0.0000000000000p+0},
    { -0x1.c0de2f022e000p+47,   0x0.0000000000000p+0},
    { -0x1.a30f6de308f20p+55,   0x0.0000000000000p+0},
    { -0x1.7d690b03a3244p+63,  -0x1.0000000000000p+8},
    { -0x1.552ef6da40ef7p+71,  0x1.ec00000000000p+14},
    { -0x1.2e88a4e387945p+79,  0x1.1444000000000p+22},
    { -0x1.078907a2331a3p+87, -0x1.37ac620000000p+31},
    { -0x1.bec64ef0faa26p+94, -0x1.dc188f8900000p+40},
    {-0x1.703d73109e93ep+102,  0x1.2f97c1b215000p+48},
    {-0x1.1e3fc471eb44fp+110,  0x1.1ff7bc8b72980p+53},
    {-0x1.b88e245754182p+117,  0x1.df666905d3cbcp+63},
    {-0x1.4c232965d6663p+125,  0x1.616c352d64acap+71},
    {-0x1.ef9c77c5f5ec7p+132, -0x1.b141114c878cep+77},
    {-0x1.5e69a0aef6e03p+140,  0x1.35acfec4e4296p+85},
    {-0x1.ea07b6b4ad3d8p+147,  0x1.087f623ab88f0p+89},
    {-0x1.4b27367819129p+155,  0x1.595f0ab0d75c5p+98}
};
}

namespace FP8 {
constexpr double2 P[19] = {
    { -0x1.1040000000000p+20,    0x0.0000000000000p+0}, // -p[0]*p[1]
    { -0x1.ff00200000000p+29,    0x0.0000000000000p+0}, // -p[0]*p[1]*p[2]
    { -0x1.a3adda4800000p+39,    0x0.0000000000000p+0}, // -p[0]*p[1]*p[2]*p[3]
    { -0x1.0026dc7a72000p+49,    0x0.0000000000000p+0},
    { -0x1.08a826cc82c90p+58,    0x0.0000000000000p+0},
    { -0x1.f143f0e641bbbp+66,   0x1.c000000000000p+12},
    { -0x1.d1370fdf6a7f1p+75,  -0x1.3700000000000p+18},
    { -0x1.aa24f00a270d6p+84,  -0x1.2d9c300000000p+30},
    { -0x1.84b0b0f1429ebp+93,  -0x1.d719f5c800000p+39},
    {-0x1.5f7dc8022bbe8p+102,  -0x1.b403f9c25c000p+48},
    {-0x1.3c7ac095f4631p+111,   0x1.ced4d73d00540p+56},
    {-0x1.1a7b90e5d8a27p+120,   0x1.dc8e7d0ef9658p+66},
    {-0x1.ef72b92320f4ep+128,   0x1.a7d5e957436b0p+74},
    {-0x1.acadc32fe503ep+137,   0x1.7ab7956500d51p+83},
    {-0x1.6f8efcdb90dcdp+146,   0x1.4170d130346d6p+91},
    {-0x1.36d86cd7b002cp+155,  0x1.0debf474a22b4p+101},
    {-0x1.05ab2f9f90aa5p+164,  0x1.767040905d06dp+109},
    {-0x1.ae52855168e81p+172, -0x1.80e5974a5c0f1p+115},
    {-0x1.6028881a1f59fp+181,  0x1.92a0839614b53p+127},
};
}

template <gemmul8::Backend backend, typename doublex_t> __forceinline__ doublex_t get_P(unsigned num_moduli);
template <> __forceinline__ double get_P<gemmul8::Backend::INT8, double>(unsigned num_moduli) { return INT8::P[num_moduli - 2].x; }
template <> __forceinline__ double2 get_P<gemmul8::Backend::INT8, double2>(unsigned num_moduli) { return INT8::P[num_moduli - 2]; }
template <> __forceinline__ double get_P<gemmul8::Backend::FP8, double>(unsigned num_moduli) { return FP8::P[num_moduli - 2].x; }
template <> __forceinline__ double2 get_P<gemmul8::Backend::FP8, double2>(unsigned num_moduli) { return FP8::P[num_moduli - 2]; }

//==========
// invP[i] = 1/P[i] in double
//==========
namespace INT8 {
constexpr double invP[19] = {
    0x1.0101010101010p-16, 0x1.040d287a7051fp-24, 0x1.093b510fbf0d4p-32, 0x1.12e5617d255d8p-40, 0x1.2401777d7fdb6p-48,
    0x1.38c6a8b145786p-56, 0x1.57a6a12c3f24ap-64, 0x1.802b2f252aa3fp-72, 0x1.b13f5ca3b64a6p-80, 0x1.f15c410568cccp-88,
    0x1.255fb5199b040p-95, 0x1.63f115f5d0b39p-103, 0x1.c9e518641aa18p-111, 0x1.2983f5dbae8acp-118, 0x1.8aa1c572fa163p-126,
    0x1.0877227a9f8e3p-133, 0x1.760ceb764616fp-141, 0x1.0b7a38d26e2fep-148, 0x1.8bce042d07acep-156};
}

namespace FP8 {
constexpr double invP[19] = {
    0x1.e1709a3611655p-21, 0x1.0080301005018p-30, 0x1.3850970ef07b9p-40, 0x1.ffb252d5b63e6p-50, 0x1.ef40ad1677488p-59,
    0x1.0795ea39ba6dep-67, 0x1.19beb4e250a04p-76, 0x1.33939a58c52fcp-85, 0x1.5136ee4a4cf32p-94, 0x1.74e70ad38bd50p-103,
    0x1.9e280794e0291p-112, 0x1.d000087e75d11p-121, 0x1.088d6ae69afa5p-129, 0x1.31c21260a09fep-138, 0x1.649a089aae815p-147,
    0x1.a5a9b895cb631p-156, 0x1.f4e880fdf9551p-165, 0x1.30971bf77899ep-173, 0x1.74323bd5cea26p-182};
}

template <gemmul8::Backend backend> __forceinline__ double get_invP(unsigned num_moduli) {
    if constexpr (backend == gemmul8::Backend::INT8) return INT8::invP[num_moduli - 2];
    else return FP8::invP[num_moduli - 2];
}

//==========
// log2P[i] = round-down( log2(P-1)/2 - 0.5 ) in float
//==========
template <gemmul8::Backend backend, int num_moduli> inline constexpr float log2P = 0.0F;

// INT8
template <> inline constexpr float log2P<gemmul8::Backend::INT8, 2>  = 0x1.dfd1ec0000000p+2F;
template <> inline constexpr float log2P<gemmul8::Backend::INT8, 3>  = 0x1.6fa3360000000p+3F;
template <> inline constexpr float log2P<gemmul8::Backend::INT8, 4>  = 0x1.ef2ea60000000p+3F;
template <> inline constexpr float log2P<gemmul8::Backend::INT8, 5>  = 0x1.372d940000000p+4F;
template <> inline constexpr float log2P<gemmul8::Backend::INT8, 6>  = 0x1.767b2e0000000p+4F;
template <> inline constexpr float log2P<gemmul8::Backend::INT8, 7>  = 0x1.b5b0280000000p+4F;
template <> inline constexpr float log2P<gemmul8::Backend::INT8, 8>  = 0x1.f49a020000000p+4F;
template <> inline constexpr float log2P<gemmul8::Backend::INT8, 9>  = 0x1.19a8580000000p+5F;
template <> inline constexpr float log2P<gemmul8::Backend::INT8, 10> = 0x1.38f6bc0000000p+5F;
template <> inline constexpr float log2P<gemmul8::Backend::INT8, 11> = 0x1.582ada0000000p+5F;
template <> inline constexpr float log2P<gemmul8::Backend::INT8, 12> = 0x1.7736ae0000000p+5F;
template <> inline constexpr float log2P<gemmul8::Backend::INT8, 13> = 0x1.9619160000000p+5F;
template <> inline constexpr float log2P<gemmul8::Backend::INT8, 14> = 0x1.b4a4fe0000000p+5F;
template <> inline constexpr float log2P<gemmul8::Backend::INT8, 15> = 0x1.d321f80000000p+5F;
template <> inline constexpr float log2P<gemmul8::Backend::INT8, 16> = 0x1.f180a60000000p+5F;
template <> inline constexpr float log2P<gemmul8::Backend::INT8, 17> = 0x1.07e7f80000000p+6F;
template <> inline constexpr float log2P<gemmul8::Backend::INT8, 18> = 0x1.16e7e20000000p+6F;
template <> inline constexpr float log2P<gemmul8::Backend::INT8, 19> = 0x1.25df9a0000000p+6F;
template <> inline constexpr float log2P<gemmul8::Backend::INT8, 20> = 0x1.34be220000000p+6F;

// FP8
template <> inline constexpr float log2P<gemmul8::Backend::FP8, 2>  = 0x1.316bae0000000p+3F;
template <> inline constexpr float log2P<gemmul8::Backend::FP8, 3>  = 0x1.cff4720000000p+3F;
template <> inline constexpr float log2P<gemmul8::Backend::FP8, 4>  = 0x1.35b4840000000p+4F;
template <> inline constexpr float log2P<gemmul8::Backend::FP8, 5>  = 0x1.8001c00000000p+4F;
template <> inline constexpr float log2P<gemmul8::Backend::FP8, 6>  = 0x1.c862420000000p+4F;
template <> inline constexpr float log2P<gemmul8::Backend::FP8, 7>  = 0x1.07d4dc0000000p+5F;
template <> inline constexpr float log2P<gemmul8::Backend::FP8, 8>  = 0x1.2b726e0000000p+5F;
template <> inline constexpr float log2P<gemmul8::Backend::FP8, 9>  = 0x1.4ef0d60000000p+5F;
template <> inline constexpr float log2P<gemmul8::Backend::FP8, 10> = 0x1.7268ee0000000p+5F;
template <> inline constexpr float log2P<gemmul8::Backend::FP8, 11> = 0x1.95d4520000000p+5F;
template <> inline constexpr float log2P<gemmul8::Backend::FP8, 12> = 0x1.b9394e0000000p+5F;
template <> inline constexpr float log2P<gemmul8::Backend::FP8, 13> = 0x1.dc916c0000000p+5F;
template <> inline constexpr float log2P<gemmul8::Backend::FP8, 14> = 0x1.ffcf720000000p+5F;
template <> inline constexpr float log2P<gemmul8::Backend::FP8, 15> = 0x1.117ccc0000000p+6F;
template <> inline constexpr float log2P<gemmul8::Backend::FP8, 16> = 0x1.230b2c0000000p+6F;
template <> inline constexpr float log2P<gemmul8::Backend::FP8, 17> = 0x1.348f620000000p+6F;
template <> inline constexpr float log2P<gemmul8::Backend::FP8, 18> = 0x1.46102c0000000p+6F;
template <> inline constexpr float log2P<gemmul8::Backend::FP8, 19> = 0x1.577fa00000000p+6F;
template <> inline constexpr float log2P<gemmul8::Backend::FP8, 20> = 0x1.68eb8e0000000p+6F;

//==========
// 2^i mod moduli (constant memory)
//==========
// INT8: mod_pow2[i][j] = mod(2^j, p[i+1])
namespace INT8 {
__constant__ __device__ int8_t mod_pow2[19][57];
constexpr int8_t mod_pow2_h[19][57] = {
    {-127,  1,   2,    4,    8,   16,  32,   64, -127,   1,    2,    4,   8,   16,   32,   64, -127,    1,   2,   4,    8,   16,   32,  64, -127,    1,    2,   4,    8,   16,  32,  64, -127,   1,   2,    4,    8,   16,   32,   64, -127,    1,    2,    4,   8,  16,  32,  64, -127,    1,    2,    4,    8,   16,   32,  64, -127},
    {-125,  3,   6,   12,   24,   48,  96,  -61, -122,   9,   18,   36,  72, -109,   35,   70, -113,   27,  54, 108,  -37,  -74,  105, -43,  -86,   81,  -91,  71, -111,   31,  62, 124,   -5, -10, -20,  -40,  -80,   93,  -67,  119,  -15,  -30,  -60, -120,  13,  26,  52, 104,  -45,  -90,   73, -107,   39,   78,  -97,  59,  118},
    {-123,  5,  10,   20,   40,   80, -91,   69, -113,  25,   50,  100, -51, -102,   47,   94,  -63,  125,  -1,  -2,   -4,   -8,  -16, -32,  -64,  123,   -5, -10,  -20,  -40, -80,  91,  -69, 113, -25,  -50, -100,   51,  102,  -47,  -94,   63, -125,    1,   2,   4,   8,  16,   32,   64, -123,    5,   10,   20,   40,  80,  -91},
    {-119,  9,  18,   36,   72, -103,  41,   82,  -83,  81,  -85,   77, -93,   61,  122,   -3,   -6,  -12, -24, -48,  -96,   55,  110, -27,  -54, -108,   31,  62, -123,    1,   2,   4,    8,  16,  32,   64, -119,    9,   18,   36,   72, -103,   41,   82, -83,  81, -85,  77,  -93,   61,  122,   -3,   -6,  -12,  -24, -48,  -96},
    {-113, 15,  30,   60,  120,   -1,  -2,   -4,   -8, -16,  -32,  -64, 113,  -15,  -30,  -60, -120,    1,   2,   4,    8,   16,   32,  64, -113,   15,   30,  60,  120,   -1,  -2,  -4,   -8, -16, -32,  -64,  113,  -15,  -30,  -60, -120,    1,    2,    4,   8,  16,  32,  64, -113,   15,   30,   60,  120,   -1,   -2,  -4,   -8},
    {-111, 17,  34,   68, -103,   33,  66, -107,   25,  50,  100,  -39, -78,   83,  -73,   93,  -53, -106,  27,  54,  108,  -23,  -46, -92,   55,  110,  -19, -38,  -76,   87, -65, 109,  -21, -42, -84,   71,  -97,   45,   90,  -59, -118,    3,    6,   12,  24,  48,  96, -47,  -94,   51,  102,  -35,  -70,   99,  -41, -82,   75},
    {-105, 23,  46,   92,  -49,  -98,  37,   74,  -85,  63, -107,   19,  38,   76,  -81,   71,  -91,   51, 102, -29,  -58, -116,    1,   2,    4,    8,   16,  32,   64, -105,  23,  46,   92, -49, -98,   37,   74,  -85,   63, -107,   19,   38,   76,  -81,  71, -91,  51, 102,  -29,  -58, -116,    1,    2,    4,    8,  16,   32},
    {-101, 27,  54,  108,  -13,  -26, -52, -104,   21,  42,   84,  -61, 107,  -15,  -30,  -60,  109,  -11, -22, -44,  -88,   53,  106, -17,  -34,  -68,   93, -43,  -86,   57, 114,  -1,   -2,  -4,  -8,  -16,  -32,  -64,  101,  -27,  -54, -108,   13,   26,  52, 104, -21, -42,  -84,   61, -107,   15,   30,   60, -109,  11,   22},
    { -99, 29,  58, -111,    5,   10,  20,   40,   80, -67,   93,  -41, -82,   63, -101,   25,   50,  100, -27, -54, -108,   11,   22,  44,   88,  -51, -102,  23,   46,   92, -43, -86,   55, 110,  -7,  -14,  -28,  -56, -112,    3,    6,   12,   24,   48,  96, -35, -70,  87,  -53, -106,   15,   30,   60, -107,   13,  26,   52},
    { -95, 33,  66,  -91,   41,   82, -59,  105,  -13, -26,  -52, -104,  15,   30,   60, -103,   17,   34,  68, -87,   49,   98,  -27, -54, -108,    7,   14,  28,   56, -111,   1,   2,    4,   8,  16,   32,   64,  -95,   33,   66,  -91,   41,   82,  -59, 105, -13, -26, -52, -104,   15,   30,   60, -103,   17,   34,  68,  -87},
    { -89, 39,  78,  -61,   95,  -27, -54, -108,    1,   2,    4,    8,  16,   32,   64,  -89,   39,   78, -61,  95,  -27,  -54, -108,   1,    2,    4,    8,  16,   32,   64, -89,  39,   78, -61,  95,  -27,  -54, -108,    1,    2,    4,    8,   16,   32,  64, -89,  39,  78,  -61,   95,  -27,  -54, -108,    1,    2,   4,    8},
    { -83, 45,  90,  -31,  -62,   87, -37,  -74,   63, -85,   41,   82, -47,  -94,   23,   46,   92,  -27, -54, 103,   -5,  -10,  -20, -40,  -80,   51,  102,  -7,  -14,  -28, -56,  99,  -13, -26, -52, -104,    3,    6,   12,   24,   48,   96,  -19,  -38, -76,  59, -93,  25,   50,  100,  -11,  -22,  -44,  -88,   35,  70,  -71},
    { -71, 57, -85,   29,   58,  -83,  33,   66,  -67,  65,  -69,   61, -77,   45,   90,  -19,  -38,  -76,  47,  94,  -11,  -22,  -44, -88,   23,   46,   92, -15,  -30,  -60,  79, -41,  -82,  35,  70,  -59,   81,  -37,  -74,   51,  -97,    5,   10,   20,  40,  80, -39, -78,   43,   86,  -27,  -54,   91,  -17,  -34, -68,   63},
    { -69, 59, -79,   39,   78,  -41, -82,   33,   66, -65,   67,  -63,  71,  -55,   87,  -23,  -46,  -92,  13,  26,   52,  -93,   11,  22,   44,   88,  -21, -42,  -84,   29,  58, -81,   35,  70, -57,   83,  -31,  -62,   73,  -51,   95,   -7,  -14,  -28, -56,  85, -27, -54,   89,  -19,  -38,  -76,   45,   90,  -17, -34,  -68},
    { -65, 63, -67,   59,  -75,   43,  86,  -21,  -42, -84,   25,   50, -93,    7,   14,   28,   56,  -81,  31,  62,  -69,   55,  -83,  27,   54,  -85,   23,  46,   92,   -9, -18, -36,  -72,  49, -95,    3,    6,   12,   24,   48,   96,   -1,   -2,   -4,  -8, -16, -32, -64,   65,  -63,   67,  -59,   75,  -43,  -86,  21,   42},
    { -63, 65, -61,   69,  -53,   85, -21,  -42,  -84,  23,   46,   92,  -7,  -14,  -28,  -56,   79,  -33, -66,  59,  -73,   45,   90, -11,  -22,  -44,  -88,  15,   30,   60, -71,  49,  -93,   5,  10,   20,   40,   80,  -31,  -62,   67,  -57,   77,  -37, -74,  43,  86, -19,  -38,  -76,   39,   78,  -35,  -70,   51, -89,   13},
    { -53, 75, -31,  -62,   57,  -67,  47,  -87,    7,  14,   28,   56, -69,   43,   86,   -9,  -18,  -36, -72,  37,   74,  -33,  -66,  49,  -83,   15,   30,  60,  -61,   59, -63,  55,  -71,  39,  78,  -25,  -50,   81,  -19,  -38,  -76,   29,   58,  -65,  51, -79,  23,  46,  -89,    3,    6,   12,   24,   48,  -85,  11,   22},
    { -51, 77, -25,  -50,   79,  -21, -42,  -84,   11,  22,   44,   88,  -3,   -6,  -12,  -24,  -48,   83, -13, -26,  -52,   75,  -29, -58,   63,  -53,   73, -33,  -66,   47, -85,   9,   18,  36,  72,  -35,  -70,   39,   78,  -23,  -46,   87,   -5,  -10, -20, -40, -80,  19,   38,   76,  -27,  -54,   71,  -37,  -74,  31,   62},
    { -45, 83,  -7,  -14,  -28,  -56,  61,  -51,   71, -31,  -62,   49, -75,   23,   46,  -81,   11,   22,  44, -85,    3,    6,   12,  24,   48,  -77,   19,  38,   76,  -21, -42, -84,    5,  10,  20,   40,   80,  -13,  -26,  -52,   69,  -35,  -70,   33,  66, -41, -82,   9,   18,   36,   72,  -29,  -58,   57,  -59,  55,  -63},
};
} // namespace INT8

// FP8: mod_pow2[i][j] = mod(2^j, p[i+1])  (mod_pow2[0][j] = mod(2^j, 1089))
namespace FP8 {
__constant__ __device__ int16_t mod_pow2[19][64];
constexpr int16_t mod_pow2_h[19][64] = {
    { 256,  512,  -65, -130, -260, -520,   49,   98,  196,  392, -305,  479, -131, -262, -524,   41,   82,  164,  328, -433,  223,  446, -197, -394,  301, -487,  115,  230,  460, -169, -338,  413, -263, -526,   37,   74,  148,  296, -497,   95,  190,  380, -329,  431, -227, -454,  181,  362, -365,  359, -371,  347, -395,  299, -491,  107,  214,  428, -233, -466,  157,  314, -461,  167},
    { 256, -449,   63,  126,  252, -457,   47,   94,  188,  376, -209, -418,  125,  250, -461,   39,   78,  156,  312, -337,  287, -387,  187,  374, -213, -426,  109,  218,  436,  -89, -178, -356,  249, -463,   35,   70,  140,  280, -401,  159,  318, -325,  311, -339,  283, -395,  171,  342, -277,  407, -147, -294,  373, -215, -430,  101,  202,  404, -153, -306,  349, -263,  435,  -91},
    { 256, -329,  183,  366, -109, -218,  405,  -31,  -62, -124, -248,  345, -151, -302,  237, -367,  107,  214, -413,   15,   30,   60,  120,  240, -361,  119,  238, -365,  111,  222, -397,   47,   94,  188,  376,  -89, -178, -356,  129,  258, -325,  191,  382,  -77, -154, -308,  225, -391,   59,  118,  236, -369,  103,  206,  412,  -17,  -34,  -68, -136, -272,  297, -247,  347, -147},
    { 256, -113, -226,  173, -279,   67,  134,  268,  -89, -178,  269,  -87, -174,  277,  -71, -142, -284,   57,  114,  228, -169,  287,  -51, -102, -204,  217, -191,  243, -139, -278,   69,  138,  276,  -73, -146, -292,   41,   82,  164, -297,   31,   62,  124,  248, -129, -258,  109,  218, -189,  247, -131, -262,  101,  202, -221,  183, -259,  107,  214, -197,  231, -163,  299,  -27},
    { 256,  -17,  -34,  -68, -136,  257,  -15,  -30,  -60, -120, -240,   49,   98,  196, -137,  255,  -19,  -38,  -76, -152,  225,  -79, -158,  213, -103, -206,  117,  234,  -61, -122, -244,   41,   82,  164, -201,  127,  254,  -21,  -42,  -84, -168,  193, -143,  243,  -43,  -86, -172,  185, -159,  211, -107, -214,  101,  202, -125, -250,   29,   58,  116,  232,  -65, -130, -260,    9},
    {-225,   31,   62,  124, -233,   15,   30,   60,  120,  240,   -1,   -2,   -4,   -8,  -16,  -32,  -64, -128,  225,  -31,  -62, -124,  233,  -15,  -30,  -60, -120, -240,    1,    2,    4,    8,   16,   32,   64,  128, -225,   31,   62,  124, -233,   15,   30,   60,  120,  240,   -1,   -2,   -4,   -8,  -16,  -32,  -64, -128,  225,  -31,  -62, -124,  233,  -15,  -30,  -60, -120, -240},
    {-223,   33,   66,  132, -215,   49,   98,  196,  -87, -174,  131, -217,   45,   90,  180, -119, -238,    3,    6,   12,   24,   48,   96,  192,  -95, -190,   99,  198,  -83, -166,  147, -185,  109,  218,  -43,  -86, -172,  135, -209,   61,  122, -235,    9,   18,   36,   72,  144, -191,   97,  194,  -91, -182,  115,  230,  -19,  -38,  -76, -152,  175, -129,  221,  -37,  -74, -148},
    {-213,   43,   86,  172, -125,  219,  -31,  -62, -124,  221,  -27,  -54, -108, -216,   37,   74,  148, -173,  123, -223,   23,   46,   92,  184, -101, -202,   65,  130, -209,   51,  102,  204,  -61, -122,  225,  -19,  -38,  -76, -152,  165, -139,  191,  -87, -174,  121, -227,   15,   30,   60,  120, -229,   11,   22,   44,   88,  176, -117, -234,    1,    2,    4,    8,   16,   32},
    {-211,   45,   90,  180, -107, -214,   39,   78,  156, -155,  157, -153,  161, -145,  177, -113, -226,   15,   30,   60,  120, -227,   13,   26,   52,  104,  208,  -51, -102, -204,   59,  118, -231,    5,   10,   20,   40,   80,  160, -147,  173, -121,  225,  -17,  -34,  -68, -136,  195,  -77, -154,  159, -149,  169, -129,  209,  -49,  -98, -196,   75,  150, -167,  133, -201,   65},
    {-207,   49,   98,  196,  -71, -142,  179, -105, -210,   43,   86,  172, -119,  225,  -13,  -26,  -52, -104, -208,   47,   94,  188,  -87, -174,  115,  230,   -3,   -6,  -12,  -24,  -48,  -96, -192,   79,  158, -147,  169, -125,  213,  -37,  -74, -148,  167, -129,  205,  -53, -106, -212,   39,   78,  156, -151,  161, -141,  181, -101, -202,   59,  118, -227,    9,   18,   36,   72},
    {-205,   51,  102,  204,  -53, -106, -212,   37,   74,  148, -165,  131, -199,   63,  126, -209,   43,   86,  172, -117,  227,   -7,  -14,  -28,  -56, -112, -224,   13,   26,   52,  104,  208,  -45,  -90, -180,  101,  202,  -57, -114, -228,    5,   10,   20,   40,   80,  160, -141,  179, -103, -206,   49,   98,  196,  -69, -138,  185,  -91, -182,   97,  194,  -73, -146,  169, -123},
    {-201,   55,  110,  220,  -17,  -34,  -68, -136,  185,  -87, -174,  109,  218,  -21,  -42,  -84, -168,  121, -215,   27,   54,  108,  216,  -25,  -50, -100, -200,   57,  114,  228,   -1,   -2,   -4,   -8,  -16,  -32,  -64, -128,  201,  -55, -110, -220,   17,   34,   68,  136, -185,   87,  174, -109, -218,   21,   42,   84,  168, -121,  215,  -27,  -54, -108, -216,   25,   50,  100},
    {-193,   63,  126, -197,   55,  110,  220,   -9,  -18,  -36,  -72, -144,  161, -127,  195,  -59, -118,  213,  -23,  -46,  -92, -184,   81,  162, -125,  199,  -51, -102, -204,   41,   82,  164, -121,  207,  -35,  -70, -140,  169, -111, -222,    5,   10,   20,   40,   80,  160, -129,  191,  -67, -134,  181,  -87, -174,  101,  202,  -45,  -90, -180,   89,  178,  -93, -186,   77,  154},
    {-187,   69,  138, -167,  109,  218,   -7,  -14,  -28,  -56, -112,  219,   -5,  -10,  -20,  -40,  -80, -160,  123, -197,   49,   98,  196,  -51, -102, -204,   35,   70,  140, -163,  117, -209,   25,   50,  100,  200,  -43,  -86, -172,   99,  198,  -47,  -94, -188,   67,  134, -175,   93,  186,  -71, -142,  159, -125,  193,  -57, -114,  215,  -13,  -26,  -52, -104, -208,   27,   54},
    {-183,   73,  146, -147,  145, -149,  141, -157,  125, -189,   61,  122, -195,   49,   98,  196,  -47,  -94, -188,   63,  126, -187,   65,  130, -179,   81,  162, -115,  209,  -21,  -42,  -84, -168,  103,  206,  -27,  -54, -108, -216,    7,   14,   28,   56,  112, -215,    9,   18,   36,   72,  144, -151,  137, -165,  109,  218,   -3,   -6,  -12,  -24,  -48,  -96, -192,   55,  110},
    {-177,   79,  158, -117,  199,  -35,  -70, -140,  153, -127,  179,  -75, -150,  133, -167,   99,  198,  -37,  -74, -148,  137, -159,  115, -203,   27,   54,  108,  216,   -1,   -2,   -4,   -8,  -16,  -32,  -64, -128,  177,  -79, -158,  117, -199,   35,   70,  140, -153,  127, -179,   75,  150, -133,  167,  -99, -198,   37,   74,  148, -137,  159, -115,  203,  -27,  -54, -108, -216},
    {-175,   81,  162, -107, -214,    3,    6,   12,   24,   48,   96,  192,  -47,  -94, -188,   55,  110, -211,    9,   18,   36,   72,  144, -143,  145, -141,  149, -133,  165, -101, -202,   27,   54,  108, -215,    1,    2,    4,    8,   16,   32,   64,  128, -175,   81,  162, -107, -214,    3,    6,   12,   24,   48,   96,  192,  -47,  -94, -188,   55,  110, -211,    9,   18,   36},
    {-165,   91,  182,  -57, -114,  193,  -35,  -70, -140,  141, -139,  143, -135,  151, -119,  183,  -55, -110,  201,  -19,  -38,  -76, -152,  117, -187,   47,   94,  188,  -45,  -90, -180,   61,  122, -177,   67,  134, -153,  115, -191,   39,   78,  156, -109,  203,  -15,  -30,  -60, -120,  181,  -59, -118,  185,  -51, -102, -204,   13,   26,   52,  104,  208,   -5,  -10,  -20,  -40},
    {-163,   93,  186,  -47,  -94, -188,   43,   86,  172,  -75, -150,  119, -181,   57,  114, -191,   37,   74,  148, -123,  173,  -73, -146,  127, -165,   89,  178,  -63, -126,  167,  -85, -170,   79,  158, -103, -206,    7,   14,   28,   56,  112, -195,   29,   58,  116, -187,   45,   90,  180,  -59, -118,  183,  -53, -106,  207,   -5,  -10,  -20,  -40,  -80, -160,   99,  198,  -23},
};
} // namespace FP8

template <gemmul8::Backend backend, int IDX> __device__ __forceinline__ int get_mod_pow2(int exp) {
    if constexpr (backend == gemmul8::Backend::INT8) {
        if constexpr (0 < IDX && IDX < 20) return (exp <= 6) ? (1 << exp) : INT8::mod_pow2[IDX - 1][exp - 7];
        return 0;
    } else {
        if constexpr (IDX == 0) return (exp <= 7) ? (1 << exp) : FP8::mod_pow2[0][exp - 8];
        if constexpr (1 < IDX && IDX < 20) return (exp <= 7) ? (1 << exp) : FP8::mod_pow2[IDX - 1][exp - 8];
        return 0;
    }
}

//==========
// q[i]*P[i]/p[i]
//==========
namespace INT8 {

// qPi_1[i] = double(q[i]*P[i]/p[i]), where q[i]*P[i]/p[i] == 1 mod p[i]
inline constexpr double qPi_1[19][20] = {
    {0x1.fc02000000000p+15, 0x1.0000000000000p+8},
    {0x1.50ac020000000p+23, 0x1.f60c000000000p+22, 0x1.a45a000000000p+23},
    {0x1.0688601000000p+28, 0x1.f01e000000000p+28, 0x1.4826900000000p+28, 0x1.6654440000000p+31},
    {0x1.99c1435808000p+37, 0x1.d553914600000p+39, 0x1.cf9d0d8400000p+38, 0x1.2ff09e4000000p+38, 0x1.dae0172c00000p+39},
    {0x1.24d0f0aa6c020p+47, 0x1.00ffb685c4000p+47, 0x1.7820600df8000p+45, 0x1.b28fb528de000p+47, 0x1.765c060a1c000p+47,
     0x1.56b441a210000p+47},
    {0x1.49071d4742060p+55, 0x1.5fae947039b40p+55, 0x1.42fdb9e1948e0p+55, 0x1.187c8ee783700p+55, 0x1.e89ef222a1c00p+52,
     0x1.0316493fe27a0p+55, 0x1.1f8e561d65780p+53},
    {0x1.4f3952ae3262ep+63, 0x1.f094cf17cf626p+61, 0x1.0f5bef8d36588p+63, 0x1.e02e9274c53aep+62, 0x1.a403bd5c1a42ep+61,
     0x1.a1cf7b99c2a51p+62, 0x1.a54e8a8f43bcfp+60, 0x1.787fdcb9fa097p+62},
    {0x1.9a7c80fe96201p+69, 0x1.43ca2f89db3d9p+71, 0x1.40f4871424cd0p+70, 0x1.2c6790ef157a1p+71, 0x1.24d66e4d76f4ep+70,
     0x1.459c5b1ee5ce3p+71, 0x1.d43c2b2519eb7p+70, 0x1.ab93da2aca3c9p+70, 0x1.dfbe1fda9333ap+70},
    {0x1.1ba01a954f1b1p+75, 0x1.b499060d20053p+76, 0x1.8d00367a835e7p+77, 0x1.348f721e1e2b2p+77, 0x1.09c9ed1acf35ap+79,
     0x1.6988bc8c2f4c5p+75, 0x1.4e2df779b91bdp+77, 0x1.54302cc6b737dp+78, 0x1.675767107d43cp+76, 0x1.1fdfa04826ca0p+77},
    {0x1.ae4dbe76d770cp+86, 0x1.258185fdee9fbp+86, 0x1.76fdabbf55de7p+85, 0x1.73ade1f8235b0p+86, 0x1.0cdeb7fb81deap+85,
     0x1.0671178918559p+87, 0x1.c416fd07412bep+86, 0x1.5350d862f82efp+86, 0x1.52567e0ff5970p+86, 0x1.d0611c1cafc1ep+85,
     0x1.814201f9bea6ep+86},
    {0x1.42dd4f0c251f6p+94, 0x1.71af2232d1654p+94, 0x1.b5f1f25063f94p+93, 0x1.0e8e8784ac2d5p+93, 0x1.0477c23ba5cfap+93,
     0x1.ac3c7c8760d50p+94, 0x1.507ba57edce57p+92, 0x1.2b20ca473f6ddp+93, 0x1.5f2d33fd22e8cp+92, 0x1.ab17cae65cfc4p+94,
     0x1.408e48b61567ep+90, 0x1.32c582e2cf7c8p+94},
    {0x1.187ecea5a8caap+102, 0x1.71af2232d1654p+94, 0x1.5a685a078a0bap+102, 0x1.48a0e93cba656p+102, 0x1.6d422253da718p+102,
     0x1.ec015f50a0e35p+101, 0x1.27d31b1922346p+99, 0x1.7b4d942fe1f7ep+100, 0x1.68332a1fe8402p+101, 0x1.7859de7afe8bdp+99,
     0x1.317d98db46b66p+102, 0x1.08b9be1306a54p+102, 0x1.411e88bd3424cp+100},
    {0x1.4af9bb23b807bp+107, 0x1.e0730f7df34a9p+109, 0x1.9e197740a03d4p+109, 0x1.11b44daf38ac0p+106, 0x1.959dba1ed526ap+109,
     0x1.d3f9c70059540p+109, 0x1.c71fc396104fap+108, 0x1.6e1a9ef49547bp+109, 0x1.067fc962e0b14p+110, 0x1.81de6aed04c27p+109,
     0x1.086d6ad9bca27p+110, 0x1.66ccfaf43fac8p+109, 0x1.d2ae54e5674d3p+109, 0x1.98842ba66fec0p+109},
    {0x1.8334edf0c0e93p+117, 0x1.d9618469e1e3bp+116, 0x1.4c97d49af8a7fp+117, 0x1.3db0f47816e79p+117, 0x1.ac11e30d56e4ap+116,
     0x1.d3f9c70059540p+109, 0x1.0210da6024681p+117, 0x1.2e86f6e52b76cp+116, 0x1.f43197eee3640p+115, 0x1.e913152bf1176p+115,
     0x1.775c686f240ffp+116, 0x1.44d556f6112cep+116, 0x1.90e26770391c0p+115, 0x1.1b5f498bca07cp+117, 0x1.9702ab51fa860p+116},
    {0x1.568442b105196p+122, 0x1.23c286bfdb74ep+125, 0x1.fffd89ae2f2d3p+124, 0x1.9f80a3facf046p+124, 0x1.6b10abb2b052bp+124,
     0x1.b90322c9142e7p+119, 0x1.ff687bb9b984bp+124, 0x1.494950989a53bp+125, 0x1.5c176f941779ap+122, 0x1.6dca3fa2e79c8p+124,
     0x1.951e4290e0c61p+122, 0x1.a671255128495p+123, 0x1.b2745cf9aee44p+124, 0x1.2c6cfd90dba85p+123, 0x1.a57e7d4e8e218p+124,
     0x1.8f40d0ef2435dp+124},
    {0x1.e01f9407c63d1p+129, 0x1.e201959d63a0bp+131, 0x1.31982160c4289p+132, 0x1.7f0fe22eef086p+132, 0x1.00d5bf9f9747cp+126,
     0x1.8ad801f1a0de6p+129, 0x1.2a9c6628027f7p+130, 0x1.d836977997e90p+131, 0x1.85903a5f3c45ap+132, 0x1.a3320451ba942p+132,
     0x1.ce462d22425e4p+132, 0x1.d67cf11ca9c0ap+132, 0x1.add7c7ba400e9p+132, 0x1.57b0afae95f50p+131, 0x1.e30840c0efae9p+128,
     0x1.5aabc9d4bfea6p+132, 0x1.82a0ee308b92fp+132},
    {0x1.06cf388339282p+134, 0x1.a1bf2dfdc2ed2p+136, 0x1.bb35a9d83d503p+137, 0x1.b0c7cfa209212p+139, 0x1.4921eae073cd6p+140,
     0x1.172ab95fd6bd4p+139, 0x1.68acfd38e8af1p+139, 0x1.f34ce4f4e91c4p+138, 0x1.01123dfc720a3p+140, 0x1.9db3f738931cfp+139,
     0x1.f6d5907a7ef5fp+138, 0x1.e7abc6d98b7c7p+139, 0x1.8e92d6501c724p+136, 0x1.1d42b11e8381cp+140, 0x1.0579b3ad70bebp+140,
     0x1.0cb5cec87da54p+138, 0x1.2009162ca2b84p+140, 0x1.3d803cbad18b8p+140},
    {0x1.b09acf4b80f05p+146, 0x1.6f0acc1cea2b1p+147, 0x1.d8992594f3fa9p+145, 0x1.4be496434b847p+146, 0x1.8cc9189a96a3dp+147,
     0x1.c776b470b2041p+143, 0x1.cd534fe2dceefp+147, 0x1.82fa017336678p+147, 0x1.946f7304e8699p+147, 0x1.551407a0b7bc5p+147,
     0x1.034c6790f7cb9p+146, 0x1.452e68b9f5e92p+145, 0x1.407e5f3ab7ac7p+147, 0x1.a514c773601f0p+147, 0x1.840b4e6816d47p+147,
     0x1.7a503c2406688p+147, 0x1.9fa0adbac081ep+147, 0x1.0c070c3e0cb90p+147, 0x1.952a21ca4d733p+145},
    {0x1.b7d01457814cap+153, 0x1.22e534ddf3e42p+150, 0x1.157cefeb36669p+153, 0x1.3ca3f6e306a29p+151, 0x1.016a241f2b53ep+152,
     0x1.e66c961dd1f94p+154, 0x1.1945b982edaa0p+155, 0x1.3e5ca23c85f9bp+152, 0x1.1ce0e513790ffp+155, 0x1.98788b5ce0e66p+154,
     0x1.b19e1bb311e81p+154, 0x1.df2e1fa0ce290p+154, 0x1.8b801f14e4ebbp+153, 0x1.38d9254eb9c6fp+153, 0x1.354ce8cdbc742p+154,
     0x1.94eef4587e294p+154, 0x1.3b8c91b979bbep+155, 0x1.0b1e3740581d2p+155, 0x1.088d7305d7f68p+155, 0x1.67ddaa2caf393p+154},
};

// idx = num_moduli - threshold<gemmul8::Backend::FP8>::P_is_double - 1
// qPi_2[idx][i][1] = first (53-ceil(log2(rho))) bits of q[i]*P[i]/p[i] for rho = sum(floor(p[:]/2)),
// qPi_2[idx][i][2] = double(q[i]*P[i]/p[i] - qPi_2[idx][i][1])
inline constexpr double2 qPi_2[14][20] = {
    {
     {0x1.49071d4742000p+55, 0x1.8080000000000p+9},
     {0x1.5fae947039800p+55, 0x1.a000000000000p+12},
     {0x1.42fdb9e194800p+55, 0x1.c000000000000p+10},
     {0x1.187c8ee783400p+55, 0x1.8000000000000p+12},
     {0x1.e89ef222a0000p+52, 0x1.c000000000000p+12},
     {0x1.0316493fe2400p+55, 0x1.d000000000000p+12},
     {0x1.1f8e561d65000p+53, 0x1.e000000000000p+11},
     },
    {
     {0x1.4f3952ae32400p+63, 0x1.16f0100000000p+20},
     {0x1.f094cf17cf000p+61, 0x1.89a0000000000p+19},
     {0x1.0f5bef8d36400p+63, 0x1.8880000000000p+19},
     {0x1.e02e9274c5000p+62, 0x1.d740000000000p+19},
     {0x1.a403bd5c1a000p+61, 0x1.0b80000000000p+19},
     {0x1.a1cf7b99c2800p+62, 0x1.2880000000000p+19},
     {0x1.a54e8a8f42000p+60, 0x1.bcf0000000000p+20},
     {0x1.787fdcb9fa000p+62, 0x1.2d80000000000p+17},
     },
    {
     {0x1.9a7c80fe96000p+69, 0x1.008cc04000000p+26},
     {0x1.43ca2f89db000p+71, 0x1.eca4600000000p+28},
     {0x1.40f4871424000p+70, 0x1.9a00780000000p+29},
     {0x1.2c6790ef15000p+71, 0x1.e855180000000p+29},
     {0x1.24d66e4d76000p+70, 0x1.e9c7f00000000p+29},
     {0x1.459c5b1ee5800p+71, 0x1.38caf00000000p+29},
     {0x1.d43c2b2519000p+70, 0x1.d6d0600000000p+29},
     {0x1.ab93da2aca000p+70, 0x1.e459400000000p+27},
     {0x1.dfbe1fda93000p+70, 0x1.9cd8200000000p+27},
     },
    {
     {0x1.1ba01a9548000p+75, 0x1.c6c29fa008000p+37},
     {0x1.b499060d20000p+76, 0x1.4ddc380000000p+30},
     {0x1.8d00367a82000p+77, 0x1.5e72640800000p+37},
     {0x1.348f721e1e000p+77, 0x1.5939d00000000p+34},
     {0x1.09c9ed1acf000p+79, 0x1.acce161000000p+36},
     {0x1.6988bc8c28000p+75, 0x1.d3148d7000000p+37},
     {0x1.4e2df779b8000p+77, 0x1.1bca621000000p+37},
     {0x1.54302cc6b7000p+78, 0x1.be65b8a000000p+35},
     {0x1.675767107c000p+76, 0x1.43b8ee6000000p+36},
     {0x1.1fdfa04826000p+77, 0x1.940b60e000000p+36},
     },
    {
     {0x1.ae4dbe76d7000p+86, 0x1.c311739de0100p+44},
     {0x1.258185fdee000p+86, 0x1.3f5c901690000p+45},
     {0x1.76fdabbf54000p+85, 0x1.de7087e210000p+45},
     {0x1.73ade1f823000p+86, 0x1.6bfc28bd30000p+44},
     {0x1.0cdeb7fb80000p+85, 0x1.de9bee2d48000p+45},
     {0x1.0671178918000p+87, 0x1.5646b56780000p+45},
     {0x1.c416fd0741000p+86, 0x1.5ee3b89260000p+43},
     {0x1.5350d862f8000p+86, 0x1.77449328c0000p+43},
     {0x1.52567e0ff5000p+86, 0x1.2e0367d338000p+45},
     {0x1.d0611c1cae000p+85, 0x1.c1e3b22c60000p+45},
     {0x1.814201f9be000p+86, 0x1.4dba603168000p+45},
     },
    {
     {0x1.42dd4f0c25000p+94, 0x1.f5cc036fee804p+50},
     {0x1.71af2232d1000p+94, 0x1.9502088c71500p+52},
     {0x1.b5f1f25063000p+93, 0x1.f27fe97ac8c00p+52},
     {0x1.0e8e8784ac000p+93, 0x1.6a7fd4fb91000p+50},
     {0x1.0477c23ba5000p+93, 0x1.9f4e1d77bb800p+52},
     {0x1.ac3c7c8760800p+94, 0x1.541aed8de8f00p+52},
     {0x1.507ba57edc000p+92, 0x1.cad9eee787600p+51},
     {0x1.2b20ca473f000p+93, 0x1.b754bf1ae1c00p+51},
     {0x1.5f2d33fd22000p+92, 0x1.d1793fc3ce200p+51},
     {0x1.ab17cae65c800p+94, 0x1.f0f2278772b00p+52},
     {0x1.408e48b610000p+90, 0x1.59f94c68de600p+52},
     {0x1.32c582e2cf000p+94, 0x1.f1f52b3aa8500p+52},
     },
    {
     {0x1.187ecea5a8800p+102, 0x1.2a800bf67755ap+60},
     {0x1.71af223280000p+94, 0x1.459502088c715p+60},
     {0x1.5a685a078a000p+102, 0x1.73141ccb58410p+57},
     {0x1.48a0e93cba000p+102, 0x1.956a7a15d56e0p+60},
     {0x1.6d422253da000p+102, 0x1.c5f9191e4aa91p+60},
     {0x1.ec015f50a0000p+101, 0x1.c69660c475d7bp+60},
     {0x1.27d31b1920000p+99, 0x1.1a2f5dd7b0278p+60},
     {0x1.7b4d942fe0000p+100, 0x1.f7e2f13271df4p+60},
     {0x1.68332a1fe8000p+101, 0x1.008e6afbfbd20p+59},
     {0x1.7859de7afc000p+99, 0x1.45eaf7cf70b15p+60},
     {0x1.317d98db46800p+102, 0x1.b2d9a1321591ap+59},
     {0x1.08b9be1306800p+102, 0x1.2a3c04a60a8a8p+59},
     {0x1.411e88bd34000p+100, 0x1.25d2c634e54f0p+57},
     },
    {
     {0x1.4af9bb23b8000p+107, 0x1.ed366131bfd87p+61},
     {0x1.e0730f7df3000p+109, 0x1.2a2b688425b37p+67},
     {0x1.9e197740a0000p+109, 0x1.ea31249d190dbp+66},
     {0x1.11b44daf38000p+106, 0x1.57f8ce0e05580p+65},
     {0x1.959dba1ed5000p+109, 0x1.34d662a4fdd1cp+66},
     {0x1.d3f9c70059000p+109, 0x1.5013290076958p+67},
     {0x1.c71fc39610000p+108, 0x1.3e7351822d438p+66},
     {0x1.6e1a9ef495000p+109, 0x1.1ebdf25941f8bp+67},
     {0x1.067fc962e0800p+110, 0x1.89ef93ae85687p+67},
     {0x1.81de6aed04000p+109, 0x1.84d8fc60d93d4p+68},
     {0x1.086d6ad9bc800p+110, 0x1.1340b8f1c34bfp+67},
     {0x1.66ccfaf43f000p+109, 0x1.590ded2a35e12p+68},
     {0x1.d2ae54e567000p+109, 0x1.34cf70f07ae33p+67},
     {0x1.98842ba66f000p+109, 0x1.d80e799d28f38p+68},
     },
    {
     {0x1.8334edf0c0800p+117, 0x1.a4b62a6fdb1e1p+75},
     {0x1.d9618469e1000p+116, 0x1.c75bd2f612d0fp+75},
     {0x1.4c97d49af8800p+117, 0x1.3f90fd4ad5142p+74},
     {0x1.3db0f47816800p+117, 0x1.9e3c7f45d92bfp+75},
     {0x1.ac11e30d56000p+116, 0x1.c9413fbd969ffp+75},
     {0x1.d3f9c70000000p+109, 0x1.6550132900769p+75},
     {0x1.0210da6024000p+117, 0x1.a05bb4379a4c1p+75},
     {0x1.2e86f6e52b000p+116, 0x1.dae820f5ffc00p+74},
     {0x1.f43197eee2000p+115, 0x1.6405781ac87d9p+75},
     {0x1.e913152bf0000p+115, 0x1.175dbeffba9cdp+75},
     {0x1.775c686f24000p+116, 0x1.fe04b43a93e73p+71},
     {0x1.44d556f611000p+116, 0x1.67335f8e813b8p+73},
     {0x1.90e2677038000p+115, 0x1.1bfe09769edb0p+75},
     {0x1.1b5f498bca000p+117, 0x1.f1ee6037f1f5dp+71},
     {0x1.9702ab51fa000p+116, 0x1.0c08e68bbfe9cp+75},
     },
    {
     {0x1.568442b104000p+122, 0x1.195bce21a4a4cp+82},
     {0x1.23c286bfdb000p+125, 0x1.d368142940c54p+83},
     {0x1.fffd89ae2f000p+124, 0x1.6961d82d67d29p+81},
     {0x1.9f80a3facf000p+124, 0x1.19861fe645aefp+78},
     {0x1.6b10abb2b0000p+124, 0x1.4aa2ee5f58c0cp+82},
     {0x1.b90322c900000p+119, 0x1.42e6d8398ebf0p+83},
     {0x1.ff687bb9b9000p+124, 0x1.09590940ec246p+83},
     {0x1.494950989a000p+125, 0x1.4ed69939f54a5p+83},
     {0x1.5c176f9414000p+122, 0x1.bccd986816af6p+83},
     {0x1.6dca3fa2e7000p+124, 0x1.38fff5b887f40p+83},
     {0x1.951e4290e0000p+122, 0x1.8c2bed86953acp+81},
     {0x1.a671255128000p+123, 0x1.2544a485cce86p+81},
     {0x1.b2745cf9ae000p+124, 0x1.c88c3ec2fb90fp+83},
     {0x1.2c6cfd90da000p+123, 0x1.a84f9682c93f7p+83},
     {0x1.a57e7d4e8e000p+124, 0x1.0beb05e6abcdfp+81},
     {0x1.8f40d0ef24000p+124, 0x1.aeb1b1661a570p+81},
     },
    {
     {0x1.e01f9407c4000p+129, 0x1.1e87e3b708c22p+90},
     {0x1.e201959d63000p+131, 0x1.4160efcbeef78p+90},
     {0x1.31982160c4000p+132, 0x1.4480003b19f81p+89},
     {0x1.7f0fe22eef000p+132, 0x1.0b25d1ed6a121p+87},
     {0x1.00d5bf9f80000p+126, 0x1.747bf6c0d8b31p+90},
     {0x1.8ad801f1a0000p+129, 0x1.bcbc193dd346cp+88},
     {0x1.2a9c662802000p+130, 0x1.fddf745f1ee5ap+88},
     {0x1.d836977997000p+131, 0x1.d1f525311dabfp+90},
     {0x1.85903a5f3c000p+132, 0x1.1660b883eb1a4p+90},
     {0x1.a3320451ba800p+132, 0x1.41dedd270b797p+88},
     {0x1.ce462d2242000p+132, 0x1.79125e4f2418ap+90},
     {0x1.d67cf11ca9800p+132, 0x1.0272e6220fc37p+90},
     {0x1.add7c7ba40000p+132, 0x1.d2e0f92de9773p+87},
     {0x1.57b0afae95000p+131, 0x1.ea083b704edc0p+90},
     {0x1.e30840c0e8000p+128, 0x1.eba48f8e2a378p+90},
     {0x1.5aabc9d4bf800p+132, 0x1.a977e531befa8p+90},
     {0x1.82a0ee308b800p+132, 0x1.2ed72602864a3p+88},
     },
    {
     {0x1.06cf388320000p+134, 0x1.928222f7c81d9p+98},
     {0x1.a1bf2dfdc0000p+136, 0x1.7691a1e475ec2p+97},
     {0x1.bb35a9d83c000p+137, 0x1.50297632195fep+97},
     {0x1.b0c7cfa209000p+139, 0x1.08e60e4fc6baep+96},
     {0x1.4921eae073800p+140, 0x1.3586f9a06cbf4p+98},
     {0x1.172ab95fd6000p+139, 0x1.7a70fdd8b6610p+98},
     {0x1.68acfd38e8000p+139, 0x1.5e172302320e2p+98},
     {0x1.f34ce4f4e8000p+138, 0x1.1c4557a753b6cp+98},
     {0x1.01123dfc72000p+140, 0x1.4622df275a365p+95},
     {0x1.9db3f73893000p+139, 0x1.cf7c96f698830p+95},
     {0x1.f6d5907a7e000p+138, 0x1.ebd8e9c0e37a5p+97},
     {0x1.e7abc6d98b000p+139, 0x1.f1a79e989b457p+97},
     {0x1.8e92d65018000p+136, 0x1.1c8ffe978c39ep+98},
     {0x1.1d42b11e83800p+140, 0x1.c0c39f95f19abp+92},
     {0x1.0579b3ad70800p+140, 0x1.f55e10b41e4a2p+97},
     {0x1.0cb5cec87c000p+138, 0x1.a546c43a54205p+98},
     {0x1.2009162ca2800p+140, 0x1.c22d132ce1471p+97},
     {0x1.3d803cbad1800p+140, 0x1.6f3d636bc541bp+95},
     },
    {
     {0x1.b09acf4b80000p+146, 0x1.e0958b3fc5a41p+105},
     {0x1.6f0acc1cea000p+147, 0x1.586ae0321d89fp+104},
     {0x1.d8992594f0000p+145, 0x1.fd46b49aa2b42p+106},
     {0x1.4be496434a000p+146, 0x1.846cf4df36408p+106},
     {0x1.8cc9189a96000p+147, 0x1.479c156f25f6cp+106},
     {0x1.c776b470b0000p+143, 0x1.020640df24636p+104},
     {0x1.cd534fe2dc000p+147, 0x1.ddd5bd56c8a86p+106},
     {0x1.82fa017336000p+147, 0x1.9e0161aac9805p+105},
     {0x1.946f7304e8000p+147, 0x1.a622eb926525ap+105},
     {0x1.551407a0b7000p+147, 0x1.78923f9483982p+106},
     {0x1.034c6790f6000p+146, 0x1.cb97659eca409p+106},
     {0x1.452e68b9f4000p+145, 0x1.e91af6550ad66p+105},
     {0x1.407e5f3ab7000p+147, 0x1.58efe3be140c5p+106},
     {0x1.a514c77360000p+147, 0x1.ef9fd35ae1d32p+103},
     {0x1.840b4e6816000p+147, 0x1.a8df1a86177b9p+106},
     {0x1.7a503c2406000p+147, 0x1.a20f1e945f461p+105},
     {0x1.9fa0adbac0000p+147, 0x1.03cca141443f9p+106},
     {0x1.0c070c3e0c000p+147, 0x1.71f2fcf80933dp+106},
     {0x1.952a21ca4c000p+145, 0x1.7334b3dff2d8bp+105},
     },
    {
     {0x1.b7d0145780000p+153, 0x1.4ca65aa6e2e69p+113},
     {0x1.22e534dde0000p+150, 0x1.3e4218a70dc2fp+114},
     {0x1.157cefeb34000p+153, 0x1.3349341ba5ba9p+114},
     {0x1.3ca3f6e300000p+151, 0x1.a8a4cf94c4963p+113},
     {0x1.016a241f28000p+152, 0x1.a9ef27a2e8284p+113},
     {0x1.e66c961dd0000p+154, 0x1.f9453ff49eeb9p+114},
     {0x1.1945b982ed000p+155, 0x1.54038a5103c5fp+114},
     {0x1.3e5ca23c80000p+152, 0x1.7e6af65bdb8e7p+114},
     {0x1.1ce0e51379000p+155, 0x1.feec500f6cd99p+110},
     {0x1.98788b5ce0000p+154, 0x1.ccccb6fa6a5aep+113},
     {0x1.b19e1bb310000p+154, 0x1.e8165d05819c5p+114},
     {0x1.df2e1fa0ce000p+154, 0x1.481a67408850bp+111},
     {0x1.8b801f14e4000p+153, 0x1.d767562c372cdp+112},
     {0x1.38d9254eb8000p+153, 0x1.c6ebbea0ef5f4p+113},
     {0x1.354ce8cdbc000p+154, 0x1.d062d71a7af94p+112},
     {0x1.94eef4587e000p+154, 0x1.4a1e8a895454cp+111},
     {0x1.3b8c91b979000p+155, 0x1.77cf77e873cd7p+114},
     {0x1.0b1e374058000p+155, 0x1.d1d5597316f21p+111},
     {0x1.088d7305d7000p+155, 0x1.ed07530f9a7fap+114},
     {0x1.67ddaa2cae000p+154, 0x1.3929cf709cf74p+114},
     },
};

} // namespace INT8

namespace FP8 {

// qPi_1[i] = double(q[i]*P[i]/p[i]), where q[i]*P[i]/p[i] == 1 mod p[i]
inline constexpr double qPi_1[19][20] = {
    {0x1.0c00000000000p+16, 0x1.ff00200000000p+19},
    {0x1.7764000000000p+24, 0x1.ff00200000000p+19, 0x1.f2c5400000000p+29},
    {0x1.1df3c0b000000p+39, 0x1.37e3f37802000p+39, 0x1.5353068800000p+39, 0x1.41ded42800000p+39},
    {0x1.6d0eb70b20000p+47, 0x1.09283a3ac0020p+47, 0x1.3db96496d0000p+47, 0x1.4c99f49408000p+48, 0x1.b412a4ced0000p+47},
    {0x1.f1b859e448000p+55, 0x1.cb87f75e19160p+57, 0x1.a78fa4a778120p+57, 0x1.8be2432b37ea0p+57, 0x1.95aa7f76ea0e0p+57,
     0x1.342ec14351280p+57},
    {0x1.6e36bcb20962fp+66, 0x1.29319af99d493p+64, 0x1.2a90f37866de3p+66, 0x1.3831f5e924181p+63, 0x1.456914d693df9p+66,
     0x1.9253098f3989dp+65, 0x1.763dc6dd30f04p+65},
    {0x1.a15e88b747187p+75, 0x1.eaa812bda2520p+72, 0x1.efb68aeb6605fp+73, 0x1.b9fb63b3fdfc8p+75, 0x1.395e6f0b3aa0fp+75,
     0x1.2235b4bacbdf0p+74, 0x1.2227aa1972dabp+74, 0x1.a9651f14fa3b9p+74},
    {0x1.96f8419fe91d8p+84, 0x1.36de721b67fd0p+84, 0x1.14b49f0f6214fp+83, 0x1.fcbcfdb0310e7p+83, 0x1.27ea22cad6b70p+84,
     0x1.3d6495bef8ca9p+84, 0x1.93fed398cb13ap+83, 0x1.721898e2df1c3p+80, 0x1.6b7304668b334p+83},
    {0x1.0a43e957ee2a7p+92, 0x1.578528613920cp+93, 0x1.46cea2174d975p+92, 0x1.e80edec944bc5p+91, 0x1.d4ea759defb47p+92,
     0x1.02a30fa105218p+93, 0x1.4000b3bc837cdp+93, 0x1.42267c34739bap+93, 0x1.31d056d060a0ep+93, 0x1.02045d5625a52p+92},
    {0x1.1973870017be7p+99, 0x1.1501605fb5f96p+102, 0x1.42f65e88c9a02p+102, 0x1.5bbad64e3f5ddp+97, 0x1.92ab20ee6aa36p+100,
     0x1.c67af41819970p+100, 0x1.4f6a30c1ac237p+102, 0x1.7d93c55a8a583p+99, 0x1.034f25d7069e0p+101, 0x1.d5a89d9cf0724p+99,
     0x1.45ae0c4226520p+102},
    {0x1.a0290e2b71c98p+110, 0x1.e2230164724eep+108, 0x1.513a13457aaf5p+103, 0x1.b788b15d80761p+108, 0x1.e4164893735b5p+109,
     0x1.8d3e9d8779688p+110, 0x1.1d8e28b6234b0p+111, 0x1.e6481a8db13abp+106, 0x1.8f7abd99cd50dp+107, 0x1.33ab698bc015cp+111,
     0x1.73d8a8c64e23ap+106, 0x1.2523664dcf7b6p+111},
    {0x1.688f9d718a858p+119, 0x1.7b960ab4db1a4p+115, 0x1.2c69d6086a3f8p+119, 0x1.57f34babcc6a7p+112, 0x1.3e3040e8af344p+119,
     0x1.0853af44712efp+120, 0x1.9d725346d3b08p+115, 0x1.4c9c2e9e6639ap+119, 0x1.5aee0002d0330p+116, 0x1.d08dd023ed4f1p+116,
     0x1.d358df22b10cbp+119, 0x1.24e64cf08b1adp+119, 0x1.19dd53858da83p+117},
    {0x1.21cebe19f47b8p+128, 0x1.c6ce4ff33f40cp+126, 0x1.b008ef3a12153p+128, 0x1.b2c4dce856c19p+127, 0x1.eddcda21c99b3p+128,
     0x1.5c67e75531cabp+126, 0x1.53e995e722da0p+128, 0x1.1f8bcd40bf381p+128, 0x1.84c0a2bc0a936p+126, 0x1.44a3fa21fb483p+127,
     0x1.2757bd6df5c9dp+127, 0x1.0b9b421fbf0f1p+128, 0x1.b1a7133379466p+125, 0x1.ddcb0014c36abp+128},
    {0x1.b685162f1d0cap+136, 0x1.130a7c7b7aafcp+137, 0x1.c0c08dfdeccf0p+136, 0x1.0991164b06b8cp+137, 0x1.660851eb61799p+136,
     0x1.0b6af1c6c5107p+136, 0x1.8e608c68a6860p+137, 0x1.fe1e6800579acp+136, 0x1.c72f801bee437p+136, 0x1.7ddd1db329320p+134,
     0x1.e8dc45fc0266dp+135, 0x1.d4a9fa996de1dp+134, 0x1.dc84a7ac8e89ap+135, 0x1.2054f1e4eaffcp+137, 0x1.75857d937bd89p+136},
    {0x1.e0a0a525b0f76p+144, 0x1.6918f969b4d0fp+145, 0x1.00a409cb1f5a8p+146, 0x1.6f1f1a635f12cp+144, 0x1.cd10c0aa90da0p+142,
     0x1.310663807e820p+146, 0x1.4d2bec1ed766cp+146, 0x1.3c2579d2c4374p+146, 0x1.f5928ab8fe72ap+145, 0x1.a124be7d4c25ep+145,
     0x1.66d3783b8871dp+145, 0x1.7eb516d3f4ae6p+143, 0x1.03c8c5c2f53b2p+146, 0x1.0d531b58798a7p+146, 0x1.d7457dfe5c1f2p+145,
     0x1.c0c5e85623c01p+145},
    {0x1.603c00276b0f9p+154, 0x1.7ad7c4a6de835p+150, 0x1.4f1add3a9186dp+153, 0x1.1183af3827337p+153, 0x1.7df782e0e563ap+147,
     0x1.3fa8d875462e6p+151, 0x1.408a04ba38f58p+151, 0x1.2d1c78660b3bdp+153, 0x1.0e6a6653eca1cp+155, 0x1.163add4f9551fp+155,
     0x1.bdca8fedff2fap+154, 0x1.08fe788dab267p+155, 0x1.0aa21d806ea4dp+154, 0x1.413ae15e9ba9ep+154, 0x1.c97f40d7bd205p+154,
     0x1.5f34abf571245p+154, 0x1.270d47fa40c94p+155},
    {0x1.18a6adc9ef8d8p+162, 0x1.8d1a49c5a70a7p+163, 0x1.eb351932ff0d1p+162, 0x1.a4a944f31ff75p+162, 0x1.c27d23b416f4ep+163,
     0x1.1500b60b276b3p+162, 0x1.395980e12474ap+158, 0x1.2ea3c8fbe3b57p+163, 0x1.88c831a6ebb5fp+161, 0x1.7d043a617dae0p+162,
     0x1.3b5bdb662e126p+163, 0x1.2bb2d9cea25bcp+160, 0x1.89ef3acfbafadp+163, 0x1.5901a29162f1ap+160, 0x1.94056cc824222p+163,
     0x1.e661c466fdbcdp+163, 0x1.a970310ce41c9p+158, 0x1.35a1946ad852cp+163},
    {0x1.a66b557fdedeap+172, 0x1.15c6c48eccf7dp+172, 0x1.373615d3c88a0p+172, 0x1.59e52f31e9c8dp+172, 0x1.34747eb123774p+168,
     0x1.bc26b71788d47p+171, 0x1.01a818792079fp+170, 0x1.1f7b00e168d2ap+168, 0x1.64eb78f2bd15ap+172, 0x1.983510484ad8ep+172,
     0x1.0e76370a4046bp+172, 0x1.162b58f4c3894p+171, 0x1.5a848cf70d558p+170, 0x1.121a560f5d9eap+171, 0x1.f91e687b4c553p+169,
     0x1.08a9b6e1034a8p+172, 0x1.898d185fbcc3bp+172, 0x1.7b670c6909ee7p+171, 0x1.38c69eecbaeb9p+172},
    {0x1.361e857d43564p+181, 0x1.a700af7b60a88p+180, 0x1.1a2a9945ae082p+181, 0x1.707d3539782afp+178, 0x1.b0bb6e4ecb2f0p+175,
     0x1.34e3303da0949p+181, 0x1.c2ff379d76633p+179, 0x1.d50de6fbc76b7p+180, 0x1.a8fe26b516018p+180, 0x1.60e993e7bef05p+180,
     0x1.17e6b6445035cp+178, 0x1.39330a89a09a6p+180, 0x1.24d2b7c093111p+180, 0x1.39ba0f3335964p+179, 0x1.14a39891579dep+181,
     0x1.4a7fd7c93554fp+180, 0x1.c912e8398dca7p+180, 0x1.6e0c6e958297bp+178, 0x1.c05a852583bc1p+180, 0x1.0711728044a1ep+181},
};

// idx = num_moduli - threshold<gemmul8::Backend::FP8>::P_is_double - 1
// qPi_2[idx][i][1] = first (53-ceil(log2(rho))) bits of q[i]*P[i]/p[i] for rho = sum(floor(p[:]/2)),
// qPi_2[idx][i][2] = double(q[i]*P[i]/p[i] - qPi_2[idx][i][1])
inline constexpr double2 qPi_2[15][20] = {
    {
     {0x1.f1b859e448000p+55, 0x0.0000000000000p+0},
     {0x1.cb87f75e19000p+57, 0x1.6008000000000p+13},
     {0x1.a78fa4a778000p+57, 0x1.2000000000000p+13},
     {0x1.8be2432b37000p+57, 0x1.d400000000000p+16},
     {0x1.95aa7f76ea000p+57, 0x1.c000000000000p+12},
     {0x1.342ec14351000p+57, 0x1.4000000000000p+14},
     },
    {
     {0x1.6e36bcb209000p+66, 0x1.8bc8000000000p+24},
     {0x1.29319af99c000p+64, 0x1.4928010000000p+24},
     {0x1.2a90f37866000p+66, 0x1.bc52000000000p+25},
     {0x1.3831f5e920000p+63, 0x1.0604000000000p+25},
     {0x1.456914d693000p+66, 0x1.bf22000000000p+25},
     {0x1.9253098f38000p+65, 0x1.89cc000000000p+25},
     {0x1.763dc6dd30000p+65, 0x1.e074000000000p+24},
     },
    {
     {0x1.a15e88b747000p+75, 0x1.86e9380000000p+31},
     {0x1.eaa812bda0000p+72, 0x1.2901c20080000p+33},
     {0x1.efb68aeb64000p+73, 0x1.02f7000000000p+34},
     {0x1.b9fb63b3fd000p+75, 0x1.f90c610000000p+34},
     {0x1.395e6f0b3a000p+75, 0x1.41ea230000000p+34},
     {0x1.2235b4baca000p+74, 0x1.defe430000000p+34},
     {0x1.2227aa1972000p+74, 0x1.b562f40000000p+33},
     {0x1.a9651f14fa000p+74, 0x1.dc4c980000000p+31},
     },
    {
     {0x1.96f8419fe9000p+84, 0x1.d7b370c000000p+40},
     {0x1.36de721b67000p+84, 0x1.fa0a702d80200p+43},
     {0x1.14b49f0f62000p+83, 0x1.4f357d4000000p+39},
     {0x1.fcbcfdb030000p+83, 0x1.0e779fa900000p+43},
     {0x1.27ea22cad6000p+84, 0x1.6df10a7b00000p+43},
     {0x1.3d6495bef8000p+84, 0x1.9529231f00000p+43},
     {0x1.93fed398ca000p+83, 0x1.139dbf9600000p+43},
     {0x1.721898e2d0000p+80, 0x1.e38663f900000p+43},
     {0x1.6b7304668a000p+83, 0x1.3344f97c00000p+43},
     },
    {
     {0x1.0a43e957ee000p+92, 0x1.537844389a000p+49},
     {0x1.5785286139000p+93, 0x1.062c6b562f004p+50},
     {0x1.46cea2174c000p+92, 0x1.974b955a3d000p+52},
     {0x1.e80edec944000p+91, 0x1.78aaaf7088000p+50},
     {0x1.d4ea759dee000p+92, 0x1.b46a84c700400p+52},
     {0x1.02a30fa105000p+93, 0x1.0c0a25f860000p+50},
     {0x1.4000b3bc83000p+93, 0x1.f341c111d6000p+51},
     {0x1.42267c3473000p+93, 0x1.373102fe6ec00p+52},
     {0x1.31d056d060000p+93, 0x1.41cd040861400p+52},
     {0x1.02045d5624000p+92, 0x1.a51bfa7644400p+52},
     },
    {
     {0x1.1973870010000p+99, 0x1.ef9a2a1eec316p+61},
     {0x1.1501605fb5000p+102, 0x1.f2cd7cf110a76p+61},
     {0x1.42f65e88c9000p+102, 0x1.404a0e331142ap+61},
     {0x1.5bbad64e20000p+97, 0x1.f5dd1c639bd4cp+61},
     {0x1.92ab20ee68000p+100, 0x1.51b0c0e7ee6cap+61},
     {0x1.c67af41818000p+100, 0x1.9704741a5e834p+60},
     {0x1.4f6a30c1ac000p+102, 0x1.1b9aba7c982e8p+59},
     {0x1.7d93c55a88000p+99, 0x1.2c1801c54aa24p+60},
     {0x1.034f25d706000p+101, 0x1.3c0d53bdf217cp+60},
     {0x1.d5a89d9cf0000p+99, 0x1.c8fb4c876d6c0p+57},
     {0x1.45ae0c4226000p+102, 0x1.47f28bb407014p+60},
     },
    {
     {0x1.a0290e2b70000p+110, 0x1.c9795d57e1ea4p+70},
     {0x1.e223016470000p+108, 0x1.27726ee7c0191p+69},
     {0x1.513a134500000p+103, 0x1.eabd57ae52df5p+69},
     {0x1.b788b15d80000p+108, 0x1.d850f932a00ddp+66},
     {0x1.e416489370000p+109, 0x1.ada5876b62cdbp+70},
     {0x1.8d3e9d8778000p+110, 0x1.687bf1fefa695p+70},
     {0x1.1d8e28b623000p+111, 0x1.2bfffd422e36ap+69},
     {0x1.e6481a8da0000p+106, 0x1.13aa9be1cf2d8p+70},
     {0x1.8f7abd99c0000p+107, 0x1.aa19f8972bf9ep+70},
     {0x1.33ab698bc0000p+111, 0x1.5bd19fb1c45f8p+67},
     {0x1.73d8a8c640000p+106, 0x1.c474d347d6e4bp+69},
     {0x1.2523664dcf000p+111, 0x1.ed95ba150cb98p+69},
     },
    {
     {0x1.688f9d7188000p+119, 0x1.42bde83b654e2p+80},
     {0x1.7b960ab4c0000p+115, 0x1.b1a3c7e821fd0p+79},
     {0x1.2c69d60868000p+119, 0x1.1fbfb856ed7b7p+80},
     {0x1.57f34baa00000p+112, 0x1.cc6a6da63b1c0p+80},
     {0x1.3e3040e8ac000p+119, 0x1.9a2320be2ec47p+80},
     {0x1.0853af4470000p+120, 0x1.2ef4f2ff7cbbep+80},
     {0x1.9d725346c0000p+115, 0x1.3b07f7ae0eea0p+79},
     {0x1.4c9c2e9e64000p+119, 0x1.1cceec515248bp+80},
     {0x1.5aee0002c0000p+116, 0x1.032fdf370c5dap+80},
     {0x1.d08dd023e0000p+116, 0x1.a9e1f6c57dc26p+79},
     {0x1.d358df22b0000p+119, 0x1.0cb1666846559p+79},
     {0x1.24e64cf088000p+119, 0x1.8d65680b0c356p+80},
     {0x1.19dd538580000p+117, 0x1.b506e0f29c137p+80},
     },
    {
     {0x1.21cebe19f4000p+128, 0x1.ede6fc03adc00p+86},
     {0x1.c6ce4ff338000p+126, 0x1.d030e2eea2ccep+88},
     {0x1.b008ef3a12000p+128, 0x1.52d7041d2e600p+84},
     {0x1.b2c4dce854000p+127, 0x1.60c6df870d490p+88},
     {0x1.eddcda21c8000p+128, 0x1.9b2eb2ad8cba3p+88},
     {0x1.5c67e75530000p+126, 0x1.caab40955c49fp+86},
     {0x1.53e995e722000p+128, 0x1.b4069d1c3a5edp+87},
     {0x1.1f8bcd40be000p+128, 0x1.380efee975457p+88},
     {0x1.84c0a2bc08000p+126, 0x1.49ae4136962f9p+87},
     {0x1.44a3fa21f8000p+127, 0x1.a416fa641b808p+88},
     {0x1.2757bd6df4000p+127, 0x1.c9cd54b630802p+87},
     {0x1.0b9b421fbe000p+128, 0x1.0f143b6d422b9p+88},
     {0x1.b1a7133370000p+125, 0x1.28cb0a03dadcap+88},
     {0x1.ddcb0014c2000p+128, 0x1.6ab2d7cbf9e6bp+88},
     },
    {
     {0x1.b685162f1c000p+136, 0x1.0ca36113a0ef5p+96},
     {0x1.130a7c7b7a000p+137, 0x1.5f71581db32f9p+96},
     {0x1.c0c08dfdec000p+136, 0x1.9e0aa7505737ep+95},
     {0x1.0991164b06000p+137, 0x1.71717053b80fbp+96},
     {0x1.660851eb60000p+136, 0x1.79976c2baeb1ep+96},
     {0x1.0b6af1c6c4000p+136, 0x1.106fd8eae632fp+96},
     {0x1.8e608c68a6000p+137, 0x1.0bf3017cd7a0dp+96},
     {0x1.fe1e680054000p+136, 0x1.cd61b0a5ef495p+97},
     {0x1.c72f801bec000p+136, 0x1.21b4316d933f5p+97},
     {0x1.7ddd1db320000p+134, 0x1.2640d93d4eba6p+97},
     {0x1.e8dc45fc00000p+135, 0x1.3367e43d37b47p+96},
     {0x1.d4a9fa9960000p+134, 0x1.bc39037171403p+97},
     {0x1.dc84a7ac88000p+135, 0x1.a2669e3d73b39p+97},
     {0x1.2054f1e4ea000p+137, 0x1.ff7ae3fb1d20fp+96},
     {0x1.75857d9378000p+136, 0x1.ec4bf0ef762a7p+97},
     },
    {
     {0x1.e0a0a525b0000p+144, 0x1.eeb48df9600b2p+103},
     {0x1.6918f969b4000p+145, 0x1.a1d7c0d6a9f49p+104},
     {0x1.00a409cb1e000p+146, 0x1.5a8117f608e57p+106},
     {0x1.6f1f1a6358000p+144, 0x1.c4ae49dbbfb60p+106},
     {0x1.cd10c0aa80000p+142, 0x1.0d9fe142473cep+106},
     {0x1.310663807e000p+146, 0x1.04024df853c4ap+105},
     {0x1.4d2bec1ed6000p+146, 0x1.66bf9622f9daap+106},
     {0x1.3c2579d2c4000p+146, 0x1.b9e62b3ce7a14p+103},
     {0x1.f5928ab8fc000p+145, 0x1.3950b56fdfcfbp+106},
     {0x1.a124be7d4c000p+145, 0x1.2f2afb37e08e2p+102},
     {0x1.66d3783b88000p+145, 0x1.c74806710c488p+103},
     {0x1.7eb516d3f0000p+143, 0x1.2b97c7b7f0635p+105},
     {0x1.03c8c5c2f4000p+146, 0x1.3b263e1a80e94p+106},
     {0x1.0d531b5878000p+146, 0x1.8a6a18f89a200p+106},
     {0x1.d7457dfe5c000p+145, 0x1.f26ad22e386b3p+101},
     {0x1.c0c5e85620000p+145, 0x1.e004270f9f348p+106},
     },
    {
     {0x1.603c002768000p+154, 0x1.87c73fa7e762ep+115},
     {0x1.7ad7c4a6c0000p+150, 0x1.e8354dc21b847p+114},
     {0x1.4f1add3a90000p+153, 0x1.86cf3b9370ea7p+113},
     {0x1.1183af3820000p+153, 0x1.ccda77dc103efp+115},
     {0x1.7df782e000000p+147, 0x1.cac73adb7f028p+114},
     {0x1.3fa8d87540000p+151, 0x1.8b975ab260137p+113},
     {0x1.408a04ba20000p+151, 0x1.8f5835e42b6cdp+115},
     {0x1.2d1c786608000p+153, 0x1.9deab799bda5dp+114},
     {0x1.0e6a6653ec000p+155, 0x1.4387207a09623p+114},
     {0x1.163add4f94000p+155, 0x1.51eecb6a4d832p+115},
     {0x1.bdca8fedfc000p+154, 0x1.97ccdf41c9c06p+115},
     {0x1.08fe788daa000p+155, 0x1.266cfe9e5075ap+115},
     {0x1.0aa21d806c000p+154, 0x1.5265d2d48614ap+115},
     {0x1.413ae15e98000p+154, 0x1.d4f0b50b36141p+115},
     {0x1.c97f40d7bc000p+154, 0x1.2056d4c249be0p+114},
     {0x1.5f34abf570000p+154, 0x1.2453df6800db5p+114},
     {0x1.270d47fa40000p+155, 0x1.927de7dfc04f3p+114},
     },
    {
     {0x1.18a6adc9ec000p+162, 0x1.c6c2223baca5dp+123},
     {0x1.8d1a49c5a6000p+163, 0x1.0a71f7854c09dp+123},
     {0x1.eb351932fc000p+162, 0x1.8689dfde294c2p+123},
     {0x1.a4a944f31c000p+162, 0x1.fba5c0d60b450p+123},
     {0x1.c27d23b416000p+163, 0x1.e9b70dd34156dp+122},
     {0x1.1500b60b24000p+162, 0x1.b59a850a61284p+123},
     {0x1.395980e100000p+158, 0x1.23a53678659cep+123},
     {0x1.2ea3c8fbe2000p+163, 0x1.b5779b31c9c91p+123},
     {0x1.88c831a6e8000p+161, 0x1.daf8e45296e5dp+122},
     {0x1.7d043a617c000p+162, 0x1.adf8b9b88ce27p+122},
     {0x1.3b5bdb662e000p+163, 0x1.25bc8886081c0p+119},
     {0x1.2bb2d9cea0000p+160, 0x1.2de1b9e00ee72p+121},
     {0x1.89ef3acfba000p+163, 0x1.f5901eb09def4p+122},
     {0x1.5901a29160000p+160, 0x1.78d2746006193p+121},
     {0x1.94056cc824000p+163, 0x1.10dd9992e83ebp+120},
     {0x1.e661c466fc000p+163, 0x1.bcd792ecfc278p+123},
     {0x1.a970310cc0000p+158, 0x1.20e494c711adfp+123},
     {0x1.35a1946ad8000p+163, 0x1.4ae4321f77fd2p+121},
     },
    {
     {0x1.a66b557fde000p+172, 0x1.bd4371d5db091p+131},
     {0x1.15c6c48ecc000p+172, 0x1.ef99a07433e8cp+131},
     {0x1.373615d3c8000p+172, 0x1.13fd974aab7b9p+131},
     {0x1.59e52f31e8000p+172, 0x1.c8d7962000887p+132},
     {0x1.34747eb120000p+168, 0x1.bb9ec8dd5de04p+129},
     {0x1.bc26b71788000p+171, 0x1.a8d889e0deef6p+130},
     {0x1.01a8187920000p+170, 0x1.e7aca4ffc6bc1p+128},
     {0x1.1f7b00e160000p+168, 0x1.1a53b42560a37p+131},
     {0x1.64eb78f2bc000p+172, 0x1.15a1d14b049c2p+132},
     {0x1.983510484a000p+172, 0x1.b1b33df909c27p+131},
     {0x1.0e76370a40000p+172, 0x1.1aa2a672fc1a0p+130},
     {0x1.162b58f4c0000p+171, 0x1.c49eae5ff7688p+132},
     {0x1.5a848cf708000p+170, 0x1.555eae4e01c8cp+132},
     {0x1.121a560f5c000p+171, 0x1.9ea6e4a210d34p+131},
     {0x1.f91e687b40000p+169, 0x1.8aa56b681a69bp+132},
     {0x1.08a9b6e102000p+172, 0x1.4a80190f345b2p+132},
     {0x1.898d185fbc000p+172, 0x1.875cb81982ad4p+131},
     {0x1.7b670c6908000p+171, 0x1.ee761eb00d3cbp+131},
     {0x1.38c69eecba000p+172, 0x1.d72041b74b4dcp+131},
     },
    {
     {0x1.361e857d42000p+181, 0x1.563b5b2ef7d58p+141},
     {0x1.a700af7b60000p+180, 0x1.510073018f89dp+139},
     {0x1.1a2a9945ae000p+181, 0x1.037eaa1ddd2e5p+136},
     {0x1.707d353970000p+178, 0x1.055e30172b7abp+141},
     {0x1.b0bb6e4e80000p+175, 0x1.2cbbfcee59047p+141},
     {0x1.34e3303da0000p+181, 0x1.291a3a070a299p+140},
     {0x1.c2ff379d70000p+179, 0x1.98cb9fc3334a7p+141},
     {0x1.d50de6fbc4000p+180, 0x1.b5bbb5cb9750bp+141},
     {0x1.a8fe26b514000p+180, 0x1.00c30fa125bdap+141},
     {0x1.60e993e7bc000p+180, 0x1.7827435e5a636p+141},
     {0x1.17e6b64450000p+178, 0x1.ae2aef7feb8f6p+135},
     {0x1.39330a89a0000p+180, 0x1.34c195038fdadp+139},
     {0x1.24d2b7c090000p+180, 0x1.888860be2a360p+141},
     {0x1.39ba0f3330000p+179, 0x1.658eb4ae0281bp+141},
     {0x1.14a3989156000p+181, 0x1.9de4165f77e8ep+141},
     {0x1.4a7fd7c934000p+180, 0x1.54f75370f31edp+140},
     {0x1.c912e8398c000p+180, 0x1.ca6d6a1592553p+140},
     {0x1.6e0c6e9580000p+178, 0x1.4bd5874f1df6fp+139},
     {0x1.c05a852580000p+180, 0x1.de067c0c119fap+141},
     {0x1.0711728044000p+181, 0x1.43bc7b4c5afcfp+140},
     },
};
} // namespace FP8

namespace INT8 {
template <int num_moduli, int IDX> inline constexpr double qPi_double_v   = INT8::qPi_1[num_moduli - 2][IDX];
template <int num_moduli, int IDX> inline constexpr double2 qPi_double2_v = INT8::qPi_2[num_moduli - threshold<gemmul8::Backend::INT8>::P_is_double - 1][IDX];
} // namespace INT8

namespace FP8 {
template <int num_moduli, int IDX> inline constexpr double qPi_double_v   = FP8::qPi_1[num_moduli - 2][IDX];
template <int num_moduli, int IDX> inline constexpr double2 qPi_double2_v = FP8::qPi_2[num_moduli - threshold<gemmul8::Backend::FP8>::P_is_double - 1][IDX];
} // namespace FP8

//==========
// Set constant memory
//==========
template <gemmul8::Backend backend> __forceinline__ void upload_constants(cudaStream_t stream);

template <> __forceinline__ void upload_constants<gemmul8::Backend::INT8>(cudaStream_t stream) {
    cudaMemcpyToSymbolAsync(INT8::mod_pow2, INT8::mod_pow2_h, sizeof(INT8::mod_pow2_h), 0, cudaMemcpyHostToDevice, stream);
}

template <> __forceinline__ void upload_constants<gemmul8::Backend::FP8>(cudaStream_t stream) {
    cudaMemcpyToSymbolAsync(FP8::mod_pow2, FP8::mod_pow2_h, sizeof(FP8::mod_pow2_h), 0, cudaMemcpyHostToDevice, stream);
}

} // namespace table
