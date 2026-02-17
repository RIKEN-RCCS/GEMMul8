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
template <> inline constexpr int32_t moduli<gemmul8::Backend::FP8, 6>  = 511;  // Karatsuba
template <> inline constexpr int32_t moduli<gemmul8::Backend::FP8, 7>  = 509;  // Karatsuba
template <> inline constexpr int32_t moduli<gemmul8::Backend::FP8, 8>  = 503;  // Karatsuba
template <> inline constexpr int32_t moduli<gemmul8::Backend::FP8, 9>  = 499;  // Karatsuba
template <> inline constexpr int32_t moduli<gemmul8::Backend::FP8, 10> = 491;  // Karatsuba
template <> inline constexpr int32_t moduli<gemmul8::Backend::FP8, 11> = 487;  // Karatsuba
template <> inline constexpr int32_t moduli<gemmul8::Backend::FP8, 12> = 481;  // Karatsuba
template <> inline constexpr int32_t moduli<gemmul8::Backend::FP8, 13> = 479;  // Karatsuba
template <> inline constexpr int32_t moduli<gemmul8::Backend::FP8, 14> = 467;  // Karatsuba
template <> inline constexpr int32_t moduli<gemmul8::Backend::FP8, 15> = 463;  // Karatsuba
template <> inline constexpr int32_t moduli<gemmul8::Backend::FP8, 16> = 461;  // Karatsuba
template <> inline constexpr int32_t moduli<gemmul8::Backend::FP8, 17> = 457;  // Karatsuba
template <> inline constexpr int32_t moduli<gemmul8::Backend::FP8, 18> = 449;  // Karatsuba
template <> inline constexpr int32_t moduli<gemmul8::Backend::FP8, 19> = 443;  // Karatsuba

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
    { -0x1.1040000000000p+20,    0x0.0000000000000p+0},
    { -0x1.ff00200000000p+29,    0x0.0000000000000p+0},
    { -0x1.a3adda4800000p+39,    0x0.0000000000000p+0},
    { -0x1.0026dc7a72000p+49,    0x0.0000000000000p+0},
    { -0x1.08a826cc82c90p+58,    0x0.0000000000000p+0},
    { -0x1.0823d2b91c87ap+67,   0x1.2000000000000p+13},
    { -0x1.06979cfd06dcdp+76,  -0x1.6c00000000000p+16},
    { -0x1.01f9f2ba943dfp+85,   0x1.91a6600000000p+29},
    { -0x1.f6da3421aef4bp+93,  -0x1.f0462cb800000p+39},
    {-0x1.e23a40fe4d47bp+102,   0x1.c852d07630000p+46},
    {-0x1.caae68d1e281bp+111,   0x1.ca0ac5486aa80p+55},
    {-0x1.aee8d9792d4adp+120,  -0x1.0a6c386a5df35p+66},
    {-0x1.9322d774dddf8p+129,  -0x1.7d00fb1e0b947p+73},
    {-0x1.6fb44785185f6p+138,   0x1.b23dcd7a0c381p+83},
    {-0x1.4c8386acdb8a4p+147,   0x1.e0aee34fde0cbp+92},
    {-0x1.2b646cc2a3abfp+156, -0x1.cf99442a4b48dp+102},
    {-0x1.0b3b2313bb170p+165,  0x1.d0cecaa0ff362p+109},
    {-0x1.d4b2b8859b235p+173, -0x1.e862a74dd0311p+118},
    {-0x1.9588a2a799bb1p+182, -0x1.0522b783a744ep+126},
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
    0x1.f038c97b34e2fp-68, 0x1.f32581bdd19d6p-77, 0x1.fc13db17bc6d5p-86, 0x1.04a832d64f39ap-94, 0x1.0fce2774976fbp-103,
    0x1.1dc2221f1c51bp-112, 0x1.302cd946e7533p-121, 0x1.4521822aa7184p-130, 0x1.6475de320d42ap-139, 0x1.8a2f679c884e2p-148,
    0x1.b5cb23a9f657ep-157, 0x1.ea7b6503e103cp-166, 0x1.17a6b5e36c568p-174, 0x1.4335687788312p-183};
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
template <> inline constexpr float log2P<gemmul8::Backend::FP8, 7>  = 0x1.082e3e0000000p+5F;
template <> inline constexpr float log2P<gemmul8::Backend::FP8, 8>  = 0x1.2c258e0000000p+5F;
template <> inline constexpr float log2P<gemmul8::Backend::FP8, 9>  = 0x1.500b5c0000000p+5F;
template <> inline constexpr float log2P<gemmul8::Backend::FP8, 10> = 0x1.73e55c0000000p+5F;
template <> inline constexpr float log2P<gemmul8::Backend::FP8, 11> = 0x1.97a77e0000000p+5F;
template <> inline constexpr float log2P<gemmul8::Backend::FP8, 12> = 0x1.bb5d8a0000000p+5F;
template <> inline constexpr float log2P<gemmul8::Backend::FP8, 13> = 0x1.df01440000000p+5F;
template <> inline constexpr float log2P<gemmul8::Backend::FP8, 14> = 0x1.014f6c0000000p+6F;
template <> inline constexpr float log2P<gemmul8::Backend::FP8, 15> = 0x1.130b780000000p+6F;
template <> inline constexpr float log2P<gemmul8::Backend::FP8, 16> = 0x1.24c1280000000p+6F;
template <> inline constexpr float log2P<gemmul8::Backend::FP8, 17> = 0x1.3673a80000000p+6F;
template <> inline constexpr float log2P<gemmul8::Backend::FP8, 18> = 0x1.481fb60000000p+6F;
template <> inline constexpr float log2P<gemmul8::Backend::FP8, 19> = 0x1.59beba0000000p+6F;
template <> inline constexpr float log2P<gemmul8::Backend::FP8, 20> = 0x1.6b53cc0000000p+6F;

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
    {-255,    1,    2,    4,    8,   16,   32,   64,  128, -255,    1,    2,    4,    8,   16,   32,   64,  128, -255,    1,    2,    4,    8,   16,   32,   64,  128, -255,    1,    2,    4,    8,   16,   32,   64,  128, -255,    1,    2,    4,    8,   16,   32,   64,  128, -255,    1,    2,    4,    8,   16,   32,   64,  128, -255,    1,    2,    4,    8,   16,   32,   64,  128, -255},
    {-253,    3,    6,   12,   24,   48,   96,  192, -125, -250,    9,   18,   36,   72,  144, -221,   67,  134, -241,   27,   54,  108,  216,  -77, -154,  201, -107, -214,   81,  162, -185,  139, -231,   47,   94,  188, -133,  243,  -23,  -46,  -92, -184,  141, -227,   55,  110,  220,  -69, -138,  233,  -43,  -86, -172,  165, -179,  151, -207,   95,  190, -129,  251,   -7,  -14,  -28},
    {-247,    9,   18,   36,   72,  144, -215,   73,  146, -211,   81,  162, -179,  145, -213,   77,  154, -195,  113,  226,  -51, -102, -204,   95,  190, -123, -246,   11,   22,   44,   88,  176, -151,  201, -101, -202,   99,  198, -107, -214,   75,  150, -203,   97,  194, -115, -230,   43,   86,  172, -159,  185, -133,  237,  -29,  -58, -116, -232,   39,   78,  156, -191,  121,  242},
    {-243,   13,   26,   52,  104,  208,  -83, -166,  167, -165,  169, -161,  177, -145,  209,  -81, -162,  175, -149,  201,  -97, -194,  111,  222,  -55, -110, -220,   59,  118,  236,  -27,  -54, -108, -216,   67,  134, -231,   37,   74,  148, -203,   93,  186, -127,  245,   -9,  -18,  -36,  -72, -144,  211,  -77, -154,  191, -117, -234,   31,   62,  124,  248,   -3,   -6,  -12,  -24},
    {-235,   21,   42,   84,  168, -155,  181, -129,  233,  -25,  -50, -100, -200,   91,  182, -127,  237,  -17,  -34,  -68, -136,  219,  -53, -106, -212,   67,  134, -223,   45,   90,  180, -131,  229,  -33,  -66, -132,  227,  -37,  -74, -148,  195, -101, -202,   87,  174, -143,  205,  -81, -162,  167, -157,  177, -137,  217,  -57, -114, -228,   35,   70,  140, -211,   69,  138, -215},
    {-231,   25,   50,  100,  200,  -87, -174,  139, -209,   69,  138, -211,   65,  130, -227,   33,   66,  132, -223,   41,   82,  164, -159,  169, -149,  189, -109, -218,   51,  102,  204,  -79, -158,  171, -145,  197,  -93, -186,  115,  230,  -27,  -54, -108, -216,   55,  110,  220,  -47,  -94, -188,  111,  222,  -43,  -86, -172,  143, -201,   85,  170, -147,  193, -101, -202,   83},
    {-225,   31,   62,  124, -233,   15,   30,   60,  120,  240,   -1,   -2,   -4,   -8,  -16,  -32,  -64, -128,  225,  -31,  -62, -124,  233,  -15,  -30,  -60, -120, -240,    1,    2,    4,    8,   16,   32,   64,  128, -225,   31,   62,  124, -233,   15,   30,   60,  120,  240,   -1,   -2,   -4,   -8,  -16,  -32,  -64, -128,  225,  -31,  -62, -124,  233,  -15,  -30,  -60, -120, -240},
    {-223,   33,   66,  132, -215,   49,   98,  196,  -87, -174,  131, -217,   45,   90,  180, -119, -238,    3,    6,   12,   24,   48,   96,  192,  -95, -190,   99,  198,  -83, -166,  147, -185,  109,  218,  -43,  -86, -172,  135, -209,   61,  122, -235,    9,   18,   36,   72,  144, -191,   97,  194,  -91, -182,  115,  230,  -19,  -38,  -76, -152,  175, -129,  221,  -37,  -74, -148},
    {-211,   45,   90,  180, -107, -214,   39,   78,  156, -155,  157, -153,  161, -145,  177, -113, -226,   15,   30,   60,  120, -227,   13,   26,   52,  104,  208,  -51, -102, -204,   59,  118, -231,    5,   10,   20,   40,   80,  160, -147,  173, -121,  225,  -17,  -34,  -68, -136,  195,  -77, -154,  159, -149,  169, -129,  209,  -49,  -98, -196,   75,  150, -167,  133, -201,   65},
    {-207,   49,   98,  196,  -71, -142,  179, -105, -210,   43,   86,  172, -119,  225,  -13,  -26,  -52, -104, -208,   47,   94,  188,  -87, -174,  115,  230,   -3,   -6,  -12,  -24,  -48,  -96, -192,   79,  158, -147,  169, -125,  213,  -37,  -74, -148,  167, -129,  205,  -53, -106, -212,   39,   78,  156, -151,  161, -141,  181, -101, -202,   59,  118, -227,    9,   18,   36,   72},
    {-205,   51,  102,  204,  -53, -106, -212,   37,   74,  148, -165,  131, -199,   63,  126, -209,   43,   86,  172, -117,  227,   -7,  -14,  -28,  -56, -112, -224,   13,   26,   52,  104,  208,  -45,  -90, -180,  101,  202,  -57, -114, -228,    5,   10,   20,   40,   80,  160, -141,  179, -103, -206,   49,   98,  196,  -69, -138,  185,  -91, -182,   97,  194,  -73, -146,  169, -123},
    {-201,   55,  110,  220,  -17,  -34,  -68, -136,  185,  -87, -174,  109,  218,  -21,  -42,  -84, -168,  121, -215,   27,   54,  108,  216,  -25,  -50, -100, -200,   57,  114,  228,   -1,   -2,   -4,   -8,  -16,  -32,  -64, -128,  201,  -55, -110, -220,   17,   34,   68,  136, -185,   87,  174, -109, -218,   21,   42,   84,  168, -121,  215,  -27,  -54, -108, -216,   25,   50,  100},
    {-193,   63,  126, -197,   55,  110,  220,   -9,  -18,  -36,  -72, -144,  161, -127,  195,  -59, -118,  213,  -23,  -46,  -92, -184,   81,  162, -125,  199,  -51, -102, -204,   41,   82,  164, -121,  207,  -35,  -70, -140,  169, -111, -222,    5,   10,   20,   40,   80,  160, -129,  191,  -67, -134,  181,  -87, -174,  101,  202,  -45,  -90, -180,   89,  178,  -93, -186,   77,  154},
    {-187,   69,  138, -167,  109,  218,   -7,  -14,  -28,  -56, -112,  219,   -5,  -10,  -20,  -40,  -80, -160,  123, -197,   49,   98,  196,  -51, -102, -204,   35,   70,  140, -163,  117, -209,   25,   50,  100,  200,  -43,  -86, -172,   99,  198,  -47,  -94, -188,   67,  134, -175,   93,  186,  -71, -142,  159, -125,  193,  -57, -114,  215,  -13,  -26,  -52, -104, -208,   27,   54},
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
    {0x1.27ea1263b6a59p+66, 0x1.4dc944c8eb8d6p+66, 0x1.6f3654a282c48p+64, 0x1.69d17d2f7c9dcp+63, 0x1.2ce880a500f64p+64,
     0x1.b3680e68c9e2fp+65, 0x1.0696d67ee9c37p+67},
    {0x1.1831e11b270eep+75, 0x1.f61ae6add09f3p+75, 0x1.039624ea25af6p+76, 0x1.882bc9e299511p+75, 0x1.37bfdfa674482p+75,
     0x1.4d939784353d5p+72, 0x1.082245bce254dp+75, 0x1.3ab2ae0e8afd9p+75},
    {0x1.fd51d91acadf6p+84, 0x1.047ee39966b09p+84, 0x1.3ffcc17414d86p+82, 0x1.12d8fb4c1ba76p+80, 0x1.55c45f77cb5aap+83,
     0x1.b2fff53c8ba11p+84, 0x1.a305b35621611p+84, 0x1.13b7249d293fbp+81, 0x1.12e6b858e32f3p+83},
    {0x1.82f39b11b12e1p+93, 0x1.0e93e98d1fe43p+93, 0x1.f6543fcc6761bp+89, 0x1.f3dcdd115998ep+93, 0x1.14c51aa52d188p+90,
     0x1.b910b15b01052p+92, 0x1.f7d61f3147988p+87, 0x1.2a5a3c1a053a3p+93, 0x1.b3df3cf7069a2p+92, 0x1.5aa7de2ab7334p+90},
    {0x1.e553c81c6185fp+99, 0x1.9b96b8790cf5bp+101, 0x1.1f073c432f548p+100, 0x1.db58c69c4f38ap+102, 0x1.1bef826d9796cp+98,
     0x1.2affd32056a19p+102, 0x1.8106ef4246b88p+101, 0x1.f8048e5dc8f9fp+101, 0x1.80709a9ef6dd9p+102, 0x1.b8ac418c8e71cp+99,
     0x1.8cc82d22940d2p+102},
    {0x1.5602a250cd907p+111, 0x1.5802ce9d69e14p+103, 0x1.5004299fe1deap+106, 0x1.8c8175b3e98c7p+111, 0x1.314c5b113a9a5p+109,
     0x1.ae1161cc0cff2p+111, 0x1.bb6bfd2112143p+110, 0x1.37cb8e2b16b3bp+110, 0x1.63a338555fb78p+110, 0x1.b9376a78d5c72p+109,
     0x1.6c544f16498d1p+111, 0x1.a5f2f8de839ecp+105},
    {0x1.a82ecb43cd9d2p+119, 0x1.011d6fc30dc66p+120, 0x1.3babe91ab3b1ap+116, 0x1.269dfbe14f0fep+120, 0x1.77104de0aa7e9p+116,
     0x1.0e7042ccc6c50p+119, 0x1.76692127a91b7p+118, 0x1.45164481a46ecp+115, 0x1.7fcaceed33072p+115, 0x1.19841d3607044p+120,
     0x1.8caebaa4bbaecp+120, 0x1.1b24d44950c02p+118, 0x1.5b982b6f0da64p+120},
    {0x1.74690f706d0d0p+129, 0x1.91f47d5346392p+129, 0x1.e84b3bd868a54p+127, 0x1.2d08aaa0eebc3p+129, 0x1.b931309bc06fep+127,
     0x1.db886ee7bca4fp+125, 0x1.a6dbe7821952ap+128, 0x1.e3213480fdd40p+128, 0x1.4f02ff63a99ecp+128, 0x1.0d06d58d9abbap+129,
     0x1.9dcf4cbcc8a2bp+127, 0x1.8adbafd009272p+129, 0x1.e114d2475679fp+128, 0x1.033812d2e53f0p+129},
    {0x1.570e37e97c0ebp+136, 0x1.57494ec541c10p+136, 0x1.25db7628d568bp+136, 0x1.d6739778ab34fp+137, 0x1.01afe6c1372c5p+138,
     0x1.9e4683dffe6c7p+137, 0x1.28774d33b3767p+138, 0x1.4231544b614d0p+137, 0x1.0d042612ff001p+136, 0x1.5c8b976aaa575p+138,
     0x1.0fd8b417d71a3p+138, 0x1.707591c92df99p+137, 0x1.505c87ad60463p+137, 0x1.9b75d2f857218p+136, 0x1.c5873263799b7p+134},
    {0x1.3c54b1045b50fp+146, 0x1.dea3535bd2048p+146, 0x1.c06d2e07bfe99p+146, 0x1.62427685a1c14p+146, 0x1.3271d7751de2fp+147,
     0x1.50e9eca198466p+146, 0x1.61fcbd7664752p+143, 0x1.6b37a7bb4f5f3p+145, 0x1.fbb1eaa327812p+143, 0x1.3a85a4d6cea29p+147,
     0x1.233419231d126p+146, 0x1.83d1a602e7dfdp+144, 0x1.608fb6b637036p+146, 0x1.9403bff5d5d67p+146, 0x1.461b07a497b86p+147,
     0x1.cba159665e773p+143},
    {0x1.cb1f870cdce74p+155, 0x1.a1830bab6e3edp+154, 0x1.b7e5de673dbc8p+155, 0x1.11c2b6dac3edap+156, 0x1.1b1af7d79714ep+156,
     0x1.1551e40423762p+155, 0x1.071119eda75c0p+156, 0x1.1b82ce9921c79p+155, 0x1.ffe22864b4154p+155, 0x1.8aca320b6ea74p+155,
     0x1.ea3f4b6120dd2p+155, 0x1.3ac2f18e13a16p+151, 0x1.b1375e0a34ae5p+155, 0x1.dec713787681ap+155, 0x1.2d50c98d9eb17p+155,
     0x1.1e75a68aae4b8p+156, 0x1.2d57320ca6f54p+153},
    {0x1.d4b279b3a62b2p+164, 0x1.8a9551cb223ffp+162, 0x1.77f56dd1b3d4bp+163, 0x1.4e46ee15a31aap+163, 0x1.f6d2653a03de4p+164,
     0x1.85fc4a6b30f87p+164, 0x1.46533c5e3b219p+161, 0x1.f2c2de6330ed1p+164, 0x1.a904fd427f920p+162, 0x1.b4feb7909e30ep+164,
     0x1.2d84f0d516e01p+164, 0x1.f23ef31370cccp+164, 0x1.f63d1d3373d83p+163, 0x1.2d4333d61c228p+163, 0x1.e0ac26f947c02p+163,
     0x1.6827e413f4f58p+164, 0x1.f4d73ba1eccd4p+163, 0x1.8a1f332c39756p+164},
    {0x1.d064ea13124efp+173, 0x1.d1f3ac70d2baap+172, 0x1.aa443d42e699ap+173, 0x1.720dde174c85dp+170, 0x1.47b6b715bcc11p+173,
     0x1.85d7eb43170e6p+172, 0x1.c345600063459p+171, 0x1.0a1e0e3cc5394p+173, 0x1.7c2d367fb654ap+173, 0x1.820abceb38fb4p+173,
     0x1.d1d59ab7e2a9ap+171, 0x1.03da6974a6f6bp+173, 0x1.6696927eb35d9p+169, 0x1.cbe445e3d49d9p+172, 0x1.210c17e4c258bp+172,
     0x1.bd6a41e0d5c80p+171, 0x1.c265c5dd2076fp+173, 0x1.ec4970954aff9p+168, 0x1.e45b2f93c319bp+169},
    {0x1.73a56a837f7c3p+181, 0x1.72498284889e9p+182, 0x1.65d944ac669f3p+182, 0x1.8e4cf8a364c35p+182, 0x1.9cabcee313387p+180,
     0x1.d3a0f1c2b713bp+181, 0x1.ef36214ce9d0ep+180, 0x1.13aad62c85abfp+181, 0x1.dbacfde5abab1p+181, 0x1.22f1ac9d2e29ep+182,
     0x1.828988c975da7p+182, 0x1.aa5a080c322e5p+177, 0x1.28c5eb96dc02cp+180, 0x1.96615f1e98fb9p+179, 0x1.9052cd84021a4p+182,
     0x1.8be62641331b4p+181, 0x1.c265c5dd2076fp+173, 0x1.3daed1344ef9ap+182, 0x1.f9c9c5bbf032ep+179, 0x1.bb10f26e4ca76p+180},
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
     {0x1.27ea1263b6000p+66, 0x1.4b1c000000000p+25},
     {0x1.4dc944c8ea000p+66, 0x1.8d5d004000000p+26},
     {0x1.6f3654a280000p+64, 0x1.6242000000000p+25},
     {0x1.69d17d2f70000p+63, 0x1.93b8000000000p+26},
     {0x1.2ce880a500000p+64, 0x1.ec78000000000p+23},
     {0x1.b3680e68c8000p+65, 0x1.e2ec000000000p+25},
     {0x1.0696d67ee9000p+67, 0x1.86dc000000000p+26},
     },
    {
     {0x1.1831e11b26000p+75, 0x1.0ee0638000000p+35},
     {0x1.f61ae6add0000p+75, 0x1.3e6fc70040000p+34},
     {0x1.039624ea25000p+76, 0x1.5ecb190000000p+35},
     {0x1.882bc9e298000p+75, 0x1.5114ae0000000p+35},
     {0x1.37bfdfa674000p+75, 0x1.2060e20000000p+33},
     {0x1.4d93978430000p+72, 0x1.4f520e0000000p+34},
     {0x1.082245bce2000p+75, 0x1.535d4a0000000p+33},
     {0x1.3ab2ae0e8a000p+75, 0x1.fb20870000000p+34},
     },
    {
     {0x1.fd51d91aca000p+84, 0x1.beb5388b80000p+43},
     {0x1.047ee39966000p+84, 0x1.611284db00200p+43},
     {0x1.3ffcc17414000p+82, 0x1.b0cfe32200000p+41},
     {0x1.12d8fb4c10000p+80, 0x1.74eb178600000p+43},
     {0x1.55c45f77ca000p+83, 0x1.5aa5f5b980000p+43},
     {0x1.b2fff53c8b000p+84, 0x1.4218f46300000p+43},
     {0x1.a305b35621000p+84, 0x1.84583fad00000p+42},
     {0x1.13b7249d28000p+81, 0x1.3facd88400000p+41},
     {0x1.12e6b858e2000p+83, 0x1.2f29d7d100000p+43},
     },
    {
     {0x1.82f39b11b1000p+93, 0x1.70bb2331ca000p+50},
     {0x1.0e93e98d1f000p+93, 0x1.c85bc04e14801p+52},
     {0x1.f6543fcc60000p+89, 0x1.d86bb16c32000p+51},
     {0x1.f3dcdd1159000p+93, 0x1.31b661270f000p+52},
     {0x1.14c51aa528000p+90, 0x1.4620a84aab400p+52},
     {0x1.b910b15b00000p+92, 0x1.051b206ad6000p+52},
     {0x1.f7d61f3140000p+87, 0x1.e62100c690000p+49},
     {0x1.2a5a3c1a05000p+93, 0x1.d17d53294a000p+50},
     {0x1.b3df3cf706000p+92, 0x1.3433de464d000p+51},
     {0x1.5aa7de2ab0000p+90, 0x1.ccce19121bc00p+52},
     },
    {
     {0x1.e553c81c60000p+99, 0x1.85f7819d23628p+59},
     {0x1.9b96b8790c000p+101, 0x1.eb5c650aa2264p+60},
     {0x1.1f073c432c000p+100, 0x1.aa3e80fc660e6p+61},
     {0x1.db58c69c4f000p+102, 0x1.c53e3208b4108p+59},
     {0x1.1bef826d90000p+98, 0x1.e5aff7cfc3d2cp+60},
     {0x1.2affd32056000p+102, 0x1.431cc2dd89d50p+61},
     {0x1.8106ef4246000p+101, 0x1.70f623cdeab90p+60},
     {0x1.f8048e5dc8000p+101, 0x1.f3e3e18630128p+60},
     {0x1.80709a9ef6000p+102, 0x1.bb11b4e44ac76p+61},
     {0x1.b8ac418c88000p+99, 0x1.9c6e9949f02eep+61},
     {0x1.8cc82d2294000p+102, 0x1.a31fcbafa4980p+57},
     },
    {
     {0x1.5602a250cd000p+111, 0x1.20d4c7e2c60cdp+70},
     {0x1.5802ce9d00000p+103, 0x1.a7850aa1dfb02p+69},
     {0x1.5004299fe0000p+106, 0x1.de99ba4ef3df3p+66},
     {0x1.8c8175b3e9000p+111, 0x1.18d27b1e6314ap+70},
     {0x1.314c5b1138000p+109, 0x1.4d278b9bbffb9p+70},
     {0x1.ae1161cc0c000p+111, 0x1.fe4fe18deb726p+70},
     {0x1.bb6bfd2112000p+110, 0x1.4292aa4ee1e30p+66},
     {0x1.37cb8e2b16000p+110, 0x1.676085b6025fbp+69},
     {0x1.63a338555e000p+110, 0x1.b78179bda10eap+70},
     {0x1.b9376a78d4000p+109, 0x1.c721ad7d97642p+69},
     {0x1.6c544f1649000p+111, 0x1.1a21c8f17f2afp+70},
     {0x1.a5f2f8de80000p+105, 0x1.cf5c385bc4cc5p+66},
     },
    {
     {0x1.a82ecb43cc000p+119, 0x1.9d1d042203844p+79},
     {0x1.011d6fc30c000p+120, 0x1.c66617e04aa5ep+80},
     {0x1.3babe91aa0000p+116, 0x1.3b19f6b872c1bp+80},
     {0x1.269dfbe14e000p+120, 0x1.0fd9cf2fe1366p+80},
     {0x1.77104de0a0000p+116, 0x1.4fd24540b5d6ap+79},
     {0x1.0e7042ccc4000p+119, 0x1.627d9dc73baeap+80},
     {0x1.76692127a8000p+118, 0x1.1b71b501e4a34p+78},
     {0x1.4516448180000p+115, 0x1.23762a397b834p+80},
     {0x1.7fcaceed00000p+115, 0x1.9838ccbd2965fp+80},
     {0x1.19841d3606000p+120, 0x1.0440bb4d14b58p+80},
     {0x1.8caebaa4ba000p+120, 0x1.aeb9ed06f0104p+80},
     {0x1.1b24d44950000p+118, 0x1.804cf621bff33p+77},
     {0x1.5b982b6f0c000p+120, 0x1.a64604e3d67f2p+80},
     },
    {
     {0x1.74690f706c000p+129, 0x1.0d006f3ebb503p+89},
     {0x1.91f47d5346000p+129, 0x1.c8cedf19d30dbp+86},
     {0x1.e84b3bd868000p+127, 0x1.4a75164fa10f0p+86},
     {0x1.2d08aaa0ee000p+129, 0x1.7869b7f2da929p+88},
     {0x1.b931309bc0000p+127, 0x1.bf93d3c11ea52p+85},
     {0x1.db886ee7a0000p+125, 0x1.ca4ebfa87b0afp+89},
     {0x1.a6dbe78218000p+128, 0x1.52a4e14dadba5p+88},
     {0x1.e3213480fc000p+128, 0x1.d40750e786788p+88},
     {0x1.4f02ff63a8000p+128, 0x1.9ec325353c942p+88},
     {0x1.0d06d58d9a000p+129, 0x1.7742a0aeb646ep+88},
     {0x1.9dcf4cbcc8000p+127, 0x1.455f572d49eb5p+86},
     {0x1.8adbafd008000p+129, 0x1.27235cffe8d1cp+89},
     {0x1.e114d24754000p+128, 0x1.3cfa3c3767881p+89},
     {0x1.033812d2e4000p+129, 0x1.3f03a11467bffp+89},
     },
    {
     {0x1.570e37e978000p+136, 0x1.03abd095d4560p+98},
     {0x1.57494ec540000p+136, 0x1.c10775309a562p+96},
     {0x1.25db7628d0000p+136, 0x1.5a2dac8612970p+98},
     {0x1.d6739778a8000p+137, 0x1.9a74782f8edf2p+98},
     {0x1.01afe6c136000p+138, 0x1.2c4d2ff4be6aap+98},
     {0x1.9e4683dffc000p+137, 0x1.36342b15d53b7p+98},
     {0x1.28774d33b2000p+138, 0x1.7677c1050bd58p+98},
     {0x1.4231544b60000p+137, 0x1.4d0533453bdf3p+97},
     {0x1.0d042612f8000p+136, 0x1.c002b009f6ffap+98},
     {0x1.5c8b976aaa000p+138, 0x1.5d42664811da0p+96},
     {0x1.0fd8b417d6000p+138, 0x1.1a2f16b09c0a6p+98},
     {0x1.707591c92c000p+137, 0x1.f993fe22e4773p+97},
     {0x1.505c87ad60000p+137, 0x1.18cedcf8215e2p+95},
     {0x1.9b75d2f850000p+136, 0x1.c85f628dc9f68p+98},
     {0x1.c587326360000p+134, 0x1.99b71aca11a82p+98},
     },
    {
     {0x1.3c54b10458000p+146, 0x1.a87ba0bb2e555p+107},
     {0x1.dea3535bd0000p+146, 0x1.023ddc1441cb9p+107},
     {0x1.c06d2e07bc000p+146, 0x1.f4c9859af09c9p+107},
     {0x1.62427685a0000p+146, 0x1.c13d71e968a85p+106},
     {0x1.3271d7751c000p+147, 0x1.e2ee33fd22722p+107},
     {0x1.50e9eca198000p+146, 0x1.1983694868c87p+104},
     {0x1.61fcbd7660000p+143, 0x1.1d46a3f511c67p+105},
     {0x1.6b37a7bb48000p+145, 0x1.d7cb0b92662cep+107},
     {0x1.fbb1eaa320000p+143, 0x1.e047dad596aa1p+105},
     {0x1.3a85a4d6ce000p+147, 0x1.452066129281dp+106},
     {0x1.233419231c000p+146, 0x1.125941cd28577p+106},
     {0x1.83d1a602e0000p+144, 0x1.f7f206b2b7e46p+106},
     {0x1.608fb6b634000p+146, 0x1.81b35d58a9dc5p+107},
     {0x1.9403bff5d4000p+146, 0x1.d6750616e0b40p+106},
     {0x1.461b07a496000p+147, 0x1.b85bf973caa26p+107},
     {0x1.cba1596640000p+143, 0x1.e7733c2657e4fp+107},
     },
    {
     {0x1.cb1f870cdc000p+155, 0x1.ce8b1bfcd8982p+114},
     {0x1.a1830bab68000p+154, 0x1.8fb3f280bc0f0p+116},
     {0x1.b7e5de673c000p+155, 0x1.bc785437cab1bp+115},
     {0x1.11c2b6dac2000p+156, 0x1.ed9ad8dcde87dp+116},
     {0x1.1b1af7d796000p+156, 0x1.14e3a650d2508p+116},
     {0x1.1551e40420000p+155, 0x1.bb13a84511652p+116},
     {0x1.071119eda6000p+156, 0x1.5bfc505fd104ep+116},
     {0x1.1b82ce9920000p+155, 0x1.c78f1861b7422p+115},
     {0x1.ffe22864b4000p+155, 0x1.538ca5c61b3d3p+111},
     {0x1.8aca320b6c000p+155, 0x1.539fef2a85992p+116},
     {0x1.ea3f4b6120000p+155, 0x1.ba37c9bb9a25ep+114},
     {0x1.3ac2f18e00000p+151, 0x1.3a1580aedab89p+115},
     {0x1.b1375e0a34000p+155, 0x1.5ca421acba983p+114},
     {0x1.dec7137874000p+155, 0x1.40d0ea461ef56p+116},
     {0x1.2d50c98d9c000p+155, 0x1.58ba5c3eeac84p+116},
     {0x1.1e75a68aae000p+156, 0x1.2e127643365d7p+114},
     {0x1.2d57320ca0000p+153, 0x1.bd51a630c0ffdp+115},
     },
    {
     {0x1.d4b279b3a6000p+164, 0x1.592fa22749f3ep+121},
     {0x1.8a9551cb20000p+162, 0x1.1ff9a8d754672p+123},
     {0x1.77f56dd1b0000p+163, 0x1.ea58c329464f1p+124},
     {0x1.4e46ee15a0000p+163, 0x1.8d50bf74ec1c6p+124},
     {0x1.f6d2653a02000p+164, 0x1.de423735adbddp+124},
     {0x1.85fc4a6b30000p+164, 0x1.f0d75fa3777afp+123},
     {0x1.46533c5e30000p+161, 0x1.6432d299cc078p+124},
     {0x1.f2c2de6330000p+164, 0x1.da1a6ac0bd2a5p+123},
     {0x1.a904fd4278000p+162, 0x1.e47f8854b60c5p+124},
     {0x1.b4feb7909e000p+164, 0x1.8716eec35e60cp+121},
     {0x1.2d84f0d516000p+164, 0x1.c023872bc940cp+123},
     {0x1.f23ef31370000p+164, 0x1.998cf1290a4dcp+123},
     {0x1.f63d1d3370000p+163, 0x1.ec17d730e71c1p+124},
     {0x1.2d4333d61c000p+163, 0x1.13c701cce738cp+120},
     {0x1.e0ac26f944000p+163, 0x1.e00db7532d424p+124},
     {0x1.6827e413f4000p+164, 0x1.eb0a08e18790fp+123},
     {0x1.f4d73ba1ec000p+163, 0x1.9a7d5bee5137bp+122},
     {0x1.8a1f332c38000p+164, 0x1.756079230aeebp+124},
     },
    {
     {0x1.d064ea1312000p+173, 0x1.3bd9b34dec159p+131},
     {0x1.d1f3ac70d0000p+172, 0x1.5d4f69861352ep+133},
     {0x1.aa443d42e6000p+173, 0x1.33348f3003652p+132},
     {0x1.720dde1740000p+170, 0x1.90ba587c2fa26p+133},
     {0x1.47b6b715bc000p+173, 0x1.8217d19d8a301p+132},
     {0x1.85d7eb4314000p+172, 0x1.8730e421a9a4ep+133},
     {0x1.c345600060000p+171, 0x1.a2caaaa41ed43p+132},
     {0x1.0a1e0e3cc4000p+173, 0x1.3946b96d83239p+133},
     {0x1.7c2d367fb6000p+173, 0x1.528eab7888f87p+131},
     {0x1.820abceb38000p+173, 0x1.f679a6e9b447bp+132},
     {0x1.d1d59ab7e0000p+171, 0x1.54d07eb0422d0p+132},
     {0x1.03da6974a6000p+173, 0x1.ed5941e7ebed4p+132},
     {0x1.6696927ea0000p+169, 0x1.35d8e5caa1d46p+133},
     {0x1.cbe445e3d4000p+172, 0x1.3b2c3b1961138p+131},
     {0x1.210c17e4c0000p+172, 0x1.2c550c4c305bdp+133},
     {0x1.bd6a41e0d0000p+171, 0x1.7200020c987b4p+133},
     {0x1.c265c5dd20000p+173, 0x1.dba5047d5a946p+131},
     {0x1.ec49709540000p+168, 0x1.5ff25b4509c2ap+131},
     {0x1.e45b2f93c0000p+169, 0x1.8cd65ac4985e2p+130},
     },
    {
     {0x1.73a56a837c000p+181, 0x1.be178c9441c0fp+142},
     {0x1.7249828488000p+182, 0x1.3d1d94e0e620fp+141},
     {0x1.65d944ac66000p+182, 0x1.3e6b4a542dc2dp+141},
     {0x1.8e4cf8a364000p+182, 0x1.86982617c0c21p+141},
     {0x1.9cabcee310000p+180, 0x1.9c35b73f21293p+141},
     {0x1.d3a0f1c2b4000p+181, 0x1.89dbd6ec557b5p+142},
     {0x1.ef36214ce8000p+180, 0x1.d0e530dabd06fp+140},
     {0x1.13aad62c84000p+181, 0x1.abec630d09925p+141},
     {0x1.dbacfde5a8000p+181, 0x1.d58937f13e54ep+142},
     {0x1.22f1ac9d2e000p+182, 0x1.4ee8195f1cc36p+139},
     {0x1.828988c974000p+182, 0x1.da7349324bdadp+142},
     {0x1.aa5a080c00000p+177, 0x1.9172a105df214p+142},
     {0x1.28c5eb96d8000p+180, 0x1.00b1aefa1c409p+142},
     {0x1.96615f1e90000p+179, 0x1.1f72e674951dep+142},
     {0x1.9052cd8402000p+182, 0x1.a3d34368db002p+138},
     {0x1.8be6264130000p+181, 0x1.8da05c144375ep+142},
     {0x1.c265c5dc00000p+173, 0x1.2076e9411f56ap+141},
     {0x1.3daed1344e000p+182, 0x1.f33a989242cf6p+141},
     {0x1.f9c9c5bbf0000p+179, 0x1.96e64210d85fcp+136},
     {0x1.bb10f26e48000p+180, 0x1.29d94ed6a113cp+142},
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
