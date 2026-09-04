#pragma once
#include "common.hpp"

namespace gemmul8::common::table {

//==========
// moduli
//==========
template <Backend BACKEND, unsigned IDX> inline constexpr int32_t moduli = 0;

// INT8: moduli
template <> inline constexpr int32_t moduli<Backend::INT8, 0U>  = 256;
template <> inline constexpr int32_t moduli<Backend::INT8, 1U>  = 255;
template <> inline constexpr int32_t moduli<Backend::INT8, 2U>  = 253;
template <> inline constexpr int32_t moduli<Backend::INT8, 3U>  = 251;
template <> inline constexpr int32_t moduli<Backend::INT8, 4U>  = 247;
template <> inline constexpr int32_t moduli<Backend::INT8, 5U>  = 241;
template <> inline constexpr int32_t moduli<Backend::INT8, 6U>  = 239;
template <> inline constexpr int32_t moduli<Backend::INT8, 7U>  = 233;
template <> inline constexpr int32_t moduli<Backend::INT8, 8U>  = 229;
template <> inline constexpr int32_t moduli<Backend::INT8, 9U>  = 227;
template <> inline constexpr int32_t moduli<Backend::INT8, 10U> = 223;
template <> inline constexpr int32_t moduli<Backend::INT8, 11U> = 217;
template <> inline constexpr int32_t moduli<Backend::INT8, 12U> = 211;
template <> inline constexpr int32_t moduli<Backend::INT8, 13U> = 199;
template <> inline constexpr int32_t moduli<Backend::INT8, 14U> = 197;
template <> inline constexpr int32_t moduli<Backend::INT8, 15U> = 193;
template <> inline constexpr int32_t moduli<Backend::INT8, 16U> = 191;
template <> inline constexpr int32_t moduli<Backend::INT8, 17U> = 181;
template <> inline constexpr int32_t moduli<Backend::INT8, 18U> = 179;
template <> inline constexpr int32_t moduli<Backend::INT8, 19U> = 173;

// FP8: moduli
template <> inline constexpr int32_t moduli<Backend::FP8, 0U>  = 2401; // base-49
template <> inline constexpr int32_t moduli<Backend::FP8, 1U>  = 2209; // base-47
template <> inline constexpr int32_t moduli<Backend::FP8, 2U>  = 2025; // base-45
template <> inline constexpr int32_t moduli<Backend::FP8, 3U>  = 1849; // base-43
template <> inline constexpr int32_t moduli<Backend::FP8, 4U>  = 1681; // base-41
template <> inline constexpr int32_t moduli<Backend::FP8, 5U>  = 1369; // base-37
template <> inline constexpr int32_t moduli<Backend::FP8, 6U>  = 1193; // Karatsuba
template <> inline constexpr int32_t moduli<Backend::FP8, 7U>  = 1097; // Karatsuba
template <> inline constexpr int32_t moduli<Backend::FP8, 8U>  = 1033; // Karatsuba
template <> inline constexpr int32_t moduli<Backend::FP8, 9U>  = 1024; // base-32
template <> inline constexpr int32_t moduli<Backend::FP8, 10U> = 1003; // Karatsuba
template <> inline constexpr int32_t moduli<Backend::FP8, 11U> = 997;  // Karatsuba
template <> inline constexpr int32_t moduli<Backend::FP8, 12U> = 961;  // base-31
template <> inline constexpr int32_t moduli<Backend::FP8, 13U> = 941;  // Karatsuba
template <> inline constexpr int32_t moduli<Backend::FP8, 14U> = 937;  // Karatsuba
template <> inline constexpr int32_t moduli<Backend::FP8, 15U> = 911;  // Karatsuba
template <> inline constexpr int32_t moduli<Backend::FP8, 16U> = 907;  // Karatsuba
template <> inline constexpr int32_t moduli<Backend::FP8, 17U> = 863;  // Karatsuba
template <> inline constexpr int32_t moduli<Backend::FP8, 18U> = 859;  // Karatsuba
template <> inline constexpr int32_t moduli<Backend::FP8, 19U> = 841;  // base-29

inline constexpr int32_t moduli_int8[20] = {
    moduli<Backend::INT8, 0U>,
    moduli<Backend::INT8, 1U>,
    moduli<Backend::INT8, 2U>,
    moduli<Backend::INT8, 3U>,
    moduli<Backend::INT8, 4U>,
    moduli<Backend::INT8, 5U>,
    moduli<Backend::INT8, 6U>,
    moduli<Backend::INT8, 7U>,
    moduli<Backend::INT8, 8U>,
    moduli<Backend::INT8, 9U>,
    moduli<Backend::INT8, 10U>,
    moduli<Backend::INT8, 11U>,
    moduli<Backend::INT8, 12U>,
    moduli<Backend::INT8, 13U>,
    moduli<Backend::INT8, 14U>,
    moduli<Backend::INT8, 15U>,
    moduli<Backend::INT8, 16U>,
    moduli<Backend::INT8, 17U>,
    moduli<Backend::INT8, 18U>,
    moduli<Backend::INT8, 19U>,
};

inline constexpr int32_t moduli_fp8[20] = {
    moduli<Backend::FP8, 0U>,
    moduli<Backend::FP8, 1U>,
    moduli<Backend::FP8, 2U>,
    moduli<Backend::FP8, 3U>,
    moduli<Backend::FP8, 4U>,
    moduli<Backend::FP8, 5U>,
    moduli<Backend::FP8, 6U>,
    moduli<Backend::FP8, 7U>,
    moduli<Backend::FP8, 8U>,
    moduli<Backend::FP8, 9U>,
    moduli<Backend::FP8, 10U>,
    moduli<Backend::FP8, 11U>,
    moduli<Backend::FP8, 12U>,
    moduli<Backend::FP8, 13U>,
    moduli<Backend::FP8, 14U>,
    moduli<Backend::FP8, 15U>,
    moduli<Backend::FP8, 16U>,
    moduli<Backend::FP8, 17U>,
    moduli<Backend::FP8, 18U>,
    moduli<Backend::FP8, 19U>,
};

inline constexpr unsigned k_block_first_fp8[20] = {
    16384, // 2401
    18432, // 2209
    21248, // 2025
    24576, // 1849
    28928, // 1681
    41728, // 1369
    16384, // 1193
    16384, // 1097
    16384, // 1033
    65536, // 1024
    16384, // 1003
    16384, // 997
    74496, // 961
    16384, // 941
    16384, // 937
    16384, // 911
    16384, // 907
    16384, // 863
    16384, // 859
    85504  // 841
};

inline constexpr unsigned k_block_next_fp8[20] = {
    16128, // 2401
    18432, // 2209
    21248, // 2025
    24576, // 1849
    28928, // 1681
    41728, // 1369
    16128, // 1193
    16128, // 1097
    16128, // 1033
    65280, // 1024
    16128, // 1003
    16128, // 997
    74496, // 961
    16128, // 941
    16128, // 937
    16128, // 911
    16128, // 907
    16128, // 863
    16128, // 859
    85504  // 841
};

// 2^32 / moduli
template <Backend BACKEND, unsigned IDX>
inline constexpr int32_t p_inv_32_v = int32_t(4294967296ULL / uint64_t(moduli<BACKEND, IDX>));

// 2^64 / moduli
inline constexpr uint64_t UINT64_MAX_V = 18446744073709551615ULL;
template <Backend BACKEND, unsigned IDX>
inline constexpr int64_t p_inv_64_v =
    int64_t(UINT64_MAX_V / uint64_t(moduli<BACKEND, IDX>)) +
    int64_t(UINT64_MAX_V % uint64_t(moduli<BACKEND, IDX>) == uint64_t(moduli<BACKEND, IDX> - 1));

// FP8: sqrt(moduli)
template <unsigned IDX> inline constexpr int32_t sqrt_moduli = 0;
template <> inline constexpr int32_t sqrt_moduli<0U>         = 49;
template <> inline constexpr int32_t sqrt_moduli<1U>         = 47;
template <> inline constexpr int32_t sqrt_moduli<2U>         = 45;
template <> inline constexpr int32_t sqrt_moduli<3U>         = 43;
template <> inline constexpr int32_t sqrt_moduli<4U>         = 41;
template <> inline constexpr int32_t sqrt_moduli<5U>         = 37;
template <> inline constexpr int32_t sqrt_moduli<9U>         = 32;
template <> inline constexpr int32_t sqrt_moduli<12U>        = 31;
template <> inline constexpr int32_t sqrt_moduli<19U>        = 29;

constexpr bool isKaratsuba[20] = {
    false, false, false, false, false,
    false, true, true, true, false,
    true, true, false, true, true,
    true, true, true, true, false};

constexpr unsigned num_mat_fp8[21] = {
    0,
    2, 4, 6, 8, 10, 12, 15, 18, 21, 23,
    26, 29, 31, 34, 37, 40, 43, 46, 49, 51};

//==========
// number of matrices for workspace of A/B
//==========
template <Backend BACKEND> inline unsigned num_mat(unsigned NUM_MODULI) {
    if constexpr (BACKEND == Backend::INT8) {
        return NUM_MODULI;
    } else {
        return num_mat_fp8[NUM_MODULI];
    }
};
template <Backend BACKEND, unsigned NUM_MODULI>
__host__ __device__ constexpr unsigned num_mat_constexpr() {
    if constexpr (BACKEND == Backend::INT8) {
        return NUM_MODULI;
    } else {
        return num_mat_fp8[NUM_MODULI];
    }
}
template <Backend BACKEND, unsigned NUM_MODULI>
inline constexpr unsigned num_mat_v = num_mat_constexpr<BACKEND, NUM_MODULI>();

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
    { -0x1.43b8040000000p+22,    0x0.0000000000000p+0},
    { -0x1.401552f480000p+33,    0x0.0000000000000p+0},
    { -0x1.20fb4084fe100p+44,    0x0.0000000000000p+0},
    { -0x1.da6474aa5211cp+54,   -0x1.0000000000000p+0},
    { -0x1.3d1c667c5a1c2p+65,   -0x1.2400000000000p+6},
    { -0x1.717256665ffb4p+75,  -0x1.1ca1880000000p+21},
    { -0x1.8bc8bd0f2c52fp+85,   0x1.ec4fd03800000p+29},
    { -0x1.8f4340b88e76bp+95,   0x1.d528e0f31f800p+41},
    {-0x1.8f4340b88e76bp+105,   0x1.d528e0f31f800p+51},
    {-0x1.87131fa4c58acp+115,   0x1.9289ca56231aap+61},
    {-0x1.7cc35e8f2d555p+125,  -0x1.504d5efe89495p+69},
    {-0x1.6556597dde4b5p+135,  -0x1.8f9c9c6660571p+79},
    {-0x1.485f99bcea86bp+145,  -0x1.acce2aae45020p+91},
    {-0x1.2c797a6d1d99cp+155, -0x1.2b5fa68df6a51p+101},
    {-0x1.0b5112aa93159p+165,  0x1.a3534f2667a54p+110},
    {-0x1.d98c1e912b8ebp+174,  0x1.9ad414b6889b5p+119},
    {-0x1.8f17d6c2d8758p+184,  0x1.e07878e9ab41cp+128},
    {-0x1.4ec93f67f3149p+194,  0x1.cc343db811a59p+136},
    {-0x1.12f4c8531f63ap+204, -0x1.72c14309f2704p+149},
};
}

template <Backend BACKEND, typename doublex_t> __forceinline__ doublex_t get_P(unsigned NUM_MODULI);
template <> __forceinline__ double get_P<Backend::INT8, double>(unsigned NUM_MODULI) { return INT8::P[NUM_MODULI - 2].x; }
template <> __forceinline__ double2 get_P<Backend::INT8, double2>(unsigned NUM_MODULI) { return INT8::P[NUM_MODULI - 2]; }
template <> __forceinline__ double get_P<Backend::FP8, double>(unsigned NUM_MODULI) { return FP8::P[NUM_MODULI - 2].x; }
template <> __forceinline__ double2 get_P<Backend::FP8, double2>(unsigned NUM_MODULI) { return FP8::P[NUM_MODULI - 2]; }

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
    0x1.94e504ced568fp-23, 0x1.997e4ff4b4f12p-34, 0x1.c590c0ad928ecp-45, 0x1.144b62a3e8951p-55, 0x1.9d54e995463f5p-66,
    0x1.62c77d2ca796bp-76, 0x1.4b2ba0f34f317p-86, 0x1.4848fcbaab304p-96, 0x1.4848fcbaab304p-106, 0x1.4f2891b7af89ep-116,
    0x1.583c27c41b41dp-126, 0x1.6ecd4901fa69dp-136, 0x1.8f27c1fb144f4p-146, 0x1.b437787735118p-156, 0x1.ea532555f0e70p-166,
    0x1.14c99bb573f30p-175, 0x1.486cb2d326ccbp-185, 0x1.878278c9a4911p-195, 0x1.dcb38fb8f1c64p-205};
}

template <Backend BACKEND>
__forceinline__ double get_invP(unsigned NUM_MODULI) {
    if constexpr (BACKEND == Backend::INT8) return INT8::invP[NUM_MODULI - 2];
    else return FP8::invP[NUM_MODULI - 2];
}

//==========
// log2P[i] = round-down( log2(P-1)/2 - 0.5 ) in float
//==========
template <Backend BACKEND, unsigned NUM_MODULI> inline constexpr float log2P = 0.0F;

// INT8
template <> inline constexpr float log2P<Backend::INT8, 2U>  = 0x1.dfd1ec0000000p+2F;
template <> inline constexpr float log2P<Backend::INT8, 3U>  = 0x1.6fa3360000000p+3F;
template <> inline constexpr float log2P<Backend::INT8, 4U>  = 0x1.ef2ea60000000p+3F;
template <> inline constexpr float log2P<Backend::INT8, 5U>  = 0x1.372d940000000p+4F;
template <> inline constexpr float log2P<Backend::INT8, 6U>  = 0x1.767b2e0000000p+4F;
template <> inline constexpr float log2P<Backend::INT8, 7U>  = 0x1.b5b0280000000p+4F;
template <> inline constexpr float log2P<Backend::INT8, 8U>  = 0x1.f49a020000000p+4F;
template <> inline constexpr float log2P<Backend::INT8, 9U>  = 0x1.19a8580000000p+5F;
template <> inline constexpr float log2P<Backend::INT8, 10U> = 0x1.38f6bc0000000p+5F;
template <> inline constexpr float log2P<Backend::INT8, 11U> = 0x1.582ada0000000p+5F;
template <> inline constexpr float log2P<Backend::INT8, 12U> = 0x1.7736ae0000000p+5F;
template <> inline constexpr float log2P<Backend::INT8, 13U> = 0x1.9619160000000p+5F;
template <> inline constexpr float log2P<Backend::INT8, 14U> = 0x1.b4a4fe0000000p+5F;
template <> inline constexpr float log2P<Backend::INT8, 15U> = 0x1.d321f80000000p+5F;
template <> inline constexpr float log2P<Backend::INT8, 16U> = 0x1.f180a60000000p+5F;
template <> inline constexpr float log2P<Backend::INT8, 17U> = 0x1.07e7f80000000p+6F;
template <> inline constexpr float log2P<Backend::INT8, 18U> = 0x1.16e7e20000000p+6F;
template <> inline constexpr float log2P<Backend::INT8, 19U> = 0x1.25df9a0000000p+6F;
template <> inline constexpr float log2P<Backend::INT8, 20U> = 0x1.34be220000000p+6F;

// FP8
template <> inline constexpr float log2P<Backend::FP8, 2U>  = 0x1.556ae40000000p+3F;
template <> inline constexpr float log2P<Backend::FP8, 3U>  = 0x1.0294120000000p+4F;
template <> inline constexpr float log2P<Backend::FP8, 4U>  = 0x1.59660e0000000p+4F;
template <> inline constexpr float log2P<Backend::FP8, 5U>  = 0x1.af1e960000000p+4F;
template <> inline constexpr float log2P<Backend::FP8, 6U>  = 0x1.013c400000000p+5F;
template <> inline constexpr float log2P<Backend::FP8, 7U>  = 0x1.2a1dec0000000p+5F;
template <> inline constexpr float log2P<Backend::FP8, 8U>  = 0x1.5283a60000000p+5F;
template <> inline constexpr float log2P<Backend::FP8, 9U>  = 0x1.7a90940000000p+5F;
template <> inline constexpr float log2P<Backend::FP8, 10U> = 0x1.a290940000000p+5F;
template <> inline constexpr float log2P<Backend::FP8, 11U> = 0x1.ca71f80000000p+5F;
template <> inline constexpr float log2P<Backend::FP8, 12U> = 0x1.f24a7e0000000p+5F;
template <> inline constexpr float log2P<Backend::FP8, 13U> = 0x1.0cf6580000000p+6F;
template <> inline constexpr float log2P<Backend::FP8, 14U> = 0x1.20b7e80000000p+6F;
template <> inline constexpr float log2P<Backend::FP8, 15U> = 0x1.3476520000000p+6F;
template <> inline constexpr float log2P<Backend::FP8, 16U> = 0x1.481ff20000000p+6F;
template <> inline constexpr float log2P<Backend::FP8, 17U> = 0x1.5bc6540000000p+6F;
template <> inline constexpr float log2P<Backend::FP8, 18U> = 0x1.6f47fa0000000p+6F;
template <> inline constexpr float log2P<Backend::FP8, 19U> = 0x1.82c6300000000p+6F;
template <> inline constexpr float log2P<Backend::FP8, 20U> = 0x1.9634c40000000p+6F;

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

// idx = NUM_MODULI - threshold<Backend::FP8>::P_is_double - 1
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
    {0x1.4716e80000000p+21, 0x1.4059280000000p+21},
    {0x1.7768718000000p+32, 0x1.2c164f4780000p+33, 0x1.30c03bc400000p+32},
    {0x1.1b90b801b3400p+42, 0x1.f54d67ed4f800p+41, 0x1.12485b535f200p+43, 0x1.24937dd37d000p+40},
    {0x1.cb92e5aa7021ep+53, 0x1.ad14489be7f4cp+53, 0x1.d015a9a0ee4f9p+54, 0x1.4a33a585eb974p+54, 0x1.7120ef69e0867p+53},
    {0x1.1669c11536b27p+63, 0x1.bde102231596bp+64, 0x1.5833cff3d2b37p+64, 0x1.4949a84f48a9ep+62, 0x1.c27b7edb9922fp+64,
     0x1.b892cb592db80p+64},
    {0x1.5ed3ffbbbceacp+73, 0x1.c63d8b06da5f0p+73, 0x1.e77cb33ce03e8p+73, 0x1.f8515024bd2bdp+73, 0x1.c21ae8e1d8afcp+69,
     0x1.cf17143943c37p+74, 0x1.98c69c1c4c284p+72},
    {0x1.6d49d9e32452fp+84, 0x1.246756f24b297p+82, 0x1.e3ee7e3cc6d1ep+83, 0x1.60f93899c8f75p+85, 0x1.5a93745d1636ap+80,
     0x1.0dbc0b69b5ad0p+85, 0x1.7ed87e3a73d58p+85, 0x1.c5263df991ba3p+84},
    {0x1.12366621fdd4fp+95, 0x1.4e03a5fa3c8abp+94, 0x1.916e7974a906fp+94, 0x1.9f06a2882f4b5p+94, 0x1.b13a33a5afc04p+90,
     0x1.8ae355843b5f2p+95, 0x1.5406a1fb60a65p+92, 0x1.e4cb024689fa4p+93, 0x1.5ef2ffa4744d9p+94},
    {0x1.75d1ffd26efeap+103, 0x1.79c11417f53b5p+104, 0x1.a0cf801ed9c00p+104, 0x1.41182e18a585cp+105, 0x1.4f5f052c9cd32p+103,
     0x1.6c8e905db7001p+105, 0x1.353c7513fc2e7p+105, 0x1.3617c32f6d034p+105, 0x1.230a5d0668580p+105, 0x1.c1f34a6ffc8ccp+104},
    {0x1.80658931421d8p+112, 0x1.79cc03737ca98p+114, 0x1.33d6ab8e37a02p+115, 0x1.82a211297a0e8p+114, 0x1.866073f4bc0d5p+115,
     0x1.cfeb4d5e4c233p+114, 0x1.5832bc5c3bb44p+115, 0x1.1e9f2b071d43dp+113, 0x1.083ff92145ca1p+115, 0x1.57560a07e86d9p+115,
     0x1.61a4d193782fap+115},
    {0x1.7d3d29afa99adp+124, 0x1.2d20f0da3a735p+125, 0x1.49ce8478fb661p+124, 0x1.b919b5ad7afd6p+124, 0x1.8a2093385f0c5p+123,
     0x1.b0c61ece7df36p+124, 0x1.29761439c724dp+125, 0x1.f3d115363a1b8p+123, 0x1.c0378362686a6p+120, 0x1.37f912b9cf64ap+125,
     0x1.84bc3ac3b0b50p+119, 0x1.51fd471927ba2p+125},
    {0x1.6614d99fee011p+134, 0x1.5eb478be48d7cp+133, 0x1.401a880b5d518p+134, 0x1.eebeb7fe28014p+130, 0x1.af86953098668p+134,
     0x1.047fa667b9cd1p+135, 0x1.3901d1ee256d2p+135, 0x1.a0f2a3b1e723dp+133, 0x1.3386490d28a18p+135, 0x1.6f1bb5f04f5f6p+133,
     0x1.7614f8de1e90bp+134, 0x1.8f4585a170093p+134, 0x1.64f728a63a800p+131},
    {0x1.43dc351b21f1fp+139, 0x1.8eafa35cf7322p+144, 0x1.e7c6c5449dbb5p+143, 0x1.94eabd4c15b3dp+142, 0x1.081af04d9fd7bp+142,
     0x1.ff63b3c527ac9p+144, 0x1.04efe91b29849p+144, 0x1.0c34e172dcf56p+144, 0x1.017c477b1a87ep+145, 0x1.4a4c292385e68p+144,
     0x1.35b64ea3e3119p+145, 0x1.da47d50fd4a8ep+144, 0x1.47592cd66eb8fp+145, 0x1.1e2625a9cb025p+143},
    {0x1.7bb1117ad9620p+154, 0x1.d8d0c4d13f9d4p+153, 0x1.08b6df6755dfdp+153, 0x1.e82b66ad1f569p+154, 0x1.ab34dca9c1eb4p+154,
     0x1.e03b8386efeacp+154, 0x1.598eec8010045p+153, 0x1.ddb105d73e3a4p+153, 0x1.173d982fe0081p+155, 0x1.1cec30d6f7d19p+155,
     0x1.7f74ff448ce3ep+153, 0x1.28dba49a88d4cp+155, 0x1.06f44c810e2efp+155, 0x1.d607c0deec322p+152, 0x1.27aa13ec992a5p+154},
    {0x1.619a0c8959899p+162, 0x1.8a044e0a9a730p+164, 0x1.eb55ad7f4dceap+164, 0x1.c27dd3dc95d8cp+163, 0x1.90bc8b7595122p+162,
     0x1.01bfadc47224ep+162, 0x1.8fbe1e0f2b584p+162, 0x1.265d7a117a874p+162, 0x1.ae9ae1d906b3fp+159, 0x1.0ded5d593d854p+164,
     0x1.9eb35f1046f44p+163, 0x1.f1a1f71913dd2p+162, 0x1.421d7f2f0a35bp+164, 0x1.ec059458b3503p+163, 0x1.7c93a977149cep+164,
     0x1.ab3cba132616ap+162},
    {0x1.71368ce752afep+174, 0x1.52b51cb364cdep+174, 0x1.b60086f88e8fcp+174, 0x1.a58e98a656673p+174, 0x1.c34ae71146e70p+173,
     0x1.95675035f8093p+172, 0x1.4fcf4cd318aa4p+174, 0x1.c7d942878b3dfp+171, 0x1.ca6b643845ec0p+174, 0x1.3db3c181e4f90p+174,
     0x1.cbdb04ecf0decp+173, 0x1.e19f32dc0e005p+173, 0x1.65bf5daf7803ep+173, 0x1.6c58382a0c4afp+172, 0x1.231a4d42045f6p+168,
     0x1.aa3e9d571ccf2p+173, 0x1.e055ad8a804acp+173},
    {0x1.b27f75466e0f9p+183, 0x1.621b68e9251d8p+181, 0x1.260c41305dbf2p+183, 0x1.0e04ea0d50221p+184, 0x1.502daa65b0b36p+183,
     0x1.334388bc28fd8p+183, 0x1.0b9f8b1425dc3p+184, 0x1.98ea73f77b12ep+183, 0x1.a206229efc012p+183, 0x1.0c87ca409a251p+184,
     0x1.84598edfcd45bp+180, 0x1.d38b1bf2e679cp+183, 0x1.0ec4c91da2936p+182, 0x1.89946188621e7p+180, 0x1.fada3eebfb22cp+183,
     0x1.6659f8c210339p+183, 0x1.fae58bd553483p+180, 0x1.d225ee16e6e07p+180},
    {0x1.32c26c35e0c1bp+193, 0x1.3662bc596b811p+186, 0x1.7848774e05f58p+192, 0x1.35fafe7b1b7f3p+193, 0x1.5e850994e3178p+191,
     0x1.a78e9925d39aep+193, 0x1.ffdc5aef2f697p+193, 0x1.23251891f0906p+193, 0x1.a16dea7140edcp+193, 0x1.454e0c5c41724p+194,
     0x1.dea5b95820935p+193, 0x1.df835c19a3082p+192, 0x1.e38a8f423986fp+192, 0x1.85ee7a46786f0p+191, 0x1.9f2d7bb62fe8ep+193,
     0x1.ddbdabb147c49p+192, 0x1.4dadc4b33f9cdp+193, 0x1.3658aa5057b1cp+192, 0x1.3e06ff23447dap+194},
    {0x1.ade5ee77f0e41p+203, 0x1.88550cb08d358p+202, 0x1.046d7958b7b5cp+203, 0x1.1b6eb0df07d9cp+203, 0x1.1fe0c31e87821p+202,
     0x1.9b546278fbe84p+200, 0x1.5b8e40d012380p+202, 0x1.06acb51bd61c7p+204, 0x1.d3ee658646457p+203, 0x1.ecfce32d0947ap+202,
     0x1.0ca6aecd42a85p+204, 0x1.8f55a561df044p+202, 0x1.e8af0f54116dbp+203, 0x1.38a64119b057cp+203, 0x1.45228c2a8e93ep+202,
     0x1.5eb6661648a28p+203, 0x1.44f9cc4f0d1bcp+202, 0x1.bc22aeaba4e8cp+203, 0x1.3490c4338f686p+202, 0x1.c9b724b01e562p+203},
};

// idx = NUM_MODULI - threshold<Backend::FP8>::P_is_double - 1
// qPi_2[idx][i][1] = first (53-ceil(log2(rho))) bits of q[i]*P[i]/p[i] for rho = sum(floor(p[:]/2)),
// qPi_2[idx][i][2] = double(q[i]*P[i]/p[i] - qPi_2[idx][i][1])
inline constexpr double2 qPi_2[15][20] = {
    {
     {0x1.1669c11534000p+63, 0x1.5938a70000000p+24},
     {0x1.bde1022314000p+64, 0x1.96af390000000p+24},
     {0x1.5833cff3d2000p+64, 0x1.66e6d60000000p+23},
     {0x1.4949a84f48000p+62, 0x1.53bf800000000p+21},
     {0x1.c27b7edb98000p+64, 0x1.22edda0000000p+24},
     {0x1.b892cb592c000p+64, 0x1.b7ffc70000000p+24},
     },
    {
     {0x1.5ed3ffbbbc000p+73, 0x1.d588495a00000p+32},
     {0x1.c63d8b06d8000p+73, 0x1.2f80a5e5c0000p+34},
     {0x1.e77cb33ce0000p+73, 0x1.f4388df000000p+30},
     {0x1.f8515024bc000p+73, 0x1.2bcc55ff80000p+33},
     {0x1.c21ae8e1c0000p+69, 0x1.8afc3a1000000p+33},
     {0x1.cf17143942000p+74, 0x1.c376e53280000p+34},
     {0x1.98c69c1c48000p+72, 0x1.0a105bc340000p+34},
     },
    {
     {0x1.6d49d9e324000p+84, 0x1.4bbff6cc8d000p+42},
     {0x1.246756f240000p+82, 0x1.652e706077600p+45},
     {0x1.e3ee7e3cc0000p+83, 0x1.b4778ce87cd80p+45},
     {0x1.60f93899c8000p+85, 0x1.ee9b384cf3100p+44},
     {0x1.5a93745d00000p+80, 0x1.636a623d1de00p+44},
     {0x1.0dbc0b69b4000p+85, 0x1.ad049d67fe280p+45},
     {0x1.7ed87e3a72000p+85, 0x1.d579edb4f1100p+45},
     {0x1.c5263df990000p+84, 0x1.ba31f47883400p+44},
     },
    {
     {0x1.12366621fc000p+95, 0x1.d4f5b5217379ap+55},
     {0x1.4e03a5fa3c000p+94, 0x1.1556efa9b118ep+53},
     {0x1.916e7974a8000p+94, 0x1.06ea90bb8cadap+54},
     {0x1.9f06a2882c000p+94, 0x1.a5a671b542a11p+55},
     {0x1.b13a33a580000p+90, 0x1.7e0210788d45dp+55},
     {0x1.8ae355843a000p+95, 0x1.5f1b2248c5a3bp+55},
     {0x1.5406a1fb60000p+92, 0x1.4c93613ee524ep+51},
     {0x1.e4cb024688000p+93, 0x1.fa422f0ca6ad2p+53},
     {0x1.5ef2ffa474000p+94, 0x1.36246dd4e9796p+52},
     },
    {
     {0x1.75d1ffd268000p+103, 0x1.bfa71b6c1e59ep+65},
     {0x1.79c11417f4000p+104, 0x1.3b4e1cdb629edp+64},
     {0x1.a0cf801ed8000p+104, 0x1.c003a55da46cep+64},
     {0x1.41182e18a4000p+105, 0x1.85c7795602231p+65},
     {0x1.4f5f052c98000p+103, 0x1.34c96a7ec9372p+65},
     {0x1.6c8e905db6000p+105, 0x1.00129f156d7a6p+65},
     {0x1.353c7513fc000p+105, 0x1.734a880f2153cp+62},
     {0x1.3617c32f6c000p+105, 0x1.033f172cfe376p+65},
     {0x1.230a5d0668000p+105, 0x1.5fee47e93c8b4p+63},
     {0x1.c1f34a6ffc000p+104, 0x1.19722a3b73f02p+63},
     },
    {
     {0x1.8065893140000p+112, 0x1.0ec2b170dbcb2p+73},
     {0x1.79cc037378000p+114, 0x1.2a5e916129105p+76},
     {0x1.33d6ab8e34000p+115, 0x1.d012bb3ed733ap+76},
     {0x1.82a2112978000p+114, 0x1.073d7ab44a8a7p+75},
     {0x1.866073f4bc000p+115, 0x1.a9727bb90374fp+70},
     {0x1.cfeb4d5e48000p+114, 0x1.08ca6e0775b3ap+76},
     {0x1.5832bc5c38000p+115, 0x1.da1ff0f6b4b9bp+76},
     {0x1.1e9f2b0710000p+113, 0x1.a87aba32c938ep+76},
     {0x1.083ff92144000p+115, 0x1.ca09022f63b0dp+75},
     {0x1.57560a07e8000p+115, 0x1.b62a299879ce1p+73},
     {0x1.61a4d19378000p+115, 0x1.7ccf2e4158150p+72},
     },
    {
     {0x1.7d3d29afa8000p+124, 0x1.9ad5d711034fbp+84},
     {0x1.2d20f0da38000p+125, 0x1.39ab0c1c46e1bp+86},
     {0x1.49ce8478f8000p+124, 0x1.b30bb3a7a3c49p+85},
     {0x1.b919b5ad78000p+124, 0x1.7eacf4161c068p+85},
     {0x1.8a20933850000p+123, 0x1.e18a376451a8cp+86},
     {0x1.b0c61ece78000p+124, 0x1.7cd9059d2f303p+86},
     {0x1.29761439c4000p+125, 0x1.92699450dd385p+86},
     {0x1.f3d1153630000p+123, 0x1.437085978bf85p+86},
     {0x1.c037836200000p+120, 0x1.a1a984d295b26p+86},
     {0x1.37f912b9cc000p+125, 0x1.b253afc5b24a8p+86},
     {0x1.84bc3ac300000p+119, 0x1.616a0f71af97ep+86},
     {0x1.51fd471924000p+125, 0x1.dd128034d3c22p+86},
     },
    {
     {0x1.6614d99fe8000p+134, 0x1.804366599cffcp+96},
     {0x1.5eb478be40000p+133, 0x1.1af89f59e7d83p+96},
     {0x1.401a880b58000p+134, 0x1.54603412b3a67p+96},
     {0x1.eebeb7fe00000p+130, 0x1.400a2cb41e94bp+95},
     {0x1.af86953098000p+134, 0x1.9a1e9cbc27b80p+92},
     {0x1.047fa667b8000p+135, 0x1.cd110267b2879p+95},
     {0x1.3901d1ee24000p+135, 0x1.6d1b3cfdb3b33p+95},
     {0x1.a0f2a3b1e0000p+133, 0x1.c8f5f437f3bf9p+95},
     {0x1.3386490d28000p+135, 0x1.42f080b9cb85dp+94},
     {0x1.6f1bb5f040000p+133, 0x1.ebec19513c95ap+96},
     {0x1.7614f8de18000p+134, 0x1.a42cb31cdd20bp+96},
     {0x1.8f4585a170000p+134, 0x1.25c6522688c74p+89},
     {0x1.64f728a600000p+131, 0x1.d3ffe1da44487p+96},
     },
    {
     {0x1.43dc351b00000p+139, 0x1.0f8f7e53a1e3dp+104},
     {0x1.8eafa35cf0000p+144, 0x1.cc869b9220950p+106},
     {0x1.e7c6c54490000p+143, 0x1.b76ade2614d3ep+106},
     {0x1.94eabd4c00000p+142, 0x1.5b3d27908d169p+106},
     {0x1.081af04d80000p+142, 0x1.fd7b66d1e31e5p+106},
     {0x1.ff63b3c520000p+144, 0x1.eb2548f4c6dcfp+106},
     {0x1.04efe91b28000p+144, 0x1.849578be16ca7p+104},
     {0x1.0c34e172d8000p+144, 0x1.3d5869abe451bp+106},
     {0x1.017c477b18000p+145, 0x1.43eee8d647651p+106},
     {0x1.4a4c292380000p+144, 0x1.79a031515fee5p+106},
     {0x1.35b64ea3e0000p+145, 0x1.88c7fe3fbb742p+106},
     {0x1.da47d50fd0000p+144, 0x1.2a393d60356dbp+106},
     {0x1.47592cd66c000p+145, 0x1.5c741afef132fp+106},
     {0x1.1e2625a9c0000p+143, 0x1.604a0a000d27cp+106},
     },
    {
     {0x1.7bb1117ad8000p+154, 0x1.620723191a4e8p+114},
     {0x1.d8d0c4d130000p+153, 0x1.f3a7f5db13356p+116},
     {0x1.08b6df6750000p+153, 0x1.77f5767cf962ap+115},
     {0x1.e82b66ad18000p+154, 0x1.d5a572a252935p+116},
     {0x1.ab34dca9c0000p+154, 0x1.eb39228787f5ep+114},
     {0x1.e03b8386e8000p+154, 0x1.faae83c76f8c5p+116},
     {0x1.598eec8010000p+153, 0x1.141143ba9bf7dp+107},
     {0x1.ddb105d730000p+153, 0x1.c74785a835fbbp+116},
     {0x1.173d982fe0000p+155, 0x1.02c6172205fe2p+110},
     {0x1.1cec30d6f4000p+155, 0x1.e8c79fc1e65e4p+116},
     {0x1.7f74ff4480000p+153, 0x1.9c7c47938ec65p+116},
     {0x1.28dba49a88000p+155, 0x1.a98b4c0740a4dp+114},
     {0x1.06f44c810c000p+155, 0x1.17766cf289021p+116},
     {0x1.d607c0dee0000p+152, 0x1.86434b8c60892p+115},
     {0x1.27aa13ec98000p+154, 0x1.2a4b805e85b7ap+114},
     },
    {
     {0x1.619a0c8950000p+162, 0x1.3132904e6eba8p+125},
     {0x1.8a044e0a98000p+164, 0x1.397c74aaa317bp+125},
     {0x1.eb55ad7f4c000p+164, 0x1.ce9c005ab5990p+124},
     {0x1.c27dd3dc90000p+163, 0x1.762e229d6b57bp+125},
     {0x1.90bc8b7590000p+162, 0x1.44891be3ab169p+124},
     {0x1.01bfadc470000p+162, 0x1.126e1ae67688ep+123},
     {0x1.8fbe1e0f20000p+162, 0x1.6b076bd1d588fp+125},
     {0x1.265d7a1170000p+162, 0x1.50e8e5f0886e5p+125},
     {0x1.ae9ae1d900000p+159, 0x1.acfbf4666943fp+121},
     {0x1.0ded5d593c000p+164, 0x1.85429928c1277p+124},
     {0x1.9eb35f1040000p+163, 0x1.bd104af2ae712p+125},
     {0x1.f1a1f71910000p+162, 0x1.ee91e0a55cbc6p+123},
     {0x1.421d7f2f08000p+164, 0x1.1ad5fc0f39cc6p+125},
     {0x1.ec059458b0000p+163, 0x1.a81995a0e1650p+124},
     {0x1.7c93a97714000p+164, 0x1.39cefce8a5a04p+123},
     {0x1.ab3cba1320000p+162, 0x1.85a8e9ac00d1ep+124},
     },
    {
     {0x1.71368ce750000p+174, 0x1.57ed359c89d02p+135},
     {0x1.52b51cb364000p+174, 0x1.9bc5b19ee9f0bp+133},
     {0x1.b60086f88c000p+174, 0x1.47dc7d2eb0cadp+135},
     {0x1.a58e98a654000p+174, 0x1.3399329ad819dp+135},
     {0x1.c34ae71140000p+173, 0x1.b9bfef5531ad0p+135},
     {0x1.95675035f0000p+172, 0x1.01264fff1d69ep+135},
     {0x1.4fcf4cd318000p+174, 0x1.548b7cf7523aap+133},
     {0x1.c7d9428780000p+171, 0x1.67bee9e87f90ep+134},
     {0x1.ca6b643844000p+174, 0x1.ec01b9b7e9bd9p+134},
     {0x1.3db3c181e4000p+174, 0x1.f1f0d980dc6a2p+133},
     {0x1.cbdb04ecf0000p+173, 0x1.bd8dddb92534ep+132},
     {0x1.e19f32dc08000p+173, 0x1.8012838f8791bp+135},
     {0x1.65bf5daf78000p+173, 0x1.f26accbbdf97cp+126},
     {0x1.6c58382a00000p+172, 0x1.895e6cb164362p+135},
     {0x1.231a4d4200000p+168, 0x1.17d6cc5029b81p+130},
     {0x1.aa3e9d5718000p+173, 0x1.33c6c6f6d26d8p+135},
     {0x1.e055ad8a80000p+173, 0x1.2ae374326e37fp+131},
     },
    {
     {0x1.b27f754668000p+183, 0x1.83e24d23c4a5ep+145},
     {0x1.621b68e920000p+181, 0x1.476000f7b40acp+143},
     {0x1.260c413058000p+183, 0x1.6fc775f12442bp+145},
     {0x1.0e04ea0d50000p+184, 0x1.10a566667eb1ap+141},
     {0x1.502daa65b0000p+183, 0x1.66c282c937a5bp+142},
     {0x1.334388bc28000p+183, 0x1.fb0a6309e8816p+142},
     {0x1.0b9f8b1424000p+184, 0x1.dc3319a633883p+144},
     {0x1.98ea73f778000p+183, 0x1.896f8511643e5p+144},
     {0x1.a206229ef8000p+183, 0x1.0048cd464edddp+145},
     {0x1.0c87ca4098000p+184, 0x1.12870e5b78526p+145},
     {0x1.84598edfc0000p+180, 0x1.a8b53e48a2ef9p+143},
     {0x1.d38b1bf2e0000p+183, 0x1.9e6e6005ccaaep+145},
     {0x1.0ec4c91da0000p+182, 0x1.49b03ee8728ddp+143},
     {0x1.8994618840000p+180, 0x1.10f39057b246ap+145},
     {0x1.fada3eebf8000p+183, 0x1.915ed7a433155p+144},
     {0x1.6659f8c210000p+183, 0x1.9c5b0dc2dc13ep+140},
     {0x1.fae58bd540000p+180, 0x1.3482d08b76d33p+144},
     {0x1.d225ee16c0000p+180, 0x1.370390b973b9cp+145},
     },
    {
     {0x1.32c26c35e0000p+193, 0x1.8366630155ea9p+152},
     {0x1.3662bc5800000p+186, 0x1.6b810a97c49bdp+154},
     {0x1.7848774e00000p+192, 0x1.7d5ed2a69e1d2p+154},
     {0x1.35fafe7b18000p+193, 0x1.bf95aae47e536p+154},
     {0x1.5e850994e0000p+191, 0x1.8bc0de6bfd67bp+152},
     {0x1.a78e9925d0000p+193, 0x1.cd6d51de4e80ap+154},
     {0x1.ffdc5aef28000p+193, 0x1.da5dec7c9e0bap+155},
     {0x1.23251891f0000p+193, 0x1.20ba33758937cp+152},
     {0x1.a16dea7140000p+193, 0x1.db8b257ec7ed6p+152},
     {0x1.454e0c5c40000p+194, 0x1.723a7c350f41ep+154},
     {0x1.dea5b95820000p+193, 0x1.26a0da3e32510p+152},
     {0x1.df835c19a0000p+192, 0x1.84118e47c94c9p+153},
     {0x1.e38a8f4230000p+192, 0x1.30dd7249b239dp+155},
     {0x1.85ee7a4660000p+191, 0x1.86f00f90bb8b5p+155},
     {0x1.9f2d7bb628000p+193, 0x1.fa376f0ba91dbp+155},
     {0x1.ddbdabb140000p+192, 0x1.f122a1d836e9ap+154},
     {0x1.4dadc4b338000p+193, 0x1.e7346bf0e42a0p+155},
     {0x1.3658aa5050000p+192, 0x1.ec70c468f33a7p+154},
     {0x1.3e06ff2344000p+194, 0x1.f682047ffe973p+152},
     },
    {
     {0x1.ade5ee77f0000p+203, 0x1.c82542543bbc2p+162},
     {0x1.88550cb080000p+202, 0x1.a6b07746e935fp+165},
     {0x1.046d7958b0000p+203, 0x1.ed6e917e8b142p+165},
     {0x1.1b6eb0df00000p+203, 0x1.f670cc6239acap+165},
     {0x1.1fe0c31e80000p+202, 0x1.e0839356186f3p+164},
     {0x1.9b546278c0000p+200, 0x1.df4225487c736p+165},
     {0x1.5b8e40d010000p+202, 0x1.1bfeaa603a21ap+163},
     {0x1.06acb51bd4000p+204, 0x1.0e38ceb55cbf8p+165},
     {0x1.d3ee658640000p+203, 0x1.915c7a2e7cb33p+165},
     {0x1.ecfce32d00000p+202, 0x1.28f4a23020ccbp+165},
     {0x1.0ca6aecd40000p+204, 0x1.542740e7b5a22p+165},
     {0x1.8f55a561d0000p+202, 0x1.e08821f692905p+165},
     {0x1.e8af0f5410000p+203, 0x1.6db04d96bbaf3p+163},
     {0x1.38a64119b0000p+203, 0x1.5ee3d2997d4f5p+161},
     {0x1.45228c2a80000p+202, 0x1.d27befed2e0bfp+165},
     {0x1.5eb6661648000p+203, 0x1.44fdb1be7fe89p+162},
     {0x1.44f9cc4f00000p+202, 0x1.a3772f2c675abp+165},
     {0x1.bc22aeaba0000p+203, 0x1.3a2ea5f398a0dp+165},
     {0x1.3490c43380000p+202, 0x1.ed0bacdc845abp+165},
     {0x1.c9b724b018000p+203, 0x1.958710ad0939ep+165},
     },
};
} // namespace FP8

namespace INT8 {
template <unsigned NUM_MODULI, unsigned IDX> inline constexpr double qPi_double_v   = INT8::qPi_1[NUM_MODULI - 2][IDX];
template <unsigned NUM_MODULI, unsigned IDX> inline constexpr double2 qPi_double2_v = INT8::qPi_2[NUM_MODULI - threshold<Backend::INT8>::P_is_double - 1][IDX];
} // namespace INT8

namespace FP8 {
template <unsigned NUM_MODULI, unsigned IDX> inline constexpr double qPi_double_v   = FP8::qPi_1[NUM_MODULI - 2][IDX];
template <unsigned NUM_MODULI, unsigned IDX> inline constexpr double2 qPi_double2_v = FP8::qPi_2[NUM_MODULI - threshold<Backend::FP8>::P_is_double - 1][IDX];
} // namespace FP8

template <Backend BACKEND, unsigned NUM_MODULI, unsigned IDX>
__device__ __forceinline__ constexpr double qPi_double() {
    if constexpr (BACKEND == Backend::INT8) {
        return INT8::qPi_double_v<NUM_MODULI, IDX>;
    } else {
        return FP8::qPi_double_v<NUM_MODULI, IDX>;
    }
}

template <Backend BACKEND, unsigned NUM_MODULI, unsigned IDX>
__device__ __forceinline__ constexpr double2 qPi_double2() {
    if constexpr (BACKEND == Backend::INT8) {
        return INT8::qPi_double2_v<NUM_MODULI, IDX>;
    } else {
        return FP8::qPi_double2_v<NUM_MODULI, IDX>;
    }
}

} // namespace gemmul8::common::table
