#pragma once
#include "common.hpp"
#include "fp8_plan.hpp"

namespace gemmul8::common::table {

//==========
// moduli
//==========
template <Backend BACKEND, unsigned IDX, bool COMPLEX = false>
inline constexpr int32_t moduli = (BACKEND == Backend::FP8 && IDX >= 20U) ? int32_t(IDX - 20U) : 0;

template <Backend BACKEND, unsigned NUM_MODULI, unsigned IDX, bool COMPLEX = false>
inline constexpr unsigned active_index = BACKEND == Backend::FP8 ? fp8_plan::index(NUM_MODULI, IDX, COMPLEX) : IDX;

// INT8 real: moduli
template <> inline constexpr int32_t moduli<Backend::INT8, 0U, false>  = 256;
template <> inline constexpr int32_t moduli<Backend::INT8, 1U, false>  = 255;
template <> inline constexpr int32_t moduli<Backend::INT8, 2U, false>  = 253;
template <> inline constexpr int32_t moduli<Backend::INT8, 3U, false>  = 251;
template <> inline constexpr int32_t moduli<Backend::INT8, 4U, false>  = 247;
template <> inline constexpr int32_t moduli<Backend::INT8, 5U, false>  = 241;
template <> inline constexpr int32_t moduli<Backend::INT8, 6U, false>  = 239;
template <> inline constexpr int32_t moduli<Backend::INT8, 7U, false>  = 233;
template <> inline constexpr int32_t moduli<Backend::INT8, 8U, false>  = 229;
template <> inline constexpr int32_t moduli<Backend::INT8, 9U, false>  = 227;
template <> inline constexpr int32_t moduli<Backend::INT8, 10U, false> = 223;
template <> inline constexpr int32_t moduli<Backend::INT8, 11U, false> = 217;
template <> inline constexpr int32_t moduli<Backend::INT8, 12U, false> = 211;
template <> inline constexpr int32_t moduli<Backend::INT8, 13U, false> = 199;
template <> inline constexpr int32_t moduli<Backend::INT8, 14U, false> = 197;
template <> inline constexpr int32_t moduli<Backend::INT8, 15U, false> = 193;
template <> inline constexpr int32_t moduli<Backend::INT8, 16U, false> = 191;
template <> inline constexpr int32_t moduli<Backend::INT8, 17U, false> = 181;
template <> inline constexpr int32_t moduli<Backend::INT8, 18U, false> = 179;
template <> inline constexpr int32_t moduli<Backend::INT8, 19U, false> = 173;

// INT8 complex: moduli
template <> inline constexpr int32_t moduli<Backend::INT8, 0U, true>  = 241;
template <> inline constexpr int32_t moduli<Backend::INT8, 1U, true>  = 233;
template <> inline constexpr int32_t moduli<Backend::INT8, 2U, true>  = 229;
template <> inline constexpr int32_t moduli<Backend::INT8, 3U, true>  = 221;
template <> inline constexpr int32_t moduli<Backend::INT8, 4U, true>  = 205;
template <> inline constexpr int32_t moduli<Backend::INT8, 5U, true>  = 197;
template <> inline constexpr int32_t moduli<Backend::INT8, 6U, true>  = 193;
template <> inline constexpr int32_t moduli<Backend::INT8, 7U, true>  = 181;
template <> inline constexpr int32_t moduli<Backend::INT8, 8U, true>  = 173;
template <> inline constexpr int32_t moduli<Backend::INT8, 9U, true>  = 157;
template <> inline constexpr int32_t moduli<Backend::INT8, 10U, true> = 149;
template <> inline constexpr int32_t moduli<Backend::INT8, 11U, true> = 137;
template <> inline constexpr int32_t moduli<Backend::INT8, 12U, true> = 113;
template <> inline constexpr int32_t moduli<Backend::INT8, 13U, true> = 109;
template <> inline constexpr int32_t moduli<Backend::INT8, 14U, true> = 101;
template <> inline constexpr int32_t moduli<Backend::INT8, 15U, true> = 97;
template <> inline constexpr int32_t moduli<Backend::INT8, 16U, true> = 89;
template <> inline constexpr int32_t moduli<Backend::INT8, 17U, true> = 73;
template <> inline constexpr int32_t moduli<Backend::INT8, 18U, true> = 61;
template <> inline constexpr int32_t moduli<Backend::INT8, 19U, true> = 53;

// FP8 real: moduli
template <> inline constexpr int32_t moduli<Backend::FP8, 0U, false>  = 2401; // base-49
template <> inline constexpr int32_t moduli<Backend::FP8, 1U, false>  = 2209; // base-47
template <> inline constexpr int32_t moduli<Backend::FP8, 2U, false>  = 2025; // base-45
template <> inline constexpr int32_t moduli<Backend::FP8, 3U, false>  = 1849; // base-43
template <> inline constexpr int32_t moduli<Backend::FP8, 4U, false>  = 1681; // base-41
template <> inline constexpr int32_t moduli<Backend::FP8, 5U, false>  = 1369; // base-37
template <> inline constexpr int32_t moduli<Backend::FP8, 6U, false>  = 1193; // BASE49
template <> inline constexpr int32_t moduli<Backend::FP8, 7U, false>  = 1097; // BASE49
template <> inline constexpr int32_t moduli<Backend::FP8, 8U, false>  = 1033; // BASE33_SUM
template <> inline constexpr int32_t moduli<Backend::FP8, 9U, false>  = 1024; // base-32
template <> inline constexpr int32_t moduli<Backend::FP8, 10U, false> = 1003; // BASE49
template <> inline constexpr int32_t moduli<Backend::FP8, 11U, false> = 997;  // BASE49
template <> inline constexpr int32_t moduli<Backend::FP8, 12U, false> = 961;  // base-31
template <> inline constexpr int32_t moduli<Backend::FP8, 13U, false> = 941;  // BASE32_DIFF
template <> inline constexpr int32_t moduli<Backend::FP8, 14U, false> = 937;  // BASE49
template <> inline constexpr int32_t moduli<Backend::FP8, 15U, false> = 911;  // BASE32_SUM
template <> inline constexpr int32_t moduli<Backend::FP8, 16U, false> = 907;  // BASE32_SUM
template <> inline constexpr int32_t moduli<Backend::FP8, 17U, false> = 863;  // BASE49
template <> inline constexpr int32_t moduli<Backend::FP8, 18U, false> = 859;  // BASE49
template <> inline constexpr int32_t moduli<Backend::FP8, 19U, false> = 841;  // base-29

// FP8 complex: moduli
template <> inline constexpr int32_t moduli<Backend::FP8, 0U, true>  = 1681; // base-41
template <> inline constexpr int32_t moduli<Backend::FP8, 1U, true>  = 1369; // base-37
template <> inline constexpr int32_t moduli<Backend::FP8, 2U, true>  = 1193; // BASE49
template <> inline constexpr int32_t moduli<Backend::FP8, 3U, true>  = 1097; // BASE49
template <> inline constexpr int32_t moduli<Backend::FP8, 4U, true>  = 1033; // BASE33_SUM
template <> inline constexpr int32_t moduli<Backend::FP8, 5U, true>  = 997;  // BASE49
template <> inline constexpr int32_t moduli<Backend::FP8, 6U, true>  = 941;  // BASE32_DIFF
template <> inline constexpr int32_t moduli<Backend::FP8, 7U, true>  = 937;  // BASE49
template <> inline constexpr int32_t moduli<Backend::FP8, 8U, true>  = 905;  // BASE32_SUM
template <> inline constexpr int32_t moduli<Backend::FP8, 9U, true>  = 901;  // BASE32_SUM
template <> inline constexpr int32_t moduli<Backend::FP8, 10U, true> = 857;  // BASE49
template <> inline constexpr int32_t moduli<Backend::FP8, 11U, true> = 853;  // BASE32_SUM
template <> inline constexpr int32_t moduli<Backend::FP8, 12U, true> = 841;  // base-29
template <> inline constexpr int32_t moduli<Backend::FP8, 13U, true> = 809;  // BASE32_DIFF
template <> inline constexpr int32_t moduli<Backend::FP8, 14U, true> = 797;  // BASE32_DIFF
template <> inline constexpr int32_t moduli<Backend::FP8, 15U, true> = 793;  // BASE32_SUM
template <> inline constexpr int32_t moduli<Backend::FP8, 16U, true> = 769;  // BASE32_SUM
template <> inline constexpr int32_t moduli<Backend::FP8, 17U, true> = 761;  // BASE32_SUM
template <> inline constexpr int32_t moduli<Backend::FP8, 18U, true> = 757;  // BASE32_SUM
template <> inline constexpr int32_t moduli<Backend::FP8, 19U, true> = 709;  // BASE32_SUM

//==========
// FP8 Karatsuba type
//==========
enum KaratsubaType : uint8_t {
    NONE,
    BASE32_SUM,
    BASE32_DIFF,
    BASE33_SUM,
    BASE31_DIFF,
    BASE49
};

template <unsigned IDX, bool COMPLEX = false> inline constexpr KaratsubaType kara_type = KaratsubaType::NONE;

template <> inline constexpr KaratsubaType kara_type<0u, false>  = KaratsubaType::NONE;        // 2401, base-49
template <> inline constexpr KaratsubaType kara_type<1u, false>  = KaratsubaType::NONE;        // 2209, base-47
template <> inline constexpr KaratsubaType kara_type<2u, false>  = KaratsubaType::NONE;        // 2025, base-45
template <> inline constexpr KaratsubaType kara_type<3u, false>  = KaratsubaType::NONE;        // 1849, base-43
template <> inline constexpr KaratsubaType kara_type<4u, false>  = KaratsubaType::NONE;        // 1681, base-41
template <> inline constexpr KaratsubaType kara_type<5u, false>  = KaratsubaType::NONE;        // 1369, base-37
template <> inline constexpr KaratsubaType kara_type<6u, false>  = KaratsubaType::BASE49;      // 1193
template <> inline constexpr KaratsubaType kara_type<7u, false>  = KaratsubaType::BASE49;      // 1097
template <> inline constexpr KaratsubaType kara_type<8u, false>  = KaratsubaType::BASE33_SUM;  // 1033
template <> inline constexpr KaratsubaType kara_type<9u, false>  = KaratsubaType::NONE;        // 1024, base-32
template <> inline constexpr KaratsubaType kara_type<10u, false> = KaratsubaType::BASE49;      // 1003
template <> inline constexpr KaratsubaType kara_type<11u, false> = KaratsubaType::BASE49;      // 997
template <> inline constexpr KaratsubaType kara_type<12u, false> = KaratsubaType::NONE;        // 961, base-31
template <> inline constexpr KaratsubaType kara_type<13u, false> = KaratsubaType::BASE32_DIFF; // 941
template <> inline constexpr KaratsubaType kara_type<14u, false> = KaratsubaType::BASE49;      // 937
template <> inline constexpr KaratsubaType kara_type<15u, false> = KaratsubaType::BASE32_SUM;  // 911
template <> inline constexpr KaratsubaType kara_type<16u, false> = KaratsubaType::BASE32_SUM;  // 907
template <> inline constexpr KaratsubaType kara_type<17u, false> = KaratsubaType::BASE49;      // 863
template <> inline constexpr KaratsubaType kara_type<18u, false> = KaratsubaType::BASE49;      // 859
template <> inline constexpr KaratsubaType kara_type<19u, false> = KaratsubaType::NONE;        // 841, base-29

template <> inline constexpr KaratsubaType kara_type<0u, true>  = KaratsubaType::NONE;        // 1681, base-41
template <> inline constexpr KaratsubaType kara_type<1u, true>  = KaratsubaType::NONE;        // 1369, base-37
template <> inline constexpr KaratsubaType kara_type<2u, true>  = KaratsubaType::BASE49;      // 1193
template <> inline constexpr KaratsubaType kara_type<3u, true>  = KaratsubaType::BASE49;      // 1097
template <> inline constexpr KaratsubaType kara_type<4u, true>  = KaratsubaType::BASE33_SUM;  // 1033
template <> inline constexpr KaratsubaType kara_type<5u, true>  = KaratsubaType::BASE49;      // 997
template <> inline constexpr KaratsubaType kara_type<6u, true>  = KaratsubaType::BASE32_DIFF; // 941
template <> inline constexpr KaratsubaType kara_type<7u, true>  = KaratsubaType::BASE49;      // 937
template <> inline constexpr KaratsubaType kara_type<8u, true>  = KaratsubaType::BASE32_SUM;  // 905
template <> inline constexpr KaratsubaType kara_type<9u, true>  = KaratsubaType::BASE32_SUM;  // 901
template <> inline constexpr KaratsubaType kara_type<10u, true> = KaratsubaType::BASE49;      // 857
template <> inline constexpr KaratsubaType kara_type<11u, true> = KaratsubaType::BASE32_SUM;  // 853
template <> inline constexpr KaratsubaType kara_type<12u, true> = KaratsubaType::NONE;        // 841, base-29
template <> inline constexpr KaratsubaType kara_type<13u, true> = KaratsubaType::BASE32_DIFF; // 809
template <> inline constexpr KaratsubaType kara_type<14u, true> = KaratsubaType::BASE32_DIFF; // 797
template <> inline constexpr KaratsubaType kara_type<15u, true> = KaratsubaType::BASE32_SUM;  // 793
template <> inline constexpr KaratsubaType kara_type<16u, true> = KaratsubaType::BASE32_SUM;  // 769
template <> inline constexpr KaratsubaType kara_type<17u, true> = KaratsubaType::BASE32_SUM;  // 761
template <> inline constexpr KaratsubaType kara_type<18u, true> = KaratsubaType::BASE32_SUM;  // 757
template <> inline constexpr KaratsubaType kara_type<19u, true> = KaratsubaType::BASE32_SUM;  // 709

inline constexpr int32_t moduli_int8[20] = {
    moduli<Backend::INT8, 0U, false>, moduli<Backend::INT8, 1U, false>,
    moduli<Backend::INT8, 2U, false>, moduli<Backend::INT8, 3U, false>,
    moduli<Backend::INT8, 4U, false>, moduli<Backend::INT8, 5U, false>,
    moduli<Backend::INT8, 6U, false>, moduli<Backend::INT8, 7U, false>,
    moduli<Backend::INT8, 8U, false>, moduli<Backend::INT8, 9U, false>,
    moduli<Backend::INT8, 10U, false>, moduli<Backend::INT8, 11U, false>,
    moduli<Backend::INT8, 12U, false>, moduli<Backend::INT8, 13U, false>,
    moduli<Backend::INT8, 14U, false>, moduli<Backend::INT8, 15U, false>,
    moduli<Backend::INT8, 16U, false>, moduli<Backend::INT8, 17U, false>,
    moduli<Backend::INT8, 18U, false>, moduli<Backend::INT8, 19U, false>};

inline constexpr int32_t moduli_int8_complex[20] = {
    moduli<Backend::INT8, 0U, true>, moduli<Backend::INT8, 1U, true>,
    moduli<Backend::INT8, 2U, true>, moduli<Backend::INT8, 3U, true>,
    moduli<Backend::INT8, 4U, true>, moduli<Backend::INT8, 5U, true>,
    moduli<Backend::INT8, 6U, true>, moduli<Backend::INT8, 7U, true>,
    moduli<Backend::INT8, 8U, true>, moduli<Backend::INT8, 9U, true>,
    moduli<Backend::INT8, 10U, true>, moduli<Backend::INT8, 11U, true>,
    moduli<Backend::INT8, 12U, true>, moduli<Backend::INT8, 13U, true>,
    moduli<Backend::INT8, 14U, true>, moduli<Backend::INT8, 15U, true>,
    moduli<Backend::INT8, 16U, true>, moduli<Backend::INT8, 17U, true>,
    moduli<Backend::INT8, 18U, true>, moduli<Backend::INT8, 19U, true>};

// 2^32 / moduli
template <Backend BACKEND, unsigned IDX, bool COMPLEX = false>
inline constexpr int32_t p_inv_32_v = int32_t(4294967296ULL / uint64_t(moduli<BACKEND, IDX, COMPLEX>));

// 2^64 / moduli
inline constexpr uint64_t UINT64_MAX_V = 18446744073709551615ULL;
template <Backend BACKEND, unsigned IDX, bool COMPLEX = false>
inline constexpr int64_t p_inv_64_v =
    int64_t(UINT64_MAX_V / uint64_t(moduli<BACKEND, IDX, COMPLEX>)) +
    int64_t(UINT64_MAX_V % uint64_t(moduli<BACKEND, IDX, COMPLEX>) == uint64_t(moduli<BACKEND, IDX, COMPLEX> - 1));

// FP8: sqrt(moduli)
template <unsigned IDX, bool COMPLEX = false> inline constexpr int32_t sqrt_moduli = 0;

template <> inline constexpr int32_t sqrt_moduli<0U, false>  = 49;
template <> inline constexpr int32_t sqrt_moduli<1U, false>  = 47;
template <> inline constexpr int32_t sqrt_moduli<2U, false>  = 45;
template <> inline constexpr int32_t sqrt_moduli<3U, false>  = 43;
template <> inline constexpr int32_t sqrt_moduli<4U, false>  = 41;
template <> inline constexpr int32_t sqrt_moduli<5U, false>  = 37;
template <> inline constexpr int32_t sqrt_moduli<9U, false>  = 32;
template <> inline constexpr int32_t sqrt_moduli<12U, false> = 31;
template <> inline constexpr int32_t sqrt_moduli<19U, false> = 29;

template <> inline constexpr int32_t sqrt_moduli<0U, true>  = 41;
template <> inline constexpr int32_t sqrt_moduli<1U, true>  = 37;
template <> inline constexpr int32_t sqrt_moduli<12U, true> = 29;

//==========
// number of matrices for workspace of A/B
//==========
template <Backend BACKEND, unsigned NUM_MODULI, bool COMPLEX = false>
__host__ __device__ constexpr unsigned num_mat_constexpr() {
    static_assert(NUM_MODULI <= 20U, "NUM_MODULI must be in [0, 20]");
    if constexpr (BACKEND == Backend::INT8) {
        return NUM_MODULI;
    } else {
        return fp8_plan::plane_count<NUM_MODULI, COMPLEX>;
    }
}
template <Backend BACKEND, unsigned NUM_MODULI, bool COMPLEX = false>
inline constexpr unsigned num_mat_v = num_mat_constexpr<BACKEND, NUM_MODULI, COMPLEX>();

inline constexpr unsigned num_mat_fp8[21] = {
    num_mat_v<Backend::FP8, 0U, false>,
    num_mat_v<Backend::FP8, 1U, false>,
    num_mat_v<Backend::FP8, 2U, false>,
    num_mat_v<Backend::FP8, 3U, false>,
    num_mat_v<Backend::FP8, 4U, false>,
    num_mat_v<Backend::FP8, 5U, false>,
    num_mat_v<Backend::FP8, 6U, false>,
    num_mat_v<Backend::FP8, 7U, false>,
    num_mat_v<Backend::FP8, 8U, false>,
    num_mat_v<Backend::FP8, 9U, false>,
    num_mat_v<Backend::FP8, 10U, false>,
    num_mat_v<Backend::FP8, 11U, false>,
    num_mat_v<Backend::FP8, 12U, false>,
    num_mat_v<Backend::FP8, 13U, false>,
    num_mat_v<Backend::FP8, 14U, false>,
    num_mat_v<Backend::FP8, 15U, false>,
    num_mat_v<Backend::FP8, 16U, false>,
    num_mat_v<Backend::FP8, 17U, false>,
    num_mat_v<Backend::FP8, 18U, false>,
    num_mat_v<Backend::FP8, 19U, false>,
    num_mat_v<Backend::FP8, 20U, false>};

inline constexpr unsigned num_mat_fp8_complex[21] = {
    num_mat_v<Backend::FP8, 0U, true>,
    num_mat_v<Backend::FP8, 1U, true>,
    num_mat_v<Backend::FP8, 2U, true>,
    num_mat_v<Backend::FP8, 3U, true>,
    num_mat_v<Backend::FP8, 4U, true>,
    num_mat_v<Backend::FP8, 5U, true>,
    num_mat_v<Backend::FP8, 6U, true>,
    num_mat_v<Backend::FP8, 7U, true>,
    num_mat_v<Backend::FP8, 8U, true>,
    num_mat_v<Backend::FP8, 9U, true>,
    num_mat_v<Backend::FP8, 10U, true>,
    num_mat_v<Backend::FP8, 11U, true>,
    num_mat_v<Backend::FP8, 12U, true>,
    num_mat_v<Backend::FP8, 13U, true>,
    num_mat_v<Backend::FP8, 14U, true>,
    num_mat_v<Backend::FP8, 15U, true>,
    num_mat_v<Backend::FP8, 16U, true>,
    num_mat_v<Backend::FP8, 17U, true>,
    num_mat_v<Backend::FP8, 18U, true>,
    num_mat_v<Backend::FP8, 19U, true>,
    num_mat_v<Backend::FP8, 20U, true>};

template <Backend BACKEND, bool COMPLEX = false>
inline unsigned num_mat(unsigned NUM_MODULI) {
    assert(NUM_MODULI <= 20U);
    if constexpr (BACKEND == Backend::INT8) {
        return NUM_MODULI;
    } else if constexpr (COMPLEX) {
        return num_mat_fp8_complex[NUM_MODULI];
    } else {
        return num_mat_fp8[NUM_MODULI];
    }
}

//==========
// log2P = round-down( log2(prod(moduli)-1)/2 - 0.5 ) in float
//==========
template <Backend BACKEND, unsigned NUM_MODULI, bool COMPLEX = false> inline constexpr float log2P = 0.0F;

// INT8 real
template <> inline constexpr float log2P<Backend::INT8, 2U, false>  = 0x1.dfd1ec0000000p+2F;
template <> inline constexpr float log2P<Backend::INT8, 3U, false>  = 0x1.6fa3360000000p+3F;
template <> inline constexpr float log2P<Backend::INT8, 4U, false>  = 0x1.ef2ea60000000p+3F;
template <> inline constexpr float log2P<Backend::INT8, 5U, false>  = 0x1.372d940000000p+4F;
template <> inline constexpr float log2P<Backend::INT8, 6U, false>  = 0x1.767b2e0000000p+4F;
template <> inline constexpr float log2P<Backend::INT8, 7U, false>  = 0x1.b5b0280000000p+4F;
template <> inline constexpr float log2P<Backend::INT8, 8U, false>  = 0x1.f49a020000000p+4F;
template <> inline constexpr float log2P<Backend::INT8, 9U, false>  = 0x1.19a8580000000p+5F;
template <> inline constexpr float log2P<Backend::INT8, 10U, false> = 0x1.38f6bc0000000p+5F;
template <> inline constexpr float log2P<Backend::INT8, 11U, false> = 0x1.582ada0000000p+5F;
template <> inline constexpr float log2P<Backend::INT8, 12U, false> = 0x1.7736ae0000000p+5F;
template <> inline constexpr float log2P<Backend::INT8, 13U, false> = 0x1.9619160000000p+5F;
template <> inline constexpr float log2P<Backend::INT8, 14U, false> = 0x1.b4a4fe0000000p+5F;
template <> inline constexpr float log2P<Backend::INT8, 15U, false> = 0x1.d321f80000000p+5F;
template <> inline constexpr float log2P<Backend::INT8, 16U, false> = 0x1.f180a60000000p+5F;
template <> inline constexpr float log2P<Backend::INT8, 17U, false> = 0x1.07e7f80000000p+6F;
template <> inline constexpr float log2P<Backend::INT8, 18U, false> = 0x1.16e7e20000000p+6F;
template <> inline constexpr float log2P<Backend::INT8, 19U, false> = 0x1.25df9a0000000p+6F;
template <> inline constexpr float log2P<Backend::INT8, 20U, false> = 0x1.34be220000000p+6F;

// INT8 complex
template <> inline constexpr float log2P<Backend::INT8, 2U, true>  = 0x1.d8de020000000p+2F;
template <> inline constexpr float log2P<Backend::INT8, 3U, true>  = 0x1.69dc460000000p+3F;
template <> inline constexpr float log2P<Backend::INT8, 4U, true>  = 0x1.e677860000000p+3F;
template <> inline constexpr float log2P<Backend::INT8, 5U, true>  = 0x1.30ab560000000p+4F;
template <> inline constexpr float log2P<Backend::INT8, 6U, true>  = 0x1.6da54c0000000p+4F;
template <> inline constexpr float log2P<Backend::INT8, 7U, true>  = 0x1.aa62a60000000p+4F;
template <> inline constexpr float log2P<Backend::INT8, 8U, true>  = 0x1.e662560000000p+4F;
template <> inline constexpr float log2P<Backend::INT8, 9U, true>  = 0x1.10ee3a0000000p+5F;
template <> inline constexpr float log2P<Backend::INT8, 10U, true> = 0x1.2e1bea0000000p+5F;
template <> inline constexpr float log2P<Backend::INT8, 11U, true> = 0x1.4afc580000000p+5F;
template <> inline constexpr float log2P<Backend::INT8, 12U, true> = 0x1.6760ba0000000p+5F;
template <> inline constexpr float log2P<Backend::INT8, 13U, true> = 0x1.82a8980000000p+5F;
template <> inline constexpr float log2P<Backend::INT8, 14U, true> = 0x1.9dbb360000000p+5F;
template <> inline constexpr float log2P<Backend::INT8, 15U, true> = 0x1.b85d380000000p+5F;
template <> inline constexpr float log2P<Backend::INT8, 16U, true> = 0x1.d2c3880000000p+5F;
template <> inline constexpr float log2P<Backend::INT8, 17U, true> = 0x1.ecaab00000000p+5F;
template <> inline constexpr float log2P<Backend::INT8, 18U, true> = 0x1.02b6880000000p+6F;
template <> inline constexpr float log2P<Backend::INT8, 19U, true> = 0x1.0e93120000000p+6F;
template <> inline constexpr float log2P<Backend::INT8, 20U, true> = 0x1.1a07c40000000p+6F;

// FP8 real
template <> inline constexpr float log2P<Backend::FP8, 2U, false>  = 0x1.a058e80000000p+3F;
template <> inline constexpr float log2P<Backend::FP8, 3U, false>  = 0x1.2395700000000p+4F;
template <> inline constexpr float log2P<Backend::FP8, 4U, false>  = 0x1.7214c80000000p+4F;
template <> inline constexpr float log2P<Backend::FP8, 5U, false>  = 0x1.bb6ba60000000p+4F;
template <> inline constexpr float log2P<Backend::FP8, 6U, false>  = 0x1.0702f00000000p+5F;
template <> inline constexpr float log2P<Backend::FP8, 7U, false>  = 0x1.304b700000000p+5F;
template <> inline constexpr float log2P<Backend::FP8, 8U, false>  = 0x1.5970e00000000p+5F;
template <> inline constexpr float log2P<Backend::FP8, 9U, false>  = 0x1.7ac8500000000p+5F;
template <> inline constexpr float log2P<Backend::FP8, 10U, false> = 0x1.a3ceac0000000p+5F;
template <> inline constexpr float log2P<Backend::FP8, 11U, false> = 0x1.ccba360000000p+5F;
template <> inline constexpr float log2P<Backend::FP8, 12U, false> = 0x1.f59be00000000p+5F;
template <> inline constexpr float log2P<Backend::FP8, 13U, false> = 0x1.0f374e0000000p+6F;
template <> inline constexpr float log2P<Backend::FP8, 14U, false> = 0x1.239e2a0000000p+6F;
template <> inline constexpr float log2P<Backend::FP8, 15U, false> = 0x1.3801400000000p+6F;
template <> inline constexpr float log2P<Backend::FP8, 16U, false> = 0x1.4c5f440000000p+6F;
template <> inline constexpr float log2P<Backend::FP8, 17U, false> = 0x1.60b6ea0000000p+6F;
template <> inline constexpr float log2P<Backend::FP8, 18U, false> = 0x1.750d460000000p+6F;
template <> inline constexpr float log2P<Backend::FP8, 19U, false> = 0x1.8955600000000p+6F;
template <> inline constexpr float log2P<Backend::FP8, 20U, false> = 0x1.9674ea0000000p+6F;

// FP8 complex
template <> inline constexpr float log2P<Backend::FP8, 2U, true>  = 0x1.72803a0000000p+3F;
template <> inline constexpr float log2P<Backend::FP8, 3U, true>  = 0x1.0c1ea00000000p+4F;
template <> inline constexpr float log2P<Backend::FP8, 4U, true>  = 0x1.5e69800000000p+4F;
template <> inline constexpr float log2P<Backend::FP8, 5U, true>  = 0x1.b040940000000p+4F;
template <> inline constexpr float log2P<Backend::FP8, 6U, true>  = 0x1.0101f40000000p+5F;
template <> inline constexpr float log2P<Backend::FP8, 7U, true>  = 0x1.197b200000000p+5F;
template <> inline constexpr float log2P<Backend::FP8, 8U, true>  = 0x1.422a6a0000000p+5F;
template <> inline constexpr float log2P<Backend::FP8, 9U, true>  = 0x1.6aba9e0000000p+5F;
template <> inline constexpr float log2P<Backend::FP8, 10U, true> = 0x1.933b0c0000000p+5F;
template <> inline constexpr float log2P<Backend::FP8, 11U, true> = 0x1.bbb0da0000000p+5F;
template <> inline constexpr float log2P<Backend::FP8, 12U, true> = 0x1.e416940000000p+5F;
template <> inline constexpr float log2P<Backend::FP8, 13U, true> = 0x1.063b740000000p+6F;
template <> inline constexpr float log2P<Backend::FP8, 14U, true> = 0x1.1a5b3a0000000p+6F;
template <> inline constexpr float log2P<Backend::FP8, 15U, true> = 0x1.2951380000000p+6F;
template <> inline constexpr float log2P<Backend::FP8, 16U, true> = 0x1.3d630a0000000p+6F;
template <> inline constexpr float log2P<Backend::FP8, 17U, true> = 0x1.5169800000000p+6F;
template <> inline constexpr float log2P<Backend::FP8, 18U, true> = 0x1.6567560000000p+6F;
template <> inline constexpr float log2P<Backend::FP8, 19U, true> = 0x1.795f5a0000000p+6F;
template <> inline constexpr float log2P<Backend::FP8, 20U, true> = 0x1.8d54740000000p+6F;

} // namespace gemmul8::common::table
