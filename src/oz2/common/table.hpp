#pragma once
#include "common.hpp"

namespace gemmul8::common::table {

//==========
// moduli
//==========
template <Backend BACKEND, unsigned IDX, bool COMPLEX = false> inline constexpr int32_t moduli = 0;

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

template <unsigned IDX, bool COMPLEX = false>
inline constexpr unsigned num_limbs = kara_type<IDX, COMPLEX> == KaratsubaType::NONE ? 2U : 3U;

inline constexpr bool isKaratsuba[20] = {
    num_limbs<0U, false> == 3U,
    num_limbs<1U, false> == 3U,
    num_limbs<2U, false> == 3U,
    num_limbs<3U, false> == 3U,
    num_limbs<4U, false> == 3U,
    num_limbs<5U, false> == 3U,
    num_limbs<6U, false> == 3U,
    num_limbs<7U, false> == 3U,
    num_limbs<8U, false> == 3U,
    num_limbs<9U, false> == 3U,
    num_limbs<10U, false> == 3U,
    num_limbs<11U, false> == 3U,
    num_limbs<12U, false> == 3U,
    num_limbs<13U, false> == 3U,
    num_limbs<14U, false> == 3U,
    num_limbs<15U, false> == 3U,
    num_limbs<16U, false> == 3U,
    num_limbs<17U, false> == 3U,
    num_limbs<18U, false> == 3U,
    num_limbs<19U, false> == 3U};

inline constexpr bool isKaratsuba_complex[20] = {
    num_limbs<0U, true> == 3U,
    num_limbs<1U, true> == 3U,
    num_limbs<2U, true> == 3U,
    num_limbs<3U, true> == 3U,
    num_limbs<4U, true> == 3U,
    num_limbs<5U, true> == 3U,
    num_limbs<6U, true> == 3U,
    num_limbs<7U, true> == 3U,
    num_limbs<8U, true> == 3U,
    num_limbs<9U, true> == 3U,
    num_limbs<10U, true> == 3U,
    num_limbs<11U, true> == 3U,
    num_limbs<12U, true> == 3U,
    num_limbs<13U, true> == 3U,
    num_limbs<14U, true> == 3U,
    num_limbs<15U, true> == 3U,
    num_limbs<16U, true> == 3U,
    num_limbs<17U, true> == 3U,
    num_limbs<18U, true> == 3U,
    num_limbs<19U, true> == 3U};

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

inline constexpr int32_t moduli_fp8[20] = {
    moduli<Backend::FP8, 0U, false>, moduli<Backend::FP8, 1U, false>,
    moduli<Backend::FP8, 2U, false>, moduli<Backend::FP8, 3U, false>,
    moduli<Backend::FP8, 4U, false>, moduli<Backend::FP8, 5U, false>,
    moduli<Backend::FP8, 6U, false>, moduli<Backend::FP8, 7U, false>,
    moduli<Backend::FP8, 8U, false>, moduli<Backend::FP8, 9U, false>,
    moduli<Backend::FP8, 10U, false>, moduli<Backend::FP8, 11U, false>,
    moduli<Backend::FP8, 12U, false>, moduli<Backend::FP8, 13U, false>,
    moduli<Backend::FP8, 14U, false>, moduli<Backend::FP8, 15U, false>,
    moduli<Backend::FP8, 16U, false>, moduli<Backend::FP8, 17U, false>,
    moduli<Backend::FP8, 18U, false>, moduli<Backend::FP8, 19U, false>};

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

inline constexpr int32_t moduli_fp8_complex[20] = {
    moduli<Backend::FP8, 0U, true>, moduli<Backend::FP8, 1U, true>,
    moduli<Backend::FP8, 2U, true>, moduli<Backend::FP8, 3U, true>,
    moduli<Backend::FP8, 4U, true>, moduli<Backend::FP8, 5U, true>,
    moduli<Backend::FP8, 6U, true>, moduli<Backend::FP8, 7U, true>,
    moduli<Backend::FP8, 8U, true>, moduli<Backend::FP8, 9U, true>,
    moduli<Backend::FP8, 10U, true>, moduli<Backend::FP8, 11U, true>,
    moduli<Backend::FP8, 12U, true>, moduli<Backend::FP8, 13U, true>,
    moduli<Backend::FP8, 14U, true>, moduli<Backend::FP8, 15U, true>,
    moduli<Backend::FP8, 16U, true>, moduli<Backend::FP8, 17U, true>,
    moduli<Backend::FP8, 18U, true>, moduli<Backend::FP8, 19U, true>};

inline constexpr unsigned k_block_first_fp8[20] = {
    16384, // 2401, B=32
    18432, // 2209, B=30
    21248, // 2025, B=28
    24576, // 1849, B=26
    28928, // 1681, B=24
    41728, // 1369, B=20
    16384, // 1193, BASE49,      B=32
    16384, // 1097, BASE49,      B=32
    18432, // 1033, BASE33_SUM,  B=30
    65536, // 1024, base-32,     B=16
    16384, // 1003, BASE49,      B=32
    28928, // 997,  BASE49,      B=24
    74496, // 961,  base-31,     B=15
    16384, // 941,  BASE32_DIFF, B=32
    21248, // 937,  BASE49,      B=28
    24576, // 911,  BASE32_SUM,  B=26
    24576, // 907,  BASE32_SUM,  B=26
    16384, // 863,  BASE49,      B=32
    16384, // 859,  BASE49,      B=32
    85504  // 841,  base-29,     B=14
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
    18432, // 1033
    65280, // 1024
    16128, // 1003
    28928, // 997
    74496, // 961
    16128, // 941
    21248, // 937
    24576, // 911
    24576, // 907
    16128, // 863
    16128, // 859
    85504  // 841
};

inline constexpr unsigned k_block_first_fp8_complex[20] = {
    28928, // 1681, base-41,     B=24
    41728, // 1369, base-37,     B=20
    16384, // 1193, BASE49,      B=32
    16384, // 1097, BASE49,      B=32
    18432, // 1033, BASE33_SUM,  B=30
    28928, // 997,  BASE49,      B=24
    16384, // 941,  BASE32_DIFF, B=32
    21248, // 937,  BASE49,      B=28
    24576, // 905,  BASE32_SUM,  B=26
    24576, // 901,  BASE32_SUM,  B=26
    16384, // 857,  BASE49,      B=32
    18432, // 853,  BASE32_SUM,  B=30
    85504, // 841,  base-29,     B=14
    16384, // 809,  BASE32_DIFF, B=32
    16384, // 797,  BASE32_DIFF, B=32
    16384, // 793,  BASE32_SUM,  B=32
    34560, // 769,  BASE32_SUM,  B=22
    34560, // 761,  BASE32_SUM,  B=22
    24576, // 757,  BASE32_SUM,  B=26
    41728  // 709,  BASE32_SUM,  B=20
};

inline constexpr unsigned k_block_next_fp8_complex[20] = {
    28928, // 1681
    41728, // 1369
    16128, // 1193
    16128, // 1097
    18432, // 1033
    28928, // 997
    16128, // 941
    21248, // 937
    24576, // 905
    24576, // 901
    16128, // 857
    18432, // 853
    85504, // 841
    16128, // 809
    16128, // 797
    16128, // 793
    34560, // 769
    34560, // 761
    24576, // 757
    41728  // 709
};

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
    } else if constexpr (NUM_MODULI == 0U) {
        return 0U;
    } else {
        return num_mat_constexpr<BACKEND, NUM_MODULI - 1U, COMPLEX>() +
               num_limbs<NUM_MODULI - 1U, COMPLEX>;
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

constexpr double2 P_complex[19] = {
    { -0x1.b6b2000000000p+15,   0x0.0000000000000p+0},
    { -0x1.886d3a0000000p+23,   0x0.0000000000000p+0},
    { -0x1.52c64b1200000p+31,   0x0.0000000000000p+0},
    { -0x1.0f48ca1d6a000p+39,   0x0.0000000000000p+0},
    { -0x1.a186071145240p+46,   0x0.0000000000000p+0},
    { -0x1.3ac60b5405202p+54,  -0x1.0000000000000p+0},
    { -0x1.bd1c0c04cf3f7p+61,  -0x1.7400000000000p+6},
    { -0x1.2ccbf41f400dep+69, -0x1.4d90000000000p+12},
    { -0x1.70f2296e54910p+76, -0x1.6324540000000p+22},
    { -0x1.ad79e43a6e70dp+83, -0x1.d6849c8000000p+25},
    { -0x1.cbac76468a34cp+90,  0x1.c433083f80000p+33},
    { -0x1.95ce406a46029p+97, -0x1.70caf2b7f1000p+40},
    {-0x1.5991a2da7f9e3p+104,  0x1.85f32d4f5cc60p+47},
    {-0x1.10acea8068b2dp+111, -0x1.46270f1fb065ep+55},
    {-0x1.9d46136a9eaf0p+117, -0x1.4b94ccbb01d6ap+63},
    {-0x1.1f5ab9802255bp+124, -0x1.9a35d9681d1cep+68},
    {-0x1.47c37b962729cp+130,  0x1.6c1a94053ecb1p+74},
    {-0x1.386651cb1d53dp+136,  0x1.96c255453ff66p+82},
    {-0x1.02b4bbbc34496p+142, -0x1.7e4e22c54e0ffp+87}
};
} // namespace INT8

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
    {-0x1.12f4c8531f63ap+204, -0x1.72c14309f2704p+149}
};

constexpr double2 P_complex[19] = {
    { -0x1.18eb480000000p+21,    0x0.0000000000000p+0},
    { -0x1.47481ca200000p+31,    0x0.0000000000000p+0},
    { -0x1.5e9d00ac8c800p+41,    0x0.0000000000000p+0},
    { -0x1.61b1e1ee10bc2p+51,    0x0.0000000000000p+0},
    { -0x1.585e713909cb3p+61,    0x1.7600000000000p+7},
    { -0x1.3c74c98baa3ffp+71,  -0x1.3ea1000000000p+16},
    { -0x1.2191dd6c0c890p+81,  -0x1.d18ed24000000p+26},
    { -0x1.ffd653e17c283p+90,   0x1.122e1abbc0000p+34},
    {-0x1.c25b554e267e6p+100,   0x1.2d4fc4416cac0p+46},
    {-0x1.78e8f024a7b74p+110,   0x1.8b007ff04cbccp+50},
    {-0x1.39f80a0a88b56p+120,  -0x1.7edbd755b4505p+66},
    {-0x1.01dbf63f26c70p+130,   0x1.3f8ff265dca91p+76},
    {-0x1.97700e96c8c56p+139,  -0x1.e7110d7e0cd6ep+85},
    {-0x1.3d1df75adbc3ap+149,  -0x1.391805c05b7e4p+95},
    {-0x1.eb28e99c39608p+158, -0x1.38edb4e86db51p+104},
    {-0x1.70d9796f9216cp+168,  0x1.0cff7ce4739cdp+114},
    {-0x1.121d9e7f2a516p+178, -0x1.4417216f3914cp+124},
    {-0x1.9548cad704115p+187, -0x1.fb2c32eef1e53p+133},
    {-0x1.189ca6715f910p+197, -0x1.1050b487e1f7dp+142}
};
} // namespace FP8

template <Backend BACKEND, typename doublex_t, bool COMPLEX = false> __forceinline__ doublex_t get_P(unsigned NUM_MODULI);
template <> __forceinline__ double get_P<Backend::INT8, double, false>(unsigned NUM_MODULI) { return INT8::P[NUM_MODULI - 2].x; }
template <> __forceinline__ double2 get_P<Backend::INT8, double2, false>(unsigned NUM_MODULI) { return INT8::P[NUM_MODULI - 2]; }
template <> __forceinline__ double get_P<Backend::FP8, double, false>(unsigned NUM_MODULI) { return FP8::P[NUM_MODULI - 2].x; }
template <> __forceinline__ double2 get_P<Backend::FP8, double2, false>(unsigned NUM_MODULI) { return FP8::P[NUM_MODULI - 2]; }
template <> __forceinline__ double get_P<Backend::INT8, double, true>(unsigned NUM_MODULI) { return INT8::P_complex[NUM_MODULI - 2].x; }
template <> __forceinline__ double2 get_P<Backend::INT8, double2, true>(unsigned NUM_MODULI) { return INT8::P_complex[NUM_MODULI - 2]; }
template <> __forceinline__ double get_P<Backend::FP8, double, true>(unsigned NUM_MODULI) { return FP8::P_complex[NUM_MODULI - 2].x; }
template <> __forceinline__ double2 get_P<Backend::FP8, double2, true>(unsigned NUM_MODULI) { return FP8::P_complex[NUM_MODULI - 2]; }

//==========
// invP[i] = 1/P[i] in double
//==========
namespace INT8 {
constexpr double invP[19] = {
    0x1.0101010101010p-16, 0x1.040d287a7051fp-24, 0x1.093b510fbf0d4p-32, 0x1.12e5617d255d8p-40, 0x1.2401777d7fdb6p-48,
    0x1.38c6a8b145786p-56, 0x1.57a6a12c3f24ap-64, 0x1.802b2f252aa3fp-72, 0x1.b13f5ca3b64a6p-80, 0x1.f15c410568cccp-88,
    0x1.255fb5199b040p-95, 0x1.63f115f5d0b39p-103, 0x1.c9e518641aa18p-111, 0x1.2983f5dbae8acp-118, 0x1.8aa1c572fa163p-126,
    0x1.0877227a9f8e3p-133, 0x1.760ceb764616fp-141, 0x1.0b7a38d26e2fep-148, 0x1.8bce042d07acep-156};

constexpr double invP_complex[19] = {
    0x1.2ac6df1b65161p-16, 0x1.4e00f969863eap-24, 0x1.82e67c6c5600cp-32, 0x1.e32751b2fec14p-40, 0x1.39ed5d5e35c51p-47,
    0x1.a0669e552b777p-55, 0x1.2678a88b6f5b8p-62, 0x1.b3bfdd4aa1d2fp-70, 0x1.6342be3310298p-77, 0x1.3130bee1fcb40p-84,
    0x1.1d243346fec9dp-91, 0x1.42fdf61e94264p-98, 0x1.7b4b1c52e1a4cp-105, 0x1.e0b04c732d2c0p-112, 0x1.3d27c8dfc959cp-118,
    0x1.c8223508527b5p-125, 0x1.8fe5e157f4260p-131, 0x1.a390a9474bb25p-137, 0x1.faa5065fc0d75p-143};
} // namespace INT8

namespace FP8 {
constexpr double invP[19] = {
    0x1.94e504ced568fp-23, 0x1.997e4ff4b4f12p-34, 0x1.c590c0ad928ecp-45, 0x1.144b62a3e8951p-55, 0x1.9d54e995463f5p-66,
    0x1.62c77d2ca796bp-76, 0x1.4b2ba0f34f317p-86, 0x1.4848fcbaab304p-96, 0x1.4848fcbaab304p-106, 0x1.4f2891b7af89ep-116,
    0x1.583c27c41b41dp-126, 0x1.6ecd4901fa69dp-136, 0x1.8f27c1fb144f4p-146, 0x1.b437787735118p-156, 0x1.ea532555f0e70p-166,
    0x1.14c99bb573f30p-175, 0x1.486cb2d326ccbp-185, 0x1.878278c9a4911p-195, 0x1.dcb38fb8f1c64p-205};

constexpr double invP_complex[19] = {
    0x1.d2953134e8c50p-22, 0x1.907c9fd42591bp-32, 0x1.75d61c4f7b0c1p-42, 0x1.72944e9e963a1p-52, 0x1.7c9d75777cc44p-62,
    0x1.9e2fd61165c3ep-72, 0x1.c4a4d75d998fdp-82, 0x1.0014d7c18aa49p-91, 0x1.230a558793a00p-101, 0x1.5bc11241cf28fp-111,
    0x1.a177d5177af0ap-121, 0x1.fc4ef096596b9p-131, 0x1.41b2c0226586ep-140, 0x1.9d52df16ac333p-150, 0x1.0adcbaa65dbf9p-159,
    0x1.635a855c08fc7p-169, 0x1.de29bde06cf98p-179, 0x1.4368526e889b4p-188, 0x1.d317ef4aca5aap-198};
} // namespace FP8

template <Backend BACKEND, bool COMPLEX = false>
__forceinline__ double get_invP(unsigned NUM_MODULI) {
    if constexpr (BACKEND == Backend::INT8) {
        if constexpr (!COMPLEX) return INT8::invP[NUM_MODULI - 2];
        else return INT8::invP_complex[NUM_MODULI - 2];
    } else {
        if constexpr (!COMPLEX) return FP8::invP[NUM_MODULI - 2];
        else return FP8::invP_complex[NUM_MODULI - 2];
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
template <> inline constexpr float log2P<Backend::FP8, 2U, false>  = 0x1.556ae40000000p+3F;
template <> inline constexpr float log2P<Backend::FP8, 3U, false>  = 0x1.0294120000000p+4F;
template <> inline constexpr float log2P<Backend::FP8, 4U, false>  = 0x1.59660e0000000p+4F;
template <> inline constexpr float log2P<Backend::FP8, 5U, false>  = 0x1.af1e960000000p+4F;
template <> inline constexpr float log2P<Backend::FP8, 6U, false>  = 0x1.013c400000000p+5F;
template <> inline constexpr float log2P<Backend::FP8, 7U, false>  = 0x1.2a1dec0000000p+5F;
template <> inline constexpr float log2P<Backend::FP8, 8U, false>  = 0x1.5283a60000000p+5F;
template <> inline constexpr float log2P<Backend::FP8, 9U, false>  = 0x1.7a90940000000p+5F;
template <> inline constexpr float log2P<Backend::FP8, 10U, false> = 0x1.a290940000000p+5F;
template <> inline constexpr float log2P<Backend::FP8, 11U, false> = 0x1.ca71f80000000p+5F;
template <> inline constexpr float log2P<Backend::FP8, 12U, false> = 0x1.f24a7e0000000p+5F;
template <> inline constexpr float log2P<Backend::FP8, 13U, false> = 0x1.0cf6580000000p+6F;
template <> inline constexpr float log2P<Backend::FP8, 14U, false> = 0x1.20b7e80000000p+6F;
template <> inline constexpr float log2P<Backend::FP8, 15U, false> = 0x1.3476520000000p+6F;
template <> inline constexpr float log2P<Backend::FP8, 16U, false> = 0x1.481ff20000000p+6F;
template <> inline constexpr float log2P<Backend::FP8, 17U, false> = 0x1.5bc6540000000p+6F;
template <> inline constexpr float log2P<Backend::FP8, 18U, false> = 0x1.6f47fa0000000p+6F;
template <> inline constexpr float log2P<Backend::FP8, 19U, false> = 0x1.82c6300000000p+6F;
template <> inline constexpr float log2P<Backend::FP8, 20U, false> = 0x1.9634c40000000p+6F;

// FP8 complex
template <> inline constexpr float log2P<Backend::FP8, 2U, true>  = 0x1.4224e80000000p+3F;
template <> inline constexpr float log2P<Backend::FP8, 3U, true>  = 0x1.e5ab920000000p+3F;
template <> inline constexpr float log2P<Backend::FP8, 4U, true>  = 0x1.43a1400000000p+4F;
template <> inline constexpr float log2P<Backend::FP8, 5U, true>  = 0x1.93bb1a0000000p+4F;
template <> inline constexpr float log2P<Backend::FP8, 6U, true>  = 0x1.e36c280000000p+4F;
template <> inline constexpr float log2P<Backend::FP8, 7U, true>  = 0x1.1939320000000p+5F;
template <> inline constexpr float log2P<Backend::FP8, 8U, true>  = 0x1.40b6080000000p+5F;
template <> inline constexpr float log2P<Backend::FP8, 9U, true>  = 0x1.67ff860000000p+5F;
template <> inline constexpr float log2P<Backend::FP8, 10U, true> = 0x1.8f427a0000000p+5F;
template <> inline constexpr float log2P<Backend::FP8, 11U, true> = 0x1.b63b780000000p+5F;
template <> inline constexpr float log2P<Backend::FP8, 12U, true> = 0x1.dd2d8a0000000p+5F;
template <> inline constexpr float log2P<Backend::FP8, 13U, true> = 0x1.0205580000000p+6F;
template <> inline constexpr float log2P<Backend::FP8, 14U, true> = 0x1.1557420000000p+6F;
template <> inline constexpr float log2P<Backend::FP8, 15U, true> = 0x1.289e240000000p+6F;
template <> inline constexpr float log2P<Backend::FP8, 16U, true> = 0x1.3be14e0000000p+6F;
template <> inline constexpr float log2P<Backend::FP8, 17U, true> = 0x1.4f0dc40000000p+6F;
template <> inline constexpr float log2P<Backend::FP8, 18U, true> = 0x1.6232800000000p+6F;
template <> inline constexpr float log2P<Backend::FP8, 19U, true> = 0x1.7553580000000p+6F;
template <> inline constexpr float log2P<Backend::FP8, 20U, true> = 0x1.8843ce0000000p+6F;

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

// Row index = NUM_MODULI - qPi_2_first_num_moduli
// qPi_2[idx][i].x = first (53-ceil(log2(rho))) bits of q[i]*P[i]/p[i] for rho = sum(floor(p[:]/2)),
// qPi_2[idx][i].y = double(q[i]*P[i]/p[i] - qPi_2[idx][i].x)
inline constexpr unsigned qPi_2_first_num_moduli = 7U;
inline constexpr double2 qPi_2[][20]             = {
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

// complex
inline constexpr double qPi_1_complex[19][20] = {
    {0x1.b4e0000000000p+12, 0x1.8018000000000p+15},
    {0x1.8048f80000000p+22, 0x1.3287dc0000000p+23, 0x1.1e2e1e0000000p+23},
    {0x1.289a7d3600000p+31, 0x1.0442a32600000p+31, 0x1.80a288a800000p+29, 0x1.886d3a0000000p+27},
    {0x1.713783acd0000p+36, 0x1.aecb7c4914000p+38, 0x1.aa7922e7d0000p+36, 0x1.0e0e8aa5f8000p+39, 0x1.c9e0017654000p+38},
    {0x1.d73ada4cda400p+45, 0x1.1b20a766e7780p+45, 0x1.cb754e348ab00p+44, 0x1.4a9e41d4ea4c0p+46, 0x1.04b297c4da000p+41,
     0x1.0f48ca1d6a000p+40},
    {0x1.029c55c9cf196p+54, 0x1.64a742ac1bc82p+53, 0x1.d359768c65839p+53, 0x1.172a6daaa908dp+52, 0x1.a1a69d643dbf4p+52,
     0x1.392cff9344c65p+53, 0x1.01b0b85ca8ac4p+54},
    {0x1.fdc084bc36a40p+59, 0x1.912bf4dbaeb44p+60, 0x1.75313a2591d09p+56, 0x1.403cb8b82a893p+61, 0x1.6cc5d067d80c1p+61,
     0x1.085aba77cf8dap+61, 0x1.5e8d778f0ffcfp+58, 0x1.3ac60b5405202p+54},
    {0x1.f5be6dc7c742bp+68, 0x1.bead37f9aaee2p+68, 0x1.ede27068f8438p+68, 0x1.49611244f887bp+68, 0x1.0dfbbe647a6dep+67,
     0x1.a56bb62e6189bp+67, 0x1.5d1c76db50f99p+66, 0x1.7794953c3c2dap+68, 0x1.abb8f38c9f270p+68},
    {0x1.601b289a28577p+76, 0x1.37f100f833b94p+76, 0x1.891cda89adf2dp+74, 0x1.1f24b9265fef2p+74, 0x1.4ec03e910f89dp+75,
     0x1.33248fa5eb6fap+75, 0x1.44fa716d1eb77p+75, 0x1.0d10c7f048fcep+74, 0x1.2209d694e24b9p+75, 0x1.3ae58390b70e8p+76},
    {0x1.a490db3e87c48p+82, 0x1.edfd5e68909c2p+81, 0x1.77169cd87b48dp+80, 0x1.1010fb262c7a7p+82, 0x1.0c29245b6e2afp+77,
     0x1.c9d134437c39ep+82, 0x1.640ada62d83dbp+80, 0x1.82c40bc378ac5p+83, 0x1.022ea3d62dacep+81, 0x1.144985c70431cp+83,
     0x1.88014c0539da1p+82},
    {0x1.d906735665dafp+89, 0x1.9a5a3b7069cc3p+87, 0x1.a1852bb573789p+88, 0x1.2b8410d4c49f8p+88, 0x1.ae860d5ec9df0p+85,
     0x1.54abff4a5ea14p+90, 0x1.869a81af56ee8p+88, 0x1.2698f58650168p+89, 0x1.16fe214fceb40p+90, 0x1.1c0090d1e6827p+90,
     0x1.0f7c313c8f7a3p+89, 0x1.2743cce82bed9p+88},
    {0x1.792e2cfffd169p+96, 0x1.14ec55b54b3e3p+97, 0x1.a235d05ba2b34p+96, 0x1.0deccf55bf5d3p+97, 0x1.22fe0b3979620p+97,
     0x1.4159414da87edp+96, 0x1.200f172cee0c7p+97, 0x1.623d3c790e709p+96, 0x1.2c3fd0f45d3bap+97, 0x1.c1bf14d79d7d2p+96,
     0x1.620f02f07aef5p+96, 0x1.d4026bec8c898p+96, 0x1.58c158b4e7a79p+94},
    {0x1.4ca9f0b058dc4p+103, 0x1.ab24151b428c9p+102, 0x1.eef69b011038ap+101, 0x1.80a90fa0f89fep+103, 0x1.c0658ea122785p+103,
     0x1.a4ff64f03383ap+103, 0x1.d1888130f6334p+102, 0x1.9892bae3381f7p+103, 0x1.5592e97abf40ap+104, 0x1.19bcd98fe5997p+100,
     0x1.3b6b286738596p+101, 0x1.0148ed1a44dabp+103, 0x1.5064fdae2f341p+104, 0x1.10a693476709cp+104},
    {0x1.262c3549bc5f3p+110, 0x1.26e9299292329p+109, 0x1.f1b8d56797fb8p+110, 0x1.1bc7a670c9a1cp+110, 0x1.fec4868f20872p+108,
     0x1.09c137b7a084bp+109, 0x1.06c921aa2d310p+111, 0x1.3c5d210bca17ap+110, 0x1.d5b1fe7b86799p+110, 0x1.850a466c338b6p+110,
     0x1.24ce4c90c2e5dp+109, 0x1.3e73fd50d3f81p+107, 0x1.9a383ecec64a3p+110, 0x1.633a83b30738dp+110, 0x1.c58f25bec77fap+109},
    {0x1.801f241addd8ep+115, 0x1.382c24a409fe3p+117, 0x1.dc70166252534p+116, 0x1.12e479cd4a461p+117, 0x1.408a0511c241dp+117,
     0x1.5a249a023c22dp+117, 0x1.7f4b9aa13a4d7p+117, 0x1.91db7a20f4d1dp+115, 0x1.5cc63cc77b8f4p+117, 0x1.a52bb1f3fea89p+115,
     0x1.79375f09b17bbp+116, 0x1.9a41d3bdede03p+115, 0x1.b6dff05f19c9cp+115, 0x1.6466b525edd66p+117, 0x1.475874c3f9e35p+116,
     0x1.10acea8068b2dp+117},
    {0x1.dcefcef8c973bp+120, 0x1.bbfb2533d8c1fp+123, 0x1.66e118dc65031p+123, 0x1.f34b65119bcdbp+120, 0x1.ed6894af0f3f8p+119,
     0x1.1524d218d1d88p+124, 0x1.0efa07cc655a2p+124, 0x1.89b9201c78976p+121, 0x1.34f28f356a76cp+123, 0x1.71b794fce9529p+123,
     0x1.edb5a5ccafd2dp+123, 0x1.4f98abcbd5e14p+119, 0x1.12a3bc9eb6567p+122, 0x1.274363027c870p+122, 0x1.6c2c01e6d2ccfp+123,
     0x1.07a7b21e793c2p+124, 0x1.35f48e8ff7034p+121},
    {0x1.7a158b5fa9753p+129, 0x1.3b1a6bd35fe1dp+127, 0x1.2b23517ddc326p+130, 0x1.125f4ba76520ap+130, 0x1.b2e2e769919f1p+126,
     0x1.a9ed228c8b41fp+125, 0x1.1d4ea349aab64p+127, 0x1.95a14dfc4cc12p+129, 0x1.093e0258f4c78p+129, 0x1.ba95d2e81acb2p+128,
     0x1.ecbecaf9c97aep+126, 0x1.28a97577165ddp+129, 0x1.3375ac9161e7cp+129, 0x1.f92d4b66410f0p+129, 0x1.0a1ae823ba68ep+129,
     0x1.b74553d1274a7p+129, 0x1.e61f0d96c9e7ap+129, 0x1.d2f36d7037cb4p+128},
    {0x1.10372a977e77cp+135, 0x1.29a6b535fd2dcp+136, 0x1.836e2020c9d22p+134, 0x1.58e970509d7bdp+134, 0x1.555a92d173b6bp+131,
     0x1.b8d9295a67c2fp+135, 0x1.50adf0b0747f7p+132, 0x1.21f64ebf5ade1p+136, 0x1.2f5eec48f1910p+133, 0x1.2085a1f64a5f5p+136,
     0x1.4f7676c20e4fbp+132, 0x1.6cd8aa4509f58p+132, 0x1.56cf6e251deabp+135, 0x1.ca917f1c0eddap+135, 0x1.41adcb5a04da1p+135,
     0x1.1b6a0595f5a86p+135, 0x1.c14b3efbd9a3ap+131, 0x1.565b1701a8e82p+133, 0x1.a3f27658622d8p+135},
    {0x1.be902320c03abp+139, 0x1.f16d15770dbc6p+141, 0x1.0a9d3143a2525p+141, 0x1.9c0e7d42ee808p+138, 0x1.43112661a8f18p+136,
     0x1.52d01b8d38c81p+141, 0x1.c263c384aa9cfp+140, 0x1.c961792acac2dp+137, 0x1.ea7e9da5efb61p+140, 0x1.ca1742a3ac7b7p+141,
     0x1.4d5db9373e38cp+136, 0x1.7d732004ae432p+141, 0x1.2077f7d8a70bap+141, 0x1.b4b6f673b88ebp+141, 0x1.a413b9b566a4dp+140,
     0x1.0009f68033bf6p+142, 0x1.cb469af212098p+141, 0x1.fe52fb1812f9fp+140, 0x1.17e9516b1b320p+141, 0x1.ad8cb07748533p+139},
};

inline constexpr unsigned qPi_2_complex_first_num_moduli = 6U;
inline constexpr double2 qPi_2_complex[][20]             = {
    {
     {0x1.d73ada4cda000p+45, 0x1.0000000000000p+3},
     {0x1.1b20a766e7000p+45, 0x1.e000000000000p+3},
     {0x1.cb754e348a000p+44, 0x1.6000000000000p+3},
     {0x1.4a9e41d4ea400p+46, 0x1.8000000000000p+1},
     {0x1.04b297c4d8000p+41, 0x1.0000000000000p+2},
     {0x1.0f48ca1d60000p+40, 0x1.4000000000000p+3},
     },
    {
     {0x1.029c55c9cf000p+54, 0x1.9580000000000p+10},
     {0x1.64a742ac1b800p+53, 0x1.2080000000000p+11},
     {0x1.d359768c65800p+53, 0x1.c800000000000p+6},
     {0x1.172a6daaa9000p+52, 0x1.1a00000000000p+7},
     {0x1.a1a69d643d000p+52, 0x1.7e80000000000p+11},
     {0x1.392cff9344800p+53, 0x1.1940000000000p+11},
     {0x1.01b0b85ca8800p+54, 0x1.61c0000000000p+11},
     },
    {
     {0x1.fdc084bc36000p+59, 0x1.4808400000000p+18},
     {0x1.912bf4dbae800p+60, 0x1.a1ee800000000p+17},
     {0x1.75313a2590000p+56, 0x1.d096000000000p+16},
     {0x1.403cb8b82a800p+61, 0x1.251f000000000p+16},
     {0x1.6cc5d067d8000p+61, 0x1.8128000000000p+16},
     {0x1.085aba77cf800p+61, 0x1.b48d000000000p+16},
     {0x1.5e8d778f0e000p+58, 0x1.fce9c00000000p+18},
     {0x1.3ac60b5400000p+54, 0x1.4809000000000p+16},
     },
    {
     {0x1.f5be6dc7c7400p+68, 0x1.59a3880000000p+21},
     {0x1.bead37f9aac00p+68, 0x1.71392e8000000p+25},
     {0x1.ede27068f8400p+68, 0x1.c041600000000p+21},
     {0x1.49611244f8800p+68, 0x1.eaad140000000p+22},
     {0x1.0dfbbe647a000p+67, 0x1.b7717b0000000p+25},
     {0x1.a56bb62e61800p+67, 0x1.3649640000000p+22},
     {0x1.5d1c76db50000p+66, 0x1.f31d5e0000000p+25},
     {0x1.7794953c3c000p+68, 0x1.6ceba28000000p+25},
     {0x1.abb8f38c9f000p+68, 0x1.37c3578000000p+25},
     },
    {
     {0x1.601b289a28400p+76, 0x1.76a356be00000p+32},
     {0x1.37f100f833800p+76, 0x1.c9ee992080000p+33},
     {0x1.891cda89ad000p+74, 0x1.e5abaf7680000p+33},
     {0x1.1f24b9265f000p+74, 0x1.de428e5980000p+33},
     {0x1.4ec03e910f800p+75, 0x1.3a62ac9400000p+30},
     {0x1.33248fa5eb000p+75, 0x1.be876eb900000p+33},
     {0x1.44fa716d1e800p+75, 0x1.bb97d93900000p+32},
     {0x1.0d10c7f048000p+74, 0x1.f9b60f8080000p+33},
     {0x1.2209d694e2000p+75, 0x1.2e203f3200000p+33},
     {0x1.3ae58390b7000p+76, 0x1.d0e5d32c00000p+31},
     },
    {
     {0x1.a490db3e87000p+82, 0x1.8907cdf253000p+41},
     {0x1.edfd5e6890000p+81, 0x1.38426e2fb3000p+40},
     {0x1.77169cd878000p+80, 0x1.a46871f81e800p+41},
     {0x1.1010fb262c000p+82, 0x1.e9bcd2078e000p+40},
     {0x1.0c29245b60000p+77, 0x1.c55ef1743a000p+40},
     {0x1.c9d134437c000p+82, 0x1.cf10e9d5da000p+39},
     {0x1.640ada62d8000p+80, 0x1.ed865d0ba0000p+37},
     {0x1.82c40bc378800p+83, 0x1.6281e55fff000p+40},
     {0x1.022ea3d62c000p+81, 0x1.ace1f0e319000p+41},
     {0x1.144985c704000p+83, 0x1.8df9658861000p+40},
     {0x1.88014c0539000p+82, 0x1.b42bcab4ca000p+41},
     },
    {
     {0x1.d906735665000p+89, 0x1.b5da10ec6a3c0p+48},
     {0x1.9a5a3b7068000p+87, 0x1.cc33f32374140p+47},
     {0x1.a1852bb572000p+88, 0x1.7892c2d80e240p+48},
     {0x1.2b8410d4c4000p+88, 0x1.3ef21b61a4e80p+47},
     {0x1.ae860d5ec0000p+85, 0x1.3bdf6da5991e0p+48},
     {0x1.54abff4a5e800p+90, 0x1.09f2f7731cd40p+47},
     {0x1.869a81af56000p+88, 0x1.dcf536cb59d20p+47},
     {0x1.2698f58650000p+89, 0x1.67d4ffaf64900p+45},
     {0x1.16fe214fce800p+90, 0x1.9fece9f7875a0p+47},
     {0x1.1c0090d1e6800p+90, 0x1.3625c591e2a00p+43},
     {0x1.0f7c313c8f000p+89, 0x1.e8d9d53db4f80p+47},
     {0x1.2743cce82a000p+88, 0x1.ed8f50decae60p+48},
     },
    {
     {0x1.792e2cfffd000p+96, 0x1.691225c311070p+52},
     {0x1.14ec55b54b000p+97, 0x1.f17027cb7b0f6p+54},
     {0x1.a235d05ba2000p+96, 0x1.668b8bb770a2cp+55},
     {0x1.0deccf55bf000p+97, 0x1.74dee02f93fa8p+55},
     {0x1.22fe0b3979000p+97, 0x1.87fb828ab4c1ep+55},
     {0x1.4159414da8000p+96, 0x1.fb217828bae26p+54},
     {0x1.200f172cee000p+97, 0x1.8dd5ffb808539p+52},
     {0x1.623d3c790e000p+96, 0x1.c23d3611b3111p+54},
     {0x1.2c3fd0f45d000p+97, 0x1.dcdf44ef91aa0p+54},
     {0x1.c1bf14d79d000p+96, 0x1.f48b5895d2415p+54},
     {0x1.620f02f07a000p+96, 0x1.dead7f9e2b646p+55},
     {0x1.d4026bec8c000p+96, 0x1.12ff535de9355p+55},
     {0x1.58c158b4e4000p+94, 0x1.d3c7ab366e742p+55},
     },
    {
     {0x1.4ca9f0b058000p+103, 0x1.b885702e55f9bp+62},
     {0x1.ab24151b42000p+102, 0x1.192d5ca0e143dp+61},
     {0x1.eef69b0110000p+101, 0x1.c530493aeefe1p+58},
     {0x1.80a90fa0f8000p+103, 0x1.3fb9177ea2b46p+62},
     {0x1.c0658ea122000p+103, 0x1.e15914c61dc56p+61},
     {0x1.a4ff64f033000p+103, 0x1.07459dc973089p+62},
     {0x1.d1888130f6000p+102, 0x1.9a322790b44c6p+59},
     {0x1.9892bae338000p+103, 0x1.f71d700823581p+59},
     {0x1.5592e97abf000p+104, 0x1.0270c7d83d5d3p+62},
     {0x1.19bcd98fe0000p+100, 0x1.665b282870c8ap+62},
     {0x1.3b6b286738000p+101, 0x1.659d012a4eca3p+59},
     {0x1.0148ed1a44000p+103, 0x1.b55ae52593856p+62},
     {0x1.5064fdae2f000p+104, 0x1.a04de44d1dbb2p+61},
     {0x1.10a6934767000p+104, 0x1.37277c85b1396p+59},
     },
    {
     {0x1.262c3549bc000p+110, 0x1.7cc183b65ce20p+68},
     {0x1.26e9299292000p+109, 0x1.948846d04e1e4p+66},
     {0x1.f1b8d56797000p+110, 0x1.f700757b963bdp+69},
     {0x1.1bc7a670c9000p+110, 0x1.438c35596874fp+69},
     {0x1.fec4868f20000p+108, 0x1.0e4d33c587550p+67},
     {0x1.09c137b7a0000p+109, 0x1.0963f200c3d0ap+68},
     {0x1.06c921aa2d000p+111, 0x1.87cb7b0aa7019p+68},
     {0x1.3c5d210bca000p+110, 0x1.799b2fbd07c75p+66},
     {0x1.d5b1fe7b86000p+110, 0x1.e647ab21d7303p+68},
     {0x1.850a466c33000p+110, 0x1.16b6ff9ebde85p+69},
     {0x1.24ce4c90c2000p+109, 0x1.cba2857b66ec7p+68},
     {0x1.3e73fd50d0000p+107, 0x1.fc0a16476464cp+68},
     {0x1.9a383ecec6000p+110, 0x1.28be09d6a8097p+68},
     {0x1.633a83b307000p+110, 0x1.c676f7509d9edp+67},
     {0x1.c58f25bec6000p+109, 0x1.7f9e00186a43fp+69},
     },
    {
     {0x1.801f241adc000p+115, 0x1.d8d9437d21855p+75},
     {0x1.382c24a409800p+117, 0x1.f8a514359d27ap+75},
     {0x1.dc70166252000p+116, 0x1.4d16c5f864abcp+74},
     {0x1.12e479cd4a000p+117, 0x1.18534971a2209p+75},
     {0x1.408a0511c2000p+117, 0x1.0743460adc6bep+75},
     {0x1.5a249a023c000p+117, 0x1.16ac76137c302p+74},
     {0x1.7f4b9aa13a000p+117, 0x1.35bbad2f596e5p+75},
     {0x1.91db7a20f4000p+115, 0x1.a39344e00882ep+74},
     {0x1.5cc63cc77b800p+117, 0x1.e88013e001f67p+72},
     {0x1.a52bb1f3fe000p+115, 0x1.51298a72234e5p+74},
     {0x1.79375f09b1000p+116, 0x1.eead8becc0ca0p+74},
     {0x1.9a41d3bdec000p+115, 0x1.e031fd6725785p+75},
     {0x1.b6dff05f18000p+115, 0x1.c9b879b0b4f19p+75},
     {0x1.6466b525ed800p+117, 0x1.596412e1ef120p+75},
     {0x1.475874c3f9000p+116, 0x1.c6af7e9d25945p+75},
     {0x1.10acea8068800p+117, 0x1.968a313878fd8p+74},
     },
    {
     {0x1.dcefcef8c8000p+120, 0x1.73ac9e398255bp+80},
     {0x1.bbfb2533d8000p+123, 0x1.83dc4c99376fep+82},
     {0x1.66e118dc65000p+123, 0x1.877f993dfa961p+76},
     {0x1.f34b651198000p+120, 0x1.e6d43734932d4p+81},
     {0x1.ed6894af00000p+119, 0x1.e7ef2c567c166p+82},
     {0x1.1524d218d1800p+124, 0x1.61fc6f8312bf2p+82},
     {0x1.0efa07cc65000p+124, 0x1.6860e7822ac71p+82},
     {0x1.89b9201c78000p+121, 0x1.2eb803ef7f9cap+80},
     {0x1.34f28f356a000p+123, 0x1.db13781e55965p+81},
     {0x1.71b794fce9000p+123, 0x1.4a24c81d49af7p+81},
     {0x1.edb5a5ccaf000p+123, 0x1.a59afb9185b16p+82},
     {0x1.4f98abcbd0000p+119, 0x1.7851c0a536222p+81},
     {0x1.12a3bc9eb6000p+122, 0x1.59bb5d67eb532p+80},
     {0x1.274363027c000p+122, 0x1.0e0caffdedc0dp+81},
     {0x1.6c2c01e6d2000p+123, 0x1.99d6700ed3683p+82},
     {0x1.07a7b21e79000p+124, 0x1.e1374ee7a2bcep+81},
     {0x1.35f48e8ff4000p+121, 0x1.81a1f15f33188p+82},
     },
    {
     {0x1.7a158b5fa9000p+129, 0x1.d4cbc2d728674p+87},
     {0x1.3b1a6bd35c000p+127, 0x1.f0e83e4d5191dp+88},
     {0x1.2b23517ddc000p+130, 0x1.932a282f695a4p+87},
     {0x1.125f4ba765000p+130, 0x1.050588b76c069p+87},
     {0x1.b2e2e76990000p+126, 0x1.9f095a162ddd9p+86},
     {0x1.a9ed228c80000p+125, 0x1.683e2504c7780p+88},
     {0x1.1d4ea349a8000p+127, 0x1.5b1f371597c8bp+88},
     {0x1.95a14dfc4c000p+129, 0x1.823239bb62f47p+88},
     {0x1.093e0258f4000p+129, 0x1.8f080e71d6b7ep+88},
     {0x1.ba95d2e81a000p+128, 0x1.963d61c2015f1p+87},
     {0x1.ecbecaf9c8000p+126, 0x1.7ae41709de774p+86},
     {0x1.28a9757716000p+129, 0x1.775de320b5a60p+87},
     {0x1.3375ac9161000p+129, 0x1.cf73230fd5e89p+88},
     {0x1.f92d4b6641000p+129, 0x1.e09b51c74c7e8p+84},
     {0x1.0a1ae823ba000p+129, 0x1.a3728283aab4bp+87},
     {0x1.b74553d127000p+129, 0x1.29ae7a7eb6fbap+87},
     {0x1.e61f0d96c9000p+129, 0x1.cf3df03b3a8bcp+88},
     {0x1.d2f36d7036000p+128, 0x1.cb409a9781493p+88},
     },
    {
     {0x1.10372a977e000p+135, 0x1.defe8fb4fb909p+93},
     {0x1.29a6b535fd000p+136, 0x1.6e2cf231b9668p+93},
     {0x1.836e2020c8000p+134, 0x1.d21cbab2784cdp+94},
     {0x1.58e970509c000p+134, 0x1.7bc89ef6dad2dp+94},
     {0x1.555a92d170000p+131, 0x1.db5b6f5b0eb4dp+92},
     {0x1.b8d9295a67000p+135, 0x1.85d1260ca54fcp+94},
     {0x1.50adf0b070000p+132, 0x1.1fdc21d5b15ddp+94},
     {0x1.21f64ebf5a800p+136, 0x1.7823ce89e7797p+94},
     {0x1.2f5eec48f0000p+133, 0x1.90ff01dc4d313p+93},
     {0x1.2085a1f64a000p+136, 0x1.7d4780667df06p+94},
     {0x1.4f7676c208000p+132, 0x1.93ea891064dd2p+94},
     {0x1.6cd8aa4508000p+132, 0x1.f57982fe71a7fp+92},
     {0x1.56cf6e251d000p+135, 0x1.d556dd956e0cdp+94},
     {0x1.ca917f1c0e000p+135, 0x1.bb443d8b5ac4bp+94},
     {0x1.41adcb5a04000p+135, 0x1.b4282391888b4p+94},
     {0x1.1b6a0595f5000p+135, 0x1.50c4f68823a52p+94},
     {0x1.c14b3efbd0000p+131, 0x1.3474ca3d9a698p+94},
     {0x1.565b1701a8000p+133, 0x1.d0313174500a1p+92},
     {0x1.a3f2765862000p+135, 0x1.6be16bef92ca4p+92},
     },
    {
     {0x1.be902320c0000p+139, 0x1.d589e88ff78b1p+96},
     {0x1.f16d15770d000p+141, 0x1.78b27d6739913p+100},
     {0x1.0a9d3143a2000p+141, 0x1.49597296149b8p+99},
     {0x1.9c0e7d42e8000p+138, 0x1.a01e2608f0d4cp+100},
     {0x1.43112661a0000p+136, 0x1.1e2feb7d8b7adp+99},
     {0x1.52d01b8d38000p+141, 0x1.9023c3245dc8ap+100},
     {0x1.c263c384aa000p+140, 0x1.39e29b53d8ee9p+99},
     {0x1.c961792ac0000p+137, 0x1.5859c3d8fcfb4p+100},
     {0x1.ea7e9da5ee000p+140, 0x1.b60d0fc523a99p+100},
     {0x1.ca1742a3ac000p+141, 0x1.edb1b08d8ab6cp+99},
     {0x1.4d5db93720000p+136, 0x1.e38c4473ac422p+100},
     {0x1.7d732004ae000p+141, 0x1.0c61f7387e1cfp+99},
     {0x1.2077f7d8a7000p+141, 0x1.732d6756b8224p+96},
     {0x1.b4b6f673b8000p+141, 0x1.1d57969c83715p+100},
     {0x1.a413b9b566000p+140, 0x1.499292240ced4p+99},
     {0x1.0009f68033800p+142, 0x1.fb24d7f79814bp+99},
     {0x1.cb469af212000p+141, 0x1.2f227bd8c7bf7p+96},
     {0x1.fe52fb1812000p+140, 0x1.f3efba84941bdp+99},
     {0x1.17e9516b1b000p+141, 0x1.9033086e6ab83p+98},
     {0x1.ad8cb07748000p+139, 0x1.4cd50b4cac0c8p+97},
     },
};

static_assert(sizeof(qPi_2) / sizeof(qPi_2[0]) == 21U - qPi_2_first_num_moduli);
static_assert(sizeof(qPi_2_complex) / sizeof(qPi_2_complex[0]) == 21U - qPi_2_complex_first_num_moduli);
static_assert(qPi_2_first_num_moduli <= threshold<Backend::INT8, false>::P_is_double + 1U);
static_assert(qPi_2_complex_first_num_moduli <= threshold<Backend::INT8, true>::P_is_double + 1U);
static_assert(P[threshold<Backend::INT8, false>::P_is_double - 2].y == 0.0);
static_assert(P_complex[threshold<Backend::INT8, true>::P_is_double - 2].y == 0.0);

} // namespace INT8

namespace FP8 {

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

inline constexpr unsigned qPi_2_first_num_moduli = 5U;
inline constexpr double2 qPi_2[][20]             = {
    {
     {0x1.cb92e5aa70000p+53, 0x1.0ec0000000000p+10},
     {0x1.ad14489be4000p+53, 0x1.fa5c000000000p+14},
     {0x1.d015a9a0ee000p+54, 0x1.3e50000000000p+12},
     {0x1.4a33a585ea000p+54, 0x1.973c000000000p+14},
     {0x1.7120ef69e0000p+53, 0x1.0ce0000000000p+12},
     },
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

inline constexpr double qPi_1_complex[19][20] = {
    {0x1.60f2000000000p+18, 0x1.d99a200000000p+20},
    {0x1.1c41116800000p+31, 0x1.1367810000000p+24, 0x1.4701e1d000000p+28},
    {0x1.5d92071176000p+40, 0x1.47d1d06280000p+33, 0x1.782e915a80000p+38, 0x1.fe1964a07e000p+39},
    {0x1.4621a8e14ccdcp+51, 0x1.1b299e3fce910p+51, 0x1.6165fc22706d0p+51, 0x1.f465604f2c720p+49, 0x1.46fcec60ec086p+51},
    {0x1.9cff32a800c7dp+59, 0x1.c2454193522e8p+60, 0x1.434c00b8040b5p+58, 0x1.91d0e8ec712c8p+56, 0x1.1fb246cc85d3ep+61,
     0x1.3b026d3806e79p+57},
    {0x1.a73256794a039p+68, 0x1.0be9b2ab89073p+71, 0x1.181d8efe85f61p+68, 0x1.02790f07e4e68p+67, 0x1.aba8fa8fa3e86p+70,
     0x1.bee93b7146b53p+67, 0x1.32b41cd6ccb8fp+67},
    {0x1.a81aed2f88b05p+80, 0x1.29d1ac310b359p+80, 0x1.01c5df12b81d4p+81, 0x1.3e578c8273fe8p+80, 0x1.ba5801b32e294p+80,
     0x1.2e0f0a67fa96ap+77, 0x1.10078875d3bbcp+79, 0x1.2142c039a99e7p+79},
    {0x1.7ac71e8cfa433p+90, 0x1.2b79a9b12f366p+90, 0x1.798cb6e645b1fp+88, 0x1.a20e093cb87c5p+84, 0x1.61c6f1dc88de1p+89,
     0x1.3cc0faa680f5bp+90, 0x1.5c1d55122631ap+85, 0x1.9bdf8f2da9080p+90, 0x1.5fc833fe433a7p+90},
    {0x1.50c34e6895180p+100, 0x1.852b35253d205p+100, 0x1.1b8085bd88e3ap+100, 0x1.81e71d75238d3p+98, 0x1.650f2a5af4a32p+100,
     0x1.eb7669ee0fd10p+99, 0x1.c1e0d02e8267fp+99, 0x1.9f4530e10f8cdp+100, 0x1.2797e16e6a2acp+99, 0x1.9ade8c55812a4p+99},
    {0x1.9ac44beba40cbp+107, 0x1.0c28d040f9799p+109, 0x1.5c28ed468d9e5p+110, 0x1.b5096e326fdc3p+107, 0x1.6c811dc4c7dd5p+110,
     0x1.52ba3b2f93e6bp+108, 0x1.6db1da94318d5p+110, 0x1.24d6dbe91d7dbp+109, 0x1.008c72cd6c37bp+110, 0x1.c0717d2fdfb29p+107,
     0x1.66e0c7fa46acbp+106},
    {0x1.6d55012cffee0p+118, 0x1.343c40760d695p+118, 0x1.d36685d406f4cp+119, 0x1.0a2c1fd952fe2p+120, 0x1.b0cf1bf7b8bcap+118,
     0x1.151fbee9b70f9p+118, 0x1.0f42d910fd982p+120, 0x1.9cd132160bff3p+118, 0x1.1160e1fd4a9fap+118, 0x1.6a67e4f2df261p+118,
     0x1.cd9c791b23352p+119, 0x1.11d93e7aa1db2p+120},
    {0x1.2c0b103735314p+129, 0x1.76d40783a0e62p+129, 0x1.6a41ad33988acp+129, 0x1.2469af65b7c22p+129, 0x1.6b72ff632bfffp+128,
     0x1.2230493e98827p+129, 0x1.db29989a75290p+129, 0x1.a6b3b61732b9cp+127, 0x1.cfdc970c0d2e8p+128, 0x1.52da04cc981b6p+128,
     0x1.5366390fe8f53p+128, 0x1.de3bcfe64450fp+129, 0x1.64e4f369f9663p+128},
    {0x1.7c0c8a7c81679p+134, 0x1.934565c48d7f2p+139, 0x1.6a03dc0983d2dp+137, 0x1.5ba3fc9441d32p+137, 0x1.fc03e98c2d22cp+138,
     0x1.7f5399aa5cdbap+138, 0x1.35268123c9947p+138, 0x1.96221b1c2cf19p+138, 0x1.0351c95b28f11p+138, 0x1.ea31027e49e9dp+137,
     0x1.61b70367b827fp+139, 0x1.1678b4700f5fap+139, 0x1.0f4d5afd363d7p+138, 0x1.cb4fcea07d127p+135},
    {0x1.56350c72cc8fdp+148, 0x1.a0f4599cad038p+146, 0x1.521dce02e3056p+146, 0x1.d4e1dd48087f5p+148, 0x1.20913b4afefeep+148,
     0x1.17e7187739a12p+149, 0x1.f2c29516ed5e2p+147, 0x1.c991b168992bcp+147, 0x1.a315c7e195d74p+147, 0x1.c82448cdb3718p+146,
     0x1.d6ae1ec7fa030p+146, 0x1.9ee4478593f7ep+147, 0x1.505902b3c615dp+147, 0x1.d7f3a364188c9p+148, 0x1.7c619d9ec5705p+147},
    {0x1.881c17da60abcp+157, 0x1.e230c46ebc598p+152, 0x1.42c62010d4cecp+158, 0x1.39dbd5ec1d328p+158, 0x1.d5c38182561a6p+158,
     0x1.6234dee80fb38p+158, 0x1.e88ccf571bf67p+157, 0x1.529f59343e81fp+157, 0x1.6a891fd84c4bcp+158, 0x1.e98641324acd5p+157,
     0x1.a4aaa86aeaf19p+158, 0x1.5e168f1a6b85fp+154, 0x1.13a829022476bp+156, 0x1.1c21cb2d1f9c2p+158, 0x1.4b8c5c2b247bap+158,
     0x1.d43e3f3c287eep+156},
    {0x1.360b52e3097e8p+168, 0x1.7be569105d3e2p+167, 0x1.00ed5195632f6p+168, 0x1.ef9c295a67e1ep+167, 0x1.5391f02afd32fp+168,
     0x1.5aa6e5ee60788p+168, 0x1.752947700f092p+166, 0x1.6af1dc4ea36bep+167, 0x1.389ae4cc078d3p+168, 0x1.18d5503da610cp+168,
     0x1.b8b9c55ff8efep+165, 0x1.44be44ccc558cp+168, 0x1.65e288b2dfef2p+168, 0x1.a54821a3d2353p+166, 0x1.18eb028c8a908p+168,
     0x1.3f149599a8057p+168, 0x1.b946c1e25b8cbp+165},
    {0x1.bb8ad7f3af3acp+176, 0x1.f22c7abe57c10p+177, 0x1.24447d378bc09p+177, 0x1.af49e7dfa63a8p+177, 0x1.73809d4fe8112p+175,
     0x1.297c532abebe8p+177, 0x1.29b6119d6d593p+177, 0x1.c159cf16f2966p+172, 0x1.84e93ec95a4ebp+177, 0x1.ea6d3ff04cd94p+177,
     0x1.6ca2789f99600p+177, 0x1.afe6a3e579993p+175, 0x1.9bfd07c92bd11p+175, 0x1.501f418ef4eb4p+176, 0x1.04b3c73e29883p+178,
     0x1.f325786b9846cp+177, 0x1.86ad606557f18p+176, 0x1.f1cd816413a9bp+177},
    {0x1.73c597315cd24p+187, 0x1.f6173bfdd888bp+186, 0x1.9e74ed78f55a8p+185, 0x1.a8dd7620a5490p+186, 0x1.473590ae95feap+187,
     0x1.e96e3ea5ae1fcp+185, 0x1.429740dfec1d1p+187, 0x1.29269e02768d2p+187, 0x1.e90736f1e80e2p+186, 0x1.271452e6fd903p+184,
     0x1.2ea9aa4c557bcp+184, 0x1.2f21bca0d21f9p+187, 0x1.6421382ca7097p+187, 0x1.73b829be903afp+187, 0x1.086d2a1570d3cp+186,
     0x1.8881e5b236f22p+182, 0x1.45b401e30eebap+187, 0x1.a7ec9a4ee8557p+186, 0x1.8a0a93d6ccd50p+182},
    {0x1.da15e7b4ba6ebp+196, 0x1.01731bfa0f7f4p+195, 0x1.4589bd7148e29p+196, 0x1.ca64799b18ba8p+189, 0x1.bd80816e5ac21p+196,
     0x1.48bdae311c81dp+196, 0x1.36226287a3915p+196, 0x1.0ca1fbf4ddcf2p+197, 0x1.670f2944f9936p+196, 0x1.3593882c620b1p+196,
     0x1.a85b1ef3bd078p+196, 0x1.4c42a2dd463d4p+196, 0x1.a069d4986949ap+192, 0x1.a03c37c1ccba0p+195, 0x1.498d401938aa5p+196,
     0x1.3fe3f82c04592p+194, 0x1.8156f28b0a63fp+196, 0x1.bfa6cefcd9b74p+196, 0x1.008460b0af496p+195, 0x1.c7f1e431e4938p+192},
};

inline constexpr unsigned qPi_2_complex_first_num_moduli = 6U;
inline constexpr double2 qPi_2_complex[15][20]           = {
    {
     {0x1.9cff32a800000p+59, 0x1.8f96000000000p+18},
     {0x1.c245419352000p+60, 0x1.73d9800000000p+17},
     {0x1.434c00b800000p+58, 0x1.02d5c00000000p+20},
     {0x1.91d0e8ec60000p+56, 0x1.12c8800000000p+20},
     {0x1.1fb246cc85000p+61, 0x1.a7b8300000000p+20},
     {0x1.3b026d3800000p+57, 0x1.b9e3200000000p+19},
     },
    {
     {0x1.a732567940000p+68, 0x1.4072b35200000p+31},
     {0x1.0be9b2ab88000p+71, 0x1.0733947e00000p+31},
     {0x1.181d8efe80000p+68, 0x1.7d85459000000p+30},
     {0x1.02790f07e0000p+67, 0x1.39a14bc000000p+29},
     {0x1.aba8fa8fa0000p+70, 0x1.f4315b5400000p+31},
     {0x1.bee93b7140000p+67, 0x1.ad4d23c000000p+29},
     {0x1.32b41cd6c0000p+67, 0x1.971e317400000p+30},
     },
    {
     {0x1.a81aed2f88000p+80, 0x1.609c3e5a2e000p+39},
     {0x1.29d1ac3108000p+80, 0x1.9acb743b60000p+41},
     {0x1.01c5df12b8000p+81, 0x1.d44afacf30000p+37},
     {0x1.3e578c8270000p+80, 0x1.ff3cce93ad800p+41},
     {0x1.ba5801b32c000p+80, 0x1.14a39d4caa800p+41},
     {0x1.2e0f0a67e0000p+77, 0x1.a96a06c86a800p+41},
     {0x1.10078875d0000p+79, 0x1.dddff95eb9000p+40},
     {0x1.2142c039a8000p+79, 0x1.9e73a67e54000p+39},
     },
    {
     {0x1.7ac71e8cfa000p+90, 0x1.0ccd50a81adc0p+48},
     {0x1.2b79a9b12e000p+90, 0x1.365da7c23e264p+50},
     {0x1.798cb6e640000p+88, 0x1.6c7ad03a72cf0p+50},
     {0x1.a20e093c80000p+84, 0x1.c3e264d0d9bf0p+49},
     {0x1.61c6f1dc88000p+89, 0x1.bc2224d2da4d0p+48},
     {0x1.3cc0faa680000p+90, 0x1.eb615ed7a0428p+49},
     {0x1.5c1d551200000p+85, 0x1.318d1a8bf1490p+50},
     {0x1.9bdf8f2da8000p+90, 0x1.0804fa62f5b08p+50},
     {0x1.5fc833fe42000p+90, 0x1.3a736b2902d78p+50},
     },
    {
     {0x1.50c34e6894000p+100, 0x1.17f9a4a23d6fdp+60},
     {0x1.852b35253c000p+100, 0x1.20553f19e93fep+60},
     {0x1.1b8085bd88000p+100, 0x1.c73635c60e9e8p+59},
     {0x1.81e71d7520000p+98, 0x1.c696a93ed322bp+59},
     {0x1.650f2a5af4000p+100, 0x1.4640f5ea4f137p+59},
     {0x1.eb7669ee0c000p+99, 0x1.e881ecb43fa28p+60},
     {0x1.c1e0d02e80000p+99, 0x1.33f502c9603cdp+60},
     {0x1.9f4530e10e000p+100, 0x1.8cd40b023a124p+60},
     {0x1.2797e16e68000p+99, 0x1.156221de020d2p+60},
     {0x1.9ade8c5580000p+99, 0x1.2a41abe7fd8a5p+59},
     },
    {
     {0x1.9ac44beba0000p+107, 0x1.032cabc207437p+69},
     {0x1.0c28d040f8000p+109, 0x1.79965f021ad94p+69},
     {0x1.5c28ed468c000p+110, 0x1.9e53ebc5415bbp+70},
     {0x1.b5096e3260000p+107, 0x1.fb85e0e284cabp+70},
     {0x1.6c811dc4c6000p+110, 0x1.dd549990be864p+70},
     {0x1.52ba3b2f90000p+108, 0x1.f3567cfeb3bdbp+69},
     {0x1.6db1da9430000p+110, 0x1.8d57932ac447bp+70},
     {0x1.24d6dbe91c000p+109, 0x1.7db59c7e9c135p+69},
     {0x1.008c72cd6c000p+110, 0x1.bdae92be7c7c2p+67},
     {0x1.c0717d2fd0000p+107, 0x1.f6511abaeaa34p+70},
     {0x1.66e0c7fa40000p+106, 0x1.ab2c2fe46f9bep+68},
     },
    {
     {0x1.6d55012cf8000p+118, 0x1.fb81df9ad6512p+80},
     {0x1.343c407608000p+118, 0x1.5a55fdc1f2a19p+80},
     {0x1.d36685d404000p+119, 0x1.7a5ca5b1f45cfp+80},
     {0x1.0a2c1fd952000p+120, 0x1.fc428013ddf7dp+79},
     {0x1.b0cf1bf7b8000p+118, 0x1.7932f0b02f2f5p+77},
     {0x1.151fbee9b0000p+118, 0x1.c3e326474e5eap+80},
     {0x1.0f42d910fc000p+120, 0x1.9819442644f1cp+80},
     {0x1.9cd1321608000p+118, 0x1.ff965fe6db479p+79},
     {0x1.1160e1fd48000p+118, 0x1.4fd271fd91042p+79},
     {0x1.6a67e4f2d8000p+118, 0x1.c9836e5864ceep+80},
     {0x1.cd9c791b20000p+119, 0x1.9a8f13a535619p+80},
     {0x1.11d93e7aa0000p+120, 0x1.db246e101a30bp+80},
     },
    {
     {0x1.2c0b103734000p+129, 0x1.313da793aa109p+89},
     {0x1.76d40783a0000p+129, 0x1.cc4cbd188db8dp+88},
     {0x1.6a41ad3398000p+129, 0x1.158b66205b127p+88},
     {0x1.2469af65b6000p+129, 0x1.c223e0618dbe3p+89},
     {0x1.6b72ff6328000p+128, 0x1.fff6888b6cf92p+89},
     {0x1.2230493e98000p+129, 0x1.04ed50b70ac3cp+88},
     {0x1.db29989a74000p+129, 0x1.290537bb6d3cbp+89},
     {0x1.a6b3b61730000p+127, 0x1.5ce0501ec6240p+88},
     {0x1.cfdc970c0c000p+128, 0x1.2e807e7b6e9abp+88},
     {0x1.52da04cc98000p+128, 0x1.b5b7682690668p+84},
     {0x1.5366390fe8000p+128, 0x1.ea59044f6badfp+87},
     {0x1.de3bcfe644000p+129, 0x1.43b0f113e242cp+87},
     {0x1.64e4f369f8000p+128, 0x1.6632eccf9719bp+88},
     },
    {
     {0x1.7c0c8a7c80000p+134, 0x1.6792bfc0eb86dp+94},
     {0x1.934565c48c000p+139, 0x1.7f21bc53f8151p+99},
     {0x1.6a03dc0980000p+137, 0x1.e966e512ee21ep+98},
     {0x1.5ba3fc9440000p+137, 0x1.d3236e3aaaaa4p+97},
     {0x1.fc03e98c2c000p+138, 0x1.22be26b6b0d45p+98},
     {0x1.7f5399aa5c000p+138, 0x1.b73ba3dafd0bcp+97},
     {0x1.35268123c8000p+138, 0x1.946c8ab39ebf9p+98},
     {0x1.96221b1c2c000p+138, 0x1.e328b7ab3ab04p+97},
     {0x1.0351c95b28000p+138, 0x1.e21b39c32af54p+97},
     {0x1.ea31027e48000p+137, 0x1.e9ce330e287f0p+97},
     {0x1.61b70367b8000p+139, 0x1.3f8667b67c80ap+96},
     {0x1.1678b4700e000p+139, 0x1.5f9ad8f24be6ep+99},
     {0x1.0f4d5afd34000p+138, 0x1.1eb7a6c91688dp+99},
     {0x1.cb4fcea060000p+135, 0x1.d126f1b1e60eap+99},
     },
    {
     {0x1.56350c72cc000p+148, 0x1.1f9f451090ae4p+107},
     {0x1.a0f4599ca0000p+146, 0x1.a06f6fc5bdd80p+109},
     {0x1.521dce02e0000p+146, 0x1.82b19a26c7cfdp+107},
     {0x1.d4e1dd4808000p+148, 0x1.fd4856c8fdad7p+106},
     {0x1.20913b4afc000p+148, 0x1.7f6d7882f3ea6p+109},
     {0x1.17e7187738000p+149, 0x1.a121554211bcfp+109},
     {0x1.f2c29516e8000p+147, 0x1.578979152ff2fp+109},
     {0x1.c991b16898000p+147, 0x1.2bc7f25e634ffp+107},
     {0x1.a315c7e190000p+147, 0x1.75d17e16210f3p+109},
     {0x1.c82448cdb0000p+146, 0x1.b8bf53d6ddafcp+107},
     {0x1.d6ae1ec7f0000p+146, 0x1.40609be208734p+109},
     {0x1.9ee4478590000p+147, 0x1.fbee6ac71c8e0p+108},
     {0x1.505902b3c0000p+147, 0x1.8574c75cc0648p+109},
     {0x1.d7f3a36418000p+148, 0x1.191627ebb2ea0p+107},
     {0x1.7c619d9ec0000p+147, 0x1.5c12eeb8eb98bp+109},
     },
    {
     {0x1.881c17da60000p+157, 0x1.578ee36fc113ep+116},
     {0x1.e230c46e80000p+152, 0x1.e2cc2a2446d1fp+117},
     {0x1.42c62010d4000p+158, 0x1.9d8b412e1883ep+117},
     {0x1.39dbd5ec1c000p+158, 0x1.3279ca38634c4p+118},
     {0x1.d5c3818256000p+158, 0x1.a64d4d6838913p+114},
     {0x1.6234dee80e000p+158, 0x1.b380f55bcab0fp+118},
     {0x1.e88ccf5718000p+157, 0x1.fb39a87a7c4c1p+118},
     {0x1.529f59343c000p+157, 0x1.40f6e4c38810bp+118},
     {0x1.6a891fd84c000p+158, 0x1.2f133ebb5566ap+116},
     {0x1.e986413248000p+157, 0x1.66a6ce1a9c9b5p+118},
     {0x1.a4aaa86aea000p+158, 0x1.e32aba56bf0d4p+117},
     {0x1.5e168f1a60000p+154, 0x1.70bd826e86b38p+117},
     {0x1.13a8290220000p+156, 0x1.1dacf9e6835d0p+118},
     {0x1.1c21cb2d1e000p+158, 0x1.9c25a30ef20aap+118},
     {0x1.4b8c5c2b24000p+158, 0x1.ee73586711ff3p+116},
     {0x1.d43e3f3c28000p+156, 0x1.fb85e4d787e07p+114},
     },
    {
     {0x1.360b52e308000p+168, 0x1.7e794e960d7e9p+128},
     {0x1.7be5691058000p+167, 0x1.4f87deca7c303p+129},
     {0x1.00ed519560000p+168, 0x1.97b08fb0756a6p+129},
     {0x1.ef9c295a60000p+167, 0x1.f877ca98d75ebp+129},
     {0x1.5391f02afc000p+168, 0x1.32ea596db03edp+128},
     {0x1.5aa6e5ee60000p+168, 0x1.e1ea4b2efca8cp+126},
     {0x1.7529477000000p+166, 0x1.e123fdc720dbap+129},
     {0x1.6af1dc4ea0000p+167, 0x1.b5f14bc82e013p+128},
     {0x1.389ae4cc04000p+168, 0x1.c697b67c33603p+129},
     {0x1.18d5503da4000p+168, 0x1.085eca655c872p+129},
     {0x1.b8b9c55fe0000p+165, 0x1.8efd985892b3bp+129},
     {0x1.44be44ccc4000p+168, 0x1.58bdf33ced3b8p+128},
     {0x1.65e288b2dc000p+168, 0x1.f78c1d9195671p+129},
     {0x1.a54821a3d0000p+166, 0x1.1a96ea383afebp+127},
     {0x1.18eb028c88000p+168, 0x1.484169265c6c2p+129},
     {0x1.3f149599a8000p+168, 0x1.5d8ebe4393792p+122},
     {0x1.b946c1e240000p+165, 0x1.b8cb764964223p+129},
     },
    {
     {0x1.bb8ad7f3a0000p+176, 0x1.e75835dea15ebp+139},
     {0x1.f22c7abe50000p+177, 0x1.f03f7eece390cp+139},
     {0x1.24447d3788000p+177, 0x1.e0499eb7a7a55p+138},
     {0x1.af49e7dfa0000p+177, 0x1.8ea173c43f570p+139},
     {0x1.73809d4fe0000p+175, 0x1.022434169cf9dp+138},
     {0x1.297c532ab8000p+177, 0x1.afa1d6dcc2f81p+139},
     {0x1.29b6119d68000p+177, 0x1.564de8035b89fp+139},
     {0x1.c159cf1600000p+172, 0x1.e52cbc2dfde34p+139},
     {0x1.84e93ec958000p+177, 0x1.275aa57d5224bp+138},
     {0x1.ea6d3ff048000p+177, 0x1.364e726f4777cp+139},
     {0x1.6ca2789f98000p+177, 0x1.600696365edbap+137},
     {0x1.afe6a3e560000p+175, 0x1.999281760dc53p+139},
     {0x1.9bfd07c920000p+175, 0x1.7a22d0d911840p+138},
     {0x1.501f418ef0000p+176, 0x1.3acf8e08be5b7p+138},
     {0x1.04b3c73e28000p+178, 0x1.883138ff05750p+138},
     {0x1.f325786b98000p+177, 0x1.1af2336105505p+135},
     {0x1.86ad606550000p+176, 0x1.fc6098f8af7c4p+138},
     {0x1.f1cd816410000p+177, 0x1.d4d739ea61e36p+138},
     },
    {
     {0x1.73c597315c000p+187, 0x1.a476aac68fd0dp+146},
     {0x1.f6173bfdd8000p+186, 0x1.115a81cdfc88ep+145},
     {0x1.9e74ed78f0000p+185, 0x1.569ebd7a3c3d8p+147},
     {0x1.a8dd7620a0000p+186, 0x1.5241660da9d7bp+148},
     {0x1.473590ae94000p+187, 0x1.fea32d2074e76p+147},
     {0x1.e96e3ea5a0000p+185, 0x1.c3f7519008289p+148},
     {0x1.429740dfec000p+187, 0x1.d0a2224cc9aa2p+143},
     {0x1.29269e0274000p+187, 0x1.46901cfc3c7b1p+148},
     {0x1.e90736f1e8000p+186, 0x1.c46bf94971537p+141},
     {0x1.271452e6e0000p+184, 0x1.d902dd4f5c854p+148},
     {0x1.2ea9aa4c40000p+184, 0x1.57bc697119f04p+148},
     {0x1.2f21bca0d0000p+187, 0x1.0fca25cbdef2dp+148},
     {0x1.6421382ca4000p+187, 0x1.84ba6a944c1e5p+148},
     {0x1.73b829be90000p+187, 0x1.d78b27a36d31cp+144},
     {0x1.086d2a1570000p+186, 0x1.a77e9de3cf11fp+145},
     {0x1.8881e5b200000p+182, 0x1.b79106d72fcf7p+147},
     {0x1.45b401e30c000p+187, 0x1.75d319dd2d940p+148},
     {0x1.a7ec9a4ee8000p+186, 0x1.55af1de2b2f16p+144},
     {0x1.8a0a93d680000p+182, 0x1.3354051e14010p+148},
     },
    {
     {0x1.da15e7b4b8000p+196, 0x1.3755678e97cc2p+157},
     {0x1.01731bfa00000p+195, 0x1.efe7bf4d0a8d8p+158},
     {0x1.4589bd7148000p+196, 0x1.c52ec7bd5152bp+155},
     {0x1.ca64799800000p+189, 0x1.8c5d3c4fdbd05p+158},
     {0x1.bd80816e58000p+196, 0x1.61046a324feddp+157},
     {0x1.48bdae3118000p+196, 0x1.2074aeeab74abp+158},
     {0x1.36226287a0000p+196, 0x1.c8a695290afe2p+157},
     {0x1.0ca1fbf4dc000p+197, 0x1.cf1eff458d031p+157},
     {0x1.670f2944f8000p+196, 0x1.935bd3c25ed48p+156},
     {0x1.3593882c60000p+196, 0x1.058bd40f29793p+157},
     {0x1.a85b1ef3b8000p+196, 0x1.41dff20ca198ap+158},
     {0x1.4c42a2dd40000p+196, 0x1.8f504f961ea64p+158},
     {0x1.a069d49800000p+192, 0x1.a5268e4d6841fp+158},
     {0x1.a03c37c1c0000p+195, 0x1.973f5db43c9e0p+158},
     {0x1.498d401938000p+196, 0x1.5498c046c404fp+155},
     {0x1.3fe3f82c00000p+194, 0x1.1647a0b44bc87p+156},
     {0x1.8156f28b08000p+196, 0x1.31f62462d7487p+157},
     {0x1.bfa6cefcd8000p+196, 0x1.b738e1a2f6fdep+156},
     {0x1.008460b0a0000p+195, 0x1.e92c097a43a2ap+158},
     {0x1.c7f1e43180000p+192, 0x1.924e0ba91b94dp+158},
     },
};

static_assert(sizeof(qPi_2) / sizeof(qPi_2[0]) == 21U - qPi_2_first_num_moduli);
static_assert(sizeof(qPi_2_complex) / sizeof(qPi_2_complex[0]) == 21U - qPi_2_complex_first_num_moduli);
static_assert(qPi_2_first_num_moduli <= threshold<Backend::FP8, false>::P_is_double + 1U);
static_assert(qPi_2_complex_first_num_moduli <= threshold<Backend::FP8, true>::P_is_double + 1U);
static_assert(P[threshold<Backend::FP8, false>::P_is_double - 2].y == 0.0);
static_assert(P_complex[threshold<Backend::FP8, true>::P_is_double - 2].y == 0.0);

} // namespace FP8

namespace INT8 {
template <unsigned NUM_MODULI, unsigned IDX, bool COMPLEX = false>
inline constexpr double qPi_double_v =
    (COMPLEX)
        ? INT8::qPi_1_complex[NUM_MODULI - 2][IDX]
        : INT8::qPi_1[NUM_MODULI - 2][IDX];

template <unsigned NUM_MODULI, unsigned IDX, bool COMPLEX = false>
inline constexpr double2 qPi_double2_v = [] {
    constexpr unsigned first = COMPLEX ? qPi_2_complex_first_num_moduli : qPi_2_first_num_moduli;
    static_assert(NUM_MODULI >= first && NUM_MODULI <= 20U, "CRT row is outside the double2 table");
    static_assert(IDX < NUM_MODULI, "CRT coefficient index is outside this row");
    if constexpr (COMPLEX) {
        return qPi_2_complex[NUM_MODULI - first][IDX];
    } else {
        return qPi_2[NUM_MODULI - first][IDX];
    }
}();
} // namespace INT8

namespace FP8 {
template <unsigned NUM_MODULI, unsigned IDX, bool COMPLEX = false>
inline constexpr double qPi_double_v =
    (COMPLEX)
        ? FP8::qPi_1_complex[NUM_MODULI - 2][IDX]
        : FP8::qPi_1[NUM_MODULI - 2][IDX];

template <unsigned NUM_MODULI, unsigned IDX, bool COMPLEX = false>
inline constexpr double2 qPi_double2_v = [] {
    constexpr unsigned first = COMPLEX ? qPi_2_complex_first_num_moduli : qPi_2_first_num_moduli;
    static_assert(NUM_MODULI >= first && NUM_MODULI <= 20U, "CRT row is outside the double2 table");
    static_assert(IDX < NUM_MODULI, "CRT coefficient index is outside this row");
    if constexpr (COMPLEX) {
        return qPi_2_complex[NUM_MODULI - first][IDX];
    } else {
        return qPi_2[NUM_MODULI - first][IDX];
    }
}();
} // namespace FP8

template <Backend BACKEND, unsigned NUM_MODULI, unsigned IDX, bool COMPLEX = false>
__device__ __forceinline__ constexpr double qPi_double() {
    static_assert(NUM_MODULI >= 2U && NUM_MODULI <= 20U, "CRT row must be in [2, 20]");
    static_assert(IDX < NUM_MODULI, "CRT coefficient index is outside this row");
    if constexpr (BACKEND == Backend::INT8) {
        return INT8::qPi_double_v<NUM_MODULI, IDX, COMPLEX>;
    } else {
        return FP8::qPi_double_v<NUM_MODULI, IDX, COMPLEX>;
    }
}

template <Backend BACKEND, unsigned NUM_MODULI, unsigned IDX, bool COMPLEX = false>
__device__ __forceinline__ constexpr double2 qPi_double2() {
    if constexpr (BACKEND == Backend::INT8) {
        return INT8::qPi_double2_v<NUM_MODULI, IDX, COMPLEX>;
    } else {
        return FP8::qPi_double2_v<NUM_MODULI, IDX, COMPLEX>;
    }
}

} // namespace gemmul8::common::table
