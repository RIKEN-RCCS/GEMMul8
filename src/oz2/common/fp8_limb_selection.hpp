#pragma once
#include "common.hpp"
#include "table.hpp"

namespace gemmul8::common::make_f8 {

struct limb_mask {
    uint32_t bits = 0;
    int32_t first = 0;
};
struct limb_selection_masks {
    limb_mask stage[3][3]{};
    bool valid = false;
};

template <int32_t P, table::KaratsubaType TYPE, uint32_t MAX_ABS>
struct limb_selection {
    static_assert(P > 0 && MAX_ABS <= 32U);
    static constexpr bool base49       = TYPE == table::KaratsubaType::BASE49;
    static constexpr bool sum32        = TYPE == table::KaratsubaType::BASE32_SUM;
    static constexpr bool base33       = TYPE == table::KaratsubaType::BASE33_SUM;
    static constexpr int32_t base      = base49 ? 48 : 32;
    static constexpr int32_t h         = (P + base / 2) / base;
    static constexpr int32_t t         = P - base * h;
    static constexpr int32_t carry     = P == 853 || P == 793 || P == 797 ? 1 : 0;
    static constexpr int32_t delta0[4] = {
        0, 1, base49 ? 1 - h : (base33 ? -31 : -h + carry), h};
    static constexpr int32_t delta1[4] = {
        0, -base, base49 ? -t - 48 : (base33 ? -41 : -t - 32 * carry), t};

    static constexpr bool representable(int32_t x) {
        const uint32_t a = uint32_t(x < 0 ? -x : x);
        return a <= MAX_ABS && (a & 0x11U) != 0x11U;
    }

    struct sample {
        int32_t f[3];
        unsigned choice;
    };
    static constexpr sample classify(int32_t u) {
        const int32_t q = u / base;
        const int32_t r = u - base * q;
        sample s{
            {q, r, sum32 ? r + q : r - q},
            4U
        };
        for (unsigned k = 0; k < (base33 ? 3U : 4U); ++k) {
            const int32_t a = q + delta0[k];
            const int32_t b = r + delta1[k];
            const int32_t c = sum32 ? a + b : b - a;
            if (representable(a) && representable(b) && representable(c)) {
                s.choice = k;
                break;
            }
        }
        return s;
    }

    static constexpr limb_selection_masks masks = [] {
        limb_selection_masks result;
        uint64_t good[3][3]{}, bad[3][3]{};
        sample samples[P / 2 + 1]{};
        for (int32_t u = 0; u <= P / 2; ++u) {
            const sample s = classify(u);
            if (s.choice > 3U) return result;
            samples[u] = s;
            for (unsigned j = 0; j < 3U; ++j) {
                if (s.f[j] < -16 || s.f[j] > 47) return result;
                const uint64_t bit = uint64_t(1) << unsigned(s.f[j] + 16);
                for (unsigned k = 0; k < 3U && k <= s.choice; ++k) {
                    (s.choice == k ? good[k][j] : bad[k][j]) |= bit;
                }
            }
        }
        constexpr unsigned subsets[] = {0U, 1U, 2U, 4U, 3U, 5U, 6U, 7U};
        for (unsigned k = 0; k < 3U; ++k) {
            uint64_t reject[3]{};
            limb_mask compressed[3]{};
            unsigned usable = 0;
            for (unsigned j = 0; j < 3U; ++j) {
                reject[j]     = bad[k][j] & ~good[k][j];
                uint64_t bits = reject[j];
                int32_t first = -16;
                if (bits != 0) {
                    while ((bits & 1U) == 0) {
                        bits >>= 1;
                        ++first;
                    }
                }
                if ((bits >> 32) == 0) {
                    usable |= 1U << j;
                    compressed[j] = {uint32_t(bits), first};
                }
            }
            bool found = false;
            for (const unsigned subset : subsets) {
                if (subset & ~usable) continue;
                bool covered = true;
                for (const sample s : samples) {
                    if (s.choice <= k) continue;
                    bool hit = false;
                    for (unsigned j = 0; j < 3U; ++j) {
                        if (subset & (1U << j)) {
                            hit |= ((reject[j] >> unsigned(s.f[j] + 16)) & 1U) != 0;
                        }
                    }
                    if (!hit) {
                        covered = false;
                        break;
                    }
                }
                if (covered) {
                    for (unsigned j = 0; j < 3U; ++j) {
                        if (subset & (1U << j)) result.stage[k][j] = compressed[j];
                    }
                    found = true;
                    break;
                }
            }
            if (!found) return limb_selection_masks{};
        }
        result.valid = true;
        return result;
    }();

    static constexpr uint32_t pack(const int32_t (&v)[4]) {
        return uint32_t(uint8_t(v[0])) | (uint32_t(uint8_t(v[1])) << 8) |
               (uint32_t(uint8_t(v[2])) << 16) | (uint32_t(uint8_t(v[3])) << 24);
    }
    static constexpr uint32_t packed0 = pack(delta0);
    static constexpr uint32_t packed1 = pack(delta1);
};

template <uint32_t MASK, int32_t FIRST>
__device__ __forceinline__ uint32_t limb_test(int32_t x) {
    if constexpr (MASK == 0U) return 0U;
    const uint32_t pos = uint32_t(x - FIRST);
#if defined(__CUDACC__) && !defined(__HIPCC__)
    uint32_t r;
    asm("shf.r.clamp.b32 %0, %1, 0, %2;" : "=r"(r) : "n"(MASK), "r"(pos));
    return r & 1U;
#else
    return pos < 32U ? (MASK >> (pos & 31U)) & 1U : 0U;
#endif
}

template <uint32_t PACKED>
__device__ __forceinline__ int32_t limb_offset(uint32_t pos) {
    return int32_t(int8_t(PACKED >> (pos & 24U)));
}

template <typename Plan, unsigned STAGE>
__device__ __forceinline__ uint32_t limb_reject(int32_t q, int32_t r, int32_t d) {
    constexpr auto a = Plan::masks.stage[STAGE][0];
    constexpr auto b = Plan::masks.stage[STAGE][1];
    constexpr auto c = Plan::masks.stage[STAGE][2];
    return limb_test<a.bits, a.first>(q) | limb_test<b.bits, b.first>(r) |
           limb_test<c.bits, c.first>(d);
}

template <typename Plan>
__device__ __forceinline__ void select_limbs(int32_t q, int32_t r,
                                             int32_t &a0, int32_t &a1, int32_t &a2) {
    static_assert(Plan::masks.valid);
    const int32_t d        = Plan::sum32 ? r + q : r - q;
    const uint32_t bad_b   = limb_reject<Plan, 0>(q, r, d);
    const uint32_t bad_bc  = bad_b & limb_reject<Plan, 1>(q, r, d);
    const uint32_t bad_bca = bad_bc & limb_reject<Plan, 2>(q, r, d);
    const uint32_t pos     = (bad_b + bad_bc + bad_bca) << 3;
    a0                     = q + limb_offset<Plan::packed0>(pos);
    const int32_t b        = r + limb_offset<Plan::packed1>(pos);
    if constexpr (Plan::base49 || Plan::base33) {
        a2 = b;
        a1 = b - a0;
    } else {
        a1 = b;
        a2 = Plan::sum32 ? a0 + b : b - a0;
    }
}

} // namespace gemmul8::common::make_f8
