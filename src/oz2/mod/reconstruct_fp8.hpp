#pragma once
#include "mod_core.hpp"

namespace gemmul8::mod {

// NONE  : C0=A0*B1, C1=A1*B0, C2=A1*B1; P=BASE^2.
// Others: C0=A0*B0, C1=A1*B1, C2=A2*B2.
template <int32_t P, common::table::KaratsubaType TYPE, int32_t SQUARE_BASE = 0>
struct fp8_reconstruction_traits {
    using KT = common::table::KaratsubaType;
    static_assert(P >= 3);
    static_assert(TYPE == KT::NONE || TYPE == KT::BASE49 ||
                  TYPE == KT::BASE33_SUM || TYPE == KT::BASE32_SUM ||
                  TYPE == KT::BASE32_DIFF || TYPE == KT::BASE31_DIFF);

    static constexpr bool square  = TYPE == KT::NONE;
    static constexpr bool diff    = TYPE == KT::BASE32_DIFF || TYPE == KT::BASE31_DIFF;
    static constexpr int32_t base = square                    ? SQUARE_BASE
                                    : TYPE == KT::BASE49      ? 49
                                    : TYPE == KT::BASE33_SUM  ? 33
                                    : TYPE == KT::BASE31_DIFF ? 31
                                                              : 32;
    static_assert(base > 0 && base <= 49);
    static_assert(!square || base * base == P);

    static constexpr int32_t s      = diff ? -1 : 1;
    static constexpr int32_t raw_w0 = square ? base : base * (base - s);
    static constexpr int32_t w1     = square ? base : 1 - s * base;
    static constexpr int32_t w2     = square ? 1 : s * base;

    static constexpr int32_t w0_mod = raw_w0 % P;
    static constexpr int32_t w0     = w0_mod > P / 2 ? w0_mod - P : w0_mod;
    static constexpr int64_t aw0    = w0 < 0 ? -int64_t(w0) : int64_t(w0);
    static constexpr int64_t aw1    = w1 < 0 ? -int64_t(w1) : int64_t(w1);
    static constexpr int64_t aw2    = w2 < 0 ? -int64_t(w2) : int64_t(w2);

    static constexpr int64_t product_bound = int64_t(1) << 24;
    static constexpr int64_t direct_bound  = (aw0 + aw1 + aw2) * product_bound;
    static constexpr int64_t safe_limit    = (int64_t(1) << 31) - P;
    static constexpr bool reduce_c0        = direct_bound >= safe_limit;

    static constexpr int64_t combined_bound = reduce_c0
                                                  ? aw0 * (int64_t(2) * P) + (aw1 + aw2) * product_bound
                                                  : direct_bound;
    static_assert(combined_bound < safe_limit, "FP8 reconstruction may overflow int32_t");
};

template <unsigned IDX, bool COMPLEX = false>
using fp8_reconstruction =
    fp8_reconstruction_traits<common::table::moduli<Backend::FP8, IDX, COMPLEX>,
                              common::table::kara_type<IDX, COMPLEX>,
                              common::table::sqrt_moduli<IDX, COMPLEX>>;

template <unsigned IDX, bool nowrap = false, bool COMPLEX = false>
__device__ __forceinline__ int32_t mod_f32x3_2_i32(const float C0, const float C1, const float C2) {
    using R          = fp8_reconstruction<IDX, COMPLEX>;
    int32_t c0       = __float2int_rn(C0);
    const int32_t c1 = __float2int_rn(C1);
    const int32_t c2 = __float2int_rn(C2);

    if constexpr (R::reduce_c0) {
        c0 = mod_small_nowrap<Backend::FP8, IDX, COMPLEX>(c0);
    }

    const int32_t t = R::square ? R::base * (c0 + c1) + c2
                                : R::w0 * c0 + R::w1 * c1 + R::w2 * c2;
    if constexpr (nowrap) {
        return mod_small_nowrap<Backend::FP8, IDX, COMPLEX>(t);
    } else {
        return mod_small<Backend::FP8, IDX, COMPLEX>(t);
    }
}

} // namespace gemmul8::mod
