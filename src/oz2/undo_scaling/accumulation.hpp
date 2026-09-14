#pragma once
#include "../common/common.hpp"
#include "../common/table.hpp"

namespace gemmul8::undo_scaling {

//------------------------------
// real accumulation, P fits in double
//------------------------------
template <Backend BACKEND, unsigned NUM_MODULI>
__device__ __forceinline__ double init_real_double(
    const common::mid_t<BACKEND, false> *const __restrict__ C_mid,
    const size_t incC_mid //
) {
    constexpr double q = common::table::qPi_double<BACKEND, NUM_MODULI, 0>();
    using Mid          = common::mid_t<BACKEND, false>;
    const double c     = common::Tcast<Mid, double>(C_mid[0]);
    return q * c;
}
template <Backend BACKEND, unsigned NUM_MODULI, unsigned IDX = 1>
__device__ __forceinline__ void accumulate_real_double(
    double &acc,
    const common::mid_t<BACKEND, false> *const __restrict__ C_mid,
    const size_t incC_mid //
) {
    if constexpr (IDX < NUM_MODULI) {
        constexpr double q = common::table::qPi_double<BACKEND, NUM_MODULI, IDX>();
        using Mid          = common::mid_t<BACKEND, false>;
        const double c     = common::Tcast<Mid, double>(C_mid[IDX * incC_mid]);

        acc = fma(q, c, acc);

        accumulate_real_double<BACKEND, NUM_MODULI, IDX + 1>(
            acc, C_mid, incC_mid);
    }
}

//------------------------------
// real accumulation, P represented by double2
// acc.x: error-free part
// acc.y: residual part
//------------------------------
template <Backend BACKEND, unsigned NUM_MODULI>
__device__ __forceinline__ double2 init_real_double2(
    const common::mid_t<BACKEND, false> *const __restrict__ C_mid,
    const size_t incC_mid //
) {
    constexpr double2 q = common::table::qPi_double2<BACKEND, NUM_MODULI, 0>();
    using Mid           = common::mid_t<BACKEND, false>;
    const double c      = common::Tcast<Mid, double>(C_mid[0]);
    double2 acc;
    acc.x = q.x * c;
    acc.y = q.y * c;
    return acc;
}
template <Backend BACKEND, unsigned NUM_MODULI, unsigned IDX = 1>
__device__ __forceinline__ void accumulate_real_double2(
    double2 &acc,
    const common::mid_t<BACKEND, false> *const __restrict__ C_mid,
    const size_t incC_mid //
) {
    if constexpr (IDX < NUM_MODULI) {
        constexpr double2 q = common::table::qPi_double2<BACKEND, NUM_MODULI, IDX>();
        using Mid           = common::mid_t<BACKEND, false>;
        const double c      = common::Tcast<Mid, double>(C_mid[IDX * incC_mid]);

        acc.x = fma(q.x, c, acc.x);
        acc.y = fma(q.y, c, acc.y);

        accumulate_real_double2<BACKEND, NUM_MODULI, IDX + 1>(
            acc, C_mid, incC_mid);
    }
}

//------------------------------
// complex accumulation, P fits in double
// acc.x: real part
// acc.y: imaginary part
//------------------------------
template <Backend BACKEND, unsigned NUM_MODULI>
__device__ __forceinline__ double2 init_complex_double(
    const common::mid_t<BACKEND, true> *const __restrict__ C_mid,
    const size_t incC_mid //
) {
    constexpr double q = common::table::qPi_double<BACKEND, NUM_MODULI, 0, true>();
    using Mid          = common::mid_t<BACKEND, true>;
    const double2 c    = common::Tcast<Mid, double2>(C_mid[0]);
    double2 acc;
    acc.x = q * c.x;
    acc.y = q * c.y;
    return acc;
}
template <Backend BACKEND, unsigned NUM_MODULI, unsigned IDX = 1>
__device__ __forceinline__ void accumulate_complex_double(
    double2 &acc,
    const common::mid_t<BACKEND, true> *const __restrict__ C_mid,
    const size_t incC_mid //
) {
    if constexpr (IDX < NUM_MODULI) {
        constexpr double q = common::table::qPi_double<BACKEND, NUM_MODULI, IDX, true>();
        using Mid          = common::mid_t<BACKEND, true>;
        const double2 c    = common::Tcast<Mid, double2>(C_mid[IDX * incC_mid]);

        acc.x = fma(q, c.x, acc.x);
        acc.y = fma(q, c.y, acc.y);

        accumulate_complex_double<BACKEND, NUM_MODULI, IDX + 1>(
            acc, C_mid, incC_mid);
    }
}

//------------------------------
// complex accumulation, P represented by double2
//
// hi.x: real error-free part
// hi.y: imag error-free part
// lo.x: real residual part
// lo.y: imag residual part
//------------------------------
template <Backend BACKEND, unsigned NUM_MODULI>
__device__ __forceinline__ common::double2x2_t init_complex_double2(
    const common::mid_t<BACKEND, true> *const __restrict__ C_mid,
    const size_t incC_mid //
) {
    constexpr double2 q = common::table::qPi_double2<BACKEND, NUM_MODULI, 0, true>();
    using Mid           = common::mid_t<BACKEND, true>;
    const double2 c     = common::Tcast<Mid, double2>(C_mid[0]);

    common::double2x2_t acc;

    acc.hi.x = q.x * c.x;
    acc.hi.y = q.x * c.y;
    acc.lo.x = q.y * c.x;
    acc.lo.y = q.y * c.y;

    return acc;
}
template <Backend BACKEND, unsigned NUM_MODULI, unsigned IDX = 1>
__device__ __forceinline__ void accumulate_complex_double2(
    common::double2x2_t &acc,
    const common::mid_t<BACKEND, true> *const __restrict__ C_mid,
    const size_t incC_mid //
) {
    if constexpr (IDX < NUM_MODULI) {
        constexpr double2 q = common::table::qPi_double2<BACKEND, NUM_MODULI, IDX, true>();
        using Mid           = common::mid_t<BACKEND, true>;
        const double2 c     = common::Tcast<Mid, double2>(C_mid[IDX * incC_mid]);

        acc.hi.x = fma(q.x, c.x, acc.hi.x);
        acc.hi.y = fma(q.x, c.y, acc.hi.y);
        acc.lo.x = fma(q.y, c.x, acc.lo.x);
        acc.lo.y = fma(q.y, c.y, acc.lo.y);

        accumulate_complex_double2<BACKEND, NUM_MODULI, IDX + 1>(
            acc, C_mid, incC_mid);
    }
}

} // namespace gemmul8::undo_scaling
