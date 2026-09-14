#pragma once
#include "accumulation.hpp"
#include "integer_accumulation.hpp"

namespace gemmul8::undo_scaling {

//------------------------------
// real, P fits in double
//------------------------------
template <typename T>
__device__ __forceinline__ T final_reduction_real(
    const double acc,
    const double P,
    const double invP //
) {
    const double quot    = rint(invP * acc);  // round(acc/P)
    const double crt_ans = fma(P, quot, acc); // mod(acc,P)
    return static_cast<T>(crt_ans);
}

//------------------------------
// real, P represented by double2
// acc.x: error-free part
// acc.y: residual part
//------------------------------
template <typename T>
__device__ __forceinline__ T final_reduction_real(
    const double2 acc,
    const double2 P,
    const double invP //
) {
    const double quot    = rint(invP * acc.x);                            // round(acc/P)
    const double crt_ans = fma(P.y, quot, fma(P.x, quot, acc.x) + acc.y); // mod(acc,P)
    return static_cast<T>(crt_ans);
}

//------------------------------
// complex, P fits in double
// acc.x: real part
// acc.y: imaginary part
//------------------------------
template <typename T>
__device__ __forceinline__ T final_reduction_complex(
    const double2 acc,
    const double P,
    const double invP //
) {
    using U = common::underlying_t<T>;

    const double qx = rint(invP * acc.x); // round(acc/P) real part
    const double qy = rint(invP * acc.y); // round(acc/P) complex part

    const double rx = fma(P, qx, acc.x); // mod(acc,P) real part
    const double ry = fma(P, qy, acc.y); // mod(acc,P) complex part

    T out;
    out.x = static_cast<U>(rx);
    out.y = static_cast<U>(ry);
    return out;
}

//------------------------------
// complex, P represented by double2
// hi.x: real error-free part
// hi.y: imag error-free part
// lo.x: real residual part
// lo.y: imag residual part
//------------------------------
template <typename T>
__device__ __forceinline__ T final_reduction_complex(
    const common::double2x2_t acc,
    const double2 P,
    const double invP //
) {
    using U = common::underlying_t<T>;

    const double qx = rint(invP * acc.hi.x); // round(acc/P) real part
    const double qy = rint(invP * acc.hi.y); // round(acc/P) complex part

    const double rx = fma(P.y, qx, fma(P.x, qx, acc.hi.x) + acc.lo.x); // mod(acc,P) real part
    const double ry = fma(P.y, qy, fma(P.x, qy, acc.hi.y) + acc.lo.y); // mod(acc,P) complex part

    T out;
    out.x = static_cast<U>(rx);
    out.y = static_cast<U>(ry);
    return out;
}

//------------------------------
// C_mid{0},...,C_mid{NUM_MODULI-1} -> AB
//------------------------------
template <typename T, Backend BACKEND, unsigned NUM_MODULI, typename TP>
__device__ __forceinline__ T reconstruct_from_crt(
    const common::mid_t<BACKEND, common::isComplex<T>> *const __restrict__ C_mid,
    const size_t incC_mid,
    const TP P,
    const double invP //
) {
    if constexpr (!common::isComplex<T>) {
        if constexpr (std::is_same_v<TP, double>) {

            // double C with small NUMMODULI or float C
            double acc = init_real_double<BACKEND, NUM_MODULI>(C_mid, incC_mid);
            accumulate_real_double<BACKEND, NUM_MODULI>(acc, C_mid, incC_mid);
            return final_reduction_real<T>(acc, P, invP);

        } else {

            if constexpr (use_integer_crt_high) {
                return final_reduction_real<T>(accumulate_real_integer_high<BACKEND, NUM_MODULI>(C_mid, incC_mid), P, invP);
            }

            // double C with large NUM_MODULI
            double2 acc = init_real_double2<BACKEND, NUM_MODULI>(C_mid, incC_mid);
            accumulate_real_double2<BACKEND, NUM_MODULI>(acc, C_mid, incC_mid);
            return final_reduction_real<T>(acc, P, invP);
        }
    } else {
        if constexpr (std::is_same_v<TP, double>) {

            // cuDoubleComplex C with small NUM_MODULI or cuFloatComplex C
            double2 acc = init_complex_double<BACKEND, NUM_MODULI>(C_mid, incC_mid);
            accumulate_complex_double<BACKEND, NUM_MODULI>(acc, C_mid, incC_mid);
            return final_reduction_complex<T>(acc, P, invP);

        } else {

            if constexpr (use_integer_crt_high) {
                return final_reduction_complex<T>(accumulate_complex_integer_high<BACKEND, NUM_MODULI>(C_mid, incC_mid), P, invP);
            }

            // cuDoubleComplex C with large NUM_MODULI
            common::double2x2_t acc = init_complex_double2<BACKEND, NUM_MODULI>(C_mid, incC_mid);
            accumulate_complex_double2<BACKEND, NUM_MODULI>(acc, C_mid, incC_mid);
            return final_reduction_complex<T>(acc, P, invP);
        }
    }
}

} // namespace gemmul8::undo_scaling
