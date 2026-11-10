#pragma once
#include "fp8_block_size.hpp"

namespace gemmul8::oz2::core {

template <typename TTri, typename TFull, typename TC, unsigned N>
inline double fp8_trmm_call_bytes(
    size_t m, size_t n, bool left, bool upper, bool unit,
    size_t inner_block, bool fastmode, bool beta_zero = false //
) {
    constexpr bool Complex       = common::isComplex<TTri>;
    constexpr size_t planes      = common::fp8_plan::planes(N, 0, N, Complex);
    constexpr size_t norm_planes = Complex ? 3U : 1U;
    const size_t s               = left ? m : n;
    const size_t mp = common::padding(m), sp = common::padding(s);
    const BlockSize3D b{m, n, s};
    double bytes = left ? fp8_gemm_block_bytes<TTri, TFull, TC, N, false>(m, n, s, b, fastmode)
                        : fp8_gemm_block_bytes<TFull, TTri, TC, N, false>(m, n, s, b, fastmode);

    const double triangle        = double(s) * (double(s) + 1.0) * 0.5;
    const double source_triangle = triangle - (unit ? double(s) : 0.0);
    bytes -= (fastmode ? 2.0 : 3.0) * sizeof(TTri) * (double(s) * s - source_triangle);
    bytes += (planes + (fastmode ? 0U : norm_planes)) * triangle;
    if (beta_zero) bytes -= double(sizeof(TC)) * m * double(n);

    if (sp <= fp8_k_safe_v<N, Complex>) return bytes;
    const auto whole    = fp8_k_cost<N, Complex, false>(s, s);
    double segmented    = 0.0;
    const size_t extent = left ? mp : n;
    for (size_t at = 0; at < extent; at += inner_block) {
        const size_t width = std::min(inner_block, extent - at);
        const size_t kk    = left ? (upper ? mp - at : at + width)
                                  : (upper ? (at + width == n ? sp : at + width) : sp - at);
        const auto cost    = fp8_k_cost<N, Complex, false>(kk, kk);
        const double area  = double(width) * (left ? n : mp);
        segmented += (double(cost.gemm_planes) + double(cost.reduction_planes)) * area;
    }
    bytes += 8.0 * (segmented - (double(whole.gemm_planes) + double(whole.reduction_planes)) * mp * double(n));
    return bytes;
}

template <typename TStructured, typename TFull, typename TC, unsigned N,
          bool Triangular, bool SharedSegments>
inline double fp8_structured_block_bytes(
    size_t s, size_t f, const BlockSize2D &b,
    bool left, bool upper, bool unit, size_t inner_block, bool fastmode //
) {
    const size_t ss[2] = {b.xB, s % b.xB}, ff[2] = {b.yB, f % b.yB};
    const size_t sc[2] = {s / b.xB, size_t(ss[1] != 0)};
    const size_t fc[2] = {f / b.yB, size_t(ff[1] != 0)};
    const auto dense   = [&](size_t sq, size_t fr, size_t kp) {
        const size_t m = left ? sq : fr, n = left ? fr : sq;
        const BlockSize3D call{m, n, kp};
        return left ? fp8_gemm_block_bytes<TStructured, TFull, TC, N, SharedSegments>(m, n, kp, call, fastmode)
                    : fp8_gemm_block_bytes<TFull, TStructured, TC, N, SharedSegments>(m, n, kp, call, fastmode);
    };
    double bytes = 0.0;
    for (unsigned q = 0; q < 2; ++q) {
        if (!sc[q]) continue;
        for (unsigned r = 0; r < 2; ++r) {
            if (!fc[r]) continue;
            double diagonal;
            if constexpr (Triangular) {
                const size_t m = left ? ss[q] : ff[r], n = left ? ff[r] : ss[q];
                diagonal = fp8_trmm_call_bytes<TStructured, TFull, TC, N>(
                    m, n, left, upper, unit, inner_block, fastmode, true);
            } else {
                diagonal                     = dense(ss[q], ff[r], ss[q]);
                const double strict_triangle = double(ss[q]) * (ss[q] - 1) * 0.5;
                diagonal -= (fastmode ? 1.0 : 2.0) * sizeof(TStructured) * strict_triangle;
            }
            bytes += double(sc[q]) * fc[r] * diagonal;
            for (unsigned p = 0; p < 2; ++p) {
                if (!sc[p]) continue;
                double pairs;
                if (q == p) {
                    pairs = double(sc[q]) * (sc[q] - 1);
                    if constexpr (Triangular) pairs *= 0.5;
                } else {
                    pairs = double(sc[q]) * sc[p];
                    if constexpr (Triangular) {
                        const bool p_after_q = left == upper;
                        if ((p > q) != p_after_q) pairs = 0.0;
                    }
                }
                if (pairs) bytes += pairs * fc[r] * dense(ss[q], ff[r], ss[p]);
            }
        }
    }
    return bytes;
}

template <typename TStructured, typename TFull, typename TC, unsigned N,
          bool Triangular, class Fits, class InnerBlock>
inline BlockSize2D find_block_size_structured_fp8(
    size_t s, size_t f, bool left, bool upper, bool unit, bool fastmode,
    Fits &&fits, InnerBlock &&inner_block //
) {
    const auto current    = find_block_size_structured(s, f, Triangular, fits);
    constexpr size_t safe = fp8_k_safe_v<N, common::isComplex<TStructured>>;
    if (!current || current.xB <= safe) return current;
    const auto candidate = find_block_size_structured(s, f, Triangular, fits, safe);
    if (!candidate) return current;
    const auto mn_class = [](const BlockSize2D &b) {
        return std::min(performance_class(b.xB, BLOCK_TARGET_MN), performance_class(b.yB, BLOCK_TARGET_MN));
    };
    if (mn_class(candidate) < mn_class(current)) return current;

    const size_t inner = Triangular ? size_t(inner_block()) : 1U;
    double old_bytes   = fp8_structured_block_bytes<TStructured, TFull, TC, N, Triangular, common::isCUDA>(
        s, f, current, left, upper, unit, inner, fastmode);
    if constexpr (common::isCUDA && common::isComplex<TStructured> && !common::fp8_plan::mixed_complex(N)) {
        old_bytes = std::min(old_bytes,
                             fp8_structured_block_bytes<TStructured, TFull, TC, N, Triangular, false>(
                                 s, f, current, left, upper, unit, inner, fastmode));
    }
    const double new_bytes = fp8_structured_block_bytes<TStructured, TFull, TC, N, Triangular, common::isCUDA>(
        s, f, candidate, left, upper, unit, inner, fastmode);
    return new_bytes < old_bytes ? candidate : current;
}

} // namespace gemmul8::oz2::core
