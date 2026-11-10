#pragma once
#include "fp8_structured_block_size.hpp"

namespace gemmul8::oz2::core {

namespace fp8_trtrmm_detail {

template <unsigned N, bool Complex>
inline double raw_bytes(size_t n, size_t inner, bool upperA, bool upperB, double &area) {
    constexpr double planes = common::fp8_plan::planes(N, 0, N, Complex);
    double bytes            = 0.0;
    area                    = 0.0;
    if (n <= fp8_k_safe_v<N, Complex>) {
        const size_t full = n / inner, tail = n % inner;
        area = upperA == upperB
                   ? double(inner) * inner * double(full) * (double(full) + 1.0) / 2.0 +
                         double(full) * inner * tail + double(tail) * tail
                   : double(n) * n;
        return 4.0 * planes * area;
    }
    const auto add = [&](double a, size_t k) {
        const auto cost = fp8_k_cost<N, Complex, false>(k, k);
        bytes += a * (8.0 * (double(cost.gemm_planes) + double(cost.reduction_planes)) - 4.0 * planes);
        area += a;
    };
    if (upperA == upperB) {
        const size_t full = n / inner, tail = n % inner;
        for (size_t d = 0; d < full; ++d)
            add(double(inner) * inner * (full - d), (d + 1) * inner);
        if (tail) {
            for (size_t i = 0; i < full; ++i)
                add(double(inner) * tail, n - i * inner);
            add(double(tail) * tail, tail);
        }
    } else {
        for (size_t s = 0; s < n; s += inner) {
            const size_t width = std::min(inner, n - s);
            if (upperA) add(double(width) * (2.0 * s + width), n - s);
            else add(double(width) * (2.0 * (n - s) - width), s + width);
        }
    }
    return bytes;
}

inline double choose(size_t n, unsigned r) {
    if (n < r) return 0.0;
    double result = 1.0;
    for (unsigned j = 0; j < r; ++j) result = result * double(n - j) / (j + 1U);
    return result;
}

// At most 26 full/tail weak-order patterns replace the outer i,j,p tile loops.
inline void counts(size_t full, bool tail, bool upperA, bool upperB, double (&out)[8][4]) {
    for (unsigned i = 0; i < 3; ++i) {
        for (unsigned j = 0; j < 3; ++j) {
            for (unsigned p = 0; p < 3; ++p) {
                const unsigned ranks = std::max({i, j, p}) + 1;
                if (((1U << i) | (1U << j) | (1U << p)) != (1U << ranks) - 1U) continue;
                if ((upperA ? i > p : i < p) || (upperB ? p > j : p < j)) continue;
                const unsigned kind = (i == p ? 1U : 0U) | (j == p ? 2U : 0U);
                out[0][kind] += choose(full, ranks);
                if (tail) {
                    const unsigned mask = (i + 1 == ranks ? 1U : 0U) |
                                          (j + 1 == ranks ? 2U : 0U) |
                                          (p + 1 == ranks ? 4U : 0U);
                    out[mask][kind] += choose(full, ranks - 1);
                }
            }
        }
    }
}

} // namespace fp8_trtrmm_detail

template <typename TA, typename TB, typename TC, unsigned N>
inline double fp8_trtrmm_call_bytes(
    size_t n, bool upperA, bool upperB, bool unitA, bool unitB,
    size_t inner, bool fastmode //
) {
    constexpr bool Complex       = common::isComplex<TA>;
    constexpr size_t planes      = common::fp8_plan::planes(N, 0, N, Complex);
    constexpr size_t norm_planes = Complex ? 3U : 1U;
    using Groups                 = common::crt_group_plan<Backend::FP8, N, Complex>;
    constexpr unsigned last      = Groups::groups.first[Groups::count - 1];
    const size_t padded          = common::padding(n);
    const double square          = double(n) * n;
    const double triangle        = double(n) * (double(n) + 1.0) / 2.0;
    const double packed          = 2.0 * double(padded) * padded;
    const bool same              = upperA == upperB;
    const bool grouped           = common::use_grouped_crt(padded * padded);
    const double valid           = same ? triangle : square;
    const double reduced         = same ? double(padded) * (double(padded) + (grouped ? 1.0 : 4.0)) / 2.0
                                        : double(padded) * padded;
    const size_t fused           = grouped ? common::fp8_plan::planes(N, last, N - last, Complex)
                                           : common::fp8_plan::planes(N, N - 1, 1, Complex);
    const size_t mid_bytes       = grouped ? sizeof(uint32_t) * (Complex ? 2U : 1U) * (Groups::count - 1)
                                           : sizeof(common::mid_t<Backend::FP8, Complex>) * (N - 1);
    const double source          = sizeof(TA) * (triangle - (unitA ? n : 0)) +
                                   sizeof(TB) * (triangle - (unitB ? n : 0));

    double bytes = (fastmode ? 2.0 : 3.0) * source + planes * (2.0 * packed + 2.0 * triangle);
    if (!fastmode) bytes += norm_planes * (2.0 * packed + 2.0 * triangle);
    double raw_area = 0.0;
    bytes += fp8_trtrmm_detail::raw_bytes<N, Complex>(padded, inner, upperA, upperB, raw_area);
    bytes += 4.0 * (double(planes - fused) * reduced + double(fused) * valid);
    bytes += double(mid_bytes) * (reduced + valid);
    bytes += 2.0 * sizeof(TC) * square;
    if (!fastmode) bytes += 4.0 * norm_planes * (raw_area + 2.0 * valid);
    return bytes;
}

template <typename TA, typename TB, typename TC, unsigned N, bool SharedSegments>
inline double fp8_trtrmm_block_bytes(
    size_t n, size_t b, bool upperA, bool upperB, bool unitA, bool unitB,
    size_t inner_trtrmm, size_t inner_trmm, bool fastmode //
) {
    const size_t full = n / b, tail = n % b;
    double counts[8][4]{};
    fp8_trtrmm_detail::counts(full, tail != 0, upperA, upperB, counts);
    double bytes = 0.0;
    for (unsigned mask = 0; mask < 8; ++mask) {
        const size_t m  = mask & 1 ? tail : b;
        const size_t nn = mask & 2 ? tail : b;
        const size_t k  = mask & 4 ? tail : b;
        for (unsigned kind = 0; kind < 4; ++kind) {
            if (!counts[mask][kind]) continue;
            double value;
            if (kind == 3) value = fp8_trtrmm_call_bytes<TA, TB, TC, N>(m, upperA, upperB, unitA, unitB, inner_trtrmm, fastmode);
            else if (kind == 1) value = fp8_trmm_call_bytes<TA, TB, TC, N>(m, nn, true, upperA, unitA, inner_trmm, fastmode);
            else if (kind == 2) value = fp8_trmm_call_bytes<TB, TA, TC, N>(m, nn, false, upperB, unitB, inner_trmm, fastmode);
            else value = fp8_gemm_block_bytes<TA, TB, TC, N, SharedSegments>(m, nn, k, {m, nn, k}, fastmode);
            bytes += counts[mask][kind] * value;
        }
    }
    if (upperA == upperB) {
        const double empty = (double(n) * n - double(full) * b * b - double(tail) * tail) / 2.0;
        bytes += 2.0 * sizeof(TC) * empty;
    }
    return bytes;
}

template <typename TA, typename TB, typename TC, unsigned N, class Fits, class FitsNoGemm>
inline size_t select_fp8_trtrmm_block(
    size_t n, size_t current, bool upperA, bool upperB, bool unitA, bool unitB,
    size_t inner_trtrmm, size_t inner_trmm, bool fastmode,
    Fits &&fits, FitsNoGemm &&fits_no_gemm //
) {
    constexpr size_t safe = fp8_k_safe_v<N, common::isComplex<TA>>;
    if (!current || current <= safe) return current;
    size_t candidate = 0;
    const size_t b2  = block_size_from_count(n, 2);
    if (upperA == upperB && b2 <= safe && b2 < n && block_count(n, b2) == 2 && fits_no_gemm(b2)) {
        candidate = b2;
    }
    if (!candidate) candidate = find_block_size_trtrmm(n, [&](size_t b) { return b <= safe && fits(b); });
    if (!candidate ||
        performance_class(candidate, BLOCK_TARGET_MN) < performance_class(current, BLOCK_TARGET_MN)) {
        return current;
    }
    double old_bytes = fp8_trtrmm_block_bytes<TA, TB, TC, N, common::isCUDA>(
        n, current, upperA, upperB, unitA, unitB, inner_trtrmm, inner_trmm, fastmode);
    if constexpr (common::isCUDA && common::isComplex<TA> && !common::fp8_plan::mixed_complex(N)) {
        old_bytes = std::min(old_bytes,
                             fp8_trtrmm_block_bytes<TA, TB, TC, N, false>(
                                 n, current, upperA, upperB, unitA, unitB, inner_trtrmm, inner_trmm, fastmode));
    }
    const double new_bytes = fp8_trtrmm_block_bytes<TA, TB, TC, N, common::isCUDA>(
        n, candidate, upperA, upperB, unitA, unitB, inner_trtrmm, inner_trmm, fastmode);
    return new_bytes < old_bytes ? candidate : current;
}

} // namespace gemmul8::oz2::core
