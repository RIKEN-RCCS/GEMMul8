#pragma once
#include "fp8_block_size.hpp"

namespace gemmul8::oz2::core {

namespace fp8_rankk_cost_detail {

inline double triangle(size_t n) { return double(n) * (double(n) + 1.0) * 0.5; }

inline double gemm_area(size_t n, size_t inner, cublasFillMode_t uplo) {
    const size_t q = (n - 1) / inner, r = n - q * inner;
    const double h = double(common::padding(n)), b = double(inner);
    return double(q) * b * (uplo == CUBLAS_FILL_MODE_UPPER ? double(n) : h) -
           b * b * double(q) * (double(q) - 1.0) * 0.5 + (h - double(q) * b) * r;
}

inline double mod_area(size_t n, bool grouped, cublasFillMode_t uplo) {
    const double h = double(common::padding(n));
    if (grouped) {
        return triangle(n) + (uplo == CUBLAS_FILL_MODE_LOWER ? double(n) * (h - n) : 0.0);
    }
    const double q = double(n / 4), r = double(n % 4);
    return uplo == CUBLAS_FILL_MODE_UPPER
               ? 8.0 * q * (q + 1.0) + 4.0 * r * (q + 1.0)
               : double(n) * h - 8.0 * q * (q - 1.0) - 4.0 * r * q;
}

inline double opposite_area(size_t n) {
    const double r = double(n % 32);
    return (double(n) * n + double(n / 32) * 1024.0 + r * r) * 0.5;
}

template <Func F, typename TA, typename TB, typename TC, unsigned N, bool SharedSegments>
inline double diagonal(size_t n, size_t k, size_t kB, bool fastmode,
                       cublasFillMode_t uplo, size_t inner) {
    constexpr bool Complex = common::isComplex<TA>;
    constexpr bool Herk    = F == Func::herk;
    constexpr bool Single  = F == Func::syrk || Herk;
    constexpr bool Rank2   = F == Func::syr2k || F == Func::her2k;
    const double leaves = double(div_up(k, kB)), valid = triangle(n);
    if constexpr (Rank2) {
        return fp8_gemm_block_bytes<TA, TB, TC, N, SharedSegments>(n, n, k, {n, n, kB}, fastmode) -
               2.0 * sizeof(TC) * (double(n) * n - valid) * leaves;
    } else {
        constexpr auto &d             = fp8_k_cost_detail::descriptor<N, Complex, Herk>;
        constexpr size_t products     = d.base_products;
        constexpr size_t inputs       = common::fp8_plan::planes(N, 0, N, Complex, Herk);
        constexpr size_t norm         = Complex ? 3U : 1U;
        using Groups                  = common::crt_group_plan<Backend::FP8, N, Complex, Herk>;
        constexpr unsigned last_group = Groups::groups.first[Groups::count - 1];
        const auto cost               = fp8_k_cost<N, Complex, false, Herk>(k, kB);
        if (cost.gemm_planes == std::numeric_limits<size_t>::max() ||
            cost.reduction_planes == std::numeric_limits<size_t>::max()) {
            return std::numeric_limits<double>::infinity();
        }
        const size_t h      = common::padding(n);
        const double packed = double(common::padding(k)) * (Single ? double(h) : double(h) + n);
        const double source = double(k) * n * (Single ? sizeof(TA) : sizeof(TA) + sizeof(TB));
        double bytes        = (fastmode ? 2.0 : 3.0) * source;
        if constexpr (Single) {
            bytes += (double(inputs) + 2.0 * products) * packed;
        } else {
            bytes += 2.0 * inputs * packed;
        }
        if (!fastmode) {
            bytes += (Single ? 3.0 : 2.0) * norm * packed;
        }

        const double raw_area = Herk ? double(h) * n : gemm_area(n, inner, uplo);
        bytes += (8.0 * (double(cost.gemm_planes) + double(cost.reduction_planes)) -
                  4.0 * products * leaves) *
                 raw_area;
        const bool grouped = common::use_grouped_crt(h * n);
        if constexpr (Herk) {
            constexpr size_t cross = inputs / 2;
            constexpr auto last    = d.moduli[N - 1];
            const double opposite  = opposite_area(n);
            if (grouped) {
                bytes += (4.0 * (double(products) * valid + double(cross) * opposite) +
                          16.0 * Groups::count * valid) *
                         leaves;
            } else {
                static_assert(last.parts == 1U, "HERK plain CRT tail requires a 2M last modulus");
                bytes += (4.0 * ((double(cross) - last.products) * double(n) * n +
                                 (double(products) - cross) * valid +
                                 double(last.products) * (valid + opposite)) +
                          2.0 * sizeof(common::mid_t<Backend::FP8, true>) * (N - 1) * valid) *
                         leaves;
            }
        } else {
            constexpr bool FusedGroup = F != Func::herkx;
            const size_t tail         = grouped
                                            ? (FusedGroup ? common::fp8_plan::planes(N, last_group, N - last_group, Complex) : 0U)
                                            : common::fp8_plan::planes(N, N - 1, 1, Complex);
            const size_t mid_bytes    = grouped
                                            ? sizeof(uint32_t) * (Complex ? 2U : 1U) * (Groups::count - unsigned(FusedGroup))
                                            : sizeof(common::mid_t<Backend::FP8, Complex>) * (N - 1);
            const double reduced      = mod_area(n, grouped, uplo);
            bytes += (4.0 * (double(products - tail) * reduced + double(tail) * valid) +
                      double(mid_bytes) * (reduced + valid)) *
                     leaves;
        }
        bytes += 2.0 * sizeof(TC) * valid * leaves;
        if (!fastmode) {
            const double scanned = Single ? double(n) * n : 2.0 * valid;
            bytes += 4.0 * norm * (gemm_area(n, inner, uplo) + scanned) * leaves;
        }
        return bytes;
    }
}

} // namespace fp8_rankk_cost_detail

template <Func F, typename TA, typename TB, typename TC, unsigned N, bool SharedSegments>
inline double fp8_rankk_block_bytes(size_t n, size_t k, const BlockSize2D &b,
                                    bool fastmode, cublasFillMode_t uplo, size_t inner_nB) {
    const size_t full = n / b.xB, tail = n % b.xB;
    double bytes = double(full) * fp8_rankk_cost_detail::diagonal<F, TA, TB, TC, N, SharedSegments>(
                                      b.xB, k, b.yB, fastmode, uplo, inner_nB);
    if (tail) {
        bytes += fp8_rankk_cost_detail::diagonal<F, TA, TB, TC, N, SharedSegments>(
            tail, k, b.yB, fastmode, uplo, inner_nB);
    }
    const auto pair = [&](size_t m, size_t columns) {
        double result = fp8_gemm_block_bytes<TA, TB, TC, N, SharedSegments>(
            m, columns, k, {m, columns, b.yB}, fastmode);
        if constexpr (F == Func::syr2k || F == Func::her2k)
            result += fp8_gemm_block_bytes<TB, TA, TC, N, SharedSegments>(
                m, columns, k, {m, columns, b.yB}, fastmode);
        return result;
    };
    if (full > 1) {
        bytes += double(full) * double(full - 1) * 0.5 * pair(b.xB, b.xB);
    }
    if (tail && full) {
        bytes += double(full) * (uplo == CUBLAS_FILL_MODE_UPPER
                                     ? pair(b.xB, tail)
                                     : pair(tail, b.xB));
    }
    return bytes;
}

template <Func F, typename TA, typename TB, typename TC, unsigned N, class Fits>
inline BlockSize2D find_block_size_rankk_fp8(size_t n, size_t k, bool fastmode,
                                             cublasFillMode_t uplo, int &arch, Fits &&fits) {
    constexpr bool Rank2  = F == Func::syr2k || F == Func::her2k;
    const auto current    = find_block_size_rankk<Backend::FP8, N>(n, k, Rank2, fits);
    constexpr size_t safe = F == Func::herk
                                ? std::min(fp8_k_safe_v<N, true, true>, fp8_k_safe_v<N, true>)
                                : fp8_k_safe_v<N, common::isComplex<TA>>;
    if (!current || current.yB <= safe) return current;
    const auto candidate = find_block_size_rankk<Backend::FP8, N>(n, k, Rank2, fits, safe);
    if (!candidate ||
        performance_class(candidate.xB, BLOCK_TARGET_MN) < performance_class(current.xB, BLOCK_TARGET_MN)) {
        return current;
    }
    const size_t inner = common::block_size<F>(arch, CUBLAS_FILL_MODE_FULL, CUBLAS_FILL_MODE_FULL);
    double old_bytes   = fp8_rankk_block_bytes<F, TA, TB, TC, N, common::isCUDA>(
        n, k, current, fastmode, uplo, inner);
    if constexpr (common::isCUDA && common::isComplex<TA> && !common::fp8_plan::mixed_complex(N)) {
        old_bytes = std::min(old_bytes,
                             fp8_rankk_block_bytes<F, TA, TB, TC, N, false>(n, k, current, fastmode, uplo, inner));
    }
    const double new_bytes = fp8_rankk_block_bytes<F, TA, TB, TC, N, common::isCUDA>(
        n, k, candidate, fastmode, uplo, inner);
    return new_bytes < old_bytes ? candidate : current;
}

} // namespace gemmul8::oz2::core
