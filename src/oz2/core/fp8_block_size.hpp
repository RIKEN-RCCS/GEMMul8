#pragma once
#include "block_size.hpp"
#include "fp8_k_cost.hpp"
#include "../common/crt_group_plan.hpp"

namespace gemmul8::oz2::core {

template <typename TA, typename TB, typename TC, unsigned N, bool SharedSegments>
inline double fp8_gemm_block_bytes(
    size_t m, size_t n, size_t k,
    const BlockSize3D &b, bool fastmode //
) {
    constexpr bool Complex        = common::isComplex<TA>;
    constexpr size_t planes       = common::fp8_plan::planes(N, 0, N, Complex);
    constexpr size_t norm_planes  = Complex ? 3U : 1U;
    using Groups                  = common::crt_group_plan<Backend::FP8, N, Complex>;
    constexpr unsigned last_group = Groups::groups.first[Groups::count - 1];
    constexpr size_t grouped_tail = common::fp8_plan::planes(N, last_group, N - last_group, Complex);
    constexpr size_t plain_tail   = common::fp8_plan::planes(N, N - 1, 1, Complex);

    const size_t mi_count = div_up(m, b.mB), nj_count = div_up(n, b.nB);
    const size_t k_count = div_up(k, b.kB);
    const auto cost      = fp8_k_cost<N, Complex, SharedSegments>(k, b.kB);
    if (cost.gemm_planes == std::numeric_limits<size_t>::max() ||
        cost.reduction_planes == std::numeric_limits<size_t>::max()) {
        return std::numeric_limits<double>::infinity();
    }

    const double packed = double(common::padding(k)) *
                          (double(common::padding(m)) * nj_count + double(n) * mi_count);
    const double source = double(sizeof(TA)) * m * double(k) * nj_count +
                          double(sizeof(TB)) * n * double(k) * mi_count;
    double bytes        = (fastmode ? 2.0 : 3.0) * source + 2.0 * planes * packed;
    if (!fastmode) bytes += 2.0 * norm_planes * packed;

    const size_t rows[2]       = {b.mB, m % b.mB};
    const size_t cols[2]       = {b.nB, n % b.nB};
    const size_t row_counts[2] = {m / b.mB, size_t(rows[1] != 0)};
    const size_t col_counts[2] = {n / b.nB, size_t(cols[1] != 0)};
    for (unsigned i = 0; i < 2; ++i) {
        if (!row_counts[i]) continue;
        const size_t padded_rows = common::padding(rows[i]);
        for (unsigned j = 0; j < 2; ++j) {
            if (!col_counts[j]) continue;
            const double copies     = double(row_counts[i]) * col_counts[j];
            const double area       = double(padded_rows) * cols[j];
            const double valid_area = double(rows[i]) * cols[j];
            const bool grouped      = common::use_grouped_crt(padded_rows * cols[j]);
            const size_t tail       = grouped ? grouped_tail : plain_tail;
            const size_t mid_bytes  = grouped
                                          ? sizeof(uint32_t) * (Complex ? 2U : 1U) * (Groups::count - 1)
                                          : sizeof(common::mid_t<Backend::FP8, Complex>) * (N - 1);

            double output = 8.0 * (double(cost.gemm_planes) + double(cost.reduction_planes)) * area;
            output -= 4.0 * tail * (area - valid_area) * k_count;
            output += double(mid_bytes) * (area + valid_area) * k_count;
            output += 2.0 * sizeof(TC) * valid_area * k_count;
            if (!fastmode) {
                output += 4.0 * norm_planes * (area + 2.0 * valid_area) * k_count;
            }
            bytes += copies * output;
        }
    }
    return bytes;
}

template <typename TA, typename TB, typename TC, unsigned N>
inline bool prefer_fp8_safe_block(
    size_t m, size_t n, size_t k,
    const BlockSize3D &current, const BlockSize3D &candidate, bool fastmode //
) {
    if (!candidate) return false;
    const auto mn_class = [](const BlockSize3D &b) {
        return std::min(performance_class(b.mB, BLOCK_TARGET_MN),
                        performance_class(b.nB, BLOCK_TARGET_MN));
    };
    if (mn_class(candidate) < mn_class(current)) return false;

    double old_bytes = fp8_gemm_block_bytes<TA, TB, TC, N, common::isCUDA>(m, n, k, current, fastmode);
    if constexpr (common::isCUDA && common::isComplex<TA> && !common::fp8_plan::mixed_complex(N)) {
        old_bytes = std::min(old_bytes,
                             fp8_gemm_block_bytes<TA, TB, TC, N, false>(m, n, k, current, fastmode));
    }
    const double new_bytes = fp8_gemm_block_bytes<TA, TB, TC, N, common::isCUDA>(m, n, k, candidate, fastmode);
    return new_bytes < old_bytes;
}

template <typename TA, typename TB, typename TC, unsigned N, class WorkSize>
inline BlockSize3D find_block_size_gemm_fp8(
    size_t m, size_t n, size_t k, size_t limit, bool fastmode, WorkSize &&worksize //
) {
    if (m == 0 || n == 0 || k == 0) return {};
    if (worksize(m, n, k) <= limit) return {m, n, k};

    const auto ms         = make_block_states(m, BLOCK_TARGET_MN);
    const auto ns         = make_block_states(n, BLOCK_TARGET_MN);
    const auto ks         = make_block_states(k, BLOCK_TARGET_K);
    const auto current    = find_block_size_gemm_states(ms, ns, ks, limit, worksize);
    constexpr size_t safe = fp8_k_safe_v<N, common::isComplex<TA>>;
    if (!current || current.kB <= safe) return current;

    const auto begin = std::find_if(ks.begin(), ks.end(),
                                    [](const BlockState &s) { return s.block <= safe; });
    const std::vector<BlockState> safe_ks(begin, ks.end());
    const auto candidate = find_block_size_gemm_states(ms, ns, safe_ks, limit, worksize);
    return prefer_fp8_safe_block<TA, TB, TC, N>(m, n, k, current, candidate, fastmode)
               ? candidate
               : current;
}

} // namespace gemmul8::oz2::core
