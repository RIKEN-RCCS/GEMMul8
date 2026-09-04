// test

#pragma once
#include "../common/common.hpp"

#include <algorithm>
#include <limits>
#include <utility>
#include <vector>

namespace gemmul8::oz2::core {

struct BlockSize2D {
    size_t xB = 0;
    size_t yB = 0;
    explicit operator bool() const noexcept { return xB != 0 && yB != 0; }
};

struct BlockSize3D {
    size_t mB = 0;
    size_t nB = 0;
    size_t kB = 0;
    explicit operator bool() const noexcept { return mB != 0 && nB != 0 && kB != 0; }
};

inline constexpr size_t BLOCK_TARGET_MN     = 8192;
inline constexpr size_t BLOCK_TARGET_K      = 32768;
inline constexpr size_t BLOCK_K_TO_MN_RATIO = 4;

inline constexpr size_t div_up(size_t x, size_t y) noexcept { return (x + y - 1) / y; }
inline size_t block_count(size_t n, size_t nB) noexcept { return div_up(n, nB); }
inline size_t block_size_from_count(size_t n, size_t nb) noexcept {
    if (n == 0) return 0;
    if (nb <= 1) return n;
    return std::min<size_t>(n, common::padding(div_up(n, nb)));
}

inline constexpr size_t saturating_mul(size_t a, size_t b) noexcept {
    constexpr size_t max = std::numeric_limits<size_t>::max();
    return (a != 0 && b > max / a) ? max : a * b;
}

inline constexpr size_t triangular_number(size_t n) noexcept {
    constexpr size_t max = std::numeric_limits<size_t>::max();
    if (n == max) return max;
    return (n & 1) ? saturating_mul(n, (n + 1) / 2) : saturating_mul(n / 2, n + 1);
}

inline size_t performance_class(size_t n, size_t target) noexcept {
    size_t units = std::max<size_t>(1, std::min<size_t>(n, target) / common::PAD_SIZE);
    return std::bit_floor(units);
    // size_t q     = 1;
    // while (units >>= 1) q <<= 1;
    // return q;
}

inline size_t gemm_shape_quality(size_t mB, size_t nB, size_t kB) noexcept {
    const size_t qmn = std::min<size_t>(performance_class(mB, BLOCK_TARGET_MN), performance_class(nB, BLOCK_TARGET_MN));
    const size_t qk  = performance_class(kB, BLOCK_TARGET_K);
    return qmn * std::min<size_t>(qk, BLOCK_K_TO_MN_RATIO * qmn);
}

struct BlockState {
    size_t block  = 0;
    size_t count  = 0;
    size_t pclass = 0;
};

inline std::vector<BlockState> make_block_states(size_t n, size_t target) {
    std::vector<BlockState> states;
    if (n == 0) return states;

    const size_t max_nb = std::max<size_t>(1, div_up(n, common::PAD_SIZE));

    auto append = [&](size_t nb) {
        const size_t block = block_size_from_count(n, nb);
        if (!states.empty() && states.back().block == block) return;
        states.push_back({block, block_count(n, block), performance_class(block, target)});
    };

    append(1);
    if (max_nb == 1) return states;

    const size_t x = n - 1;
    for (size_t nb = 2; nb <= max_nb;) {
        append(nb);

        const size_t q = x / nb;
        size_t last    = max_nb;
        if (q != 0) last = std::min<size_t>(max_nb, x / q);

        if (last >= max_nb) break;
        nb = last + 1;
    }

    append(max_nb);
    return states;
}

inline size_t last_state_with_pclass_at_least(const std::vector<BlockState> &states, size_t q) noexcept {
    size_t lo = 0;
    size_t hi = states.size();
    while (lo < hi) {
        const size_t mid = lo + (hi - lo) / 2;
        if (states[mid].pclass >= q) {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    return (lo == 0) ? std::numeric_limits<size_t>::max() : lo - 1;
}

template <class WorkSize>
inline size_t find_max_gemm_quality(
    const std::vector<BlockState> &m_states,
    const std::vector<BlockState> &n_states,
    const std::vector<BlockState> &k_states,
    size_t limit,
    WorkSize &worksize //
) {
    auto collect_pclasses = [](const std::vector<BlockState> &states) {
        std::vector<size_t> q;
        q.reserve(states.size());
        for (const auto &s : states) {
            if (q.empty() || q.back() != s.pclass) q.push_back(s.pclass);
        }
        return q;
    };

    const auto qm_levels = collect_pclasses(m_states);
    const auto qn_levels = collect_pclasses(n_states);
    const auto qk_levels = collect_pclasses(k_states);

    std::vector<size_t> qmn_levels;
    qmn_levels.reserve(qm_levels.size() * qn_levels.size());
    for (const size_t qm : qm_levels) {
        for (const size_t qn : qn_levels) {
            qmn_levels.push_back(std::min<size_t>(qm, qn));
        }
    }
    std::sort(qmn_levels.begin(), qmn_levels.end());
    qmn_levels.erase(std::unique(qmn_levels.begin(), qmn_levels.end()), qmn_levels.end());

    struct Probe {
        size_t quality;
        size_t mi;
        size_t ni;
        size_t ki;
    };
    std::vector<Probe> probes;
    probes.reserve(qmn_levels.size() * qk_levels.size());

    for (const size_t qmn : qmn_levels) {
        const size_t mi = last_state_with_pclass_at_least(m_states, qmn);
        const size_t ni = last_state_with_pclass_at_least(n_states, qmn);
        if (mi == std::numeric_limits<size_t>::max() || ni == std::numeric_limits<size_t>::max()) continue;

        for (const size_t qk : qk_levels) {
            const size_t ki = last_state_with_pclass_at_least(k_states, qk);
            if (ki == std::numeric_limits<size_t>::max()) continue;

            const size_t quality = gemm_shape_quality(m_states[mi].block, n_states[ni].block, k_states[ki].block);
            probes.push_back({quality, mi, ni, ki});
        }
    }

    std::sort(probes.begin(), probes.end(), [](const Probe &a, const Probe &b) {
        if (a.quality != b.quality) return a.quality > b.quality;
        if (a.mi != b.mi) return a.mi > b.mi;
        if (a.ni != b.ni) return a.ni > b.ni;
        return a.ki > b.ki;
    });

    probes.erase(std::unique(probes.begin(), probes.end(),
                             [](const Probe &a, const Probe &b) {
                                 return a.mi == b.mi && a.ni == b.ni && a.ki == b.ki;
                             }),
                 probes.end());

    for (const auto &p : probes) {
        const auto &ms = m_states[p.mi];
        const auto &ns = n_states[p.ni];
        const auto &ks = k_states[p.ki];
        if (worksize(ms.block, ns.block, ks.block) <= limit) return p.quality;
    }

    return 0;
}

inline bool better_gemm_shape(const BlockSize3D &a, const BlockSize3D &b) noexcept {
    if (!b) return true;

    const size_t a_min = std::min<size_t>(a.mB, a.nB);
    const size_t a_max = std::max<size_t>(a.mB, a.nB);
    const size_t b_min = std::min<size_t>(b.mB, b.nB);
    const size_t b_max = std::max<size_t>(b.mB, b.nB);

    if (a_min != b_min) return a_min > b_min;

    const size_t a_aspect = div_up(a_max, a_min);
    const size_t b_aspect = div_up(b_max, b_min);
    if (a_aspect != b_aspect) return a_aspect < b_aspect;
    if (a_max != b_max) return a_max > b_max;
    return a.kB > b.kB;
}

template <class Fits>
inline size_t find_block_size_1d(size_t n, Fits &&fits) {
    if (n == 0) return 0;
    if (fits(n)) return n;

    const size_t max_nb = div_up(n, common::PAD_SIZE);
    if (max_nb <= 1) return 0;

    const size_t smallest = block_size_from_count(n, max_nb);
    if (!fits(smallest)) return 0;

    size_t lo = 2;
    size_t hi = max_nb;
    while (lo < hi) {
        const size_t mid = lo + (hi - lo) / 2;
        if (fits(block_size_from_count(n, mid))) {
            hi = mid;
        } else {
            lo = mid + 1;
        }
    }
    return block_size_from_count(n, lo);
}

template <class F>
inline void for_each_coarse_block_count(size_t n, size_t target, F &&f) {
    if (n == 0) return;

    const size_t max_nb = std::max<size_t>(1, div_up(n, common::PAD_SIZE));

    f(size_t(1));

    for (size_t nb = 2; nb < max_nb;) {
        f(nb);
        if (nb > max_nb / 2) break;
        nb *= 2;
    }

    if (max_nb > 1) f(max_nb);

    const size_t target_nb = std::min<size_t>(max_nb, std::max<size_t>(1, div_up(n, target)));

    if (target_nb > 1 && target_nb < max_nb && (target_nb & (target_nb - 1)) != 0) {
        f(target_nb);
    }
}

template <class F>
inline void for_each_nearby_block_count(size_t n, size_t center, F &&f) {
    const size_t max_nb = std::max<size_t>(1, div_up(n, common::PAD_SIZE));
    const size_t begin  = (center > 2) ? center - 2 : 1;
    const size_t end    = std::min<size_t>(max_nb, center + 2);

    for (size_t nb = begin; nb <= end; ++nb) f(nb);
}

template <class WorkSize>
inline BlockSize3D find_block_size_gemm_impl(
    size_t m, size_t n, size_t k,
    size_t limit,
    WorkSize &&worksize //
) {
    if (m == 0 || n == 0 || k == 0) return {};

    if (worksize(m, n, k) <= limit) return {m, n, k};

    const auto m_states = make_block_states(m, BLOCK_TARGET_MN);
    const auto n_states = make_block_states(n, BLOCK_TARGET_MN);
    const auto k_states = make_block_states(k, BLOCK_TARGET_K);
    if (m_states.empty() || n_states.empty() || k_states.empty()) return {};

    if (worksize(m_states.back().block, n_states.back().block, k_states.back().block) > limit) return {};

    auto &ws          = worksize;
    const size_t qmax = find_max_gemm_quality(m_states, n_states, k_states, limit, ws);
    if (qmax == 0) return {};

    struct PairCandidate {
        size_t mi;
        size_t ni;
        size_t k_first;
        size_t k_last;
        size_t calls_lb;
    };

    std::vector<PairCandidate> pairs;
    pairs.reserve(m_states.size() * n_states.size());

    auto quality_for = [](size_t qmn, const BlockState &ks) noexcept {
        return qmn * std::min<size_t>(ks.pclass, BLOCK_K_TO_MN_RATIO * qmn);
    };

    for (size_t mi = 0; mi < m_states.size(); ++mi) {
        const auto &ms = m_states[mi];
        for (size_t ni = 0; ni < n_states.size(); ++ni) {
            const auto &ns   = n_states[ni];
            const size_t qmn = std::min<size_t>(ms.pclass, ns.pclass);

            auto q_at = [&](size_t ki) noexcept { return quality_for(qmn, k_states[ki]); };

            size_t lo = 0;
            size_t hi = k_states.size();
            while (lo < hi) {
                const size_t mid = lo + (hi - lo) / 2;
                if (q_at(mid) > qmax) {
                    lo = mid + 1;
                } else {
                    hi = mid;
                }
            }
            if (lo == k_states.size() || q_at(lo) != qmax) continue;
            const size_t k_first = lo;

            lo = k_first;
            hi = k_states.size();
            while (lo < hi) {
                const size_t mid = lo + (hi - lo) / 2;
                if (q_at(mid) >= qmax) {
                    lo = mid + 1;
                } else {
                    hi = mid;
                }
            }
            const size_t k_last = lo - 1;

            const size_t calls_lb = saturating_mul(saturating_mul(ms.count, ns.count), k_states[k_first].count);

            pairs.push_back({mi, ni, k_first, k_last, calls_lb});
        }
    }

    std::sort(pairs.begin(), pairs.end(), [&](const PairCandidate &a, const PairCandidate &b) {
        if (a.calls_lb != b.calls_lb) return a.calls_lb < b.calls_lb;

        const BlockSize3D ba{m_states[a.mi].block, n_states[a.ni].block, k_states[a.k_first].block};
        const BlockSize3D bb{m_states[b.mi].block, n_states[b.ni].block, k_states[b.k_first].block};
        if (better_gemm_shape(ba, bb)) return true;
        if (better_gemm_shape(bb, ba)) return false;

        if (a.mi != b.mi) return a.mi < b.mi;
        if (a.ni != b.ni) return a.ni < b.ni;
        return a.k_first < b.k_first;
    });

    BlockSize3D best{};
    size_t best_calls = std::numeric_limits<size_t>::max();

    for (const auto &p : pairs) {
        if (best && p.calls_lb > best_calls) break;

        const auto &ms = m_states[p.mi];
        const auto &ns = n_states[p.ni];

        if (best && p.calls_lb == best_calls) {
            const BlockSize3D ideal{ms.block, ns.block, k_states[p.k_first].block};
            if (!better_gemm_shape(ideal, best)) continue;
        }

        auto fits_k = [&](size_t ki) {
            return worksize(ms.block, ns.block, k_states[ki].block) <= limit;
        };

        if (!fits_k(p.k_last)) continue;

        size_t ki = p.k_first;
        if (p.k_first != p.k_last && !fits_k(p.k_first)) {
            size_t lo = p.k_first + 1; // p.k_first is known infeasible
            size_t hi = p.k_last;      // p.k_last is known feasible
            while (lo < hi) {
                const size_t mid = lo + (hi - lo) / 2;
                if (fits_k(mid)) {
                    hi = mid;
                } else {
                    lo = mid + 1;
                }
            }
            ki = lo;
        }

        const auto &ks     = k_states[ki];
        const size_t calls = saturating_mul(saturating_mul(ms.count, ns.count), ks.count);
        const BlockSize3D block{ms.block, ns.block, ks.block};

        if (!best || calls < best_calls || (calls == best_calls && better_gemm_shape(block, best))) {
            best       = block;
            best_calls = calls;
        }
    }

    return best;
}

template <Backend BACKEND, unsigned NUM_MODULI, class WorkSize>
inline BlockSize3D find_block_size_gemm(
    size_t m, size_t n, size_t k, size_t limit,
    WorkSize &&worksize //
) {
    return find_block_size_gemm_impl(m, n, k, limit, std::forward<WorkSize>(worksize));
}

template <class WorkSize>
inline BlockSize3D find_block_size_gemm_runtime(
    Backend, unsigned,
    size_t m, size_t n, size_t k, size_t limit,
    WorkSize &&worksize //
) {
    return find_block_size_gemm_impl(m, n, k, limit, std::forward<WorkSize>(worksize));
}

template <Backend BACKEND, unsigned NUM_MODULI, class Fits>
inline BlockSize2D find_block_size_rankk(
    size_t n, size_t k, bool rank2,
    Fits &&fits //
) {
    if (n == 0 || k == 0) return {};

    struct Candidate {
        BlockSize2D block{};
        size_t quality = 0;
        size_t calls   = 0;
    };
    Candidate best{};

    auto consider = [&](size_t nn, size_t nk) {
        const size_t nB = block_size_from_count(n, nn);
        const size_t kB = block_size_from_count(k, nk);

        const size_t nn_actual = block_count(n, nB);
        if (!fits(nB, kB, nn_actual > 1)) return;

        const size_t nk_actual = block_count(k, kB);
        const size_t leaves    = rank2 ? saturating_mul(nn_actual, nn_actual) : triangular_number(nn_actual);
        const size_t calls     = saturating_mul(leaves, nk_actual);
        const size_t quality   = gemm_shape_quality(nB, nB, kB);

        if (!best.block ||
            quality > best.quality ||
            (quality == best.quality && calls < best.calls) ||
            (quality == best.quality && calls == best.calls && nB > best.block.xB) ||
            (quality == best.quality && calls == best.calls && nB == best.block.xB && kB > best.block.yB)) {
            best = {
                {nB, kB},
                quality,
                calls
            };
        }
    };

    for_each_coarse_block_count(n, BLOCK_TARGET_MN, [&](size_t nn) {
        for_each_coarse_block_count(k, BLOCK_TARGET_K, [&](size_t nk) {
            consider(nn, nk);
        });
    });

    if (!best.block) return {};

    const size_t nn0 = block_count(n, best.block.xB);
    const size_t nk0 = block_count(k, best.block.yB);

    for_each_nearby_block_count(n, nn0, [&](size_t nn) {
        for_each_nearby_block_count(k, nk0, [&](size_t nk) {
            consider(nn, nk);
        });
    });

    return best.block;
}

template <class Fits>
inline BlockSize2D find_block_size_structured(
    size_t s, size_t f, bool triangular,
    Fits &&fits //
) {
    if (s == 0 || f == 0) return {};

    struct Candidate {
        BlockSize2D block{};
        size_t quality = 0;
        size_t calls   = 0;
    };
    Candidate best{};

    auto consider = [&](size_t ns, size_t nf) {
        const size_t sB = block_size_from_count(s, ns);
        const size_t fB = block_size_from_count(f, nf);

        const size_t ns_actual = block_count(s, sB);
        if (!fits(sB, fB, ns_actual > 1)) return;

        const size_t nf_actual = block_count(f, fB);
        const size_t leaves    = triangular ? triangular_number(ns_actual) : saturating_mul(ns_actual, ns_actual);
        const size_t calls     = saturating_mul(leaves, nf_actual);
        const size_t quality   = gemm_shape_quality(sB, fB, sB);

        const size_t min_sf      = std::min<size_t>(sB, fB);
        const size_t max_sf      = std::max<size_t>(sB, fB);
        const size_t best_min_sf = std::min<size_t>(best.block.xB, best.block.yB);
        const size_t best_max_sf = std::max<size_t>(best.block.xB, best.block.yB);
        const size_t aspect      = div_up(max_sf, min_sf);
        const size_t best_aspect = best.block ? div_up(best_max_sf, best_min_sf) : 0;

        if (!best.block ||
            quality > best.quality ||
            (quality == best.quality && calls < best.calls) ||
            (quality == best.quality && calls == best.calls && min_sf > best_min_sf) ||
            (quality == best.quality && calls == best.calls && min_sf == best_min_sf && aspect < best_aspect) ||
            (quality == best.quality && calls == best.calls && min_sf == best_min_sf && aspect == best_aspect && sB > best.block.xB) ||
            (quality == best.quality && calls == best.calls && min_sf == best_min_sf && aspect == best_aspect && sB == best.block.xB && fB > best.block.yB)) {
            best = {
                {sB, fB},
                quality,
                calls
            };
        }
    };

    for_each_coarse_block_count(s, BLOCK_TARGET_MN, [&](size_t ns) {
        for_each_coarse_block_count(f, BLOCK_TARGET_MN, [&](size_t nf) {
            consider(ns, nf);
        });
    });

    if (!best.block) return {};

    const size_t ns0 = block_count(s, best.block.xB);
    const size_t nf0 = block_count(f, best.block.yB);

    for_each_nearby_block_count(s, ns0, [&](size_t ns) {
        for_each_nearby_block_count(f, nf0, [&](size_t nf) {
            consider(ns, nf);
        });
    });

    return best.block;
}

template <class Fits>
inline size_t find_block_size_trtrmm(size_t n, Fits &&fits) {
    return find_block_size_1d(n, std::forward<Fits>(fits));
}

} // namespace gemmul8::oz2::core
