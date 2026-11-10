#pragma once
#include "../common/fp8_plan.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>

namespace gemmul8::oz2::core {

struct Fp8KCost {
    size_t gemm_planes      = 0;
    size_t reduction_planes = 0;
};

namespace fp8_k_cost_detail {

inline constexpr size_t alignment = 256;
inline constexpr uint64_t budget  = uint64_t(1) << 24;
inline constexpr size_t maximum   = std::numeric_limits<size_t>::max();

constexpr size_t add(size_t a, size_t b) noexcept { return b > maximum - a ? maximum : a + b; }
constexpr size_t multiply(size_t a, size_t b) noexcept { return a && b > maximum / a ? maximum : a * b; }

inline void accumulate(Fp8KCost &a, const Fp8KCost &b, size_t count = 1) noexcept {
    a.gemm_planes      = add(a.gemm_planes, multiply(b.gemm_planes, count));
    a.reduction_planes = add(a.reduction_planes, multiply(b.reduction_planes, count));
}

struct Plane {
    unsigned first = 0, next = 0;
    unsigned square = 0, remainder = 0;
    unsigned initial = 0, increment = 0, free_steps = 0, period = 0;
};

struct Modulus {
    unsigned first = 0, next = 0, products = 0, parts = 0;
    std::array<Plane, 3> planes{};
};

template <unsigned N>
struct Descriptor {
    size_t safe = maximum, base_products = 0;
    std::array<Modulus, N> moduli{};
};

template <unsigned N, bool Complex, bool Herk>
constexpr Descriptor<N> prepare() noexcept {
    static_assert(N >= 2 && N <= 20);
    static_assert(!Herk || Complex);
    Descriptor<N> result{};
    for (unsigned i = 0; i < N; ++i) {
        const auto p = common::fp8_plan::modulus(N, i, Complex, Herk);
        const auto s = common::fp8_plan::get(p);
        auto &m      = result.moduli[i];
        m.first      = common::fp8_plan::k_block(p, true);
        m.next       = common::fp8_plan::k_block(p, false);
        m.products   = s.products;
        m.parts      = Complex ? common::fp8_plan::complex_parts(p, Herk) : 1U;
        result.safe  = std::min(result.safe, size_t(m.first));
        result.base_products += m.products * m.parts;
        for (unsigned j = 0; j < m.products; ++j) {
            auto &v                  = m.planes[j];
            const unsigned magnitude = common::fp8_plan::plane_max_abs(s, j);
            v.first                  = common::fp8_plan::plane_k_block(p, j, true);
            v.next                   = common::fp8_plan::plane_k_block(p, j, false);
            v.square                 = magnitude * magnitude;
            v.remainder              = unsigned(s.q[j] / 2) + 1U;
            v.initial                = m.first * v.square;
            v.increment              = m.next * v.square;
            v.free_steps             = unsigned((budget - v.initial) / v.increment);
            v.period                 = unsigned((budget - v.remainder) / v.increment);
        }
    }
    return result;
}

template <unsigned N, bool Complex, bool Herk>
inline constexpr auto descriptor = prepare<N, Complex, Herk>();

struct Segments {
    size_t full, tail, count;
};

inline Segments split(size_t k, unsigned first, unsigned next) noexcept {
    if (k <= first) return {0, 0, 1};
    const size_t remaining = k - first;
    const size_t full = remaining / next, tail = remaining % next;
    return {full, tail, 1 + full + size_t(tail != 0)};
}

inline size_t reductions(const Plane &p, const Segments &s) noexcept {
    uint64_t used = p.initial;
    size_t count  = 0;
    if (s.full <= p.free_steps) {
        used += uint64_t(s.full) * p.increment;
    } else {
        const size_t after = s.full - p.free_steps - 1;
        count              = 1 + after / p.period;
        used               = p.remainder + uint64_t(1 + after % p.period) * p.increment;
    }
    if (s.tail && used + uint64_t(s.tail) * p.square > budget) ++count;
    return count;
}

template <bool SharedSegments>
inline Fp8KCost modulus_cost(const Modulus &m, size_t k) noexcept {
    Fp8KCost result{};
    if constexpr (SharedSegments) {
        const auto s          = split(k, m.first, m.next);
        const size_t products = size_t(m.products) * m.parts;
        result.gemm_planes    = multiply(products, s.count);
        if (s.count > 1) {
            for (unsigned j = 0; j < m.products; ++j)
                result.reduction_planes = add(result.reduction_planes, multiply(m.parts, reductions(m.planes[j], s)));
        }
    } else {
        for (unsigned j = 0; j < m.products; ++j) {
            const auto &p = m.planes[j];
            const auto s  = split(k, p.first, p.next);
            const Fp8KCost cost{multiply(m.parts, s.count), multiply(m.parts, s.count - 1)};
            accumulate(result, cost);
        }
    }
    return result;
}

template <unsigned N, bool Complex, bool SharedSegments, bool Herk>
inline Fp8KCost single(size_t k) noexcept {
    if (!k) return {};
    if (k > maximum - (alignment - 1)) return {maximum, maximum};
    k                 = (k + alignment - 1) / alignment * alignment;
    constexpr auto &d = descriptor<N, Complex, Herk>;
    if (k <= d.safe)
        return {d.base_products, 0};
    Fp8KCost result{};
    for (const auto &m : d.moduli) accumulate(result, modulus_cost<SharedSegments>(m, k));
    return result;
}

} // namespace fp8_k_cost_detail

template <unsigned N, bool Complex, bool Herk = false>
inline constexpr size_t fp8_k_safe_v = fp8_k_cost_detail::descriptor<N, Complex, Herk>.safe;

template <unsigned N, bool Complex, bool SharedSegments, bool Herk = false>
inline Fp8KCost fp8_k_cost(size_t k, size_t k_block) noexcept {
    using namespace fp8_k_cost_detail;
    if (!k) return {};
    if (!k_block) return {maximum, maximum};
    const size_t full = k / k_block, tail = k % k_block;
    Fp8KCost result{};
    if (full) accumulate(result, single<N, Complex, SharedSegments, Herk>(k_block), full);
    if (tail) accumulate(result, single<N, Complex, SharedSegments, Herk>(tail));
    return result;
}

} // namespace gemmul8::oz2::core
