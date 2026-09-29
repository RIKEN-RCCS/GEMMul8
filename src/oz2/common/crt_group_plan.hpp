#pragma once
#include "table.hpp"

namespace gemmul8::common {

inline constexpr bool use_grouped_crt(size_t sizeC) {
    return sizeC <= size_t(8192) * 8192;
}

template <Backend BACKEND, bool Complex>
constexpr uint32_t crt_modulus(unsigned n, unsigned i) {
    if constexpr (BACKEND == Backend::FP8) {
        return uint32_t(fp8_plan::modulus(n, i, Complex));
    } else {
        return uint32_t(Complex ? table::moduli_int8_complex[i] : table::moduli_int8[i]);
    }
}

template <Backend BACKEND, bool Complex>
constexpr unsigned crt_group_end(unsigned n, unsigned first) {
    unsigned end = first;
    uint64_t p   = 1;
    while (end < n && end < first + 4 && p * crt_modulus<BACKEND, Complex>(n, end) <= UINT32_MAX) {
        p *= crt_modulus<BACKEND, Complex>(n, end++);
    }
    return end;
}

template <Backend BACKEND, unsigned N, bool Complex> struct crt_moduli {
    static constexpr uint32_t get(unsigned i) {
        return crt_modulus<BACKEND, Complex>(N, i);
    }
};

template <Backend BACKEND, unsigned N, bool Complex> struct crt_group_plan {
    struct data {
        unsigned first[N + 1]{};
        unsigned count = 0;
    };
    static constexpr data groups = [] {
        data r{};
        while (r.first[r.count] < N) {
            const unsigned end = crt_group_end<BACKEND, Complex>(N, r.first[r.count]);
            r.first[++r.count] = end;
        }
        return r;
    }();
    static constexpr unsigned count = groups.count;
    static constexpr uint32_t get(unsigned g) {
        uint32_t p = 1;
        for (unsigned i = groups.first[g]; i < groups.first[g + 1]; ++i) {
            p *= crt_moduli<BACKEND, N, Complex>::get(i);
        }
        return p;
    }
};

} // namespace gemmul8::common
