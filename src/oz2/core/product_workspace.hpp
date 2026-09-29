#pragma once
#include "../common/common.hpp"
#include "../common/table.hpp"
#include "../common/crt_group_plan.hpp"

namespace gemmul8::oz2::core {

template <Backend BACKEND, bool Complex, bool Herk = false, bool Pointers = false>
struct product_workspace {
    static_assert(!Herk || Complex);
    static constexpr unsigned parts = Complex && !Herk ? 2U : 1U;

    static constexpr unsigned products(unsigned n, unsigned i) {
        if constexpr (BACKEND == Backend::FP8) {
            return common::fp8_plan::get(common::fp8_plan::modulus(n, i, Complex)).products;
        } else {
            return 1U;
        }
    }

    static constexpr unsigned planes(unsigned n, unsigned i, unsigned count) {
        unsigned sum = 0;
        for (unsigned j = 0; j < count; ++j) sum += products(n, i + j);
        return parts * sum;
    }

    static size_t pointer_count(unsigned n, unsigned i, unsigned count) {
        return Pointers ? planes(n, i, count) : 0U;
    }

    static size_t pointer_bytes(unsigned n, unsigned i, unsigned count) {
        return common::padding(3 * pointer_count(n, i, count) * sizeof(void *));
    }

    static size_t required(
        unsigned n, unsigned i, unsigned count,
        size_t sizeC, size_t blas //
    ) {
        const size_t mid = sizeof(common::mid_t<BACKEND, Complex>) * sizeC;
        const size_t gap = std::max(pointer_bytes(n, i, count) + blas, i + 1U < n ? mid : 0U);
        return i * mid + gap + sizeof(common::hi_t<BACKEND>) * sizeC * planes(n, i, count);
    }

    static size_t norm_bytes(size_t sizeC) {
        constexpr unsigned planes = Complex ? (BACKEND == Backend::INT8 ? 2U : 3U) : 1U;
        return sizeof(common::hi_t<BACKEND>) * sizeC * planes;
    }

    static size_t bytes(unsigned n, size_t sizeC, size_t blas, bool fastmode) {
        size_t result = fastmode ? 0U : norm_bytes(sizeC) + blas;
        {
            const size_t grouped_size = std::min(sizeC, size_t(8192) * 8192);
            const size_t group_bytes  = sizeof(uint32_t) * (Complex ? 2U : 1U) * grouped_size;
            for (unsigned i = 0, g = 0; i < n; ++g) {
                const unsigned end = common::crt_group_end<BACKEND, Complex>(n, i);
                const size_t gap   = std::max(group_bytes, pointer_bytes(n, i, end - i) + blas);
                result             = std::max(result, g * group_bytes + gap + sizeof(common::hi_t<BACKEND>) * grouped_size * planes(n, i, end - i));
                i                  = end;
            }
            if (common::use_grouped_crt(sizeC)) return result;
        }
        for (unsigned i = 0; i < n; ++i) {
            result = std::max(result, required(n, i, 1U, sizeC, blas));
        }
        return result;
    }
};

} // namespace gemmul8::oz2::core
