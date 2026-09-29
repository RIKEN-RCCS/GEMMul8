#pragma once
#include "undo_scaling_declaration.hpp"
#include "../common/table.hpp"
#include "../mod/reconstruct_fp8.hpp"
#include "../mod/complex_2m.hpp"
#include "crt_tree.hpp"
#include "../common/crt_group_plan.hpp"
#include "../mod/crt_group.hpp"

namespace gemmul8::undo_scaling {

template <Backend BACKEND, unsigned NUM_MODULI, bool COMPLEX> struct crt_moduli {
    static constexpr uint32_t get(unsigned i) {
        using namespace common::table;
        if constexpr (BACKEND == Backend::INT8) {
            return COMPLEX ? moduli_int8_complex[i] : moduli_int8[i];
        } else {
            return common::fp8_plan::modulus(NUM_MODULI, i, COMPLEX);
        }
    }
};

template <Backend BACKEND, unsigned IDX, bool COMPLEX>
__device__ __forceinline__ int32_t raw_product_residue(
    const common::hi_t<BACKEND> *ptr, size_t index, size_t plane_size) {
    if constexpr (BACKEND == Backend::INT8) {
        return mod::mod_small<BACKEND, IDX, COMPLEX>(ptr[index]);
    } else {
        return mod::load_fp8_product<IDX, COMPLEX>(ptr, index, plane_size);
    }
}

template <Backend BACKEND, unsigned NUM_MODULI, bool COMPLEX> struct crt_input {
    static constexpr bool complex = COMPLEX;
    const common::mid_t<BACKEND, COMPLEX> *mid;
    size_t plane_size;
    crt_tail<BACKEND> tail;
    size_t index;
    int32_t herk_opposite;

    template <unsigned IDX>
    __device__ __forceinline__ int2 load() const {
        if constexpr (IDX + 1 == NUM_MODULI) {
            if (tail.ptr0 != nullptr) {
                constexpr unsigned MOD = common::table::active_index<BACKEND, NUM_MODULI, IDX, COMPLEX>;
                const int32_t plus     = raw_product_residue<BACKEND, MOD, COMPLEX>(tail.ptr0, index, plane_size);
                if constexpr (!COMPLEX) {
                    return int2{plus, 0};
                } else {
                    const int32_t minus =
                        tail.ptr1 != nullptr
                            ? raw_product_residue<BACKEND, MOD, true>(tail.ptr1, index, plane_size)
                            : herk_opposite;
                    auto residue = mod::reconstruct_complex_2m<BACKEND, MOD, false>(plus, minus);
                    if (tail.flip_imag) {
                        residue.y = -residue.y;
                    }
                    return residue;
                }
            }
        }
        const auto r = mid[size_t(IDX) * plane_size];
        if constexpr (COMPLEX) {
            return int2{int32_t(r.x), int32_t(r.y)};
        } else {
            return int2{int32_t(r), 0};
        }
    }
};

template <Backend B, unsigned N, bool Complex, bool Fused> struct crt_group_input {
    static constexpr bool complex = Complex;

    using Group = std::conditional_t<Complex, uint2, uint32_t>;
    const Group *mid;
    size_t plane_size;
    const common::hi_t<B> *last;
    uint32_t index;

    template <unsigned I> __device__ __forceinline__ uint2 load() const {
        using Plan = common::crt_group_plan<B, N, Complex>;
        if constexpr (Fused && I + 1 == Plan::count) {
            constexpr unsigned first = Plan::groups.first[I];
            const mod::crt_group_input<B, N, Complex, first, N> input{last, index, uint32_t(plane_size)};
            const auto value = mod::crt_group_value<B, N, Complex, first, N>(input);
            if constexpr (Complex) return value;
            else return {value, 0};
        } else if constexpr (Complex) {
            return mid[size_t(I) * plane_size];
        } else {
            return {mid[size_t(I) * plane_size], 0};
        }
    }
};

template <typename T, Backend BACKEND, unsigned NUM_MODULI,
          bool GROUPED = false, bool FUSED = GROUPED>
__device__ __forceinline__ T reconstruct_from_crt(
    const common::mid_t<BACKEND, common::isComplex<T>> *const __restrict__ C_mid,
    const size_t incC_mid,
    const crt_tail<BACKEND> tail,
    const size_t index,
    const int32_t herk_opposite = 0 //
) {
    constexpr bool COMPLEX = common::isComplex<T>;
    using Plan             = common::crt_group_plan<BACKEND, NUM_MODULI, COMPLEX>;
    using Moduli           = std::conditional_t<GROUPED, Plan, crt_moduli<BACKEND, NUM_MODULI, COMPLEX>>;
    using Root             = crt::node<Moduli, 0, GROUPED ? Plan::count : NUM_MODULI>;
    constexpr auto modulus = Root::modulus;
    const auto input       = [&] {
        if constexpr (GROUPED) {
            using Input = crt_group_input<BACKEND, NUM_MODULI, COMPLEX, FUSED>;
            return Input{reinterpret_cast<const typename Input::Group *>(C_mid - index) + index, incC_mid, tail.ptr0, uint32_t(index)};
        } else {
            return crt_input<BACKEND, NUM_MODULI, COMPLEX>{C_mid, incC_mid, tail, index, herk_opposite};
        }
    }();
    const auto exact = Root::reconstruct(input);
    if constexpr (COMPLEX) {
        return T{crt::centered_double(exact.re, modulus),
                 crt::centered_double(exact.im, modulus)};
    } else {
        return crt::centered_double(exact.re, modulus);
    }
}

} // namespace gemmul8::undo_scaling
