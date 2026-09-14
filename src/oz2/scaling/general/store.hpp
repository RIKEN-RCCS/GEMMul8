#pragma once
#include "../../common/common.hpp"
#include "../../common/table.hpp"
#include "../../mod/mod.hpp"
#include "roundup.hpp"

namespace gemmul8::scaling::general {

template <typename T, Backend BACKEND, unsigned NUM_MODULI>
inline void memset_low_mats_async(
    const cudaStream_t stream,
    common::matptr_t<common::low_t<BACKEND>, common::isComplex<T>> A_lo,
    const size_t incA_lo //
) {
    constexpr size_t cmul = common::isComplex<T> ? (NUM_MODULI == 0U ? 3ULL : 2ULL) : 1ULL;
    size_t bytes          = incA_lo * sizeof(common::low_t<BACKEND>) * cmul;
    if constexpr (NUM_MODULI != 0U) {
        bytes *= common::table::num_mat_v<BACKEND, NUM_MODULI, common::isComplex<T>>;
    }
    cudaMemsetAsync(A_lo.ptr0, 0, bytes, stream);
}

template <typename T, Backend BACKEND, unsigned NUM_MODULI>
inline void memset_padding_low_mats_2d_async(
    const cudaStream_t stream,
    common::matptr_t<common::low_t<BACKEND>, common::isComplex<T>> A_lo,
    const size_t n,
    const size_t lda_lo,
    const size_t cols_lo //
) {
    if (lda_lo <= n) return;
    using LowT = common::low_t<BACKEND>;

    constexpr size_t cmul       = common::isComplex<T> ? (NUM_MODULI == 0U ? 3ULL : 2ULL) : 1ULL;
    const size_t width_in_bytes = (lda_lo - n) * sizeof(LowT);
    size_t height               = cols_lo * cmul;
    if constexpr (NUM_MODULI != 0U) {
        height *= common::table::num_mat_v<BACKEND, NUM_MODULI, common::isComplex<T>>;
    }
    cudaMemset2DAsync(A_lo.ptr0 + n, lda_lo * sizeof(LowT),
                      0, width_in_bytes, height, stream);
}

template <typename T, Backend BACKEND, unsigned NUM_MODULI, bool CONJ>
__device__ __forceinline__ void scaling_colwise_store4_real(
    common::lowx4_t<BACKEND> *__restrict__ out,
    const size_t incA_lo4,
    const unsigned i,
    const T a0, const T a1, const T a2, const T a3,
    const int32_t sft //
) {
    using ValT = decltype(trunc_scalbn<true, T, BACKEND, NUM_MODULI>::run(T{}, int32_t{}));

    const ValT v0 = trunc_scalbn<true, T, BACKEND, NUM_MODULI>::run(a0, sft);
    const ValT v1 = trunc_scalbn<true, T, BACKEND, NUM_MODULI>::run(a1, sft);
    const ValT v2 = trunc_scalbn<true, T, BACKEND, NUM_MODULI>::run(a2, sft);
    const ValT v3 = trunc_scalbn<true, T, BACKEND, NUM_MODULI>::run(a3, sft);

    mod::ModUnroll<NUM_MODULI, ValT>::run(out + i, incA_lo4, v0, v1, v2, v3);
}

template <typename T, Backend BACKEND, unsigned NUM_MODULI, bool CONJ>
__device__ __forceinline__ void scaling_colwise_store4_complex(
    common::lowx4_t<BACKEND> *__restrict__ out_1,
    common::lowx4_t<BACKEND> *__restrict__ out_2,
    const size_t incA_lo4,
    const unsigned i,
    const T a0, const T a1, const T a2, const T a3,
    const int32_t sft //
) {
    using ValT = decltype(trunc_scalbn<true, T, BACKEND, NUM_MODULI>::run(T{}, int32_t{}));

    const ValT v0 = trunc_scalbn<true, T, BACKEND, NUM_MODULI>::run(common::conj<T, CONJ>(a0), sft);
    const ValT v1 = trunc_scalbn<true, T, BACKEND, NUM_MODULI>::run(common::conj<T, CONJ>(a1), sft);
    const ValT v2 = trunc_scalbn<true, T, BACKEND, NUM_MODULI>::run(common::conj<T, CONJ>(a2), sft);
    const ValT v3 = trunc_scalbn<true, T, BACKEND, NUM_MODULI>::run(common::conj<T, CONJ>(a3), sft);

    mod::ModUnroll<NUM_MODULI, ValT>::run(out_1 + i, out_2 + i, incA_lo4, v0, v1, v2, v3);
}

template <typename T, Backend BACKEND, unsigned NUM_MODULI>
__device__ __forceinline__ void scaling_store_one_real(
    common::matptr_t<common::low_t<BACKEND>, false> A_lo,
    const size_t idx,
    const size_t incA_lo,
    const T a,
    const int32_t sft //
) {
    using ValT = decltype(trunc_scalbn<true, T, BACKEND, NUM_MODULI>::run(T{}, int32_t{}));

    common::low_t<BACKEND> *__restrict__ out = A_lo.ptr0 + idx;

    const ValT v = trunc_scalbn<true, T, BACKEND, NUM_MODULI>::run(a, sft);
    mod::ModUnroll<NUM_MODULI, ValT>::run(out, incA_lo, v);
}

template <typename T, Backend BACKEND, unsigned NUM_MODULI>
__device__ __forceinline__ void scaling_store_one_complex(
    common::matptr_t<common::low_t<BACKEND>, true> A_lo,
    const size_t idx,
    const size_t incA_lo,
    const T a,
    const int32_t sft //
) {
    using ValT = decltype(trunc_scalbn<true, T, BACKEND, NUM_MODULI>::run(T{}, int32_t{}));

    common::low_t<BACKEND> *__restrict__ out_1 = A_lo.ptr0 + idx;
    common::low_t<BACKEND> *__restrict__ out_2 = A_lo.ptr1 + idx;

    const ValT v = trunc_scalbn<true, T, BACKEND, NUM_MODULI>::run(a, sft);
    mod::ModUnroll<NUM_MODULI, ValT>::run(out_1, out_2, incA_lo, v);
}

template <typename T, Backend BACKEND, unsigned NUM_MODULI>
__device__ __forceinline__ void scaling_store_one(
    common::matptr_t<common::low_t<BACKEND>, common::isComplex<T>> A_lo,
    const size_t idx,
    const size_t incA_lo,
    const T a,
    const int32_t sft //
) {
    if constexpr (common::isComplex<T>) {
        scaling_store_one_complex<T, BACKEND, NUM_MODULI>(
            A_lo, idx, incA_lo, a, sft);
    } else {
        scaling_store_one_real<T, BACKEND, NUM_MODULI>(
            A_lo, idx, incA_lo, a, sft);
    }
}

} // namespace gemmul8::scaling::general
