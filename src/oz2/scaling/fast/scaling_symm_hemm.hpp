#pragma once
#include "config.hpp"
#include "calc_sft.hpp"

#include "../general/scaling_general_declaration.hpp"
#include "../general/helper_temp.hpp"

namespace gemmul8::scaling::fast {

template <typename T, Backend BACKEND, unsigned NUM_MODULI,
          cublasFillMode_t UPLO, bool HERM,
          bool STORE_TRANSPOSE = false>
void scaling_symm_hemm_launch(
    const cudaStream_t stream,
    const unsigned n,
    const T *const A, const size_t lda,
    common::matptr_t<common::low_t<BACKEND>, common::isComplex<T>> &A_lo,
    const size_t lda_lo, const size_t incA_lo,
    int16_t *const sftA //
) {
    const unsigned num_col_blocks =
        (n < rowwise_sft_split_threshold)
            ? 1U
            : (n + rowwise_sft_col_tile - 1U) / rowwise_sft_col_tile;

    using U                    = common::underlying_t<T>;
    const size_t partial_elems = size_t(n) * num_col_blocks;
    U *const partial_sum       = general::temporary_memory<common::low_t<BACKEND>, U>(A_lo.ptr0);
    U *const sumA              = partial_sum + partial_elems;

    calc_stat_sym_rowwise_launch<T, UPLO, HERM>(
        stream, n, A, lda, partial_sum, sumA);

    calc_sft_sym_colwise<T, BACKEND, NUM_MODULI, UPLO>
        <<<n, threads_fast, 0, stream>>>(
            n, A, lda, sftA, sumA);

    if constexpr (HERM) {

        general::scaling_hemm<T, BACKEND, NUM_MODULI, UPLO, STORE_TRANSPOSE>(
            stream, n, A, lda, A_lo, lda_lo, incA_lo, sftA);

    } else {

        general::scaling_symm<T, BACKEND, NUM_MODULI, UPLO>(
            stream, n, A, lda, A_lo, lda_lo, incA_lo, sftA);
    }
}

template <typename T, Backend BACKEND, unsigned NUM_MODULI, cublasFillMode_t UPLO_A>
void scaling_symm(
    const cudaStream_t stream,
    const unsigned n,
    const T *const A, const size_t lda,
    common::matptr_t<common::low_t<BACKEND>, common::isComplex<T>> &A_lo,
    const size_t lda_lo, const size_t incA_lo,
    int16_t *const sftA //
) {
    static_assert(UPLO_A == CUBLAS_FILL_MODE_UPPER || UPLO_A == CUBLAS_FILL_MODE_LOWER,
                  "scaling_symm requires UPLO = UPPER or LOWER.");

    scaling_symm_hemm_launch<T, BACKEND, NUM_MODULI, UPLO_A, false>(
        stream, n, A, lda, A_lo, lda_lo, incA_lo, sftA);
}

template <typename T, Backend BACKEND, unsigned NUM_MODULI, cublasFillMode_t UPLO_A,
          bool STORE_TRANSPOSE>
void scaling_hemm(
    const cudaStream_t stream,
    const unsigned n,
    const T *const A, const size_t lda,
    common::matptr_t<common::low_t<BACKEND>, common::isComplex<T>> &A_lo,
    const size_t lda_lo, const size_t incA_lo,
    int16_t *const sftA //
) {
    static_assert(UPLO_A == CUBLAS_FILL_MODE_UPPER || UPLO_A == CUBLAS_FILL_MODE_LOWER,
                  "scaling_hemm requires UPLO = UPPER or LOWER.");

    scaling_symm_hemm_launch<T, BACKEND, NUM_MODULI, UPLO_A, true, STORE_TRANSPOSE>(
        stream, n, A, lda, A_lo, lda_lo, incA_lo, sftA);
}

} // namespace gemmul8::scaling::fast
