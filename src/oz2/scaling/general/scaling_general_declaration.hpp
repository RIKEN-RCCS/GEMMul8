#pragma once
#include "../../common/common.hpp"

#include "scaling_colwise.hpp"

namespace gemmul8::scaling::general {

template <typename T, Backend BACKEND, unsigned NUM_MODULI,
          cublasFillMode_t UPLO, cublasDiagType_t DIAG, bool CONJ, bool HERK = false>
void scaling_rowwise(
    const cudaStream_t stream,
    const unsigned rows_A, const unsigned cols_A,
    const T *const A, const size_t lda,
    common::matptr_t<common::low_t<BACKEND>, common::isComplex<T>> A_lo,
    const size_t lda_lo, const size_t incA_lo,
    int16_t *const sftA);

template <typename T, Backend BACKEND, unsigned NUM_MODULI, cublasFillMode_t UPLO>
void scaling_symm(
    const cudaStream_t stream,
    const unsigned n,
    const T *const A, const size_t lda,
    common::matptr_t<common::low_t<BACKEND>, common::isComplex<T>> A_lo,
    const size_t lda_lo, const size_t incA_lo,
    const int16_t *const sftA);

template <typename T, Backend BACKEND, unsigned NUM_MODULI, cublasFillMode_t UPLO,
          bool STORE_TRANSPOSE = false>
void scaling_hemm(
    const cudaStream_t stream,
    const unsigned n,
    const T *const A, const size_t lda,
    common::matptr_t<common::low_t<BACKEND>, common::isComplex<T>> A_lo,
    const size_t lda_lo, const size_t incA_lo,
    const int16_t *const sftA);

template <typename T, Backend BACKEND, unsigned NUM_MODULI,
          cublasFillMode_t UPLO, cublasDiagType_t DIAG, bool CONJ, bool HERK = false>
__device__ __forceinline__ void scaling_colwise_device(
    const unsigned rows_A,
    const T *const __restrict__ in,
    common::matptr_t<common::low_t<BACKEND>, common::isComplex<T>> A_lo,
    const size_t lda_lo4, const size_t incA_lo4,
    const int32_t sft //
) {
    if constexpr (UPLO == CUBLAS_FILL_MODE_FULL) {

        scaling_colwise_full_device<T, BACKEND, NUM_MODULI, CONJ, false, false, HERK>(
            rows_A, in, A_lo, lda_lo4, incA_lo4, sft);

    } else if constexpr (UPLO == CUBLAS_FILL_MODE_UPPER) {

        scaling_colwise_upper_device<T, BACKEND, NUM_MODULI, DIAG, CONJ, HERK>(
            rows_A, in, A_lo, lda_lo4, incA_lo4, sft);

    } else {

        scaling_colwise_lower_device<T, BACKEND, NUM_MODULI, DIAG, CONJ, HERK>(
            rows_A, in, A_lo, lda_lo4, incA_lo4, sft);
    }
}

} // namespace gemmul8::scaling::general
