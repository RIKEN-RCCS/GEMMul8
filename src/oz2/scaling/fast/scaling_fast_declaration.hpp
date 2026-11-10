#pragma once
#include "../../common/common.hpp"

namespace gemmul8::scaling::fast {

template <typename T, Backend BACKEND, unsigned NUM_MODULI,
          cublasFillMode_t UPLO = CUBLAS_FILL_MODE_FULL,
          cublasDiagType_t DIAG = CUBLAS_DIAG_NON_UNIT, bool HERK = false>
void scaling(
    const cudaStream_t stream,
    const cublasOperation_t op_A, const cublasSideMode_t side,
    const unsigned rows_A, const unsigned cols_A,
    const T *const A, const size_t lda,
    common::matptr_t<common::low_t<BACKEND>, common::isComplex<T>> &A_lo,
    const size_t lda_lo, const size_t incA_lo,
    int16_t *const sftA);

template <typename T, Backend BACKEND, unsigned NUM_MODULI,
          cublasFillMode_t UPLO_A>
void scaling_symm(
    const cudaStream_t stream,
    const unsigned n,
    const T *const A, const size_t lda,
    common::matptr_t<common::low_t<BACKEND>, common::isComplex<T>> &A_lo,
    const size_t lda_lo, const size_t incA_lo,
    int16_t *const sftA);

template <typename T, Backend BACKEND, unsigned NUM_MODULI,
          cublasFillMode_t UPLO_A,
          bool STORE_TRANSPOSE = false>
void scaling_hemm(
    const cudaStream_t stream,
    const unsigned n,
    const T *const A, const size_t lda,
    common::matptr_t<common::low_t<BACKEND>, common::isComplex<T>> &A_lo,
    const size_t lda_lo, const size_t incA_lo,
    int16_t *const sftA);

} // namespace gemmul8::scaling::fast
