#pragma once
#include "config.hpp"
#include "calc_sft.hpp"
#include "scaling_colwise.hpp"

#include "../general/scaling_general_declaration.hpp"
#include "../general/helper_temp.hpp"

namespace gemmul8::scaling::accu {

template <typename T, Backend BACKEND, unsigned NUM_MODULI,
          cublasFillMode_t UPLO_C>
void scaling_syrk(
    const cudaStream_t stream,
    const cublasOperation_t op_A,
    const unsigned n, const unsigned k,
    const T *const A, const size_t lda,
    common::matptr_t<common::low_t<BACKEND>, common::isComplex<T>> &A_lo,
    const size_t lda_lo, const size_t incA_lo,
    int16_t *const sftA,
    common::matptr_t<common::hi_t<BACKEND>, common::isComplex<T>> &C_hi,
    const size_t ldc_hi //
) {
    using Hi                      = common::hi_t<BACKEND>;
    const unsigned num_col_blocks = rowwise_sft_num_blocks<T, BACKEND, NUM_MODULI>(n, n, incA_lo, true);
    Hi *const partial_amax        = general::temporary_memory<common::low_t<BACKEND>, Hi>(A_lo.ptr0);
    Hi *const maxC                = partial_amax + size_t(n) * num_col_blocks;

    calc_sft_after_rowwise_sym_launch<T, BACKEND, NUM_MODULI, UPLO_C>(
        stream, n, n, k, C_hi, ldc_hi, partial_amax, maxC, num_col_blocks);

    calc_sft_after_colwise_sym<T, BACKEND, NUM_MODULI, UPLO_C>
        <<<n, threads_accu, 0, stream>>>(
            n, k, C_hi, ldc_hi, maxC, sftA);

    if (op_A == CUBLAS_OP_N) {

        // SYRK: A*A^T -> A_lo^T*A_lo
        // A: n x k -> A_lo: k x n
        general::scaling_rowwise<T, BACKEND, NUM_MODULI,
                                 CUBLAS_FILL_MODE_FULL, CUBLAS_DIAG_NON_UNIT, false>(
            stream, n, k, A, lda, A_lo, lda_lo, incA_lo, sftA);

    } else if (op_A == CUBLAS_OP_T) {

        // SYRK: A^T*A -> A_lo^T*A_lo
        // A: k x n -> A_lo: k x n
        scaling_colwise_launch<T, BACKEND, NUM_MODULI,
                               CUBLAS_FILL_MODE_FULL, CUBLAS_DIAG_NON_UNIT, false>(
            stream, k, n, A, lda, A_lo, lda_lo >> 2, incA_lo >> 2, sftA);

    } else {

        // SYRK: conj(A)^T*conj(A) -> A_lo^T*A_lo
        // A: k x n -> A_lo: k x n
        constexpr bool CONJ = (common::isComplex<T>) ? true : false;
        scaling_colwise_launch<T, BACKEND, NUM_MODULI,
                               CUBLAS_FILL_MODE_FULL, CUBLAS_DIAG_NON_UNIT, CONJ>(
            stream, k, n, A, lda, A_lo, lda_lo >> 2, incA_lo >> 2, sftA);
    }
}

template <typename T, Backend BACKEND, unsigned NUM_MODULI,
          cublasFillMode_t UPLO_C>
void scaling_herk(
    const cudaStream_t stream,
    const cublasOperation_t op_A,
    const unsigned n, const unsigned k,
    const T *const A, const size_t lda,
    common::matptr_t<common::low_t<BACKEND>, common::isComplex<T>> &A_lo,
    const size_t lda_lo, const size_t incA_lo,
    int16_t *const sftA,
    common::matptr_t<common::hi_t<BACKEND>, common::isComplex<T>> &C_hi,
    const size_t ldc_hi //
) {
    using Hi                      = common::hi_t<BACKEND>;
    const unsigned num_col_blocks = rowwise_sft_num_blocks<T, BACKEND, NUM_MODULI, BACKEND == Backend::FP8>(n, n, incA_lo, true);
    Hi *const partial_amax        = general::temporary_memory<common::low_t<BACKEND>, Hi>(A_lo.ptr0);
    Hi *const maxC                = partial_amax + size_t(n) * num_col_blocks;

    calc_sft_after_rowwise_sym_launch<T, BACKEND, NUM_MODULI, UPLO_C>(
        stream, n, n, k, C_hi, ldc_hi, partial_amax, maxC, num_col_blocks);

    calc_sft_after_colwise_sym<T, BACKEND, NUM_MODULI, UPLO_C, BACKEND == Backend::FP8>
        <<<n, threads_accu, 0, stream>>>(
            n, k, C_hi, ldc_hi, maxC, sftA);

    if (op_A == CUBLAS_OP_N) {

        // HERK: A*A^H -> A_lo^H*A_lo
        // A: n x k -> A_lo: k x n
        general::scaling_rowwise<T, BACKEND, NUM_MODULI, CUBLAS_FILL_MODE_FULL, CUBLAS_DIAG_NON_UNIT, false, BACKEND == Backend::FP8>(
            stream, n, k, A, lda, A_lo, lda_lo, incA_lo, sftA);

    } else {

        // HERK: A^H*A -> A_lo^H*A_lo
        // A: k x n -> A_lo: k x n
        scaling_colwise_launch<T, BACKEND, NUM_MODULI, CUBLAS_FILL_MODE_FULL, CUBLAS_DIAG_NON_UNIT, false, BACKEND == Backend::FP8>(
            stream, k, n, A, lda, A_lo, lda_lo >> 2, incA_lo >> 2, sftA);
    }
}

} // namespace gemmul8::scaling::accu
