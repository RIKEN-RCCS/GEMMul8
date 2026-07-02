#pragma once
#include "config.hpp"
#include "find_max.hpp"
#include "calc_sft.hpp"

#include "../../mod/mod.hpp"
#include "../general/scaling_general_declaration.hpp"
#include "../general/helper_temp.hpp"

namespace gemmul8::scaling::fast {

template <typename T, Backend BACKEND, unsigned NUM_MODULI,
          cublasFillMode_t UPLO, cublasDiagType_t DIAG, bool CONJ>
__global__ void scaling_colwise(
    const unsigned rows_A,
    const T *const __restrict__ A, const size_t lda,
    common::matptr_t<common::low_t<BACKEND>, common::isComplex<T>> A_lo,
    const size_t lda_lo, const size_t incA_lo,
    int16_t *const __restrict__ sftA //
) {
    using U = common::underlying_t<T>;
    __shared__ U ssum[32];

    const unsigned col_idx         = blockIdx.x;
    const T *const __restrict__ in = A + col_idx * lda;

    const int32_t sft = calc_sft_colwise<T, BACKEND, NUM_MODULI, UPLO, DIAG>(
        rows_A, in, sftA, ssum);

    general::scaling_colwise_device<T, BACKEND, NUM_MODULI, UPLO, DIAG, CONJ>(
        rows_A, in, A_lo, lda_lo, incA_lo, sft);
}

template <typename T, Backend BACKEND, unsigned NUM_MODULI,
          cublasFillMode_t UPLO, cublasDiagType_t DIAG>
void scaling(
    const cudaStream_t stream,
    const cublasOperation_t op_A, const cublasSideMode_t side,
    const unsigned rows_A, const unsigned cols_A,
    const T *const A, const size_t lda,
    common::matptr_t<common::low_t<BACKEND>, common::isComplex<T>> &A_lo,
    const size_t lda_lo, const size_t incA_lo,
    int16_t *const sftA //
) {
    if (side == CUBLAS_SIDE_LEFT) {
        if (op_A == CUBLAS_OP_N) {

            // A: rows_A x cols_A -> A_lo: cols_A x rows_A
            using U              = common::underlying_t<T>;
            U *const partial_sum = general::temporary_memory<common::low_t<BACKEND>, U>(A_lo.ptr0);

            calc_sft_rowwise_launch<T, BACKEND, NUM_MODULI, UPLO, DIAG>(
                stream, rows_A, cols_A, A, lda, sftA, partial_sum);

            if constexpr (UPLO != CUBLAS_FILL_MODE_FULL) {
                general::memset_low_mats_async<T, BACKEND, NUM_MODULI>(stream, A_lo, incA_lo);
            }

            general::scaling_rowwise<T, BACKEND, NUM_MODULI, UPLO, DIAG, false>(
                stream, rows_A, cols_A, A, lda, A_lo, lda_lo, incA_lo, sftA);

        } else {
            if constexpr (UPLO != CUBLAS_FILL_MODE_FULL) {
                general::memset_low_mats_async<T, BACKEND, NUM_MODULI>(stream, A_lo, incA_lo);
            }

            // A: cols_A x rows_A -> A_lo: cols_A x rows_A
            if (op_A == CUBLAS_OP_T) {

                scaling_colwise<T, BACKEND, NUM_MODULI, UPLO, DIAG, false>
                    <<<rows_A, threads_fast, 0, stream>>>(
                        cols_A, A, lda, A_lo, lda_lo >> 2, incA_lo >> 2, sftA);

            } else {

                constexpr bool CONJ = (common::isComplex<T>) ? true : false;
                scaling_colwise<T, BACKEND, NUM_MODULI, UPLO, DIAG, CONJ>
                    <<<rows_A, threads_fast, 0, stream>>>(
                        cols_A, A, lda, A_lo, lda_lo >> 2, incA_lo >> 2, sftA);
            }
        }

    } else {
        if (op_A == CUBLAS_OP_N) {

            if constexpr (UPLO != CUBLAS_FILL_MODE_FULL) {
                general::memset_low_mats_async<T, BACKEND, NUM_MODULI>(stream, A_lo, incA_lo);
            }

            // A: rows_A x cols_A -> A_lo: rows_A x cols_A
            scaling_colwise<T, BACKEND, NUM_MODULI, UPLO, DIAG, false>
                <<<cols_A, threads_fast, 0, stream>>>(
                    rows_A, A, lda, A_lo, lda_lo >> 2, incA_lo >> 2, sftA);

        } else {

            // A: cols_A x rows_A -> A_lo: rows_A x cols_A
            using U              = common::underlying_t<T>;
            U *const partial_sum = general::temporary_memory<common::low_t<BACKEND>, U>(A_lo.ptr0);

            calc_sft_rowwise_launch<T, BACKEND, NUM_MODULI, UPLO, DIAG>(
                stream, cols_A, rows_A, A, lda, sftA, partial_sum);

            if constexpr (UPLO != CUBLAS_FILL_MODE_FULL) {
                general::memset_low_mats_async<T, BACKEND, NUM_MODULI>(stream, A_lo, incA_lo);
            }

            if (op_A == CUBLAS_OP_T) {

                general::scaling_rowwise<T, BACKEND, NUM_MODULI, UPLO, DIAG, false>(
                    stream, cols_A, rows_A, A, lda, A_lo, lda_lo, incA_lo, sftA);

            } else {

                constexpr bool CONJ = (common::isComplex<T>) ? true : false;
                general::scaling_rowwise<T, BACKEND, NUM_MODULI, UPLO, DIAG, CONJ>(
                    stream, cols_A, rows_A, A, lda, A_lo, lda_lo, incA_lo, sftA);
            }
        }
    }
}

} // namespace gemmul8::scaling::fast
