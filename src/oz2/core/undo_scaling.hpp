#pragma once
#include "../common/common.hpp"
#include "../undo_scaling/undo_scaling_declaration.hpp"
#include "helper_triangular.hpp"

namespace gemmul8::oz2::core {

template <Func FUNC,
          typename TC,
          Backend BACKEND, unsigned NUM_MODULI,
          typename TAlpha         = TC,
          typename TBeta          = TC,
          cublasFillMode_t UPLO_A = CUBLAS_FILL_MODE_FULL,
          cublasFillMode_t UPLO_B = CUBLAS_FILL_MODE_FULL,
          cublasFillMode_t UPLO_C = CUBLAS_FILL_MODE_FULL>
inline void undo_scaling(
    cudaStream_t stream,
    cublasOperation_t op_A, cublasOperation_t op_B,
    size_t m, size_t n,
    common::mid_t<BACKEND, common::isComplex<TC>> *C_mid,
    const size_t ldc_mid, const size_t incC_mid,
    TC *const C, const size_t ldc,
    const int16_t *const sftA, const int16_t *const sftB,
    const TAlpha *alpha, const TBeta *beta //
) {
    if constexpr (FUNC == Func::syr2k) {

        undo_scaling::undo_scaling_syr2k<TC, TAlpha, TBeta, BACKEND, NUM_MODULI, UPLO_C>(
            stream, m, C_mid, ldc_mid, incC_mid, C, ldc, sftA, sftB, alpha, beta);

    } else if constexpr (FUNC == Func::her2k) {

        undo_scaling::undo_scaling_her2k<TC, TAlpha, TBeta, BACKEND, NUM_MODULI, UPLO_C>(
            stream, m, C_mid, ldc_mid, incC_mid, C, ldc, sftA, sftB, alpha, beta);

    } else if constexpr (FUNC == Func::herk || FUNC == Func::herkx) {

        undo_scaling::undo_scaling<TC, TAlpha, TBeta, BACKEND, NUM_MODULI, UPLO_C, false, true>(
            stream, m, n, C_mid, ldc_mid, incC_mid, C, ldc, sftA, sftB, alpha, beta);

    } else if constexpr (FUNC == Func::trtrmm) {

        if (op_A == CUBLAS_OP_N) {
            if (op_B == CUBLAS_OP_N) {
                constexpr cublasFillMode_t UPLO = ((UPLO_A == UPLO_B)) ? UPLO_A : UPLO_C;
                undo_scaling::undo_scaling<TC, TAlpha, TBeta, BACKEND, NUM_MODULI, UPLO, true>(
                    stream, m, n, C_mid, ldc_mid, incC_mid, C, ldc, sftA, sftB, alpha, beta);
            } else {
                constexpr cublasFillMode_t UPLO = ((UPLO_A == flip_uplo<UPLO_B>)) ? UPLO_A : UPLO_C;
                undo_scaling::undo_scaling<TC, TAlpha, TBeta, BACKEND, NUM_MODULI, UPLO, true>(
                    stream, m, n, C_mid, ldc_mid, incC_mid, C, ldc, sftA, sftB, alpha, beta);
            }
        } else {
            if (op_B == CUBLAS_OP_N) {
                constexpr cublasFillMode_t UPLO = ((flip_uplo<UPLO_A> == UPLO_B)) ? flip_uplo<UPLO_A> : UPLO_C;
                undo_scaling::undo_scaling<TC, TAlpha, TBeta, BACKEND, NUM_MODULI, UPLO, true>(
                    stream, m, n, C_mid, ldc_mid, incC_mid, C, ldc, sftA, sftB, alpha, beta);
            } else {
                constexpr cublasFillMode_t UPLO = ((flip_uplo<UPLO_A> == flip_uplo<UPLO_B>)) ? flip_uplo<UPLO_A> : UPLO_C;
                undo_scaling::undo_scaling<TC, TAlpha, TBeta, BACKEND, NUM_MODULI, UPLO, true>(
                    stream, m, n, C_mid, ldc_mid, incC_mid, C, ldc, sftA, sftB, alpha, beta);
            }
        }

    } else {

        undo_scaling::undo_scaling<TC, TAlpha, TBeta, BACKEND, NUM_MODULI, UPLO_C>(
            stream, m, n, C_mid, ldc_mid, incC_mid, C, ldc, sftA, sftB, alpha, beta);
    }
}

} // namespace gemmul8::oz2::core
