#pragma once
#include "../common/common.hpp"

namespace gemmul8::undo_scaling {

template <Backend BACKEND> struct crt_tail {
    const common::hi_t<BACKEND> *ptr0 = nullptr;
    const common::hi_t<BACKEND> *ptr1 = nullptr;
    bool flip_imag                    = false;
    bool grouped                      = false;
};

template <typename T, typename TAlpha, typename TBeta,
          Backend BACKEND, unsigned NUM_MODULI,
          cublasFillMode_t UPLO = CUBLAS_FILL_MODE_FULL,
          bool isTRTRMM         = false,
          bool HERMITIAN        = false>
void undo_scaling(
    const cudaStream_t stream,
    const unsigned m, const unsigned n,
    common::mid_t<BACKEND, common::isComplex<T>> *C_mid,
    const size_t ldc_mid, const size_t incC_mid,
    T *const C, const size_t ldc,
    const int16_t *const sftA, const int16_t *const sftB,
    const TAlpha *alpha, const TBeta *beta,
    const crt_tail<BACKEND> tail = {});

template <typename T, typename TAlpha, typename TBeta,
          Backend BACKEND, unsigned NUM_MODULI,
          cublasFillMode_t UPLO>
void undo_scaling_syr2k(
    const cudaStream_t stream,
    const unsigned n,
    common::mid_t<BACKEND, common::isComplex<T>> *C_mid,
    const size_t ldc_mid, const size_t incC_mid,
    T *const C, const size_t ldc,
    const int16_t *const sftA, const int16_t *const sftB,
    const TAlpha *alpha, const TBeta *beta,
    const crt_tail<BACKEND> tail = {});

template <typename T, typename TAlpha, typename TBeta,
          Backend BACKEND, unsigned NUM_MODULI,
          cublasFillMode_t UPLO>
void undo_scaling_her2k(
    const cudaStream_t stream,
    const unsigned n,
    common::mid_t<BACKEND, common::isComplex<T>> *C_mid,
    const size_t ldc_mid, const size_t incC_mid,
    T *const C, const size_t ldc,
    const int16_t *const sftA, const int16_t *const sftB,
    const TAlpha *alpha, const TBeta *beta,
    const crt_tail<BACKEND> tail = {});

} // namespace gemmul8::undo_scaling
