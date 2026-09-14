#pragma once
#include "../common/common.hpp"

namespace gemmul8::mod {

template <Backend BACKEND, bool COMPLEX,
          cublasFillMode_t UPLO = CUBLAS_FILL_MODE_FULL>
void mod_hi2mid(
    const cudaStream_t stream,
    const unsigned idx,
    const size_t ldc, const unsigned n,
    common::matptr_t<common::hi_t<BACKEND>, COMPLEX> &C_hi,
    common::mid_t<BACKEND, COMPLEX> *C_mid);

template <Backend BACKEND, bool COMPLEX,
          cublasFillMode_t UPLO = CUBLAS_FILL_MODE_FULL>
void mod_hi2mid_batched(
    const cudaStream_t stream,
    const unsigned idx, const unsigned bcnt,
    const size_t ldc, const unsigned n,
    common::matptr_t<common::hi_t<BACKEND>, COMPLEX> C_hi,
    common::mid_t<BACKEND, COMPLEX> *C_mid,
    const size_t incC_hi);

template <Backend BACKEND, cublasFillMode_t UPLO, bool FLIP_IMAG = false>
void mod_hi2mid_AHA(
    const cudaStream_t stream,
    const unsigned idx,
    const size_t ldc, const unsigned n,
    common::matptr_t<common::hi_t<BACKEND>, true> &C_hi,
    common::mid_t<BACKEND, true> *C_mid);

} // namespace gemmul8::mod
