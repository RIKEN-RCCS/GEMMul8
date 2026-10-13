#pragma once
#include "../common/common.hpp"

namespace gemmul8::mod {

template <Backend BACKEND>
void mod_reduce_matprod_fp8_batch(
    cudaStream_t stream,
    float *C,
    int m, int n,
    size_t ldc,
    int64_t strideC,
    int batchCount,
    unsigned p,
    unsigned mask //
);

template <Backend BACKEND>
void mod_reduce_matprod_fp8_pointer_batch(
    cudaStream_t stream,
    void **Aarray,
    void **Barray,
    void **Carray,
    int m, int n,
    size_t ldc,
    int batchCount,
    int k_advance,
    unsigned p,
    unsigned mask //
);

template <Backend BACKEND>
void mod_reduce_matprod(
    const cudaStream_t stream,
    common::hi_t<BACKEND> *C,
    const int m,
    const int n,
    const size_t ldc,
    const unsigned modulus_idx,
    const bool complex_moduli = false //
);

template <Backend BACKEND>
void mod_reduce_matprod_strided(
    const cudaStream_t stream,
    common::hi_t<BACKEND> *C,
    const int m,
    const int n,
    const size_t ldc,
    const int64_t strideC,
    const int batchCount,
    const unsigned modulus_idx,
    const bool complex_moduli = false //
);

template <Backend BACKEND>
void mod_reduce_matprod_pointer_and_advance(
    const cudaStream_t stream,
    void **Aarray,
    void **Barray,
    void **Carray,
    const int m,
    const int n,
    const size_t ldc,
    const int batchCount,
    const int k_advance,
    const unsigned modulus_idx,
    const bool complex_moduli = false //
);

} // namespace gemmul8::mod
