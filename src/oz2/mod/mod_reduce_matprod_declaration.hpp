#pragma once
#include "../common/common.hpp"

namespace gemmul8::mod {

template <Backend BACKEND>
void mod_reduce_matprod(
    const cudaStream_t stream,
    common::hi_t<BACKEND> *C,
    const int m,
    const int n,
    const size_t ldc,
    const unsigned modulus_idx //
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
    const unsigned modulus_idx //
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
    const unsigned modulus_idx //
);

} // namespace gemmul8::mod
