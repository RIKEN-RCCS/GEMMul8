#pragma once
#include "../common/common.hpp"
#include "../mod/mod_hi2mid_declaration.hpp"
#include "helper_triangular.hpp"
#include "product_workspace.hpp"

namespace gemmul8::oz2::core {

template <Func FUNC,
          Backend BACKEND,
          bool COMPLEX,
          cublasFillMode_t UPLO_A = CUBLAS_FILL_MODE_FULL,
          cublasFillMode_t UPLO_B = CUBLAS_FILL_MODE_FULL,
          cublasFillMode_t UPLO_C = CUBLAS_FILL_MODE_FULL>
inline void mod_hi2mid(
    const cudaStream_t stream,
    cublasOperation_t op_A, cublasOperation_t op_B,
    const unsigned idx, const unsigned bcnt,
    const size_t ldc, const unsigned n,
    common::matptr_t<common::hi_t<BACKEND>, COMPLEX> &C_hi,
    common::mid_t<BACKEND, COMPLEX> *C_mid,
    const unsigned num_moduli //
) {
    using Layout       = product_workspace<BACKEND, COMPLEX, FUNC == Func::herk>;
    const size_t sizeC = ldc * size_t(n);
    const auto advance = [&](unsigned i) {
        C_hi.ptr0 += sizeC * Layout::planes(num_moduli, i, 1U);
        if constexpr (COMPLEX && FUNC != Func::herk) {
            C_hi.ptr1 = i + 1U < idx + bcnt
                            ? C_hi.ptr0 + sizeC * Layout::products(num_moduli, i + 1U)
                            : nullptr;
        }
    };

    if constexpr (FUNC == Func::herk) {

        if (op_A == CUBLAS_OP_N) {
            for (unsigned b = 0; b < bcnt; ++b) {
                mod::mod_hi2mid_AHA<BACKEND, UPLO_C, true>(stream, idx + b, ldc, n, C_hi, C_mid, num_moduli);
                advance(idx + b);
            }
        } else {
            for (unsigned b = 0; b < bcnt; ++b) {
                mod::mod_hi2mid_AHA<BACKEND, UPLO_C, false>(stream, idx + b, ldc, n, C_hi, C_mid, num_moduli);
                advance(idx + b);
            }
        }

    } else if constexpr (FUNC == Func::syr2k || FUNC == Func::her2k) {

        for (unsigned b = 0; b < bcnt; ++b) {
            mod::mod_hi2mid<BACKEND, COMPLEX, CUBLAS_FILL_MODE_FULL>(stream, idx + b, ldc, n, C_hi, C_mid, num_moduli);
            advance(idx + b);
        }

    } else if constexpr (FUNC == Func::trtrmm) {

        if (op_A == CUBLAS_OP_N) {
            if (op_B == CUBLAS_OP_N) {
                constexpr cublasFillMode_t UPLO = (UPLO_A == UPLO_B) ? UPLO_A : UPLO_C;
                for (unsigned b = 0; b < bcnt; ++b) {
                    mod::mod_hi2mid<BACKEND, COMPLEX, UPLO>(stream, idx + b, ldc, n, C_hi, C_mid, num_moduli);
                    advance(idx + b);
                }
            } else {
                constexpr cublasFillMode_t UPLO = ((UPLO_A == flip_uplo<UPLO_B>)) ? UPLO_A : UPLO_C;
                for (unsigned b = 0; b < bcnt; ++b) {
                    mod::mod_hi2mid<BACKEND, COMPLEX, UPLO>(stream, idx + b, ldc, n, C_hi, C_mid, num_moduli);
                    advance(idx + b);
                }
            }
        } else {
            if (op_B == CUBLAS_OP_N) {
                constexpr cublasFillMode_t UPLO = ((flip_uplo<UPLO_A> == UPLO_B)) ? flip_uplo<UPLO_A> : UPLO_C;
                for (unsigned b = 0; b < bcnt; ++b) {
                    mod::mod_hi2mid<BACKEND, COMPLEX, UPLO>(stream, idx + b, ldc, n, C_hi, C_mid, num_moduli);
                    advance(idx + b);
                }
            } else {
                constexpr cublasFillMode_t UPLO = ((flip_uplo<UPLO_A> == flip_uplo<UPLO_B>)) ? flip_uplo<UPLO_A> : UPLO_C;
                for (unsigned b = 0; b < bcnt; ++b) {
                    mod::mod_hi2mid<BACKEND, COMPLEX, UPLO>(stream, idx + b, ldc, n, C_hi, C_mid, num_moduli);
                    advance(idx + b);
                }
            }
        }

    } else {

        for (unsigned b = 0; b < bcnt; ++b) {
            mod::mod_hi2mid<BACKEND, COMPLEX, UPLO_C>(stream, idx + b, ldc, n, C_hi, C_mid, num_moduli);
            advance(idx + b);
        }
    }
}

} // namespace gemmul8::oz2::core
