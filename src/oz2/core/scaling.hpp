#pragma once
#include "../common/common.hpp"
#include "../scaling/fast/scaling_fast_declaration.hpp"
#include "../scaling/accu/scaling_accu_declaration.hpp"
#include "../scaling/accu/calc_sft.hpp"
#include "helper_triangular.hpp"
#include "matmult.hpp"

namespace gemmul8::oz2::core {

namespace {
template <Func FUNC,
          typename T,
          Backend BACKEND,
          common::MatStruct STRUCT_A,
          common::MatStruct STRUCT_B,
          cublasFillMode_t UPLO_A,
          cublasFillMode_t UPLO_B,
          cublasFillMode_t UPLO_C>
inline void scaling_accu_make_C_hi(
    const cudaStream_t stream,
    common::Handle_t &handle,
    const size_t ldc_hi,
    const size_t n,
    const size_t lda_lo,
    const size_t ldb_lo,
    common::matptr_t<common::low_t<BACKEND>, common::isComplex<T>> &A_lo_high,
    common::matptr_t<common::low_t<BACKEND>, common::isComplex<T>> &B_lo_high,
    common::matptr_t<common::hi_t<BACKEND>, common::isComplex<T>> &C_hi //
) {
    constexpr common::MatMulKind KIND =
        (FUNC == Func::syrk || FUNC == Func::herk)
            ? common::MatMulKind::ATxB
            : matmul_kind<FUNC, STRUCT_A, STRUCT_B, UPLO_C>();

    using HiT = common::hi_t<BACKEND>;

    constexpr HiT one  = HiT(1);
    constexpr HiT zero = HiT(0);

    if constexpr (!common::isComplex<T>) {

        matmul_block_1<KIND, BACKEND, UPLO_A, UPLO_B, UPLO_C>(
            stream, handle, ldc_hi, n, lda_lo, ldb_lo,
            &one, A_lo_high.ptr0, B_lo_high.ptr0,
            &zero, C_hi.ptr0);

    } else if constexpr (BACKEND == Backend::INT8) {

        matmul_block_3<KIND, BACKEND, UPLO_A, UPLO_B, UPLO_C>(
            stream, handle, ldc_hi, n, lda_lo, ldb_lo,
            &one, &one, &one,
            A_lo_high.ptr0, A_lo_high.ptr1, A_lo_high.ptr2,
            B_lo_high.ptr1, B_lo_high.ptr0, B_lo_high.ptr2,
            &zero, &one, &zero,
            C_hi.ptr1, C_hi.ptr1, C_hi.ptr0);

    } else {

        matmul_block_3<KIND, BACKEND, UPLO_A, UPLO_B, UPLO_C>(
            stream, handle, ldc_hi, n, lda_lo, ldb_lo,
            &one, &one, &one,
            A_lo_high.ptr0, A_lo_high.ptr1, A_lo_high.ptr2,
            B_lo_high.ptr1, B_lo_high.ptr0, B_lo_high.ptr2,
            &zero, &zero, &zero,
            C_hi.ptr1, C_hi.ptr2, C_hi.ptr0);
    }
}

template <Func FUNC,
          typename T,
          Backend BACKEND,
          common::MatStruct STRUCT_A,
          common::MatStruct STRUCT_B,
          cublasFillMode_t UPLO_A,
          cublasFillMode_t UPLO_B,
          cublasFillMode_t UPLO_C>
inline void scaling_accu_make_C_hi_launch(
    const cudaStream_t stream,
    common::Handle_t &handle,
    const size_t ldc_hi,
    const size_t n,
    const size_t lda_lo,
    const size_t ldb_lo,
    cublasOperation_t op_A, cublasOperation_t op_B,
    common::matptr_t<common::low_t<BACKEND>, common::isComplex<T>> &A_lo_high,
    common::matptr_t<common::low_t<BACKEND>, common::isComplex<T>> &B_lo_high,
    common::matptr_t<common::hi_t<BACKEND>, common::isComplex<T>> &C_hi //
) {
    constexpr common::MatMulKind KIND = matmul_kind<FUNC, STRUCT_A, STRUCT_B, UPLO_C>();

    if constexpr (KIND == common::MatMulKind::TrmmLeft) {
        if (op_A == CUBLAS_OP_N) {
            scaling_accu_make_C_hi<FUNC, T, BACKEND, STRUCT_A, STRUCT_B, UPLO_A, UPLO_B, UPLO_C>(
                stream, handle, ldc_hi, n, lda_lo, ldb_lo,
                A_lo_high, B_lo_high, C_hi);
        } else {
            scaling_accu_make_C_hi<FUNC, T, BACKEND, STRUCT_A, STRUCT_B, flip_uplo<UPLO_A>, UPLO_B, UPLO_C>(
                stream, handle, ldc_hi, n, lda_lo, ldb_lo,
                A_lo_high, B_lo_high, C_hi);
        }
    } else if constexpr (KIND == common::MatMulKind::TrmmRight) {
        if (op_B == CUBLAS_OP_N) {
            scaling_accu_make_C_hi<FUNC, T, BACKEND, STRUCT_A, STRUCT_B, UPLO_A, UPLO_B, UPLO_C>(
                stream, handle, ldc_hi, n, lda_lo, ldb_lo,
                A_lo_high, B_lo_high, C_hi);
        } else {
            scaling_accu_make_C_hi<FUNC, T, BACKEND, STRUCT_A, STRUCT_B, UPLO_A, flip_uplo<UPLO_B>, UPLO_C>(
                stream, handle, ldc_hi, n, lda_lo, ldb_lo,
                A_lo_high, B_lo_high, C_hi);
        }
    } else if constexpr (KIND == common::MatMulKind::Trtrmm) {
        if (op_A == CUBLAS_OP_N) {
            if (op_B == CUBLAS_OP_N) {
                scaling_accu_make_C_hi<FUNC, T, BACKEND, STRUCT_A, STRUCT_B, UPLO_A, UPLO_B, UPLO_C>(
                    stream, handle, ldc_hi, n, lda_lo, ldb_lo,
                    A_lo_high, B_lo_high, C_hi);
            } else {
                scaling_accu_make_C_hi<FUNC, T, BACKEND, STRUCT_A, STRUCT_B, UPLO_A, flip_uplo<UPLO_B>, UPLO_C>(
                    stream, handle, ldc_hi, n, lda_lo, ldb_lo,
                    A_lo_high, B_lo_high, C_hi);
            }
        } else {
            if (op_B == CUBLAS_OP_N) {
                scaling_accu_make_C_hi<FUNC, T, BACKEND, STRUCT_A, STRUCT_B, flip_uplo<UPLO_A>, UPLO_B, UPLO_C>(
                    stream, handle, ldc_hi, n, lda_lo, ldb_lo,
                    A_lo_high, B_lo_high, C_hi);
            } else {
                scaling_accu_make_C_hi<FUNC, T, BACKEND, STRUCT_A, STRUCT_B, flip_uplo<UPLO_A>, flip_uplo<UPLO_B>, UPLO_C>(
                    stream, handle, ldc_hi, n, lda_lo, ldb_lo,
                    A_lo_high, B_lo_high, C_hi);
            }
        }
    } else {
        scaling_accu_make_C_hi<FUNC, T, BACKEND, STRUCT_A, STRUCT_B, UPLO_A, UPLO_B, UPLO_C>(
            stream, handle, ldc_hi, n, lda_lo, ldb_lo,
            A_lo_high, B_lo_high, C_hi);
    }
}

} // namespace

// gemm, symm, syr2k, syrkx, trmm, hemm, her2k, herkx, trtrmm
template <typename TA, typename TB,
          Backend BACKEND, unsigned NUM_MODULI,
          common::MatStruct STRUCT_A = common::MatStruct::Full,
          common::MatStruct STRUCT_B = common::MatStruct::Full,
          cublasFillMode_t UPLO_A    = CUBLAS_FILL_MODE_FULL,
          cublasFillMode_t UPLO_B    = CUBLAS_FILL_MODE_FULL,
          cublasDiagType_t DIAG_A    = CUBLAS_DIAG_NON_UNIT,
          cublasDiagType_t DIAG_B    = CUBLAS_DIAG_NON_UNIT>
inline void scaling_fast(
    const cudaStream_t stream,
    const cublasOperation_t op_A, const cublasOperation_t op_B,
    const unsigned m, const unsigned n, const unsigned k,
    const TA *const A, const size_t lda,
    common::matptr_t<common::low_t<BACKEND>, common::isComplex<TA>> &A_lo,
    const size_t lda_lo, const size_t incA_lo,
    int16_t *const sftA, const bool skip_scalA,
    const TB *const B, const size_t ldb,
    common::matptr_t<common::low_t<BACKEND>, common::isComplex<TA>> &B_lo,
    const size_t ldb_lo, const size_t incB_lo,
    int16_t *const sftB, const bool skip_scalB //
) {
    static_assert(common::isComplex<TA> == common::isComplex<TB>,
                  "scaling_fast currently assumes TA and TB have the same real/complex category.");

    // A -> A_lo
    if (!skip_scalA) {
        if constexpr (common::is_symmetric<STRUCT_A>) {

            scaling::fast::scaling_symm<TA, BACKEND, NUM_MODULI, UPLO_A>(
                stream, m, A, lda, A_lo, lda_lo, incA_lo, sftA);

        } else if constexpr (common::is_hermitian<STRUCT_A>) {

            static_assert(common::isComplex<TA>,
                          "Hermitian scaling requires complex input type.");

            scaling::fast::scaling_hemm<TA, BACKEND, NUM_MODULI, UPLO_A, true>(
                stream, m, A, lda, A_lo, lda_lo, incA_lo, sftA);

        } else {

            scaling::fast::scaling<TA, BACKEND, NUM_MODULI, UPLO_A, DIAG_A>(
                stream, op_A, CUBLAS_SIDE_LEFT, m, k, A, lda, A_lo, lda_lo, incA_lo, sftA);
        }
    }

    // B -> B_lo
    if (!skip_scalB) {
        if constexpr (common::is_symmetric<STRUCT_B>) {

            scaling::fast::scaling_symm<TB, BACKEND, NUM_MODULI, UPLO_B>(
                stream, n, B, ldb, B_lo, ldb_lo, incB_lo, sftB);

        } else if constexpr (common::is_hermitian<STRUCT_B>) {

            static_assert(common::isComplex<TB>,
                          "Hermitian scaling requires complex input type.");

            scaling::fast::scaling_hemm<TB, BACKEND, NUM_MODULI, UPLO_B, false>(
                stream, n, B, ldb, B_lo, ldb_lo, incB_lo, sftB);

        } else {

            scaling::fast::scaling<TB, BACKEND, NUM_MODULI, UPLO_B, DIAG_B>(
                stream, op_B, CUBLAS_SIDE_RIGHT, k, n, B, ldb, B_lo, ldb_lo, incB_lo, sftB);
        }
    }
}

// gemm, symm, syr2k, syrkx, trmm, hemm, her2k, herkx, trtrmm
template <Func FUNC,
          typename TA, typename TB,
          Backend BACKEND, unsigned NUM_MODULI,
          common::MatStruct STRUCT_A = common::MatStruct::Full,
          common::MatStruct STRUCT_B = common::MatStruct::Full,
          cublasFillMode_t UPLO_A    = CUBLAS_FILL_MODE_FULL,
          cublasFillMode_t UPLO_B    = CUBLAS_FILL_MODE_FULL,
          cublasDiagType_t DIAG_A    = CUBLAS_DIAG_NON_UNIT,
          cublasDiagType_t DIAG_B    = CUBLAS_DIAG_NON_UNIT,
          cublasFillMode_t UPLO_C    = CUBLAS_FILL_MODE_FULL>
inline void scaling_accu(
    const cudaStream_t stream, common::Handle_t &handle,
    const cublasOperation_t op_A, const cublasOperation_t op_B,
    const unsigned m, const unsigned n, const unsigned k,
    const TA *const A, const size_t lda,
    common::matptr_t<common::low_t<BACKEND>, common::isComplex<TA>> &A_lo,
    common::matptr_t<common::low_t<BACKEND>, common::isComplex<TA>> &A_lo_high,
    const size_t lda_lo, const size_t incA_lo,
    int16_t *const sftA, const bool skip_scalA, const bool enable_skip_scalA,
    const TB *const B, const size_t ldb,
    common::matptr_t<common::low_t<BACKEND>, common::isComplex<TA>> &B_lo,
    common::matptr_t<common::low_t<BACKEND>, common::isComplex<TA>> &B_lo_high,
    const size_t ldb_lo, const size_t incB_lo,
    int16_t *const sftB, const bool skip_scalB, const bool enable_skip_scalB,
    common::matptr_t<common::hi_t<BACKEND>, common::isComplex<TA>> &C_hi, const size_t ldc_hi //
) {
    static_assert(common::isComplex<TA> == common::isComplex<TB>,
                  "scaling_accu currently assumes TA and TB have the same real/complex category.");

    constexpr cublasFillMode_t SCAL_UPLO_C =
        (FUNC == Func::syr2k || FUNC == Func::her2k) ? CUBLAS_FILL_MODE_FULL : UPLO_C;

    if (skip_scalA && skip_scalB) return;

    const size_t size_vecA = common::padding(m);
    const size_t size_vecB = common::padding(n);

    int16_t *const sftA_delta = (enable_skip_scalA) ? (sftA + size_vecA) : nullptr;
    int16_t *const sftB_delta = (enable_skip_scalB) ? (sftB + size_vecB) : nullptr;

    // A -> A_lo
    if (!skip_scalA) {

        if constexpr (common::is_symmetric<STRUCT_A>) {

            scaling::accu::extract_symm<TA, BACKEND, UPLO_A>(
                stream, m, A, lda, A_lo_high, lda_lo, incA_lo, sftA);

        } else if constexpr (common::is_hermitian<STRUCT_A>) {

            scaling::accu::extract_hemm<TA, BACKEND, UPLO_A, true>(
                stream, m, A, lda, A_lo_high, lda_lo, incA_lo, sftA);

        } else {
            scaling::accu::extract<TA, BACKEND, UPLO_A, DIAG_A>(
                stream, op_A, CUBLAS_SIDE_LEFT, m, k, A, lda, A_lo_high, lda_lo, incA_lo, sftA);
        }

        if (enable_skip_scalA) {
            cudaMemcpyAsync(sftA_delta, sftA, size_vecA * sizeof(int16_t),
                            cudaMemcpyDeviceToDevice, stream);
        }
    }

    // B -> B_lo
    if (!skip_scalB) {

        if constexpr (common::is_symmetric<STRUCT_B>) {

            scaling::accu::extract_symm<TB, BACKEND, UPLO_B>(
                stream, n, B, ldb, B_lo_high, ldb_lo, incB_lo, sftB);

        } else if constexpr (common::is_hermitian<STRUCT_B>) {

            scaling::accu::extract_hemm<TB, BACKEND, UPLO_B, false>(
                stream, n, B, ldb, B_lo_high, ldb_lo, incB_lo, sftB);

        } else {
            scaling::accu::extract<TB, BACKEND, UPLO_B, DIAG_B>(
                stream, op_B, CUBLAS_SIDE_RIGHT, k, n, B, ldb, B_lo_high, ldb_lo, incB_lo, sftB);
        }

        if (enable_skip_scalB) {
            cudaMemcpyAsync(sftB_delta, sftB, size_vecB * sizeof(int16_t),
                            cudaMemcpyDeviceToDevice, stream);
        }
    }

    // C_hi <- A_lo_high * B_lo_high
    const size_t n_mat = (FUNC == Func::trtrmm) ? common::padding(n) : n;
    scaling_accu_make_C_hi_launch<FUNC, TA, BACKEND, STRUCT_A, STRUCT_B, UPLO_A, UPLO_B, UPLO_C>(
        stream, handle, ldc_hi, n_mat, lda_lo, ldb_lo, op_A, op_B, A_lo_high, B_lo_high, C_hi);

    // A -> A_lo
    if (!skip_scalA) {

        const int16_t *const fixed_deltaB = (skip_scalB) ? sftB_delta : nullptr;

        if constexpr (common::is_symmetric<STRUCT_A>) {

            scaling::accu::scaling_symm<TA, BACKEND, NUM_MODULI, UPLO_A>(
                stream, CUBLAS_SIDE_LEFT, m, n, k, A, lda, A_lo, lda_lo, incA_lo, sftA, C_hi, ldc_hi, fixed_deltaB);

        } else if constexpr (common::is_hermitian<STRUCT_A>) {

            scaling::accu::scaling_hemm<TA, BACKEND, NUM_MODULI, UPLO_A, true>(
                stream, CUBLAS_SIDE_LEFT, m, n, k, A, lda, A_lo, lda_lo, incA_lo, sftA, C_hi, ldc_hi, fixed_deltaB);

        } else {

            if constexpr (FUNC == Func::trtrmm) {
                if (op_A == CUBLAS_OP_N) {
                    if (op_B == CUBLAS_OP_N) {
                        constexpr cublasFillMode_t UPLO = (UPLO_A == UPLO_B) ? UPLO_A : UPLO_C;
                        scaling::accu::scaling<TA, BACKEND, NUM_MODULI, UPLO_A, DIAG_A, UPLO>(
                            stream, op_A, CUBLAS_SIDE_LEFT, m, n, k, A, lda, A_lo, lda_lo, incA_lo, sftA, C_hi, ldc_hi, fixed_deltaB);
                    } else {
                        constexpr cublasFillMode_t UPLO = ((UPLO_A == flip_uplo<UPLO_B>)) ? UPLO_A : UPLO_C;
                        scaling::accu::scaling<TA, BACKEND, NUM_MODULI, UPLO_A, DIAG_A, UPLO>(
                            stream, op_A, CUBLAS_SIDE_LEFT, m, n, k, A, lda, A_lo, lda_lo, incA_lo, sftA, C_hi, ldc_hi, fixed_deltaB);
                    }
                } else {
                    if (op_B == CUBLAS_OP_N) {
                        constexpr cublasFillMode_t UPLO = ((flip_uplo<UPLO_A> == UPLO_B)) ? flip_uplo<UPLO_A> : UPLO_C;
                        scaling::accu::scaling<TA, BACKEND, NUM_MODULI, UPLO_A, DIAG_A, UPLO>(
                            stream, op_A, CUBLAS_SIDE_LEFT, m, n, k, A, lda, A_lo, lda_lo, incA_lo, sftA, C_hi, ldc_hi, fixed_deltaB);
                    } else {
                        constexpr cublasFillMode_t UPLO = ((flip_uplo<UPLO_A> == flip_uplo<UPLO_B>)) ? flip_uplo<UPLO_A> : UPLO_C;
                        scaling::accu::scaling<TA, BACKEND, NUM_MODULI, UPLO_A, DIAG_A, UPLO>(
                            stream, op_A, CUBLAS_SIDE_LEFT, m, n, k, A, lda, A_lo, lda_lo, incA_lo, sftA, C_hi, ldc_hi, fixed_deltaB);
                    }
                }
            } else {
                scaling::accu::scaling<TA, BACKEND, NUM_MODULI, UPLO_A, DIAG_A, SCAL_UPLO_C>(
                    stream, op_A, CUBLAS_SIDE_LEFT, m, n, k, A, lda, A_lo, lda_lo, incA_lo, sftA, C_hi, ldc_hi, fixed_deltaB);
            }
        }

        if (enable_skip_scalA) {
            constexpr unsigned threads_delta = 256;
            const unsigned grid_delta        = (m + threads_delta - 1) / threads_delta;
            scaling::accu::calc_sft_delta<BACKEND, NUM_MODULI, common::isComplex<TA>>
                <<<grid_delta, threads_delta, 0, stream>>>(m, sftA, sftA_delta);
        }
    }

    // B -> B_lo
    if (!skip_scalB) {

        const int16_t *const fixed_deltaA = (skip_scalA) ? sftA_delta : nullptr;

        if constexpr (common::is_symmetric<STRUCT_B>) {

            scaling::accu::scaling_symm<TB, BACKEND, NUM_MODULI, UPLO_B>(
                stream, CUBLAS_SIDE_RIGHT, m, n, k, B, ldb, B_lo, ldb_lo, incB_lo, sftB, C_hi, ldc_hi, fixed_deltaA);

        } else if constexpr (common::is_hermitian<STRUCT_B>) {

            scaling::accu::scaling_hemm<TB, BACKEND, NUM_MODULI, UPLO_B, false>(
                stream, CUBLAS_SIDE_RIGHT, m, n, k, B, ldb, B_lo, ldb_lo, incB_lo, sftB, C_hi, ldc_hi, fixed_deltaA);

        } else {

            if constexpr (FUNC == Func::trtrmm) {
                if (op_A == CUBLAS_OP_N) {
                    if (op_B == CUBLAS_OP_N) {
                        constexpr cublasFillMode_t UPLO = (UPLO_A == UPLO_B) ? UPLO_A : UPLO_C;
                        scaling::accu::scaling<TB, BACKEND, NUM_MODULI, UPLO_B, DIAG_B, UPLO>(
                            stream, op_B, CUBLAS_SIDE_RIGHT, m, n, k, B, ldb, B_lo, ldb_lo, incB_lo, sftB, C_hi, ldc_hi, fixed_deltaA);
                    } else {
                        constexpr cublasFillMode_t UPLO = ((UPLO_A == flip_uplo<UPLO_B>)) ? UPLO_A : UPLO_C;
                        scaling::accu::scaling<TB, BACKEND, NUM_MODULI, UPLO_B, DIAG_B, UPLO>(
                            stream, op_B, CUBLAS_SIDE_RIGHT, m, n, k, B, ldb, B_lo, ldb_lo, incB_lo, sftB, C_hi, ldc_hi, fixed_deltaA);
                    }
                } else {
                    if (op_B == CUBLAS_OP_N) {
                        constexpr cublasFillMode_t UPLO = ((flip_uplo<UPLO_A> == UPLO_B)) ? flip_uplo<UPLO_A> : UPLO_C;
                        scaling::accu::scaling<TB, BACKEND, NUM_MODULI, UPLO_B, DIAG_B, UPLO>(
                            stream, op_B, CUBLAS_SIDE_RIGHT, m, n, k, B, ldb, B_lo, ldb_lo, incB_lo, sftB, C_hi, ldc_hi, fixed_deltaA);
                    } else {
                        constexpr cublasFillMode_t UPLO = ((flip_uplo<UPLO_A> == flip_uplo<UPLO_B>)) ? flip_uplo<UPLO_A> : UPLO_C;
                        scaling::accu::scaling<TB, BACKEND, NUM_MODULI, UPLO_B, DIAG_B, UPLO>(
                            stream, op_B, CUBLAS_SIDE_RIGHT, m, n, k, B, ldb, B_lo, ldb_lo, incB_lo, sftB, C_hi, ldc_hi, fixed_deltaA);
                    }
                }
            } else {
                scaling::accu::scaling<TB, BACKEND, NUM_MODULI, UPLO_B, DIAG_B, SCAL_UPLO_C>(
                    stream, op_B, CUBLAS_SIDE_RIGHT, m, n, k, B, ldb, B_lo, ldb_lo, incB_lo, sftB, C_hi, ldc_hi, fixed_deltaA);
            }
        }

        if (enable_skip_scalB) {
            constexpr unsigned threads_delta = 256;
            const unsigned grid_delta        = (n + threads_delta - 1) / threads_delta;
            scaling::accu::calc_sft_delta<BACKEND, NUM_MODULI, common::isComplex<TA>>
                <<<grid_delta, threads_delta, 0, stream>>>(n, sftB, sftB_delta);
        }
    }
}

// syrk, herk
template <Func FUNC,
          typename TA,
          Backend BACKEND, unsigned NUM_MODULI>
inline void scaling_fast_rk(
    const cudaStream_t stream,
    const cublasOperation_t op_A,
    const unsigned n, const unsigned k,
    const TA *const A, const size_t lda,
    common::matptr_t<common::low_t<BACKEND>, common::isComplex<TA>> &A_lo,
    const size_t lda_lo, const size_t incA_lo,
    int16_t *const sftA,
    const bool skip_scalA //
) {
    static_assert(FUNC == Func::syrk || FUNC == Func::herk,
                  "scaling_fast_rk supports only syrk/herk.");
    if constexpr (FUNC == Func::herk) {
        static_assert(common::isComplex<TA>, "herk requires complex input type.");
    }

    if (skip_scalA) return;

    const cublasOperation_t op_scale =
        (FUNC == Func::herk && op_A != CUBLAS_OP_N) ? CUBLAS_OP_T : op_A;

    scaling::fast::scaling<TA, BACKEND, NUM_MODULI, CUBLAS_FILL_MODE_FULL, CUBLAS_DIAG_NON_UNIT>(
        stream, op_scale, CUBLAS_SIDE_LEFT, n, k, A, lda, A_lo, lda_lo, incA_lo, sftA);
}

// syrk, herk
template <Func FUNC,
          typename TA,
          Backend BACKEND, unsigned NUM_MODULI,
          cublasFillMode_t UPLO_C = CUBLAS_FILL_MODE_FULL>
inline void scaling_accu_rk(
    const cudaStream_t stream, common::Handle_t &handle,
    const cublasOperation_t op_A,
    const unsigned n, const unsigned k,
    const TA *const A, const size_t lda,
    common::matptr_t<common::low_t<BACKEND>, common::isComplex<TA>> &A_lo,
    common::matptr_t<common::low_t<BACKEND>, common::isComplex<TA>> &A_lo_high,
    const size_t lda_lo, const size_t incA_lo,
    int16_t *const sftA,
    const bool skip_scalA,
    common::matptr_t<common::hi_t<BACKEND>, common::isComplex<TA>> &C_hi,
    const size_t ldc_hi //
) {
    static_assert(FUNC == Func::syrk || FUNC == Func::herk,
                  "scaling_accu_rk supports only syrk/herk.");
    if constexpr (FUNC == Func::herk) {
        static_assert(common::isComplex<TA>, "herk requires complex input type.");
    }
    static_assert(UPLO_C == CUBLAS_FILL_MODE_UPPER || UPLO_C == CUBLAS_FILL_MODE_LOWER,
                  "scaling_accu_rk requires UPLO_C = UPPER or LOWER.");

    if (skip_scalA) return;

    scaling::accu::extract<TA, BACKEND, CUBLAS_FILL_MODE_FULL, CUBLAS_DIAG_NON_UNIT>(
        stream, op_A, CUBLAS_SIDE_LEFT, n, k, A, lda, A_lo_high, lda_lo, incA_lo, sftA);

    // C_hi <- A_lo_high^T * A_lo_high
    scaling_accu_make_C_hi<FUNC, TA, BACKEND,
                           common::MatStruct::Full, common::MatStruct::Full,
                           CUBLAS_FILL_MODE_FULL, CUBLAS_FILL_MODE_FULL, UPLO_C>(
        stream, handle, ldc_hi, n, lda_lo, lda_lo, A_lo_high, A_lo_high, C_hi);

    if constexpr (FUNC == Func::syrk) {
        scaling::accu::scaling_syrk<TA, BACKEND, NUM_MODULI, UPLO_C>(
            stream, op_A, n, k, A, lda, A_lo, lda_lo, incA_lo, sftA, C_hi, ldc_hi);

    } else if constexpr (FUNC == Func::herk) {
        scaling::accu::scaling_herk<TA, BACKEND, NUM_MODULI, UPLO_C>(
            stream, op_A, n, k, A, lda, A_lo, lda_lo, incA_lo, sftA, C_hi, ldc_hi);
    }
}

} // namespace gemmul8::oz2::core
