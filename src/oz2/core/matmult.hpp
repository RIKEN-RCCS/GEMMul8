#pragma once
#include "../common/common.hpp"
#include "../common/matmult.hpp"
#include "../common/table.hpp"
#include "helper_triangular.hpp"
#include "helper_matmult.hpp"
#include "product_workspace.hpp"

namespace gemmul8::oz2::core {

template <unsigned NUM_MODULI, class Layout>
inline unsigned batch_count(
    const int arch, const size_t n, const unsigned i,
    const size_t sizeC, const size_t lwork_blas, const size_t worksizeC //
) {
    if (arch == 121 || (arch == 90 && n > 2048)) return 1U;
    for (unsigned cnt = NUM_MODULI - i; cnt > 1U; --cnt) {
        if (Layout::required(NUM_MODULI, i, cnt, sizeC, lwork_blas) <= worksizeC)
            return cnt;
    }
    return 1U;
}

inline constexpr size_t K_BLOCK_INT8 = size_t(1 << 17);

template <Backend BACKEND, bool COMPLEX = false>
inline void configure_matprod_k_blocking(
    common::Handle_t &handle,
    const unsigned idx,
    const size_t k //
) {
    handle.modulus_idx     = BACKEND == Backend::FP8 ? common::fp8_plan::index(handle.fp8_num_moduli, idx, COMPLEX) : idx;
    handle.modulus_complex = COMPLEX;
    if constexpr (BACKEND == Backend::INT8) {
        constexpr int KB             = int(K_BLOCK_INT8);
        const bool is_mod256         = (!COMPLEX && common::table::moduli_int8[idx] == 256);
        handle.matprod_k_blocking    = !is_mod256 && k > size_t(KB);
        handle.matprod_k_block_first = KB;
        handle.matprod_k_block_next  = KB;
    } else {
        const int32_t p              = common::fp8_plan::modulus(handle.fp8_num_moduli, idx, COMPLEX);
        const int KB0                = int(common::fp8_plan::k_block(p, true));
        const int KB                 = int(common::fp8_plan::k_block(p, false));
        handle.matprod_k_blocking    = k > size_t(KB0);
        handle.matprod_k_block_first = KB0;
        handle.matprod_k_block_next  = KB;
    }
}

template <Backend BACKEND, bool COMPLEX = false>
inline bool needs_k_blocking(const unsigned idx, const size_t k, const unsigned num_moduli) {
    if constexpr (BACKEND == Backend::INT8) {
        if (!COMPLEX && common::table::moduli_int8[idx] == 256) { return false; }
        return k > K_BLOCK_INT8;
    } else {
        return k > common::fp8_plan::k_block(common::fp8_plan::modulus(num_moduli, idx, COMPLEX), true);
    }
}

template <Backend BACKEND, bool COMPLEX = false>
inline unsigned limit_batch_for_k(const unsigned i, const unsigned bcnt, const size_t k, const unsigned num_moduli) {
    unsigned safe = 0;

    while (safe < bcnt) {
        const unsigned idx = i + safe;
        if (needs_k_blocking<BACKEND, COMPLEX>(idx, k, num_moduli)) { break; }
        ++safe;
    }

    return (safe == 0) ? 1 : safe;
}

inline void upload_pointer_arrays(
    const cudaStream_t stream,
    common::Handle_t &handle,
    void *const *hA,
    void *const *hB,
    void *const *hC,
    const unsigned count //
) {
    cudaMemcpyAsync(handle.Aarray, hA, count * sizeof(void *), cudaMemcpyHostToDevice, stream);
    cudaMemcpyAsync(handle.Barray, hB, count * sizeof(void *), cudaMemcpyHostToDevice, stream);
    cudaMemcpyAsync(handle.Carray, hC, count * sizeof(void *), cudaMemcpyHostToDevice, stream);
}

template <common::MatMulKind KIND,
          Backend BACKEND,
          cublasFillMode_t UPLO_A,
          cublasFillMode_t UPLO_B,
          cublasFillMode_t UPLO_C>
inline void matmul_block_1(
    const cudaStream_t stream,
    common::Handle_t &handle,
    const size_t ldc_hi,
    const size_t n,
    const size_t lda_lo,
    const size_t ldb_lo,
    const common::hi_t<BACKEND> *alpha,
    const common::low_t<BACKEND> *A,
    const common::low_t<BACKEND> *B,
    const common::hi_t<BACKEND> *beta,
    common::hi_t<BACKEND> *C //
) {
    if constexpr (KIND == common::MatMulKind::TrmmRight) {
        common::block_matmul_1<KIND, BACKEND, UPLO_A, UPLO_B, UPLO_C>(
            stream, handle,
            static_cast<int>(ldc_hi),
            static_cast<int>(n),
            static_cast<int>(lda_lo),
            alpha,
            B, ldb_lo,
            A, lda_lo,
            beta,
            C, ldc_hi);
    } else {
        common::block_matmul_1<KIND, BACKEND, UPLO_A, UPLO_B, UPLO_C>(
            stream, handle,
            static_cast<int>(ldc_hi),
            static_cast<int>(n),
            static_cast<int>(lda_lo),
            alpha,
            A, lda_lo,
            B, ldb_lo,
            beta,
            C, ldc_hi);
    }
}

template <common::MatMulKind KIND,
          Backend BACKEND,
          cublasFillMode_t UPLO_A,
          cublasFillMode_t UPLO_B,
          cublasFillMode_t UPLO_C>
inline void matmul_block_1_strided_batched(
    const cudaStream_t stream,
    common::Handle_t &handle,
    const size_t ldc_hi,
    const size_t n,
    const size_t lda_lo,
    const size_t ldb_lo,
    const unsigned bcnt,
    const common::hi_t<BACKEND> *alpha,
    const common::low_t<BACKEND> *A,
    const int64_t strideA,
    const common::low_t<BACKEND> *B,
    const int64_t strideB,
    const common::hi_t<BACKEND> *beta,
    common::hi_t<BACKEND> *C,
    const int64_t strideC //
) {
    if constexpr (KIND == common::MatMulKind::TrmmRight) {
        common::block_matmul_1_strided_batched<KIND, BACKEND, UPLO_A, UPLO_B, UPLO_C>(
            stream, handle,
            static_cast<int>(ldc_hi),
            static_cast<int>(n),
            static_cast<int>(lda_lo),
            static_cast<int>(bcnt),
            alpha,
            B, ldb_lo, strideB,
            A, lda_lo, strideA,
            beta,
            C, ldc_hi, strideC);
    } else {
        common::block_matmul_1_strided_batched<KIND, BACKEND, UPLO_A, UPLO_B, UPLO_C>(
            stream, handle,
            static_cast<int>(ldc_hi),
            static_cast<int>(n),
            static_cast<int>(lda_lo),
            static_cast<int>(bcnt),
            alpha,
            A, lda_lo, strideA,
            B, ldb_lo, strideB,
            beta,
            C, ldc_hi, strideC);
    }
}

template <common::MatMulKind KIND,
          Backend BACKEND,
          cublasFillMode_t UPLO_A,
          cublasFillMode_t UPLO_B,
          cublasFillMode_t UPLO_C>
inline void matmul_block_3(
    const cudaStream_t stream,
    common::Handle_t &handle,
    const size_t ldc_hi,
    const size_t n,
    const size_t lda_lo,
    const size_t ldb_lo,
    const common::hi_t<BACKEND> *alpha1,
    const common::hi_t<BACKEND> *alpha2,
    const common::hi_t<BACKEND> *alpha3,
    const common::low_t<BACKEND> *A1,
    const common::low_t<BACKEND> *A2,
    const common::low_t<BACKEND> *A3,
    const common::low_t<BACKEND> *B1,
    const common::low_t<BACKEND> *B2,
    const common::low_t<BACKEND> *B3,
    const common::hi_t<BACKEND> *beta1,
    const common::hi_t<BACKEND> *beta2,
    const common::hi_t<BACKEND> *beta3,
    common::hi_t<BACKEND> *C1,
    common::hi_t<BACKEND> *C2,
    common::hi_t<BACKEND> *C3 //
) {
    if constexpr (KIND == common::MatMulKind::TrmmRight) {
        common::block_matmul_3<KIND, BACKEND, UPLO_A, UPLO_B, UPLO_C>(
            stream, handle,
            static_cast<int>(ldc_hi),
            static_cast<int>(n),
            static_cast<int>(lda_lo),
            alpha1, alpha2, alpha3,
            B1, B2, B3, ldb_lo,
            A1, A2, A3, lda_lo,
            beta1, beta2, beta3,
            C1, C2, C3, ldc_hi);
    } else {
        common::block_matmul_3<KIND, BACKEND, UPLO_A, UPLO_B, UPLO_C>(
            stream, handle,
            static_cast<int>(ldc_hi),
            static_cast<int>(n),
            static_cast<int>(lda_lo),
            alpha1, alpha2, alpha3,
            A1, A2, A3, lda_lo,
            B1, B2, B3, ldb_lo,
            beta1, beta2, beta3,
            C1, C2, C3, ldc_hi);
    }
}

template <common::MatMulKind KIND,
          Backend BACKEND,
          cublasFillMode_t UPLO_A,
          cublasFillMode_t UPLO_B,
          cublasFillMode_t UPLO_C>
inline void matmul_block_3_strided_batched(
    const cudaStream_t stream,
    common::Handle_t &handle,
    const size_t ldc_hi,
    const size_t n,
    const size_t lda_lo,
    const size_t ldb_lo,
    const unsigned bcnt,
    const common::hi_t<BACKEND> *alpha1,
    const common::hi_t<BACKEND> *alpha2,
    const common::hi_t<BACKEND> *alpha3,
    const common::low_t<BACKEND> *A1,
    const common::low_t<BACKEND> *A2,
    const common::low_t<BACKEND> *A3,
    const int64_t strideA,
    const common::low_t<BACKEND> *B1,
    const common::low_t<BACKEND> *B2,
    const common::low_t<BACKEND> *B3,
    const int64_t strideB,
    const common::hi_t<BACKEND> *beta1,
    const common::hi_t<BACKEND> *beta2,
    const common::hi_t<BACKEND> *beta3,
    common::hi_t<BACKEND> *C1,
    common::hi_t<BACKEND> *C2,
    common::hi_t<BACKEND> *C3,
    const int64_t strideC //
) {
    if constexpr (KIND == common::MatMulKind::TrmmRight) {
        common::block_matmul_3_strided_batched<KIND, BACKEND, UPLO_A, UPLO_B, UPLO_C>(
            stream, handle,
            static_cast<int>(ldc_hi),
            static_cast<int>(n),
            static_cast<int>(lda_lo),
            static_cast<int>(bcnt),
            alpha1, alpha2, alpha3,
            B1, B2, B3, ldb_lo, strideB,
            A1, A2, A3, lda_lo, strideA,
            beta1, beta2, beta3,
            C1, C2, C3, ldc_hi, strideC);
    } else {
        common::block_matmul_3_strided_batched<KIND, BACKEND, UPLO_A, UPLO_B, UPLO_C>(
            stream, handle,
            static_cast<int>(ldc_hi),
            static_cast<int>(n),
            static_cast<int>(lda_lo),
            static_cast<int>(bcnt),
            alpha1, alpha2, alpha3,
            A1, A2, A3, lda_lo, strideA,
            B1, B2, B3, ldb_lo, strideB,
            beta1, beta2, beta3,
            C1, C2, C3, ldc_hi, strideC);
    }
}

template <common::MatMulKind KIND,
          cublasFillMode_t UPLO_A,
          cublasFillMode_t UPLO_B,
          cublasFillMode_t UPLO_C>
inline void error_free_matmult_i8_real(
    const cudaStream_t stream,
    common::Handle_t &handle,
    const unsigned bcnt,
    const size_t ldc_hi,
    const size_t n,
    const size_t lda_lo,
    const size_t ldb_lo,
    const size_t sizeA,
    const size_t sizeB,
    const size_t sizeC,
    common::matptr_t<common::low_t<Backend::INT8>, false> &A_lo,
    common::matptr_t<common::low_t<Backend::INT8>, false> &B_lo,
    common::matptr_t<common::hi_t<Backend::INT8>, false> &C_hi //
) {
    using HiT = common::hi_t<Backend::INT8>;

    constexpr HiT one  = 1;
    constexpr HiT zero = 0;

    if (bcnt == 1) {

        matmul_block_1<KIND, Backend::INT8, UPLO_A, UPLO_B, UPLO_C>(
            stream, handle, ldc_hi, n, lda_lo, ldb_lo,
            &one, A_lo.ptr0, B_lo.ptr0, &zero, C_hi.ptr0);

    } else {

        matmul_block_1_strided_batched<KIND, Backend::INT8, UPLO_A, UPLO_B, UPLO_C>(
            stream, handle, ldc_hi, n, lda_lo, ldb_lo, bcnt,
            &one,
            A_lo.ptr0, static_cast<int64_t>(sizeA),
            B_lo.ptr0, static_cast<int64_t>(sizeB),
            &zero,
            C_hi.ptr0, static_cast<int64_t>(sizeC));
    }

    A_lo.shift(bcnt * sizeA);
    if constexpr (KIND != common::MatMulKind::ATxA && KIND != common::MatMulKind::AHxA) {
        B_lo.shift(bcnt * sizeB);
    }
}

template <common::MatMulKind KIND,
          cublasFillMode_t UPLO_A,
          cublasFillMode_t UPLO_B,
          cublasFillMode_t UPLO_C>
inline void error_free_matmult_i8_complex(
    const cudaStream_t stream,
    common::Handle_t &handle,
    const unsigned bcnt,
    const size_t ldc_hi,
    const size_t n,
    const size_t lda_lo,
    const size_t ldb_lo,
    const size_t sizeA,
    const size_t sizeB,
    const size_t sizeC,
    common::matptr_t<common::low_t<Backend::INT8>, true> &A_lo,
    common::matptr_t<common::low_t<Backend::INT8>, true> &B_lo,
    common::matptr_t<common::hi_t<Backend::INT8>, true> &C_hi //
) {
    using HiT         = common::hi_t<Backend::INT8>;
    constexpr HiT one = 1, zero = 0;

    constexpr auto MK = KIND;
    auto *A0          = KIND == common::MatMulKind::AHxA ? A_lo.ptr1 : A_lo.ptr0;
    auto *A1          = KIND == common::MatMulKind::AHxA ? A_lo.ptr0 : A_lo.ptr1;

    if (bcnt == 1) {
        matmul_block_1<MK, Backend::INT8, UPLO_A, UPLO_B, UPLO_C>(
            stream, handle, ldc_hi, n, lda_lo, ldb_lo,
            &one, A0, B_lo.ptr0, &zero, C_hi.ptr0);
        if constexpr (KIND != common::MatMulKind::AHxA) {
            matmul_block_1<MK, Backend::INT8, UPLO_A, UPLO_B, UPLO_C>(
                stream, handle, ldc_hi, n, lda_lo, ldb_lo,
                &one, A1, B_lo.ptr1, &zero, C_hi.ptr1);
        }
    } else {
        constexpr unsigned products = KIND == common::MatMulKind::AHxA ? 1u : 2u;
        matmul_block_1_strided_batched<MK, Backend::INT8, UPLO_A, UPLO_B, UPLO_C>(
            stream, handle, ldc_hi, n, lda_lo, ldb_lo, bcnt,
            &one, A0, int64_t(sizeA), B_lo.ptr0, int64_t(sizeB),
            &zero, C_hi.ptr0, int64_t(products * sizeC));
        if constexpr (KIND != common::MatMulKind::AHxA) {
            matmul_block_1_strided_batched<MK, Backend::INT8, UPLO_A, UPLO_B, UPLO_C>(
                stream, handle, ldc_hi, n, lda_lo, ldb_lo, bcnt,
                &one, A1, int64_t(sizeA), B_lo.ptr1, int64_t(sizeB),
                &zero, C_hi.ptr1, int64_t(2 * sizeC));
        }
    }
    A_lo.shift(bcnt * sizeA);
    if constexpr (KIND != common::MatMulKind::ATxA && KIND != common::MatMulKind::AHxA) {
        B_lo.shift(bcnt * sizeB);
    }
}

inline void error_free_matmult_f8_real(
    const cudaStream_t stream,
    common::Handle_t &handle,
    const unsigned idx,
    const unsigned bcnt,
    const size_t ldc_hi,
    const size_t n,
    const size_t lda_lo,
    const size_t ldb_lo,
    const size_t sizeA,
    const size_t sizeB,
    const size_t sizeC,
    common::matptr_t<common::low_t<Backend::FP8>, false> &A_lo,
    common::matptr_t<common::low_t<Backend::FP8>, false> &B_lo,
    common::matptr_t<common::hi_t<Backend::FP8>, false> &C_hi //
) {
    using LowT          = common::low_t<Backend::FP8>;
    constexpr float one = 1.0f, zero = 0.0f;
    constexpr unsigned parts = 1U;
    std::array<void *, 60U * parts> hA{}, hB{}, hC{};
    unsigned pcnt  = 0;
    size_t offsetA = 0, offsetB = 0;
    for (unsigned b = 0; b < bcnt; ++b) {
        const auto scheme = common::fp8_plan::get(common::fp8_plan::modulus(handle.fp8_num_moduli, idx + b, false));
        for (unsigned part = 0; part < parts; ++part) {
            LowT *a    = A_lo.ptr0 + offsetA;
            LowT *bptr = B_lo.ptr0 + offsetB;
            for (unsigned j = 0; j < scheme.products; ++j) {
                hA[pcnt] = a + j * sizeA;
                hB[pcnt] = bptr + j * sizeB;
                hC[pcnt] = C_hi.ptr0 + size_t(pcnt) * sizeC;
                ++pcnt;
            }
        }
        offsetA += scheme.products * sizeA;
        offsetB += scheme.products * sizeB;
    }
    upload_pointer_arrays(stream, handle, hA.data(), hB.data(), hC.data(), pcnt);
    common::call_gemm_tn_pointer_batched<Backend::FP8>(
        stream, handle, ldc_hi, n, lda_lo, pcnt, &one,
        handle.Aarray, lda_lo, handle.Barray, ldb_lo, &zero, handle.Carray, ldc_hi);
    A_lo.shift(offsetA);
    B_lo.shift(offsetB);
}

template <common::MatMulKind KIND,
          cublasFillMode_t UPLO_A,
          cublasFillMode_t UPLO_B,
          cublasFillMode_t UPLO_C>
inline void error_free_matmult_f8_real_strided(
    const cudaStream_t stream,
    common::Handle_t &handle,
    const unsigned idx,
    const unsigned bcnt,
    const size_t ldc_hi,
    const size_t n,
    const size_t lda_lo,
    const size_t ldb_lo,
    const size_t sizeA,
    const size_t sizeB,
    const size_t sizeC,
    common::matptr_t<common::low_t<Backend::FP8>, false> &A_lo,
    common::matptr_t<common::low_t<Backend::FP8>, false> &B_lo,
    common::matptr_t<common::hi_t<Backend::FP8>, false> &C_hi //
) {
    using LowT          = common::low_t<Backend::FP8>;
    constexpr float one = 1.0f, zero = 0.0f;
    constexpr unsigned parts = 1U;
    const unsigned products  = common::fp8_plan::get(common::fp8_plan::modulus(handle.fp8_num_moduli, idx, false)).products;
    const int64_t strideA    = products * int64_t(sizeA);
    const int64_t strideB    = products * int64_t(sizeB);
    const int64_t strideC    = products * parts * int64_t(sizeC);
    for (unsigned part = 0; part < parts; ++part) {
        LowT *a = A_lo.ptr0;
        LowT *b = B_lo.ptr0;

        for (unsigned j = 0; j < products; ++j) {
            matmul_block_1_strided_batched<KIND, Backend::FP8, UPLO_A, UPLO_B, UPLO_C>(
                stream, handle, ldc_hi, n, lda_lo, ldb_lo, bcnt, &one,
                a + j * sizeA, strideA, b + j * sizeB, strideB, &zero,
                C_hi.ptr0 + (products * part + j) * sizeC, strideC);
        }
    }
    A_lo.shift(bcnt * products * sizeA);
    if constexpr (KIND != common::MatMulKind::ATxA && KIND != common::MatMulKind::AHxA) {
        B_lo.shift(bcnt * products * sizeB);
    }
}

inline void error_free_matmult_f8_complex(
    const cudaStream_t stream,
    common::Handle_t &handle,
    const unsigned idx,
    const unsigned bcnt,
    const size_t ldc_hi,
    const size_t n,
    const size_t lda_lo,
    const size_t ldb_lo,
    const size_t sizeA,
    const size_t sizeB,
    const size_t sizeC,
    common::matptr_t<common::low_t<Backend::FP8>, true> &A_lo,
    common::matptr_t<common::low_t<Backend::FP8>, true> &B_lo,
    common::matptr_t<common::hi_t<Backend::FP8>, true> &C_hi //
) {
    using LowT          = common::low_t<Backend::FP8>;
    constexpr float one = 1.0f, zero = 0.0f;
    constexpr unsigned parts = 2U;
    std::array<void *, 60U * parts> hA{}, hB{}, hC{};
    unsigned pcnt  = 0;
    size_t offsetA = 0, offsetB = 0;
    for (unsigned b = 0; b < bcnt; ++b) {
        const auto scheme = common::fp8_plan::get(common::fp8_plan::modulus(handle.fp8_num_moduli, idx + b, true));
        for (unsigned part = 0; part < parts; ++part) {
            LowT *a    = (part == 0 ? A_lo.ptr0 : A_lo.ptr1) + offsetA;
            LowT *bptr = (part == 0 ? B_lo.ptr0 : B_lo.ptr1) + offsetB;
            for (unsigned j = 0; j < scheme.products; ++j) {
                hA[pcnt] = a + j * sizeA;
                hB[pcnt] = bptr + j * sizeB;
                hC[pcnt] = C_hi.ptr0 + size_t(pcnt) * sizeC;
                ++pcnt;
            }
        }
        offsetA += scheme.products * sizeA;
        offsetB += scheme.products * sizeB;
    }
    upload_pointer_arrays(stream, handle, hA.data(), hB.data(), hC.data(), pcnt);
    common::call_gemm_tn_pointer_batched<Backend::FP8>(
        stream, handle, ldc_hi, n, lda_lo, pcnt, &one,
        handle.Aarray, lda_lo, handle.Barray, ldb_lo, &zero, handle.Carray, ldc_hi);
    A_lo.shift(offsetA);
    B_lo.shift(offsetB);
}

template <common::MatMulKind KIND,
          cublasFillMode_t UPLO_A,
          cublasFillMode_t UPLO_B,
          cublasFillMode_t UPLO_C>
inline void error_free_matmult_f8_complex_strided(
    const cudaStream_t stream,
    common::Handle_t &handle,
    const unsigned idx,
    const unsigned bcnt,
    const size_t ldc_hi,
    const size_t n,
    const size_t lda_lo,
    const size_t ldb_lo,
    const size_t sizeA,
    const size_t sizeB,
    const size_t sizeC,
    common::matptr_t<common::low_t<Backend::FP8>, true> &A_lo,
    common::matptr_t<common::low_t<Backend::FP8>, true> &B_lo,
    common::matptr_t<common::hi_t<Backend::FP8>, true> &C_hi //
) {
    using LowT          = common::low_t<Backend::FP8>;
    constexpr float one = 1.0f, zero = 0.0f;
    constexpr unsigned parts = KIND == common::MatMulKind::AHxA ? 1U : 2U;
    const unsigned products  = common::fp8_plan::get(common::fp8_plan::modulus(handle.fp8_num_moduli, idx, true)).products;
    const int64_t strideA    = products * int64_t(sizeA);
    const int64_t strideB    = products * int64_t(sizeB);
    const int64_t strideC    = products * parts * int64_t(sizeC);
    for (unsigned part = 0; part < parts; ++part) {
        LowT *a = part == 0 ? A_lo.ptr0 : A_lo.ptr1;
        LowT *b = part == 0 ? B_lo.ptr0 : B_lo.ptr1;
        if constexpr (KIND == common::MatMulKind::AHxA) a = A_lo.ptr1;
        for (unsigned j = 0; j < products; ++j) {
            matmul_block_1_strided_batched<KIND, Backend::FP8, UPLO_A, UPLO_B, UPLO_C>(
                stream, handle, ldc_hi, n, lda_lo, ldb_lo, bcnt, &one,
                a + j * sizeA, strideA, b + j * sizeB, strideB, &zero,
                C_hi.ptr0 + (products * part + j) * sizeC, strideC);
        }
    }
    A_lo.shift(bcnt * products * sizeA);
    if constexpr (KIND != common::MatMulKind::ATxA && KIND != common::MatMulKind::AHxA) {
        B_lo.shift(bcnt * products * sizeB);
    }
}

template <common::MatMulKind KIND,
          cublasFillMode_t UPLO_A,
          cublasFillMode_t UPLO_B,
          cublasFillMode_t UPLO_C>
inline void error_free_matmult_f8_real_strided_split(
    const cudaStream_t stream,
    common::Handle_t &handle,
    unsigned idx,
    unsigned bcnt,
    const size_t ldc_hi,
    const size_t n,
    const size_t lda_lo,
    const size_t ldb_lo,
    const size_t sizeA,
    const size_t sizeB,
    const size_t sizeC,
    common::matptr_t<common::low_t<Backend::FP8>, false> &A_lo,
    common::matptr_t<common::low_t<Backend::FP8>, false> &B_lo,
    common::matptr_t<common::hi_t<Backend::FP8>, false> &C_hi //
) {
    unsigned done  = 0;
    size_t offsetC = 0;

    while (done < bcnt) {
        const unsigned cur      = idx + done;
        const unsigned products = common::fp8_plan::get(common::fp8_plan::modulus(handle.fp8_num_moduli, cur, false)).products;

        unsigned cnt = 1;
        while (done + cnt < bcnt && common::fp8_plan::get(common::fp8_plan::modulus(handle.fp8_num_moduli, cur + cnt, false)).products == products) {
            ++cnt;
        }

        auto C_part = C_hi;
        C_part.ptr0 += offsetC;

        error_free_matmult_f8_real_strided<KIND, UPLO_A, UPLO_B, UPLO_C>(
            stream, handle,
            cur, cnt,
            ldc_hi, n, lda_lo, ldb_lo,
            sizeA, sizeB, sizeC,
            A_lo, B_lo, C_part);

        offsetC += size_t(cnt) * products * sizeC;
        done += cnt;
    }
}

template <common::MatMulKind KIND,
          cublasFillMode_t UPLO_A,
          cublasFillMode_t UPLO_B,
          cublasFillMode_t UPLO_C>
inline void error_free_matmult_f8_complex_strided_split(
    const cudaStream_t stream,
    common::Handle_t &handle,
    unsigned idx,
    unsigned bcnt,
    const size_t ldc_hi,
    const size_t n,
    const size_t lda_lo,
    const size_t ldb_lo,
    const size_t sizeA,
    const size_t sizeB,
    const size_t sizeC,
    common::matptr_t<common::low_t<Backend::FP8>, true> &A_lo,
    common::matptr_t<common::low_t<Backend::FP8>, true> &B_lo,
    common::matptr_t<common::hi_t<Backend::FP8>, true> &C_hi //
) {
    unsigned done  = 0;
    size_t offsetC = 0;

    while (done < bcnt) {
        const unsigned cur      = idx + done;
        const unsigned products = common::fp8_plan::get(common::fp8_plan::modulus(handle.fp8_num_moduli, cur, true)).products;

        unsigned cnt = 1;
        while (done + cnt < bcnt && common::fp8_plan::get(common::fp8_plan::modulus(handle.fp8_num_moduli, cur + cnt, true)).products == products) {
            ++cnt;
        }

        auto C_part              = C_hi;
        constexpr unsigned parts = KIND == common::MatMulKind::AHxA ? 1U : 2U;
        C_part.ptr0 += offsetC;
        C_part.ptr1 = parts == 1U ? nullptr : C_part.ptr0 + products * sizeC;

        error_free_matmult_f8_complex_strided<KIND, UPLO_A, UPLO_B, UPLO_C>(
            stream, handle,
            cur, cnt,
            ldc_hi, n, lda_lo, ldb_lo,
            sizeA, sizeB, sizeC,
            A_lo, B_lo, C_part);

        offsetC += size_t(cnt) * parts * products * sizeC;
        done += cnt;
    }
}

template <common::MatMulKind KIND,
          cublasFillMode_t UPLO_A,
          cublasFillMode_t UPLO_B,
          cublasFillMode_t UPLO_C>
inline void error_free_matmult_i8_real_launch(
    const cudaStream_t stream,
    common::Handle_t &handle,
    const unsigned bcnt,
    const size_t ldc_hi,
    const size_t n,
    const size_t lda_lo,
    const size_t ldb_lo,
    const size_t sizeA,
    const size_t sizeB,
    const size_t sizeC,
    cublasOperation_t op_A, cublasOperation_t op_B,
    common::matptr_t<common::low_t<Backend::INT8>, false> &A_lo,
    common::matptr_t<common::low_t<Backend::INT8>, false> &B_lo,
    common::matptr_t<common::hi_t<Backend::INT8>, false> &C_hi //
) {
    if constexpr (KIND == common::MatMulKind::TrmmLeft) {
        if (op_A == CUBLAS_OP_N) {
            error_free_matmult_i8_real<KIND, UPLO_A, UPLO_B, UPLO_C>(
                stream, handle, bcnt, ldc_hi, n, lda_lo, ldb_lo,
                sizeA, sizeB, sizeC, A_lo, B_lo, C_hi);
        } else {
            error_free_matmult_i8_real<KIND, flip_uplo<UPLO_A>, UPLO_B, UPLO_C>(
                stream, handle, bcnt, ldc_hi, n, lda_lo, ldb_lo,
                sizeA, sizeB, sizeC, A_lo, B_lo, C_hi);
        }
    } else if constexpr (KIND == common::MatMulKind::TrmmRight) {
        if (op_B == CUBLAS_OP_N) {
            error_free_matmult_i8_real<KIND, UPLO_A, UPLO_B, UPLO_C>(
                stream, handle, bcnt, ldc_hi, n, lda_lo, ldb_lo,
                sizeA, sizeB, sizeC, A_lo, B_lo, C_hi);
        } else {
            error_free_matmult_i8_real<KIND, UPLO_A, flip_uplo<UPLO_B>, UPLO_C>(
                stream, handle, bcnt, ldc_hi, n, lda_lo, ldb_lo,
                sizeA, sizeB, sizeC, A_lo, B_lo, C_hi);
        }
    } else if constexpr (KIND == common::MatMulKind::Trtrmm) {
        if (op_A == CUBLAS_OP_N) {
            if (op_B == CUBLAS_OP_N) {
                error_free_matmult_i8_real<KIND, UPLO_A, UPLO_B, UPLO_C>(
                    stream, handle, bcnt, ldc_hi, n, lda_lo, ldb_lo,
                    sizeA, sizeB, sizeC, A_lo, B_lo, C_hi);
            } else {
                error_free_matmult_i8_real<KIND, UPLO_A, flip_uplo<UPLO_B>, UPLO_C>(
                    stream, handle, bcnt, ldc_hi, n, lda_lo, ldb_lo,
                    sizeA, sizeB, sizeC, A_lo, B_lo, C_hi);
            }
        } else {
            if (op_B == CUBLAS_OP_N) {
                error_free_matmult_i8_real<KIND, flip_uplo<UPLO_A>, UPLO_B, UPLO_C>(
                    stream, handle, bcnt, ldc_hi, n, lda_lo, ldb_lo,
                    sizeA, sizeB, sizeC, A_lo, B_lo, C_hi);
            } else {
                error_free_matmult_i8_real<KIND, flip_uplo<UPLO_A>, flip_uplo<UPLO_B>, UPLO_C>(
                    stream, handle, bcnt, ldc_hi, n, lda_lo, ldb_lo,
                    sizeA, sizeB, sizeC, A_lo, B_lo, C_hi);
            }
        }
    } else {
        error_free_matmult_i8_real<KIND, UPLO_A, UPLO_B, UPLO_C>(
            stream, handle, bcnt, ldc_hi, n, lda_lo, ldb_lo,
            sizeA, sizeB, sizeC, A_lo, B_lo, C_hi);
    }
}

template <common::MatMulKind KIND,
          cublasFillMode_t UPLO_A,
          cublasFillMode_t UPLO_B,
          cublasFillMode_t UPLO_C>
inline void error_free_matmult_i8_complex_launch(
    const cudaStream_t stream,
    common::Handle_t &handle,
    const unsigned bcnt,
    const size_t ldc_hi,
    const size_t n,
    const size_t lda_lo,
    const size_t ldb_lo,
    const size_t sizeA,
    const size_t sizeB,
    const size_t sizeC,
    cublasOperation_t op_A, cublasOperation_t op_B,
    common::matptr_t<common::low_t<Backend::INT8>, true> &A_lo,
    common::matptr_t<common::low_t<Backend::INT8>, true> &B_lo,
    common::matptr_t<common::hi_t<Backend::INT8>, true> &C_hi //
) {
    if constexpr (KIND == common::MatMulKind::TrmmLeft) {
        if (op_A == CUBLAS_OP_N) {
            error_free_matmult_i8_complex<KIND, UPLO_A, UPLO_B, UPLO_C>(
                stream, handle, bcnt, ldc_hi, n, lda_lo, ldb_lo,
                sizeA, sizeB, sizeC, A_lo, B_lo, C_hi);
        } else {
            error_free_matmult_i8_complex<KIND, flip_uplo<UPLO_A>, UPLO_B, UPLO_C>(
                stream, handle, bcnt, ldc_hi, n, lda_lo, ldb_lo,
                sizeA, sizeB, sizeC, A_lo, B_lo, C_hi);
        }
    } else if constexpr (KIND == common::MatMulKind::TrmmRight) {
        if (op_B == CUBLAS_OP_N) {
            error_free_matmult_i8_complex<KIND, UPLO_A, UPLO_B, UPLO_C>(
                stream, handle, bcnt, ldc_hi, n, lda_lo, ldb_lo,
                sizeA, sizeB, sizeC, A_lo, B_lo, C_hi);
        } else {
            error_free_matmult_i8_complex<KIND, UPLO_A, flip_uplo<UPLO_B>, UPLO_C>(
                stream, handle, bcnt, ldc_hi, n, lda_lo, ldb_lo,
                sizeA, sizeB, sizeC, A_lo, B_lo, C_hi);
        }
    } else if constexpr (KIND == common::MatMulKind::Trtrmm) {
        if (op_A == CUBLAS_OP_N) {
            if (op_B == CUBLAS_OP_N) {
                error_free_matmult_i8_complex<KIND, UPLO_A, UPLO_B, UPLO_C>(
                    stream, handle, bcnt, ldc_hi, n, lda_lo, ldb_lo,
                    sizeA, sizeB, sizeC, A_lo, B_lo, C_hi);
            } else {
                error_free_matmult_i8_complex<KIND, UPLO_A, flip_uplo<UPLO_B>, UPLO_C>(
                    stream, handle, bcnt, ldc_hi, n, lda_lo, ldb_lo,
                    sizeA, sizeB, sizeC, A_lo, B_lo, C_hi);
            }
        } else {
            if (op_B == CUBLAS_OP_N) {
                error_free_matmult_i8_complex<KIND, flip_uplo<UPLO_A>, UPLO_B, UPLO_C>(
                    stream, handle, bcnt, ldc_hi, n, lda_lo, ldb_lo,
                    sizeA, sizeB, sizeC, A_lo, B_lo, C_hi);
            } else {
                error_free_matmult_i8_complex<KIND, flip_uplo<UPLO_A>, flip_uplo<UPLO_B>, UPLO_C>(
                    stream, handle, bcnt, ldc_hi, n, lda_lo, ldb_lo,
                    sizeA, sizeB, sizeC, A_lo, B_lo, C_hi);
            }
        }
    } else {
        error_free_matmult_i8_complex<KIND, UPLO_A, UPLO_B, UPLO_C>(
            stream, handle, bcnt, ldc_hi, n, lda_lo, ldb_lo,
            sizeA, sizeB, sizeC, A_lo, B_lo, C_hi);
    }
}

template <common::MatMulKind KIND,
          cublasFillMode_t UPLO_A,
          cublasFillMode_t UPLO_B,
          cublasFillMode_t UPLO_C>
inline void error_free_matmult_f8_real_strided_split_launch(
    const cudaStream_t stream,
    common::Handle_t &handle,
    unsigned idx,
    unsigned bcnt,
    const size_t ldc_hi,
    const size_t n,
    const size_t lda_lo,
    const size_t ldb_lo,
    const size_t sizeA,
    const size_t sizeB,
    const size_t sizeC,
    cublasOperation_t op_A, cublasOperation_t op_B,
    common::matptr_t<common::low_t<Backend::FP8>, false> &A_lo,
    common::matptr_t<common::low_t<Backend::FP8>, false> &B_lo,
    common::matptr_t<common::hi_t<Backend::FP8>, false> &C_hi //
) {
    if constexpr (KIND == common::MatMulKind::TrmmLeft) {
        if (op_A == CUBLAS_OP_N) {
            error_free_matmult_f8_real_strided_split<KIND, UPLO_A, UPLO_B, UPLO_C>(
                stream, handle, idx, bcnt, ldc_hi, n, lda_lo, ldb_lo,
                sizeA, sizeB, sizeC, A_lo, B_lo, C_hi);
        } else {
            error_free_matmult_f8_real_strided_split<KIND, flip_uplo<UPLO_A>, UPLO_B, UPLO_C>(
                stream, handle, idx, bcnt, ldc_hi, n, lda_lo, ldb_lo,
                sizeA, sizeB, sizeC, A_lo, B_lo, C_hi);
        }
    } else if constexpr (KIND == common::MatMulKind::TrmmRight) {
        if (op_B == CUBLAS_OP_N) {
            error_free_matmult_f8_real_strided_split<KIND, UPLO_A, UPLO_B, UPLO_C>(
                stream, handle, idx, bcnt, ldc_hi, n, lda_lo, ldb_lo,
                sizeA, sizeB, sizeC, A_lo, B_lo, C_hi);
        } else {
            error_free_matmult_f8_real_strided_split<KIND, UPLO_A, flip_uplo<UPLO_B>, UPLO_C>(
                stream, handle, idx, bcnt, ldc_hi, n, lda_lo, ldb_lo,
                sizeA, sizeB, sizeC, A_lo, B_lo, C_hi);
        }
    } else if constexpr (KIND == common::MatMulKind::Trtrmm) {
        if (op_A == CUBLAS_OP_N) {
            if (op_B == CUBLAS_OP_N) {
                error_free_matmult_f8_real_strided_split<KIND, UPLO_A, UPLO_B, UPLO_C>(
                    stream, handle, idx, bcnt, ldc_hi, n, lda_lo, ldb_lo,
                    sizeA, sizeB, sizeC, A_lo, B_lo, C_hi);
            } else {
                error_free_matmult_f8_real_strided_split<KIND, UPLO_A, flip_uplo<UPLO_B>, UPLO_C>(
                    stream, handle, idx, bcnt, ldc_hi, n, lda_lo, ldb_lo,
                    sizeA, sizeB, sizeC, A_lo, B_lo, C_hi);
            }
        } else {
            if (op_B == CUBLAS_OP_N) {
                error_free_matmult_f8_real_strided_split<KIND, flip_uplo<UPLO_A>, UPLO_B, UPLO_C>(
                    stream, handle, idx, bcnt, ldc_hi, n, lda_lo, ldb_lo,
                    sizeA, sizeB, sizeC, A_lo, B_lo, C_hi);
            } else {
                error_free_matmult_f8_real_strided_split<KIND, flip_uplo<UPLO_A>, flip_uplo<UPLO_B>, UPLO_C>(
                    stream, handle, idx, bcnt, ldc_hi, n, lda_lo, ldb_lo,
                    sizeA, sizeB, sizeC, A_lo, B_lo, C_hi);
            }
        }
    } else {
        error_free_matmult_f8_real_strided_split<KIND, UPLO_A, UPLO_B, UPLO_C>(
            stream, handle, idx, bcnt, ldc_hi, n, lda_lo, ldb_lo,
            sizeA, sizeB, sizeC, A_lo, B_lo, C_hi);
    }
}

template <common::MatMulKind KIND,
          cublasFillMode_t UPLO_A,
          cublasFillMode_t UPLO_B,
          cublasFillMode_t UPLO_C>
inline void error_free_matmult_f8_complex_strided_split_launch(
    const cudaStream_t stream,
    common::Handle_t &handle,
    unsigned idx,
    unsigned bcnt,
    const size_t ldc_hi,
    const size_t n,
    const size_t lda_lo,
    const size_t ldb_lo,
    const size_t sizeA,
    const size_t sizeB,
    const size_t sizeC,
    cublasOperation_t op_A, cublasOperation_t op_B,
    common::matptr_t<common::low_t<Backend::FP8>, true> &A_lo,
    common::matptr_t<common::low_t<Backend::FP8>, true> &B_lo,
    common::matptr_t<common::hi_t<Backend::FP8>, true> &C_hi //
) {
    if constexpr (KIND == common::MatMulKind::TrmmLeft) {
        if (op_A == CUBLAS_OP_N) {
            error_free_matmult_f8_complex_strided_split<KIND, UPLO_A, UPLO_B, UPLO_C>(
                stream, handle, idx, bcnt, ldc_hi, n, lda_lo, ldb_lo,
                sizeA, sizeB, sizeC, A_lo, B_lo, C_hi);
        } else {
            error_free_matmult_f8_complex_strided_split<KIND, flip_uplo<UPLO_A>, UPLO_B, UPLO_C>(
                stream, handle, idx, bcnt, ldc_hi, n, lda_lo, ldb_lo,
                sizeA, sizeB, sizeC, A_lo, B_lo, C_hi);
        }
    } else if constexpr (KIND == common::MatMulKind::TrmmRight) {
        if (op_B == CUBLAS_OP_N) {
            error_free_matmult_f8_complex_strided_split<KIND, UPLO_A, UPLO_B, UPLO_C>(
                stream, handle, idx, bcnt, ldc_hi, n, lda_lo, ldb_lo,
                sizeA, sizeB, sizeC, A_lo, B_lo, C_hi);
        } else {
            error_free_matmult_f8_complex_strided_split<KIND, UPLO_A, flip_uplo<UPLO_B>, UPLO_C>(
                stream, handle, idx, bcnt, ldc_hi, n, lda_lo, ldb_lo,
                sizeA, sizeB, sizeC, A_lo, B_lo, C_hi);
        }
    } else if constexpr (KIND == common::MatMulKind::Trtrmm) {
        if (op_A == CUBLAS_OP_N) {
            if (op_B == CUBLAS_OP_N) {
                error_free_matmult_f8_complex_strided_split<KIND, UPLO_A, UPLO_B, UPLO_C>(
                    stream, handle, idx, bcnt, ldc_hi, n, lda_lo, ldb_lo,
                    sizeA, sizeB, sizeC, A_lo, B_lo, C_hi);
            } else {
                error_free_matmult_f8_complex_strided_split<KIND, UPLO_A, flip_uplo<UPLO_B>, UPLO_C>(
                    stream, handle, idx, bcnt, ldc_hi, n, lda_lo, ldb_lo,
                    sizeA, sizeB, sizeC, A_lo, B_lo, C_hi);
            }
        } else {
            if (op_B == CUBLAS_OP_N) {
                error_free_matmult_f8_complex_strided_split<KIND, flip_uplo<UPLO_A>, UPLO_B, UPLO_C>(
                    stream, handle, idx, bcnt, ldc_hi, n, lda_lo, ldb_lo,
                    sizeA, sizeB, sizeC, A_lo, B_lo, C_hi);
            } else {
                error_free_matmult_f8_complex_strided_split<KIND, flip_uplo<UPLO_A>, flip_uplo<UPLO_B>, UPLO_C>(
                    stream, handle, idx, bcnt, ldc_hi, n, lda_lo, ldb_lo,
                    sizeA, sizeB, sizeC, A_lo, B_lo, C_hi);
            }
        }
    } else {
        error_free_matmult_f8_complex_strided_split<KIND, UPLO_A, UPLO_B, UPLO_C>(
            stream, handle, idx, bcnt, ldc_hi, n, lda_lo, ldb_lo,
            sizeA, sizeB, sizeC, A_lo, B_lo, C_hi);
    }
}

template <Func FUNC,
          Backend BACKEND,
          bool COMPLEX,
          common::MatStruct STRUCT_A = common::MatStruct::Full,
          common::MatStruct STRUCT_B = common::MatStruct::Full,
          cublasFillMode_t UPLO_A    = CUBLAS_FILL_MODE_FULL,
          cublasFillMode_t UPLO_B    = CUBLAS_FILL_MODE_FULL,
          cublasFillMode_t UPLO_C    = CUBLAS_FILL_MODE_FULL>
inline void error_free_matmult(
    const cudaStream_t stream,
    common::Handle_t &handle,
    const unsigned idx,
    const unsigned bcnt,
    const size_t ldc_hi,
    const size_t n,
    const size_t lda_lo,
    const size_t ldb_lo,
    const size_t sizeA,
    const size_t sizeB,
    const size_t sizeC,
    cublasOperation_t op_A, cublasOperation_t op_B,
    common::matptr_t<common::low_t<BACKEND>, COMPLEX> &A_lo,
    common::matptr_t<common::low_t<BACKEND>, COMPLEX> &B_lo,
    common::matptr_t<common::hi_t<BACKEND>, COMPLEX> &C_hi //
) {
    constexpr common::MatMulKind KIND = matmul_kind<FUNC, STRUCT_A, STRUCT_B, UPLO_C>();

    if constexpr (BACKEND == Backend::INT8 && !COMPLEX) {

        error_free_matmult_i8_real_launch<KIND, UPLO_A, UPLO_B, UPLO_C>(
            stream, handle, bcnt, ldc_hi, n, lda_lo, ldb_lo,
            sizeA, sizeB, sizeC, op_A, op_B, A_lo, B_lo, C_hi);

    } else if constexpr (BACKEND == Backend::INT8 && COMPLEX) {

        error_free_matmult_i8_complex_launch<KIND, UPLO_A, UPLO_B, UPLO_C>(
            stream, handle, bcnt, ldc_hi, n, lda_lo, ldb_lo,
            sizeA, sizeB, sizeC, op_A, op_B, A_lo, B_lo, C_hi);

    } else if constexpr (BACKEND == Backend::FP8 && !COMPLEX) {

        if constexpr (common::isCUDA && (KIND == common::MatMulKind::Gemm)) {
            error_free_matmult_f8_real(
                stream, handle, idx, bcnt, ldc_hi, n, lda_lo, ldb_lo,
                sizeA, sizeB, sizeC, A_lo, B_lo, C_hi);
        } else {
            error_free_matmult_f8_real_strided_split_launch<KIND, UPLO_A, UPLO_B, UPLO_C>(
                stream, handle, idx, bcnt, ldc_hi, n, lda_lo, ldb_lo,
                sizeA, sizeB, sizeC, op_A, op_B, A_lo, B_lo, C_hi);
        }

    } else {

        if constexpr (common::isCUDA && (KIND == common::MatMulKind::Gemm)) {
            error_free_matmult_f8_complex(
                stream, handle, idx, bcnt, ldc_hi, n, lda_lo, ldb_lo,
                sizeA, sizeB, sizeC, A_lo, B_lo, C_hi);
        } else {
            error_free_matmult_f8_complex_strided_split_launch<KIND, UPLO_A, UPLO_B, UPLO_C>(
                stream, handle, idx, bcnt, ldc_hi, n, lda_lo, ldb_lo,
                sizeA, sizeB, sizeC, op_A, op_B, A_lo, B_lo, C_hi);
        }
    }
}

template <Func FUNC,
          Backend BACKEND,
          bool COMPLEX,
          cublasFillMode_t UPLO_C>
inline void error_free_matmult_rk(
    const cudaStream_t stream,
    common::Handle_t &handle,
    const unsigned idx,
    const unsigned bcnt,
    const size_t ldc_hi,
    const size_t n,
    const size_t lda_lo,
    const size_t sizeA,
    const size_t sizeC,
    common::matptr_t<common::low_t<BACKEND>, COMPLEX> &A_lo,
    common::matptr_t<common::hi_t<BACKEND>, COMPLEX> &C_hi //
) {
    static_assert(FUNC == Func::syrk || FUNC == Func::herk,
                  "error_free_matmult_rk supports only syrk/herk.");
    if constexpr (FUNC == Func::herk) {
        static_assert(COMPLEX, "herk requires complex input type.");
    }

    constexpr common::MatMulKind KIND = (FUNC == Func::syrk) ? common::MatMulKind::ATxA : common::MatMulKind::AHxA;

    if constexpr (BACKEND == Backend::INT8 && !COMPLEX) {

        error_free_matmult_i8_real<KIND, CUBLAS_FILL_MODE_FULL, CUBLAS_FILL_MODE_FULL, UPLO_C>(
            stream, handle, bcnt, ldc_hi, n, lda_lo, lda_lo,
            sizeA, sizeA, sizeC, A_lo, A_lo, C_hi);

    } else if constexpr (BACKEND == Backend::INT8 && COMPLEX) {

        error_free_matmult_i8_complex<KIND, CUBLAS_FILL_MODE_FULL, CUBLAS_FILL_MODE_FULL, UPLO_C>(
            stream, handle, bcnt, ldc_hi, n, lda_lo, lda_lo,
            sizeA, sizeA, sizeC, A_lo, A_lo, C_hi);

    } else if constexpr (BACKEND == Backend::FP8 && !COMPLEX) {

        error_free_matmult_f8_real_strided_split<KIND, CUBLAS_FILL_MODE_FULL, CUBLAS_FILL_MODE_FULL, UPLO_C>(
            stream, handle, idx, bcnt, ldc_hi, n, lda_lo, lda_lo,
            sizeA, sizeA, sizeC, A_lo, A_lo, C_hi);

    } else {

        error_free_matmult_f8_complex_strided_split<KIND, CUBLAS_FILL_MODE_FULL, CUBLAS_FILL_MODE_FULL, UPLO_C>(
            stream, handle, idx, bcnt, ldc_hi, n, lda_lo, lda_lo,
            sizeA, sizeA, sizeC, A_lo, A_lo, C_hi);
    }
}

} // namespace gemmul8::oz2::core
