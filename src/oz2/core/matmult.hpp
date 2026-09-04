#pragma once
#include "../common/common.hpp"
#include "../common/matmult.hpp"
#include "../common/table.hpp"
#include "helper_triangular.hpp"
#include "helper_matmult.hpp"

namespace gemmul8::oz2::core {

template <unsigned NUM_MODULI>
inline unsigned batch_count(
    const int arch,
    const size_t n,
    const unsigned i,
    const size_t sizeC_Mid,
    const size_t sizeC_Hi,
    const size_t lwork_blas,
    const size_t worksizeC,
    const unsigned pointer_products_per_modulus = 0u //
) {
    if (arch == 121 && sizeC_Hi < (size_t(64) << 20)) return 1u;
    if (arch == 90 && n > 2048) return 1u;
    const size_t available  = worksizeC - i * sizeC_Mid;
    const unsigned rem_calc = NUM_MODULI - i;
    for (unsigned cnt = rem_calc; cnt >= 1; --cnt) {
        const size_t mid_bytes             = cnt * sizeC_Mid;
        const size_t hi_bytes              = cnt * sizeC_Hi;
        const unsigned pointer_batch_count = pointer_products_per_modulus * cnt;
        const size_t ptr_bytes             = size_t(3 * pointer_batch_count) * sizeof(void *);
        const size_t ptr_gap               = common::padding(ptr_bytes);
        const size_t gap                   = std::max<size_t>(ptr_gap + lwork_blas, mid_bytes);
        const size_t req                   = gap + hi_bytes;
        if (req <= available) return cnt;
    }
    return 1u;
}

inline constexpr size_t K_BLOCK_INT8 = size_t(1 << 17);

template <Backend BACKEND>
inline void configure_matprod_k_blocking(
    common::Handle_t &handle,
    const unsigned idx,
    const size_t k //
) {
    handle.modulus_idx = idx;
    if constexpr (BACKEND == Backend::INT8) {
        constexpr int KB             = int(K_BLOCK_INT8);
        const bool is_mod256         = (common::table::moduli_int8[idx] == 256);
        handle.matprod_k_blocking    = !is_mod256 && k > size_t(KB);
        handle.matprod_k_block_first = KB;
        handle.matprod_k_block_next  = KB;
    } else {
        const int KB0                = int(common::table::k_block_first_fp8[idx]);
        const int KB                 = int(common::table::k_block_next_fp8[idx]);
        handle.matprod_k_blocking    = k > size_t(KB0);
        handle.matprod_k_block_first = KB0;
        handle.matprod_k_block_next  = KB;
    }
}

template <Backend BACKEND>
inline bool needs_k_blocking(const unsigned idx, const size_t k) {
    if constexpr (BACKEND == Backend::INT8) {
        if (common::table::moduli_int8[idx] == 256) { return false; }
        return k > K_BLOCK_INT8;
    } else {
        return k > common::table::k_block_first_fp8[idx];
    }
}

template <Backend BACKEND>
inline unsigned limit_batch_for_k(const unsigned i, const unsigned bcnt, const size_t k) {
    unsigned safe = 0;

    while (safe < bcnt) {
        const unsigned idx = i + safe;
        if (needs_k_blocking<BACKEND>(idx, k)) { break; }
        ++safe;
    }

    return (safe == 0) ? 1 : safe;
}

inline unsigned fp8_planes_consumed(const unsigned i0, const unsigned bcnt) {
    return common::table::num_mat_fp8[i0 + bcnt] - common::table::num_mat_fp8[i0];
}

inline size_t fp8_plane_offset_from_group_start(
    const unsigned i0,
    const unsigned b,
    const size_t sizeX //
) {
    return size_t(common::table::num_mat_fp8[i0 + b] - common::table::num_mat_fp8[i0]) * sizeX;
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
    using HiT = common::hi_t<Backend::INT8>;

    constexpr HiT one  = 1;
    constexpr HiT zero = 0;

    if (bcnt == 1) {
        if constexpr (KIND == common::MatMulKind::AHxA) {
            matmul_block_3<KIND, Backend::INT8, UPLO_A, UPLO_B, UPLO_C>(
                stream, handle, ldc_hi, n, lda_lo, ldb_lo,
                &one, &one, &one,
                A_lo.ptr0, A_lo.ptr1, A_lo.ptr2,
                A_lo.ptr1, A_lo.ptr0, A_lo.ptr2,
                &zero, &zero, &zero,
                C_hi.ptr0, C_hi.ptr1, C_hi.ptr2);
        } else {
            matmul_block_3<KIND, Backend::INT8, UPLO_A, UPLO_B, UPLO_C>(
                stream, handle, ldc_hi, n, lda_lo, ldb_lo,
                &one, &one, &one,
                A_lo.ptr0, A_lo.ptr1, A_lo.ptr2,
                B_lo.ptr0, B_lo.ptr1, B_lo.ptr2,
                &zero, &zero, &zero,
                C_hi.ptr0, C_hi.ptr1, C_hi.ptr2);
        }

    } else {
        if constexpr (KIND == common::MatMulKind::AHxA) {
            matmul_block_3_strided_batched<KIND, Backend::INT8, UPLO_A, UPLO_B, UPLO_C>(
                stream, handle, ldc_hi, n, lda_lo, ldb_lo, bcnt,
                &one, &one, &one,
                A_lo.ptr0, A_lo.ptr1, A_lo.ptr2, static_cast<int64_t>(sizeA),
                A_lo.ptr1, A_lo.ptr0, A_lo.ptr2, static_cast<int64_t>(sizeA),
                &zero, &zero, &zero,
                C_hi.ptr0, C_hi.ptr1, C_hi.ptr2, static_cast<int64_t>(3 * sizeC));
        } else {
            matmul_block_3_strided_batched<KIND, Backend::INT8, UPLO_A, UPLO_B, UPLO_C>(
                stream, handle, ldc_hi, n, lda_lo, ldb_lo, bcnt,
                &one, &one, &one,
                A_lo.ptr0, A_lo.ptr1, A_lo.ptr2, static_cast<int64_t>(sizeA),
                B_lo.ptr0, B_lo.ptr1, B_lo.ptr2, static_cast<int64_t>(sizeB),
                &zero, &zero, &zero,
                C_hi.ptr0, C_hi.ptr1, C_hi.ptr2, static_cast<int64_t>(3 * sizeC));
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
    using LowT = common::low_t<Backend::FP8>;
    using HiT  = common::hi_t<Backend::FP8>;

    constexpr HiT one                = 1.0f;
    constexpr HiT zero               = 0.0f;
    constexpr unsigned max_ptr_batch = 3u * 20u;

    std::array<void *, max_ptr_batch> hA{};
    std::array<void *, max_ptr_batch> hB{};
    std::array<void *, max_ptr_batch> hC{};

    unsigned pcnt = 0;

    for (unsigned b = 0; b < bcnt; ++b) {
        const unsigned mod_idx = idx + b;

        const size_t offA = fp8_plane_offset_from_group_start(idx, b, sizeA);
        const size_t offB = fp8_plane_offset_from_group_start(idx, b, sizeB);

        LowT *A0 = A_lo.ptr0 + offA;
        LowT *B0 = B_lo.ptr0 + offB;
        HiT *C0  = C_hi.ptr0 + b * 3 * sizeC;

        if (!common::table::isKaratsuba[mod_idx]) {
            LowT *Ahi = A0;
            LowT *Alo = A0 + sizeA;
            LowT *Bhi = B0;
            LowT *Blo = B0 + sizeB;

            hA[pcnt] = Ahi;
            hB[pcnt] = Blo;
            hC[pcnt] = C0 + 0 * sizeC;
            ++pcnt;

            hA[pcnt] = Alo;
            hB[pcnt] = Bhi;
            hC[pcnt] = C0 + 1 * sizeC;
            ++pcnt;

            hA[pcnt] = Alo;
            hB[pcnt] = Blo;
            hC[pcnt] = C0 + 2 * sizeC;
            ++pcnt;
        } else {
            LowT *Ahi = A0;
            LowT *Alo = A0 + sizeA;
            LowT *As  = A0 + 2 * sizeA;
            LowT *Bhi = B0;
            LowT *Blo = B0 + sizeB;
            LowT *Bs  = B0 + 2 * sizeB;

            hA[pcnt] = Ahi;
            hB[pcnt] = Bhi;
            hC[pcnt] = C0 + 0 * sizeC;
            ++pcnt;

            hA[pcnt] = Alo;
            hB[pcnt] = Blo;
            hC[pcnt] = C0 + 1 * sizeC;
            ++pcnt;

            hA[pcnt] = As;
            hB[pcnt] = Bs;
            hC[pcnt] = C0 + 2 * sizeC;
            ++pcnt;
        }
    }

    upload_pointer_arrays(stream, handle, hA.data(), hB.data(), hC.data(), pcnt);

    common::call_gemm_tn_pointer_batched<Backend::FP8>(
        stream, handle, ldc_hi, n, lda_lo, pcnt,
        &one,
        handle.Aarray, lda_lo,
        handle.Barray, ldb_lo,
        &zero,
        handle.Carray, ldc_hi);

    const unsigned consumed = fp8_planes_consumed(idx, bcnt);

    A_lo.shift(consumed * sizeA);
    B_lo.shift(consumed * sizeB);
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
    using LowT = common::low_t<Backend::FP8>;
    using HiT  = common::hi_t<Backend::FP8>;

    constexpr HiT one  = 1.0f;
    constexpr HiT zero = 0.0f;

    if (!common::table::isKaratsuba[idx]) {
        const int64_t strideA = 2 * int64_t(sizeA);
        const int64_t strideB = 2 * int64_t(sizeB);
        const int64_t strideC = 3 * int64_t(sizeC);

        LowT *Ahi = A_lo.ptr0;
        LowT *Alo = A_lo.ptr0 + sizeA;
        LowT *Bhi = B_lo.ptr0;
        LowT *Blo = B_lo.ptr0 + sizeB;

        matmul_block_3_strided_batched<KIND, Backend::FP8, UPLO_A, UPLO_B, UPLO_C>(
            stream, handle, ldc_hi, n, lda_lo, ldb_lo, bcnt,
            &one, &one, &one,
            Ahi, Alo, Alo, strideA,
            Blo, Bhi, Blo, strideB,
            &zero, &zero, &zero,
            C_hi.ptr0 + 0 * sizeC,
            C_hi.ptr0 + 1 * sizeC,
            C_hi.ptr0 + 2 * sizeC,
            strideC);

        A_lo.shift(bcnt * 2 * sizeA);
        if constexpr (KIND != common::MatMulKind::ATxA && KIND != common::MatMulKind::AHxA) {
            B_lo.shift(bcnt * 2 * sizeB);
        }

    } else {
        const int64_t strideA = 3 * int64_t(sizeA);
        const int64_t strideB = 3 * int64_t(sizeB);
        const int64_t strideC = 3 * int64_t(sizeC);

        LowT *Ahi = A_lo.ptr0;
        LowT *Alo = A_lo.ptr0 + sizeA;
        LowT *As  = A_lo.ptr0 + 2 * sizeA;
        LowT *Bhi = B_lo.ptr0;
        LowT *Blo = B_lo.ptr0 + sizeB;
        LowT *Bs  = B_lo.ptr0 + 2 * sizeB;

        matmul_block_3_strided_batched<KIND, Backend::FP8, UPLO_A, UPLO_B, UPLO_C>(
            stream, handle, ldc_hi, n, lda_lo, ldb_lo, bcnt,
            &one, &one, &one,
            Ahi, Alo, As, strideA,
            Bhi, Blo, Bs, strideB,
            &zero, &zero, &zero,
            C_hi.ptr0 + 0 * sizeC,
            C_hi.ptr0 + 1 * sizeC,
            C_hi.ptr0 + 2 * sizeC,
            strideC);

        A_lo.shift(bcnt * 3 * sizeA);
        if constexpr (KIND != common::MatMulKind::ATxA && KIND != common::MatMulKind::AHxA) {
            B_lo.shift(bcnt * 3 * sizeB);
        }
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
    using LowT = common::low_t<Backend::FP8>;
    using HiT  = common::hi_t<Backend::FP8>;

    constexpr HiT one                = 1.0f;
    constexpr HiT zero               = 0.0f;
    constexpr unsigned max_ptr_batch = 9u * 20u;

    std::array<void *, max_ptr_batch> hA{};
    std::array<void *, max_ptr_batch> hB{};
    std::array<void *, max_ptr_batch> hC{};

    LowT *Aptr[3] = {A_lo.ptr0, A_lo.ptr1, A_lo.ptr2};
    LowT *Bptr[3] = {B_lo.ptr0, B_lo.ptr1, B_lo.ptr2};

    unsigned pcnt = 0;

    for (unsigned b = 0; b < bcnt; ++b) {
        const unsigned mod_idx = idx + b;

        const size_t offA = fp8_plane_offset_from_group_start(idx, b, sizeA);
        const size_t offB = fp8_plane_offset_from_group_start(idx, b, sizeB);

        HiT *Cbase = C_hi.ptr0 + b * 9 * sizeC;

#pragma unroll
        for (unsigned p = 0; p < 3; ++p) {
            LowT *A0 = Aptr[p] + offA;
            LowT *B0 = Bptr[p] + offB;
            HiT *C0  = Cbase + 3 * p * sizeC;

            if (!common::table::isKaratsuba[mod_idx]) {
                LowT *Ahi = A0;
                LowT *Alo = A0 + sizeA;
                LowT *Bhi = B0;
                LowT *Blo = B0 + sizeB;

                hA[pcnt] = Ahi;
                hB[pcnt] = Blo;
                hC[pcnt] = C0 + 0 * sizeC;
                ++pcnt;

                hA[pcnt] = Alo;
                hB[pcnt] = Bhi;
                hC[pcnt] = C0 + 1 * sizeC;
                ++pcnt;

                hA[pcnt] = Alo;
                hB[pcnt] = Blo;
                hC[pcnt] = C0 + 2 * sizeC;
                ++pcnt;
            } else {
                LowT *Ahi = A0;
                LowT *Alo = A0 + sizeA;
                LowT *As  = A0 + 2 * sizeA;
                LowT *Bhi = B0;
                LowT *Blo = B0 + sizeB;
                LowT *Bs  = B0 + 2 * sizeB;

                hA[pcnt] = Ahi;
                hB[pcnt] = Bhi;
                hC[pcnt] = C0 + 0 * sizeC;
                ++pcnt;

                hA[pcnt] = Alo;
                hB[pcnt] = Blo;
                hC[pcnt] = C0 + 1 * sizeC;
                ++pcnt;

                hA[pcnt] = As;
                hB[pcnt] = Bs;
                hC[pcnt] = C0 + 2 * sizeC;
                ++pcnt;
            }
        }
    }

    upload_pointer_arrays(stream, handle, hA.data(), hB.data(), hC.data(), pcnt);

    common::call_gemm_tn_pointer_batched<Backend::FP8>(
        stream, handle, ldc_hi, n, lda_lo, pcnt,
        &one,
        handle.Aarray, lda_lo,
        handle.Barray, ldb_lo,
        &zero,
        handle.Carray, ldc_hi);

    const unsigned consumed = fp8_planes_consumed(idx, bcnt);

    A_lo.shift(consumed * sizeA);
    B_lo.shift(consumed * sizeB);
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
    using LowT = common::low_t<Backend::FP8>;
    using HiT  = common::hi_t<Backend::FP8>;

    constexpr HiT one  = 1.0f;
    constexpr HiT zero = 0.0f;

    LowT *Aptr[3] = {A_lo.ptr0, A_lo.ptr1, A_lo.ptr2};
    LowT *Bptr[3];
    if constexpr (KIND == common::MatMulKind::AHxA) {
        Bptr[0] = B_lo.ptr1;
        Bptr[1] = B_lo.ptr0;
        Bptr[2] = B_lo.ptr2;
    } else {
        Bptr[0] = B_lo.ptr0;
        Bptr[1] = B_lo.ptr1;
        Bptr[2] = B_lo.ptr2;
    }

    const int64_t strideC = 9 * int64_t(sizeC);

    if (!common::table::isKaratsuba[idx]) {
        const int64_t strideA = 2 * int64_t(sizeA);
        const int64_t strideB = 2 * int64_t(sizeB);

#pragma unroll
        for (unsigned p = 0; p < 3; ++p) {
            LowT *Ahi = Aptr[p];
            LowT *Alo = Aptr[p] + sizeA;
            LowT *Bhi = Bptr[p];
            LowT *Blo = Bptr[p] + sizeB;

            HiT *C0 = C_hi.ptr0 + 3 * p * sizeC;

            matmul_block_3_strided_batched<KIND, Backend::FP8, UPLO_A, UPLO_B, UPLO_C>(
                stream, handle, ldc_hi, n, lda_lo, ldb_lo, bcnt,
                &one, &one, &one,
                Ahi, Alo, Alo, strideA,
                Blo, Bhi, Blo, strideB,
                &zero, &zero, &zero,
                C0 + 0 * sizeC,
                C0 + 1 * sizeC,
                C0 + 2 * sizeC,
                strideC);
        }

        A_lo.shift(bcnt * 2 * sizeA);
        if constexpr (KIND != common::MatMulKind::ATxA && KIND != common::MatMulKind::AHxA) {
            B_lo.shift(bcnt * 2 * sizeB);
        }

    } else {
        const int64_t strideA = 3 * int64_t(sizeA);
        const int64_t strideB = 3 * int64_t(sizeB);

#pragma unroll
        for (unsigned p = 0; p < 3; ++p) {
            LowT *Ahi = Aptr[p];
            LowT *Alo = Aptr[p] + sizeA;
            LowT *As  = Aptr[p] + 2 * sizeA;
            LowT *Bhi = Bptr[p];
            LowT *Blo = Bptr[p] + sizeB;
            LowT *Bs  = Bptr[p] + 2 * sizeB;

            HiT *C0 = C_hi.ptr0 + 3 * p * sizeC;

            matmul_block_3_strided_batched<KIND, Backend::FP8, UPLO_A, UPLO_B, UPLO_C>(
                stream, handle, ldc_hi, n, lda_lo, ldb_lo, bcnt,
                &one, &one, &one,
                Ahi, Alo, As, strideA,
                Bhi, Blo, Bs, strideB,
                &zero, &zero, &zero,
                C0 + 0 * sizeC,
                C0 + 1 * sizeC,
                C0 + 2 * sizeC,
                strideC);
        }

        A_lo.shift(bcnt * 3 * sizeA);
        if constexpr (KIND != common::MatMulKind::ATxA && KIND != common::MatMulKind::AHxA) {
            B_lo.shift(bcnt * 3 * sizeB);
        }
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
    unsigned done = 0;

    while (done < bcnt) {
        const unsigned cur   = idx + done;
        const bool karatsuba = common::table::isKaratsuba[cur];

        unsigned cnt = 1;
        while (done + cnt < bcnt && common::table::isKaratsuba[cur + cnt] == karatsuba) {
            ++cnt;
        }

        auto C_part = C_hi;
        C_part.shift(size_t(done) * 3 * sizeC);

        error_free_matmult_f8_real_strided<KIND, UPLO_A, UPLO_B, UPLO_C>(
            stream, handle,
            cur, cnt,
            ldc_hi, n, lda_lo, ldb_lo,
            sizeA, sizeB, sizeC,
            A_lo, B_lo, C_part);

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
    unsigned done = 0;

    while (done < bcnt) {
        const unsigned cur   = idx + done;
        const bool karatsuba = common::table::isKaratsuba[cur];

        unsigned cnt = 1;
        while (done + cnt < bcnt && common::table::isKaratsuba[cur + cnt] == karatsuba) {
            ++cnt;
        }

        auto C_part = C_hi;
        C_part.shift(size_t(done) * 9 * sizeC);

        error_free_matmult_f8_complex_strided<KIND, UPLO_A, UPLO_B, UPLO_C>(
            stream, handle,
            cur, cnt,
            ldc_hi, n, lda_lo, ldb_lo,
            sizeA, sizeB, sizeC,
            A_lo, B_lo, C_part);

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
