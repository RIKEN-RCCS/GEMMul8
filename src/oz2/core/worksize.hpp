#pragma once
#include "../common/common.hpp"
#include "../common/table.hpp"
#include "helper_matmult.hpp"

namespace gemmul8::oz2::core {

// gemm, symm, syr2k, syrkx, trmm, hemm, her2k, herkx, trtrmm
template <bool is_Complex, Backend BACKEND,
          common::MatMulKind KIND>
inline size_t workSize(
    size_t m, size_t n, size_t k, unsigned NUM_MODULI,
    bool enable_skip_scalA, bool enable_skip_scalB,
    size_t *workSizeA, size_t *workSizeB, bool fastmode //
) {
    using LowT = common::low_t<BACKEND>;
    using MidT = common::mid_t<BACKEND, is_Complex>;
    using HiT  = common::hi_t<BACKEND>;

    // sizes
    const size_t m_pad          = common::padding(m);
    const size_t n_pad          = common::padding(n);
    const size_t k_pad          = common::padding(k);
    const size_t n_work         = (KIND == common::MatMulKind::Trtrmm) ? n_pad : n;
    const size_t sizeA          = k_pad * m_pad;
    const size_t sizeB          = k_pad * n_work;
    const size_t sizeC          = m_pad * n_work;
    const size_t size_vecA      = m_pad;
    const size_t size_vecB      = n_pad;
    const unsigned num_mat      = common::table::num_mat<BACKEND>(NUM_MODULI);
    constexpr size_t lwork_blas = size_t(32) << 20; // 32 MiB

    unsigned num_A_lo  = num_mat + ((enable_skip_scalA && !fastmode) ? 1 : 0); // +1 for skip_scalA in accurate mode
    unsigned num_B_lo  = num_mat + ((enable_skip_scalB && !fastmode) ? 1 : 0); // +1 for skip_scalB in accurate mode
    unsigned num_C_mid = NUM_MODULI;
    unsigned num_C_hi  = (BACKEND == Backend::INT8) ? 1 : 3;

    if constexpr (is_Complex) {
        num_A_lo *= 3;
        num_B_lo *= 3;
        num_C_hi *= 3;
    }

    constexpr bool use_pointer_arrays               = common::isCUDA && (BACKEND == Backend::FP8) && (KIND == common::MatMulKind::Gemm);
    constexpr unsigned pointer_products_per_modulus = use_pointer_arrays ? ((is_Complex) ? 9u : 3u) : 0u;
    const unsigned pointer_batch_count_max          = pointer_products_per_modulus * NUM_MODULI;
    const size_t pointer_array_bytes                = use_pointer_arrays ? common::padding(3 * pointer_batch_count_max * sizeof(void *)) : 0u;

    const size_t sizeC_Mid = sizeof(MidT) * sizeC;
    const size_t sizeC_Hi  = sizeof(HiT) * sizeC * num_C_hi;

    size_t total_size_A = common::PAD_SIZE - 1;
    size_t total_size_B = common::PAD_SIZE - 1;
    size_t total_size_C = common::PAD_SIZE - 1;

    total_size_A += sizeof(LowT) * sizeA * num_A_lo;
    total_size_A += sizeof(int16_t) * size_vecA * ((enable_skip_scalA && !fastmode) ? 2u : 1u);

    total_size_B += sizeof(LowT) * sizeB * num_B_lo;
    total_size_B += sizeof(int16_t) * size_vecB * ((enable_skip_scalB && !fastmode) ? 2u : 1u);

    total_size_C += sizeC_Mid * (num_C_mid - 1);
    total_size_C += std::max<size_t>(pointer_array_bytes + lwork_blas, sizeC_Mid);
    total_size_C += sizeC_Hi;

    if (workSizeA != nullptr) *workSizeA = total_size_A;
    if (workSizeB != nullptr) *workSizeB = total_size_B;

    return total_size_A + total_size_B + total_size_C;
}

// syrk, herk
template <bool is_Complex, Backend BACKEND>
inline size_t workSize_rk(
    size_t n, size_t k, unsigned NUM_MODULI,
    size_t *workSizeA //
) {
    using LowT = common::low_t<BACKEND>;
    using MidT = common::mid_t<BACKEND, is_Complex>;
    using HiT  = common::hi_t<BACKEND>;

    // sizes
    const size_t n_pad          = common::padding(n);
    const size_t k_pad          = common::padding(k);
    const size_t lda_lo         = k_pad;
    const size_t ldc_hi         = n_pad;
    const size_t sizeA          = k_pad * n_pad;
    const size_t sizeC          = n_pad * n;
    const size_t size_vecA      = n_pad;
    const unsigned num_mat      = common::table::num_mat<BACKEND>(NUM_MODULI);
    constexpr size_t lwork_blas = size_t(32) << 20; // 32 MiB

    unsigned num_A_lo  = num_mat;
    unsigned num_C_mid = NUM_MODULI;
    unsigned num_C_hi  = (BACKEND == Backend::INT8) ? 1 : 3;

    if constexpr (is_Complex) {
        num_A_lo *= 3;
        num_C_hi *= 3;
    }

    const size_t sizeC_Mid = sizeof(MidT) * sizeC;
    const size_t sizeC_Hi  = sizeof(HiT) * sizeC * num_C_hi;

    size_t total_size_A = common::PAD_SIZE - 1;
    size_t total_size_C = common::PAD_SIZE - 1;

    total_size_A += sizeof(LowT) * sizeA * num_A_lo;
    total_size_A += sizeof(int16_t) * size_vecA;

    total_size_C += sizeC_Mid * (num_C_mid - 1);
    total_size_C += std::max<size_t>(lwork_blas, sizeC_Mid);
    total_size_C += sizeC_Hi;

    if (workSizeA != nullptr) *workSizeA = total_size_A;

    return total_size_A + total_size_C;
}

} // namespace gemmul8::oz2::core
