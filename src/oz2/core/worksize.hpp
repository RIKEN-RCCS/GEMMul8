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
    const unsigned num_mat      = common::table::num_mat<BACKEND, is_Complex>(NUM_MODULI);
    constexpr size_t lwork_blas = size_t(32) << 20; // 32 MiB

    unsigned num_A_lo  = num_mat + ((enable_skip_scalA && !fastmode) ? 1 : 0); // +1 for skip_scalA in accurate mode
    unsigned num_B_lo  = num_mat + ((enable_skip_scalB && !fastmode) ? 1 : 0); // +1 for skip_scalB in accurate mode
    unsigned num_C_mid = NUM_MODULI;
    unsigned num_C_hi  = (BACKEND == Backend::INT8) ? 1 : 3;

    if constexpr (is_Complex) {
        num_A_lo *= 2;
        num_B_lo *= 2;
        num_C_hi *= 2;
        num_A_lo += (enable_skip_scalA && !fastmode) ? 1u : 0u;
        num_B_lo += (enable_skip_scalB && !fastmode) ? 1u : 0u;
    }

    constexpr bool use_pointer_arrays               = common::isCUDA && (BACKEND == Backend::FP8) && (KIND == common::MatMulKind::Gemm);
    constexpr unsigned pointer_products_per_modulus = use_pointer_arrays ? ((is_Complex) ? 6u : 3u) : 0u;
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

template <bool COMPLEX, Backend BACKEND, bool HERK = false>
inline size_t rankk_worksize_C(size_t sizeC, unsigned num_moduli, bool fastmode) {
    static_assert(!HERK || COMPLEX);
    constexpr size_t lwork_blas    = size_t(32) << 20;
    constexpr unsigned products    = (COMPLEX && !HERK) ? 2u : 1u;
    constexpr unsigned limbs       = (BACKEND == Backend::INT8) ? 1u : 3u;
    const size_t mid               = sizeof(common::mid_t<BACKEND, COMPLEX>) * sizeC;
    const size_t hi                = sizeof(common::hi_t<BACKEND>) * sizeC * products * limbs;
    const size_t multiplication    = (num_moduli - 1u) * mid + std::max(lwork_blas, mid) + hi;
    constexpr unsigned norm_planes = COMPLEX ? ((BACKEND == Backend::INT8) ? 2u : 3u) : 1u;
    const size_t norm              = fastmode ? 0u : sizeof(common::hi_t<BACKEND>) * sizeC * norm_planes + lwork_blas;
    return std::max(multiplication, norm);
}

// syrk, herk
template <bool is_Complex, Backend BACKEND, bool HERK = false>
inline size_t workSize_rk(
    size_t n, size_t k, unsigned NUM_MODULI,
    size_t *workSizeA, bool fastmode = false //
) {
    using LowT = common::low_t<BACKEND>;

    // sizes
    const size_t n_pad     = common::padding(n);
    const size_t k_pad     = common::padding(k);
    const size_t sizeA     = k_pad * n_pad;
    const size_t sizeC     = n_pad * n;
    const size_t size_vecA = n_pad;
    const unsigned num_mat = common::table::num_mat<BACKEND, is_Complex>(NUM_MODULI);
    unsigned num_A_lo      = num_mat;
    if constexpr (is_Complex) num_A_lo *= 2;

    size_t total_size_A = common::PAD_SIZE - 1;
    size_t total_size_C = common::PAD_SIZE - 1;

    total_size_A += sizeof(LowT) * sizeA * num_A_lo;
    total_size_A += sizeof(int16_t) * size_vecA;

    total_size_C += rankk_worksize_C<is_Complex, BACKEND, HERK>(sizeC, NUM_MODULI, fastmode);

    if (workSizeA != nullptr) *workSizeA = total_size_A;

    return total_size_A + total_size_C;
}

} // namespace gemmul8::oz2::core
