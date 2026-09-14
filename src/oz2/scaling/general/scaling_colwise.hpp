#pragma once
#include "../../common/common.hpp"
#include "../../common/table.hpp"
#include "../../mod/mod.hpp"

#include "config.hpp"
#include "roundup.hpp"
#include "helper_triangular.hpp"
#include "store.hpp"

namespace gemmul8::scaling::general {

template <typename T, Backend BACKEND, unsigned NUM_MODULI, bool CONJ>
__device__ __forceinline__ void scaling_colwise_full_device(
    const unsigned rows_A,
    const T *const __restrict__ in,
    common::matptr_t<common::low_t<BACKEND>, common::isComplex<T>> A_lo,
    const size_t lda_lo4, const size_t incA_lo4,
    const int32_t sft //
) {
    using Low4              = common::lowx4_t<BACKEND>;
    const unsigned col_idx  = blockIdx.x;
    Low4 *__restrict__ out0 = reinterpret_cast<Low4 *>(A_lo.ptr0) + col_idx * lda_lo4;
    Low4 *__restrict__ out1 = nullptr;
    if constexpr (common::isComplex<T>) {
        out1 = reinterpret_cast<Low4 *>(A_lo.ptr1) + col_idx * lda_lo4;
    }

    for (unsigned i = threadIdx.x; i < lda_lo4; i += blockDim.x) {
        const unsigned row0 = i << 2;
        if (row0 >= rows_A) {
            if constexpr (common::isComplex<T>) {
                mod::ModUnrollFillZero<BACKEND, NUM_MODULI, Low4>::run(out0 + i, out1 + i, incA_lo4);
            } else {
                mod::ModUnrollFillZero<BACKEND, NUM_MODULI, Low4>::run(out0 + i, incA_lo4);
            }
            continue;
        }
        const T z  = common::Tconst<T>::zero();
        const T a0 = in[row0];
        const T a1 = row0 + 1U < rows_A ? in[row0 + 1U] : z;
        const T a2 = row0 + 2U < rows_A ? in[row0 + 2U] : z;
        const T a3 = row0 + 3U < rows_A ? in[row0 + 3U] : z;
        if constexpr (common::isComplex<T>) {
            scaling_colwise_store4_complex<T, BACKEND, NUM_MODULI, CONJ>(
                out0, out1, incA_lo4, i, a0, a1, a2, a3, sft);
        } else {
            scaling_colwise_store4_real<T, BACKEND, NUM_MODULI, CONJ>(
                out0, incA_lo4, i, a0, a1, a2, a3, sft);
        }
    }
}

template <typename T, Backend BACKEND, unsigned NUM_MODULI,
          cublasDiagType_t DIAG, bool CONJ>
__device__ __forceinline__ void scaling_colwise_upper_device(
    const unsigned rows_A,
    const T *const __restrict__ in,
    common::matptr_t<common::low_t<BACKEND>, common::isComplex<T>> A_lo,
    const size_t lda_lo4, const size_t incA_lo4,
    const int32_t sft //
) {
    using Low4             = common::lowx4_t<BACKEND>;
    const unsigned col_idx = blockIdx.x;

    unsigned boundary_i;
    bool has_boundary;
    if constexpr (DIAG == CUBLAS_DIAG_UNIT) {
        const bool diag_flag = (col_idx < rows_A);
        boundary_i           = (diag_flag ? col_idx : rows_A) >> 2;
        has_boundary         = diag_flag ? true : ((rows_A & 3U) != 0U);
    } else {
        const unsigned active_end = min(rows_A, col_idx + 1U);
        boundary_i                = active_end >> 2;
        has_boundary              = ((active_end & 3U) != 0U);
    }

    if constexpr (common::isComplex<T>) {
        Low4 *__restrict__ out_1 = reinterpret_cast<Low4 *>(A_lo.ptr0) + col_idx * lda_lo4;
        Low4 *__restrict__ out_2 = reinterpret_cast<Low4 *>(A_lo.ptr1) + col_idx * lda_lo4;

        for (unsigned i = threadIdx.x; i < boundary_i; i += blockDim.x) {
            const unsigned row0 = i << 2;

            scaling_colwise_store4_complex<T, BACKEND, NUM_MODULI, CONJ>(
                out_1, out_2, incA_lo4, i,
                in[row0], in[row0 + 1U], in[row0 + 2U], in[row0 + 3U], sft);
        }

        if (threadIdx.x == 0 && has_boundary) {
            const unsigned row0 = boundary_i << 2;

            const T a0 = tri_col_value<true, T, DIAG>(in, row0, col_idx, rows_A);
            const T a1 = tri_col_value<true, T, DIAG>(in, row0 + 1U, col_idx, rows_A);
            const T a2 = tri_col_value<true, T, DIAG>(in, row0 + 2U, col_idx, rows_A);
            const T a3 = tri_col_value<true, T, DIAG>(in, row0 + 3U, col_idx, rows_A);

            scaling_colwise_store4_complex<T, BACKEND, NUM_MODULI, CONJ>(
                out_1, out_2, incA_lo4, boundary_i,
                a0, a1, a2, a3, sft);
        }

    } else {
        Low4 *__restrict__ out = reinterpret_cast<Low4 *>(A_lo.ptr0) + col_idx * lda_lo4;

        for (unsigned i = threadIdx.x; i < boundary_i; i += blockDim.x) {
            const unsigned row0 = i << 2;

            scaling_colwise_store4_real<T, BACKEND, NUM_MODULI, CONJ>(
                out, incA_lo4, i,
                in[row0], in[row0 + 1U], in[row0 + 2U], in[row0 + 3U], sft);
        }

        if (threadIdx.x == 0 && has_boundary) {
            const unsigned row0 = boundary_i << 2;

            const T a0 = tri_col_value<true, T, DIAG>(in, row0, col_idx, rows_A);
            const T a1 = tri_col_value<true, T, DIAG>(in, row0 + 1U, col_idx, rows_A);
            const T a2 = tri_col_value<true, T, DIAG>(in, row0 + 2U, col_idx, rows_A);
            const T a3 = tri_col_value<true, T, DIAG>(in, row0 + 3U, col_idx, rows_A);

            scaling_colwise_store4_real<T, BACKEND, NUM_MODULI, CONJ>(
                out, incA_lo4, boundary_i,
                a0, a1, a2, a3, sft);
        }
    }
}

template <typename T, Backend BACKEND, unsigned NUM_MODULI,
          cublasDiagType_t DIAG, bool CONJ>
__device__ __forceinline__ void scaling_colwise_lower_device(
    const unsigned rows_A,
    const T *const __restrict__ in,
    common::matptr_t<common::low_t<BACKEND>, common::isComplex<T>> A_lo,
    const size_t lda_lo4, const size_t incA_lo4,
    const int32_t sft //
) {
    using Low4             = common::lowx4_t<BACKEND>;
    const unsigned col_idx = blockIdx.x;
    if (col_idx >= rows_A) return;

    const unsigned first4 = col_idx >> 2;
    const unsigned tail4  = rows_A >> 2;
    bool has_head_boundary;
    if constexpr (DIAG == CUBLAS_DIAG_UNIT) {
        has_head_boundary = true;
    } else {
        has_head_boundary = ((col_idx & 3U) != 0U);
    }
    const bool has_tail_boundary = ((rows_A & 3U) != 0U);
    const unsigned dense_begin4  = first4 + (has_head_boundary ? 1U : 0U);

    if constexpr (common::isComplex<T>) {
        Low4 *__restrict__ out_1 = reinterpret_cast<Low4 *>(A_lo.ptr0) + col_idx * lda_lo4;
        Low4 *__restrict__ out_2 = reinterpret_cast<Low4 *>(A_lo.ptr1) + col_idx * lda_lo4;

        for (unsigned i = dense_begin4 + threadIdx.x; i < tail4; i += blockDim.x) {
            const unsigned row0 = i << 2;

            scaling_colwise_store4_complex<T, BACKEND, NUM_MODULI, CONJ>(
                out_1, out_2, incA_lo4, i,
                in[row0], in[row0 + 1U], in[row0 + 2U], in[row0 + 3U], sft);
        }

        if (threadIdx.x == 0) {
            if (has_head_boundary) {
                const unsigned row0 = first4 << 2;

                const T a0 = tri_col_value<false, T, DIAG>(in, row0, col_idx, rows_A);
                const T a1 = tri_col_value<false, T, DIAG>(in, row0 + 1U, col_idx, rows_A);
                const T a2 = tri_col_value<false, T, DIAG>(in, row0 + 2U, col_idx, rows_A);
                const T a3 = tri_col_value<false, T, DIAG>(in, row0 + 3U, col_idx, rows_A);

                scaling_colwise_store4_complex<T, BACKEND, NUM_MODULI, CONJ>(
                    out_1, out_2, incA_lo4, first4,
                    a0, a1, a2, a3, sft);
            }

            if (has_tail_boundary && (!has_head_boundary || tail4 != first4)) {
                const unsigned row0 = tail4 << 2;

                const T a0 = tri_col_value<false, T, DIAG>(in, row0, col_idx, rows_A);
                const T a1 = tri_col_value<false, T, DIAG>(in, row0 + 1U, col_idx, rows_A);
                const T a2 = tri_col_value<false, T, DIAG>(in, row0 + 2U, col_idx, rows_A);
                const T a3 = tri_col_value<false, T, DIAG>(in, row0 + 3U, col_idx, rows_A);

                scaling_colwise_store4_complex<T, BACKEND, NUM_MODULI, CONJ>(
                    out_1, out_2, incA_lo4, tail4,
                    a0, a1, a2, a3, sft);
            }
        }

    } else {
        Low4 *__restrict__ out = reinterpret_cast<Low4 *>(A_lo.ptr0) + col_idx * lda_lo4;

        for (unsigned i = dense_begin4 + threadIdx.x; i < tail4; i += blockDim.x) {
            const unsigned row0 = i << 2;

            scaling_colwise_store4_real<T, BACKEND, NUM_MODULI, CONJ>(
                out, incA_lo4, i,
                in[row0], in[row0 + 1U], in[row0 + 2U], in[row0 + 3U], sft);
        }

        if (threadIdx.x == 0) {
            if (has_head_boundary) {
                const unsigned row0 = first4 << 2;

                const T a0 = tri_col_value<false, T, DIAG>(in, row0, col_idx, rows_A);
                const T a1 = tri_col_value<false, T, DIAG>(in, row0 + 1U, col_idx, rows_A);
                const T a2 = tri_col_value<false, T, DIAG>(in, row0 + 2U, col_idx, rows_A);
                const T a3 = tri_col_value<false, T, DIAG>(in, row0 + 3U, col_idx, rows_A);

                scaling_colwise_store4_real<T, BACKEND, NUM_MODULI, CONJ>(
                    out, incA_lo4, first4,
                    a0, a1, a2, a3, sft);
            }

            if (has_tail_boundary && (!has_head_boundary || tail4 != first4)) {
                const unsigned row0 = tail4 << 2;

                const T a0 = tri_col_value<false, T, DIAG>(in, row0, col_idx, rows_A);
                const T a1 = tri_col_value<false, T, DIAG>(in, row0 + 1U, col_idx, rows_A);
                const T a2 = tri_col_value<false, T, DIAG>(in, row0 + 2U, col_idx, rows_A);
                const T a3 = tri_col_value<false, T, DIAG>(in, row0 + 3U, col_idx, rows_A);

                scaling_colwise_store4_real<T, BACKEND, NUM_MODULI, CONJ>(
                    out, incA_lo4, tail4,
                    a0, a1, a2, a3, sft);
            }
        }
    }
}

template <typename T, Backend BACKEND, unsigned NUM_MODULI, bool CONJ>
__global__ void scaling_colwise_full_tiled_kernel(
    const unsigned rows_A,
    const unsigned cols_A,
    const T *const __restrict__ A,
    const size_t lda,
    common::matptr_t<common::low_t<BACKEND>, common::isComplex<T>> A_lo,
    const size_t lda_lo4,
    const size_t incA_lo4,
    const int16_t *const __restrict__ sftA //
) {
    using Low4 = common::lowx4_t<BACKEND>;

    const unsigned r4      = blockIdx.x * threads_x_colwise_full_tiled + threadIdx.x;
    const unsigned col_idx = blockIdx.y * threads_y_colwise_full_tiled + threadIdx.y;

    if (col_idx >= cols_A || r4 >= lda_lo4) return;

    const unsigned row0 = r4 << 2;
    const int32_t sft   = -sftA[col_idx];

    const T *const __restrict__ in = A + col_idx * lda;

    if (row0 >= rows_A) {
        Low4 *__restrict__ out0 = reinterpret_cast<Low4 *>(A_lo.ptr0) + col_idx * lda_lo4 + r4;
        if constexpr (common::isComplex<T>) {
            Low4 *__restrict__ out1 = reinterpret_cast<Low4 *>(A_lo.ptr1) + col_idx * lda_lo4 + r4;
            mod::ModUnrollFillZero<BACKEND, NUM_MODULI, Low4>::run(out0, out1, incA_lo4);
        } else {
            mod::ModUnrollFillZero<BACKEND, NUM_MODULI, Low4>::run(out0, incA_lo4);
        }
        return;
    }

    // Share the expensive modular conversion between dense and tail groups.
    const T z               = common::Tconst<T>::zero();
    const T a0              = in[row0];
    const T a1              = row0 + 1U < rows_A ? in[row0 + 1U] : z;
    const T a2              = row0 + 2U < rows_A ? in[row0 + 2U] : z;
    const T a3              = row0 + 3U < rows_A ? in[row0 + 3U] : z;
    Low4 *__restrict__ out0 = reinterpret_cast<Low4 *>(A_lo.ptr0) + col_idx * lda_lo4;
    if constexpr (common::isComplex<T>) {
        Low4 *__restrict__ out1 = reinterpret_cast<Low4 *>(A_lo.ptr1) + col_idx * lda_lo4;
        scaling_colwise_store4_complex<T, BACKEND, NUM_MODULI, CONJ>(
            out0, out1, incA_lo4, r4, a0, a1, a2, a3, sft);
    } else {
        scaling_colwise_store4_real<T, BACKEND, NUM_MODULI, CONJ>(
            out0, incA_lo4, r4, a0, a1, a2, a3, sft);
    }
}

template <typename T>
__device__ __forceinline__ common::vec4_t<T> load_vec4_aligned(const T *ptr) {
    return *reinterpret_cast<const common::vec4_t<T> *>(ptr);
}

template <typename T, Backend BACKEND, unsigned NUM_MODULI, bool CONJ>
__global__ void scaling_colwise_full_tiled_aligned_kernel(
    const unsigned rows_A,
    const unsigned cols_A,
    const T *const __restrict__ A,
    const size_t lda,
    common::matptr_t<common::low_t<BACKEND>, common::isComplex<T>> A_lo,
    const size_t lda_lo4,
    const size_t incA_lo4,
    const int16_t *const __restrict__ sftA //
) {
    using Low4 = common::lowx4_t<BACKEND>;

    const unsigned r4      = blockIdx.x * threads_x_colwise_full_tiled + threadIdx.x;
    const unsigned col_idx = blockIdx.y * threads_y_colwise_full_tiled + threadIdx.y;

    if (col_idx >= cols_A) return;

    const unsigned row0 = r4 << 2;
    const int32_t sft   = -sftA[col_idx];

    const T *const __restrict__ in = A + col_idx * lda;

    if constexpr (common::isComplex<T>) {
        Low4 *__restrict__ out_1 = reinterpret_cast<Low4 *>(A_lo.ptr0) + col_idx * lda_lo4;
        Low4 *__restrict__ out_2 = reinterpret_cast<Low4 *>(A_lo.ptr1) + col_idx * lda_lo4;

        scaling_colwise_store4_complex<T, BACKEND, NUM_MODULI, CONJ>(
            out_1, out_2, incA_lo4, r4,
            in[row0 + 0U], in[row0 + 1U], in[row0 + 2U], in[row0 + 3U], sft);

    } else {
        Low4 *__restrict__ out     = reinterpret_cast<Low4 *>(A_lo.ptr0) + col_idx * lda_lo4;
        const common::vec4_t<T> a4 = *reinterpret_cast<const common::vec4_t<T> *>(in + row0);

        scaling_colwise_store4_real<T, BACKEND, NUM_MODULI, CONJ>(
            out, incA_lo4, r4,
            a4.x, a4.y, a4.z, a4.w, sft);
    }
}

template <typename T, Backend BACKEND, unsigned NUM_MODULI,
          cublasFillMode_t UPLO, cublasDiagType_t DIAG, bool CONJ>
__global__ void scaling_colwise_tri_tiled_kernel(
    const unsigned rows_A,
    const unsigned cols_A,
    const T *const __restrict__ A,
    const size_t lda,
    common::matptr_t<common::low_t<BACKEND>, common::isComplex<T>> A_lo,
    const size_t lda_lo4,
    const size_t incA_lo4,
    const int16_t *const __restrict__ sftA //
) {
    static_assert(UPLO == CUBLAS_FILL_MODE_UPPER || UPLO == CUBLAS_FILL_MODE_LOWER);

    using Low4 = common::lowx4_t<BACKEND>;

    const unsigned r4      = blockIdx.x * threads_x_colwise_full_tiled + threadIdx.x;
    const unsigned col_idx = blockIdx.y * threads_y_colwise_full_tiled + threadIdx.y;

    const unsigned rows4 = (rows_A + 3U) >> 2;

    if (col_idx >= cols_A || r4 >= rows4) return;

    bool dense    = false;
    bool boundary = false;

    if constexpr (UPLO == CUBLAS_FILL_MODE_UPPER) {
        unsigned boundary_i;
        bool has_boundary;

        if constexpr (DIAG == CUBLAS_DIAG_UNIT) {
            const bool diag_flag = (col_idx < rows_A);
            boundary_i           = (diag_flag ? col_idx : rows_A) >> 2;
            has_boundary         = diag_flag ? true : ((rows_A & 3U) != 0U);
        } else {
            const unsigned active_end = min(rows_A, col_idx + 1U);
            boundary_i                = active_end >> 2;
            has_boundary              = ((active_end & 3U) != 0U);
        }

        dense    = (r4 < boundary_i);
        boundary = (r4 == boundary_i && has_boundary);

    } else {
        if (col_idx >= rows_A) return;

        const unsigned first4 = col_idx >> 2;
        const unsigned tail4  = rows_A >> 2;

        bool has_head_boundary;
        if constexpr (DIAG == CUBLAS_DIAG_UNIT) {
            has_head_boundary = true;
        } else {
            has_head_boundary = ((col_idx & 3U) != 0U);
        }

        const bool has_tail_boundary = ((rows_A & 3U) != 0U);
        const unsigned dense_begin4  = first4 + (has_head_boundary ? 1U : 0U);

        dense = (r4 >= dense_begin4 && r4 < tail4);

        const bool head_boundary = (r4 == first4 && has_head_boundary);
        const bool tail_boundary =
            (r4 == tail4 && has_tail_boundary && (!has_head_boundary || tail4 != first4));

        boundary = head_boundary || tail_boundary;
    }

    if (!dense && !boundary) return;

    const unsigned row0 = r4 << 2;
    const int32_t sft   = -sftA[col_idx];

    const T *const __restrict__ in = A + size_t(col_idx) * lda;

    T a0, a1, a2, a3;

    if (dense) {
        a0 = in[row0 + 0U];
        a1 = in[row0 + 1U];
        a2 = in[row0 + 2U];
        a3 = in[row0 + 3U];
    } else {
        if constexpr (UPLO == CUBLAS_FILL_MODE_UPPER) {
            a0 = tri_col_value<true, T, DIAG>(in, row0 + 0U, col_idx, rows_A);
            a1 = tri_col_value<true, T, DIAG>(in, row0 + 1U, col_idx, rows_A);
            a2 = tri_col_value<true, T, DIAG>(in, row0 + 2U, col_idx, rows_A);
            a3 = tri_col_value<true, T, DIAG>(in, row0 + 3U, col_idx, rows_A);
        } else {
            a0 = tri_col_value<false, T, DIAG>(in, row0 + 0U, col_idx, rows_A);
            a1 = tri_col_value<false, T, DIAG>(in, row0 + 1U, col_idx, rows_A);
            a2 = tri_col_value<false, T, DIAG>(in, row0 + 2U, col_idx, rows_A);
            a3 = tri_col_value<false, T, DIAG>(in, row0 + 3U, col_idx, rows_A);
        }
    }

    if constexpr (common::isComplex<T>) {
        Low4 *__restrict__ out_1 = reinterpret_cast<Low4 *>(A_lo.ptr0) + size_t(col_idx) * lda_lo4;
        Low4 *__restrict__ out_2 = reinterpret_cast<Low4 *>(A_lo.ptr1) + size_t(col_idx) * lda_lo4;

        scaling_colwise_store4_complex<T, BACKEND, NUM_MODULI, CONJ>(
            out_1, out_2, incA_lo4, r4,
            a0, a1, a2, a3, sft);

    } else {
        Low4 *__restrict__ out = reinterpret_cast<Low4 *>(A_lo.ptr0) + size_t(col_idx) * lda_lo4;

        scaling_colwise_store4_real<T, BACKEND, NUM_MODULI, CONJ>(
            out, incA_lo4, r4,
            a0, a1, a2, a3, sft);
    }
}

} // namespace gemmul8::scaling::general
