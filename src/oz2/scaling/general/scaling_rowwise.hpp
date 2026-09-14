#pragma once
#include "../../common/common.hpp"
#include "../../mod/mod.hpp"

#include "config.hpp"
#include "roundup.hpp"
#include "helper_triangular.hpp"
#include "store.hpp"

namespace gemmul8::scaling::general {

template <typename T, Backend BACKEND, unsigned NUM_MODULI, bool CONJ>
__global__ void scaling_rowwise_full_kernel(
    const unsigned rows_A, const unsigned cols_A,
    const T *const __restrict__ A, const size_t lda,
    common::matptr_t<common::low_t<BACKEND>, common::isComplex<T>> A_lo,
    const size_t lda_lo, const size_t incA_lo,
    int16_t *const __restrict__ sftA //
) {
    using ValT               = decltype(trunc_scalbn<true, T, BACKEND, NUM_MODULI>::run(T{}, int32_t{}));
    constexpr bool cast_flag = (sizeof(ValT) <= sizeof(T)) && common::isCUDA;
    using shm_t              = std::conditional_t<(cast_flag), ValT, T>;

    __shared__ shm_t tile[common::TILE_DIM][common::TILE_DIM + 1];

    const unsigned rowBase = blockIdx.x * common::TILE_DIM;
    const unsigned colBase = blockIdx.y * common::TILE_DIM;

    const unsigned in_row = rowBase + threadIdx.x;
    const int32_t sft     = (in_row < rows_A) ? -sftA[in_row] : 0;

#pragma unroll
    for (unsigned j = 0; j < common::TILE_DIM; j += threads_y_rowwise<BACKEND>) {
        const unsigned yy = threadIdx.y + j;
        if (yy >= common::TILE_DIM) continue;

        const unsigned in_col = colBase + yy;
        const T Atmp          = (in_row < rows_A && in_col < cols_A) ? A[in_col * lda + in_row] : common::Tconst<T>::zero();
        const shm_t A_scaled  = trunc_scalbn<cast_flag, T, BACKEND, NUM_MODULI>::run(common::conj<T, CONJ>(Atmp), sft);
        tile[yy][threadIdx.x] = A_scaled;
    }
    __syncthreads();

    constexpr unsigned mod_unroll = sizeof(common::underlying_t<T>) == sizeof(float)
                                        ? common::TILE_DIM / threads_y_rowwise<BACKEND>
                                        : 1U;
#pragma unroll mod_unroll
    for (unsigned j = 0; j < common::TILE_DIM; j += threads_y_rowwise<BACKEND>) {
        const unsigned yy = threadIdx.y + j;
        if (yy >= common::TILE_DIM) continue;

        const unsigned out_col = rowBase + yy;
        const unsigned out_row = colBase + threadIdx.x;
        if (out_col >= rows_A || out_row >= cols_A) continue;

        const size_t idx = out_col * lda_lo + out_row;
        const ValT in    = trunc_scalbn<cast_flag, T, BACKEND, NUM_MODULI>::cast(tile[threadIdx.x][yy]);

        if constexpr (common::isComplex<T>) {
            common::low_t<BACKEND> *__restrict__ out_1 = A_lo.ptr0 + idx;
            common::low_t<BACKEND> *__restrict__ out_2 = A_lo.ptr1 + idx;
            mod::ModUnroll<NUM_MODULI, ValT>::run(out_1, out_2, incA_lo, in);
        } else {
            common::low_t<BACKEND> *__restrict__ out = A_lo.ptr0 + idx;
            mod::ModUnroll<NUM_MODULI, ValT>::run(out, incA_lo, in);
        }
    }
}

template <typename T, Backend BACKEND, unsigned NUM_MODULI, bool CONJ>
__global__ void scaling_rowwise_full_vec4_kernel(
    const unsigned rows_A, const unsigned cols_A,
    const T *const __restrict__ A, const size_t lda,
    common::matptr_t<common::low_t<BACKEND>, common::isComplex<T>> A_lo,
    const size_t lda_lo, const size_t incA_lo,
    int16_t *const __restrict__ sftA //
) {
    using ValT               = decltype(trunc_scalbn<true, T, BACKEND, NUM_MODULI>::run(T{}, int32_t{}));
    constexpr bool cast_flag = sizeof(ValT) <= sizeof(T) && common::isCUDA;
    using shm_t              = std::conditional_t<cast_flag, ValT, T>;
    using Low4               = common::lowx4_t<BACKEND>;
    constexpr unsigned ty    = BACKEND == Backend::FP8 ? 8U : threads_y_rowwise<BACKEND>;
    static_assert(common::TILE_DIM == 32U && sizeof(common::low_t<BACKEND>) == 1U);
    static_assert(ty <= 8U && common::TILE_DIM % (4U * ty) == 0U);

    __shared__ shm_t tile[common::TILE_DIM][common::TILE_DIM + 1];
    const unsigned rowBase = blockIdx.x * common::TILE_DIM;
    const unsigned colBase = blockIdx.y * common::TILE_DIM;
    const unsigned in_row  = rowBase + threadIdx.x;
    const int32_t sft      = in_row < rows_A ? -sftA[in_row] : 0;

#pragma unroll
    for (unsigned j = 0; j < common::TILE_DIM; j += ty) {
        const unsigned yy     = threadIdx.y + j;
        const unsigned in_col = colBase + yy;
        const T a             = in_row < rows_A && in_col < cols_A ? A[in_col * lda + in_row] : common::Tconst<T>::zero();
        tile[yy][threadIdx.x] = trunc_scalbn<cast_flag, T, BACKEND, NUM_MODULI>::run(common::conj<T, CONJ>(a), sft);
    }
    __syncthreads();

    const unsigned xx      = threadIdx.x & 28U;
    const unsigned out_row = colBase + xx;
#pragma unroll 1
    for (unsigned j = 0; j < common::TILE_DIM; j += 4U * ty) {
        const unsigned yy      = 4U * threadIdx.y + (threadIdx.x & 3U) + j;
        const unsigned out_col = rowBase + yy;
        if (out_row >= cols_A || out_col >= rows_A) continue;

        const ValT v0 = trunc_scalbn<cast_flag, T, BACKEND, NUM_MODULI>::cast(tile[xx][yy]);
        const ValT v1 = trunc_scalbn<cast_flag, T, BACKEND, NUM_MODULI>::cast(tile[xx + 1U][yy]);
        const ValT v2 = trunc_scalbn<cast_flag, T, BACKEND, NUM_MODULI>::cast(tile[xx + 2U][yy]);
        const ValT v3 = trunc_scalbn<cast_flag, T, BACKEND, NUM_MODULI>::cast(tile[xx + 3U][yy]);

        const size_t idx        = out_col * (lda_lo >> 2) + (out_row >> 2);
        Low4 *__restrict__ out0 = reinterpret_cast<Low4 *>(A_lo.ptr0) + idx;
        if constexpr (common::isComplex<T>) {
            Low4 *__restrict__ out1 = reinterpret_cast<Low4 *>(A_lo.ptr1) + idx;
            mod::ModUnroll<NUM_MODULI, ValT>::run(out0, out1, incA_lo >> 2, v0, v1, v2, v3);
        } else {
            mod::ModUnroll<NUM_MODULI, ValT>::run(out0, incA_lo >> 2, v0, v1, v2, v3);
        }
    }
}

template <bool UPPER,
          typename T, Backend BACKEND, unsigned NUM_MODULI,
          cublasDiagType_t DIAG, bool CONJ>
__global__ void scaling_rowwise_tri_kernel(
    const unsigned rows_A, const unsigned cols_A,
    const T *const __restrict__ A, const size_t lda,
    common::matptr_t<common::low_t<BACKEND>, common::isComplex<T>> A_lo,
    const size_t lda_lo, const size_t incA_lo,
    int16_t *const __restrict__ sftA //
) {
    using ValT               = decltype(trunc_scalbn<true, T, BACKEND, NUM_MODULI>::run(T{}, int32_t{}));
    constexpr bool cast_flag = (sizeof(ValT) <= sizeof(T)) && common::isCUDA;
    using shm_t              = std::conditional_t<(cast_flag), ValT, T>;

    __shared__ shm_t tile[common::TILE_DIM][common::TILE_DIM + 1];

    const unsigned rowBase = blockIdx.x * common::TILE_DIM;
    const unsigned colBase = blockIdx.y * common::TILE_DIM;
    if (tri_tile_zero<UPPER>(rowBase, colBase)) return;
    const bool full_active = tri_tile_full_active<UPPER>(rowBase, colBase);

    const unsigned in_row = rowBase + threadIdx.x;
    const int32_t sft     = (in_row < rows_A) ? -sftA[in_row] : 0;

#pragma unroll
    for (unsigned j = 0; j < common::TILE_DIM; j += threads_y_rowwise<BACKEND>) {
        const unsigned yy = threadIdx.y + j;
        if (yy >= common::TILE_DIM) continue;

        const unsigned in_col = colBase + yy;
        T Atmp;
        if (full_active) {
            Atmp = (in_row < rows_A && in_col < cols_A) ? A[in_col * lda + in_row] : common::Tconst<T>::zero();
        } else {
            Atmp = tri_mat_value<UPPER, T, DIAG>(A, lda, in_row, in_col, rows_A, cols_A);
        }
        const shm_t A_scaled  = trunc_scalbn<cast_flag, T, BACKEND, NUM_MODULI>::run(common::conj<T, CONJ>(Atmp), sft);
        tile[yy][threadIdx.x] = A_scaled;
    }
    __syncthreads();

    constexpr unsigned mod_unroll = sizeof(common::underlying_t<T>) == sizeof(float)
                                        ? common::TILE_DIM / threads_y_rowwise<BACKEND>
                                        : 1U;
#pragma unroll mod_unroll
    for (unsigned j = 0; j < common::TILE_DIM; j += threads_y_rowwise<BACKEND>) {
        const unsigned yy = threadIdx.y + j;
        if (yy >= common::TILE_DIM) continue;

        const unsigned out_col = rowBase + yy;
        const unsigned out_row = colBase + threadIdx.x;
        if (out_col >= rows_A || out_row >= cols_A) continue;
        if (!full_active) {
            if (!tri_elem_active<UPPER>(out_col, out_row)) continue;
        }

        const size_t idx = out_col * lda_lo + out_row;
        const ValT in    = trunc_scalbn<cast_flag, T, BACKEND, NUM_MODULI>::cast(tile[threadIdx.x][yy]);

        if constexpr (common::isComplex<T>) {
            common::low_t<BACKEND> *__restrict__ out_1 = A_lo.ptr0 + idx;
            common::low_t<BACKEND> *__restrict__ out_2 = A_lo.ptr1 + idx;
            mod::ModUnroll<NUM_MODULI, ValT>::run(out_1, out_2, incA_lo, in);
        } else {
            common::low_t<BACKEND> *__restrict__ out = A_lo.ptr0 + idx;
            mod::ModUnroll<NUM_MODULI, ValT>::run(out, incA_lo, in);
        }
    }
}

template <typename T, Backend BACKEND, unsigned NUM_MODULI,
          cublasFillMode_t UPLO, cublasDiagType_t DIAG, bool CONJ>
void scaling_rowwise(
    const cudaStream_t stream,
    const unsigned rows_A, const unsigned cols_A,
    const T *const A, const size_t lda,
    common::matptr_t<common::low_t<BACKEND>, common::isComplex<T>> A_lo,
    const size_t lda_lo, const size_t incA_lo,
    int16_t *const sftA //
) {

    constexpr dim3 threads(threads_x_general, threads_y_rowwise<BACKEND>);
    dim3 grid((rows_A + threads_x_general - 1) / threads_x_general,
              (cols_A + common::TILE_DIM - 1) / common::TILE_DIM);

    if constexpr (UPLO == CUBLAS_FILL_MODE_FULL) {

        memset_padding_low_mats_2d_async<T, BACKEND, NUM_MODULI>(
            stream, A_lo, cols_A, lda_lo, incA_lo / lda_lo);

        if constexpr (common::isCUDA && (!common::isComplex<T> || NUM_MODULI <= common::threshold<BACKEND, true>::S)) {
            size_t alignment = lda_lo | incA_lo | reinterpret_cast<uintptr_t>(A_lo.ptr0);
            if constexpr (common::isComplex<T>) alignment |= reinterpret_cast<uintptr_t>(A_lo.ptr1);
            if ((alignment & 3U) == 0U) {
                constexpr unsigned ty = BACKEND == Backend::FP8 ? 8U : threads_y_rowwise<BACKEND>;
                constexpr dim3 threads_vec4(threads_x_general, ty);
                scaling_rowwise_full_vec4_kernel<T, BACKEND, NUM_MODULI, CONJ>
                    <<<grid, threads_vec4, 0, stream>>>(
                        rows_A, cols_A, A, lda, A_lo, lda_lo, incA_lo, sftA);
                return;
            }
        }

        scaling_rowwise_full_kernel<T, BACKEND, NUM_MODULI, CONJ>
            <<<grid, threads, 0, stream>>>(
                rows_A, cols_A, A, lda, A_lo, lda_lo, incA_lo, sftA);

    } else {

        constexpr bool isUPPER = UPLO == CUBLAS_FILL_MODE_UPPER;
        scaling_rowwise_tri_kernel<isUPPER, T, BACKEND, NUM_MODULI, DIAG, CONJ>
            <<<grid, threads, 0, stream>>>(
                rows_A, cols_A, A, lda, A_lo, lda_lo, incA_lo, sftA);
    }
}

} // namespace gemmul8::scaling::general
