#pragma once
#include "../../common/launch_grid.hpp"
#include "config.hpp"
#include "store.hpp"

#include "../general/helper_triangular.hpp"
#include "../general/roundup.hpp"

namespace gemmul8::scaling::accu {

template <typename T, Backend BACKEND, bool EXTENDED = false>
__global__ void extract_rowwise_full_kernel(
    const unsigned rows_A, const unsigned cols_A,
    const T *const __restrict__ A, const size_t lda,
    common::matptr_t<common::low_t<BACKEND>, common::isComplex<T>> A_lo,
    const size_t lda_lo,
    int16_t *const __restrict__ sftA //
) {
    using ValT = upperBound_t<T, BACKEND>;
    __shared__ ValT tile[common::TILE_DIM][common::TILE_DIM + 1];

    const unsigned rowBase = blockIdx.x * common::TILE_DIM;
    const unsigned colBase = common::block_y<EXTENDED>() * common::TILE_DIM;

    const unsigned in_row = rowBase + threadIdx.x;
    const int32_t sft     = (in_row < rows_A) ? sftA[in_row] : 0;

#pragma unroll
    for (unsigned j = 0; j < common::TILE_DIM; j += threads_y_extract_rowwise_full) {
        const unsigned yy = threadIdx.y + j;
        if (yy >= common::TILE_DIM) continue;

        const unsigned in_col = colBase + yy;
        const T Atmp          = (in_row < rows_A && in_col < cols_A) ? A[in_col * lda + in_row] : common::Tconst<T>::zero();
        const ValT A_upper    = upperBound_lo<T, BACKEND>(Atmp, sft);
        tile[yy][threadIdx.x] = A_upper;
    }
    __syncthreads();

#pragma unroll
    for (unsigned j = 0; j < common::TILE_DIM; j += threads_y_extract_rowwise_full) {
        const unsigned yy = threadIdx.y + j;
        if (yy >= common::TILE_DIM) continue;

        const unsigned out_col = rowBase + yy;
        const unsigned out_row = colBase + threadIdx.x;
        if (out_col >= rows_A || out_row >= lda_lo) continue;

        const ValT out   = tile[threadIdx.x][yy];
        const size_t idx = out_col * lda_lo + out_row;
        if constexpr (common::isComplex<T>) {
            A_lo.ptr0[idx] = out.x;
            A_lo.ptr1[idx] = out.y;
            A_lo.ptr2[idx] = sub_ru_8bit(out.x, out.y);
        } else {
            A_lo.ptr0[idx] = out;
        }
    }
}

template <bool UPPER, typename T, Backend BACKEND, cublasDiagType_t DIAG, bool EXTENDED = false>
__global__ void extract_rowwise_tri_kernel(
    const unsigned rows_A, const unsigned cols_A,
    const T *const __restrict__ A, const size_t lda,
    common::matptr_t<common::low_t<BACKEND>, common::isComplex<T>> A_lo,
    const size_t lda_lo,
    int16_t *const __restrict__ sftA //
) {
    using ValT = upperBound_t<T, BACKEND>;
    __shared__ ValT tile[common::TILE_DIM][common::TILE_DIM + 1];

    const unsigned rowBase = blockIdx.x * common::TILE_DIM;
    const unsigned colBase = common::block_y<EXTENDED>() * common::TILE_DIM;
    if (general::tri_tile_zero<UPPER>(rowBase, colBase)) return;

    const bool full_active = general::tri_tile_full_active<UPPER>(rowBase, colBase);
    const unsigned in_row  = rowBase + threadIdx.x;
    const int32_t sft      = (in_row < rows_A) ? sftA[in_row] : 0;

#pragma unroll
    for (unsigned j = 0; j < common::TILE_DIM; j += threads_y_extract_rowwise_tri) {
        const unsigned yy = threadIdx.y + j;
        if (yy >= common::TILE_DIM) continue;

        const unsigned in_col = colBase + yy;
        T Atmp;
        if (full_active) {
            Atmp = (in_row < rows_A && in_col < cols_A) ? A[in_col * lda + in_row] : common::Tconst<T>::zero();
        } else {
            Atmp = general::tri_mat_value<UPPER, T, DIAG>(A, lda, in_row, in_col, rows_A, cols_A);
        }
        const ValT A_upper    = upperBound_lo<T, BACKEND>(Atmp, sft);
        tile[yy][threadIdx.x] = A_upper;
    }
    __syncthreads();

#pragma unroll
    for (unsigned j = 0; j < common::TILE_DIM; j += threads_y_extract_rowwise_tri) {
        const unsigned yy = threadIdx.y + j;
        if (yy >= common::TILE_DIM) continue;

        const unsigned out_col = rowBase + yy;          // original row
        const unsigned out_row = colBase + threadIdx.x; // original col
        if (out_col >= rows_A || out_row >= cols_A) continue;
        if (!full_active) {
            if (!general::tri_elem_active<UPPER>(out_col, out_row)) continue;
        }

        const ValT out   = tile[threadIdx.x][yy];
        const size_t idx = out_col * lda_lo + out_row;
        if constexpr (common::isComplex<T>) {
            A_lo.ptr0[idx] = out.x;
            A_lo.ptr1[idx] = out.y;
            A_lo.ptr2[idx] = sub_ru_8bit(out.x, out.y);
        } else {
            A_lo.ptr0[idx] = out;
        }
    }
}

} // namespace gemmul8::scaling::accu
