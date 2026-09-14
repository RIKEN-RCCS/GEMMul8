#pragma once
#include "../../common/common.hpp"
#include "../../mod/mod.hpp"

#include "config.hpp"
#include "roundup.hpp"
#include "helper_triangular.hpp"
#include "helper_hermitian.hpp"
#include "store.hpp"

namespace gemmul8::scaling::general {

template <bool UPPER, bool HERM, bool STORE_TRANSPOSE,
          typename T, Backend BACKEND, unsigned NUM_MODULI>
__global__ void scaling_symm_hemm_offdiag_kernel(
    const unsigned n,
    const T *const __restrict__ A, const size_t lda,
    common::matptr_t<common::low_t<BACKEND>, common::isComplex<T>> A_lo,
    const size_t lda_lo, const size_t incA_lo,
    const int16_t *const __restrict__ sftA //
) {
    __shared__ T tile[common::TILE_DIM][common::TILE_DIM + 1];
    __shared__ int16_t sft_row[common::TILE_DIM];
    __shared__ int16_t sft_col[common::TILE_DIM];

    const unsigned tile_i = blockIdx.x;
    const unsigned tile_j = blockIdx.y;
    if constexpr (UPPER) {
        if (tile_i >= tile_j) return;
    } else {
        if (tile_i <= tile_j) return;
    }

    const unsigned rowBase = tile_i * common::TILE_DIM;
    const unsigned colBase = tile_j * common::TILE_DIM;
    if (threadIdx.y == 0) {
        const unsigned r = rowBase + threadIdx.x;
        const unsigned c = colBase + threadIdx.x;

        sft_row[threadIdx.x] = (r < n) ? -sftA[r] : 0;
        sft_col[threadIdx.x] = (c < n) ? -sftA[c] : 0;
    }

    const unsigned row = rowBase + threadIdx.x;

#pragma unroll(common::isComplex<T> ? 2 : 1)
    for (unsigned j = 0; j < common::TILE_DIM; j += threads_y_symm_hemm<T>) {
        const unsigned yy = threadIdx.y + j;
        if (yy >= common::TILE_DIM) continue;

        const unsigned col    = colBase + yy;
        const T a             = (row < n && col < n) ? A[col * lda + row] : common::Tconst<T>::zero();
        tile[yy][threadIdx.x] = a;
    }
    __syncthreads();

#pragma unroll(common::isComplex<T> ? 2 : 1)
    for (unsigned j = 0; j < common::TILE_DIM; j += threads_y_symm_hemm<T>) {
        const unsigned yy = threadIdx.y + j;
        if (yy >= common::TILE_DIM) continue;

        const unsigned col = colBase + yy;
        if (row < n && col < n) {
            const size_t idx = col * lda_lo + row;
            const T a        = general::hemm_active_store_value<HERM, STORE_TRANSPOSE, T>(tile[yy][threadIdx.x], row, col);
            scaling_store_one<T, BACKEND, NUM_MODULI>(A_lo, idx, incA_lo, a, int32_t(sft_col[yy]));
        }

        const unsigned mrow = colBase + threadIdx.x;
        const unsigned mcol = rowBase + yy;
        if (mrow < n && mcol < n) {
            const size_t midx = mcol * lda_lo + mrow;
            const T b         = general::hemm_mirror_store_value<HERM, STORE_TRANSPOSE, T>(tile[threadIdx.x][yy]);
            scaling_store_one<T, BACKEND, NUM_MODULI>(A_lo, midx, incA_lo, b, int32_t(sft_row[yy]));
        }
    }
}

template <bool UPPER, bool HERM, bool STORE_TRANSPOSE,
          typename T, Backend BACKEND, unsigned NUM_MODULI>
__global__ void scaling_symm_hemm_diag_kernel(
    const unsigned n,
    const T *const __restrict__ A, const size_t lda,
    common::matptr_t<common::low_t<BACKEND>, common::isComplex<T>> A_lo,
    const size_t lda_lo, const size_t incA_lo,
    const int16_t *const __restrict__ sftA //
) {
    __shared__ T tile[common::TILE_DIM][common::TILE_DIM + 1];
    __shared__ int16_t sft_tile[common::TILE_DIM];

    const unsigned tile_id = blockIdx.x;
    const unsigned base    = tile_id * common::TILE_DIM;

    if (threadIdx.y == 0) {
        const unsigned idx    = base + threadIdx.x;
        sft_tile[threadIdx.x] = (idx < n) ? -sftA[idx] : 0;
    }

#pragma unroll(common::isComplex<T> ? 2 : 1)
    for (unsigned j = 0; j < common::TILE_DIM; j += threads_y_symm_hemm<T>) {
        const unsigned yy = threadIdx.y + j;
        if (yy >= common::TILE_DIM) continue;

        const unsigned row  = base + threadIdx.x;
        const unsigned col  = base + yy;
        const bool in_range = (row < n && col < n);
        const bool active   = in_range && tri_elem_active<UPPER>(row, col);

        T a = common::Tconst<T>::zero();
        if (active) {
            a                     = A[col * lda + row];
            a                     = hemm_diag_value<HERM, T>(a, row, col);
            tile[yy][threadIdx.x] = a;
        }
    }
    __syncthreads();

#pragma unroll(common::isComplex<T> ? 2 : 1)
    for (unsigned j = 0; j < common::TILE_DIM; j += threads_y_symm_hemm<T>) {
        const unsigned yy = threadIdx.y + j;
        if (yy >= common::TILE_DIM) continue;

        const unsigned row = base + threadIdx.x;
        const unsigned col = base + yy;
        if (row >= n || col >= n) continue;

        const bool active = tri_elem_active<UPPER>(row, col);
        const size_t idx  = col * lda_lo + row;

        if (active) {
            const T a = general::hemm_active_store_value<HERM, STORE_TRANSPOSE, T>(tile[yy][threadIdx.x], row, col);
            scaling_store_one<T, BACKEND, NUM_MODULI>(A_lo, idx, incA_lo, a, int32_t(sft_tile[yy]));
        } else {
            const T b = general::hemm_mirror_store_value<HERM, STORE_TRANSPOSE, T>(tile[threadIdx.x][yy]);
            scaling_store_one<T, BACKEND, NUM_MODULI>(A_lo, idx, incA_lo, b, int32_t(sft_tile[yy]));
        }
    }
}

template <typename T, Backend BACKEND, unsigned NUM_MODULI, cublasFillMode_t UPLO>
void scaling_symm(
    const cudaStream_t stream,
    const unsigned n,
    const T *const A, const size_t lda,
    common::matptr_t<common::low_t<BACKEND>, common::isComplex<T>> A_lo,
    const size_t lda_lo, const size_t incA_lo,
    const int16_t *const sftA //
) {
    static_assert(UPLO == CUBLAS_FILL_MODE_UPPER || UPLO == CUBLAS_FILL_MODE_LOWER,
                  "scaling_symm requires UPLO = UPPER or LOWER.");

    constexpr dim3 threads(threads_x_general, threads_y_symm_hemm<T>);
    const unsigned nt      = (n + common::TILE_DIM - 1) / common::TILE_DIM;
    constexpr bool isUPPER = UPLO == CUBLAS_FILL_MODE_UPPER;

    scaling_symm_hemm_offdiag_kernel<isUPPER, false, false, T, BACKEND, NUM_MODULI>
        <<<dim3(nt, nt), threads, 0, stream>>>(
            n, A, lda, A_lo, lda_lo, incA_lo, sftA);

    scaling_symm_hemm_diag_kernel<isUPPER, false, false, T, BACKEND, NUM_MODULI>
        <<<nt, threads, 0, stream>>>(
            n, A, lda, A_lo, lda_lo, incA_lo, sftA);

    memset_padding_low_mats_2d_async<T, BACKEND, NUM_MODULI>(
        stream, A_lo, n, lda_lo, incA_lo / lda_lo);
}

template <typename T, Backend BACKEND, unsigned NUM_MODULI, cublasFillMode_t UPLO,
          bool STORE_TRANSPOSE>
void scaling_hemm(
    const cudaStream_t stream,
    const unsigned n,
    const T *const A, const size_t lda,
    common::matptr_t<common::low_t<BACKEND>, common::isComplex<T>> A_lo,
    const size_t lda_lo, const size_t incA_lo,
    const int16_t *const sftA //
) {
    static_assert(UPLO == CUBLAS_FILL_MODE_UPPER || UPLO == CUBLAS_FILL_MODE_LOWER,
                  "scaling_hemm requires UPLO = UPPER or LOWER.");

    constexpr dim3 threads(threads_x_general, threads_y_symm_hemm<T>);
    const unsigned nt      = (n + common::TILE_DIM - 1) / common::TILE_DIM;
    constexpr bool isUPPER = UPLO == CUBLAS_FILL_MODE_UPPER;

    scaling_symm_hemm_offdiag_kernel<isUPPER, true, STORE_TRANSPOSE, T, BACKEND, NUM_MODULI>
        <<<dim3(nt, nt), threads, 0, stream>>>(
            n, A, lda, A_lo, lda_lo, incA_lo, sftA);

    scaling_symm_hemm_diag_kernel<isUPPER, true, STORE_TRANSPOSE, T, BACKEND, NUM_MODULI>
        <<<nt, threads, 0, stream>>>(
            n, A, lda, A_lo, lda_lo, incA_lo, sftA);

    memset_padding_low_mats_2d_async<T, BACKEND, NUM_MODULI>(
        stream, A_lo, n, lda_lo, incA_lo / lda_lo);
}

} // namespace gemmul8::scaling::general
