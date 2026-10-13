#pragma once

#include "../general/scaling_general_declaration.hpp"

namespace gemmul8::scaling::accu {

template <typename T, Backend BACKEND, unsigned NUM_MODULI,
          cublasFillMode_t UPLO, cublasDiagType_t DIAG, bool CONJ>
__device__ __forceinline__ void scaling_colwise_body(
    const unsigned rows_A,
    const T *const __restrict__ A, const size_t lda,
    common::matptr_t<common::low_t<BACKEND>, common::isComplex<T>> A_lo,
    const size_t lda_lo, const size_t incA_lo,
    int16_t *const __restrict__ sftA //
) {
    const unsigned col_idx         = blockIdx.x;
    const T *const __restrict__ in = A + col_idx * lda;
    const int32_t sft              = -sftA[col_idx];

#if defined(__CUDA_ARCH__) && __CUDA_ARCH__ >= 1000 && __CUDA_ARCH__ < 1200
    if constexpr (BACKEND == Backend::FP8 &&
                  std::is_same_v<T, float> &&
                  UPLO == CUBLAS_FILL_MODE_FULL &&
                  (NUM_MODULI >= 8U && NUM_MODULI <= 11U)) {
        general::scaling_colwise_full_device<T, BACKEND, NUM_MODULI, CONJ, true>(
            rows_A, in, A_lo, lda_lo, incA_lo, sft);
    } else
#endif
        if constexpr (BACKEND == Backend::FP8 && common::isComplex<T> &&
                      general::fp8_complex_column_registers<T, NUM_MODULI> != 0U && UPLO == CUBLAS_FILL_MODE_FULL) {
        general::scaling_colwise_full_device<T, BACKEND, NUM_MODULI, CONJ, false, true>(
            rows_A, in, A_lo, lda_lo, incA_lo, sft);
    } else general::scaling_colwise_device<T, BACKEND, NUM_MODULI, UPLO, DIAG, CONJ>(
        rows_A, in, A_lo, lda_lo, incA_lo, sft);
}

GEMMUL8_FP8_KERNEL_EX(scaling_colwise,
                      (typename T, Backend BACKEND, unsigned NUM_MODULI, cublasFillMode_t UPLO, cublasDiagType_t DIAG, bool CONJ),
                      (T, BACKEND, NUM_MODULI, UPLO, DIAG, CONJ),
                      (BACKEND == Backend::FP8 && ((!common::isComplex<T> && NUM_MODULI >= 6U) || general::fp8_complex_column_registers<T, NUM_MODULI> != 0U) && UPLO == CUBLAS_FILL_MODE_FULL),
                      (common::isComplex<T> ? general::fp8_complex_column_registers<T, NUM_MODULI> : general::fp8_column_registers<T, NUM_MODULI>),
                      (BACKEND == Backend::FP8 && UPLO == CUBLAS_FILL_MODE_FULL ? general::fp8_refined_registers<T, NUM_MODULI, general::fp8_kernel_kind::accurate> : 0U))

template <typename T, Backend BACKEND, unsigned NUM_MODULI,
          cublasFillMode_t UPLO, cublasDiagType_t DIAG, bool CONJ>
inline void scaling_colwise_launch(
    const cudaStream_t stream,
    const unsigned rows_A,
    const unsigned cols_A,
    const T *const A,
    const size_t lda,
    common::matptr_t<common::low_t<BACKEND>, common::isComplex<T>> A_lo,
    const size_t lda_lo4,
    const size_t incA_lo4,
    int16_t *const sftA //
) {
    constexpr dim3 threads(general::threads_x_colwise_full_tiled,
                           general::threads_y_colwise_full_tiled);

    if constexpr (UPLO == CUBLAS_FILL_MODE_FULL) {

        constexpr unsigned max_columns = 65535U * general::threads_y_colwise_full_tiled;
        if (cols_A > max_columns) {
            for (unsigned first = 0; first < cols_A;) {
                const unsigned count = std::min(max_columns, cols_A - first);
                auto part_lo         = A_lo;
                part_lo.shift(first * lda_lo4 * 4U);
                scaling_colwise_launch<T, BACKEND, NUM_MODULI, UPLO, DIAG, CONJ>(
                    stream, rows_A, count, A + first * lda, lda,
                    part_lo, lda_lo4, incA_lo4, sftA + first);
                first += count;
            }
            return;
        }

        bool vector_load_aligned;
        if constexpr (common::isComplex<T>) {
            vector_load_aligned = false;
        } else {
            vector_load_aligned = common::is_aligned_as<T>(A) &&
                                  ((lda * sizeof(T)) % common::vec4_alignment_v<T> == 0ULL);
        }

        const bool aligned_full =
            ((rows_A & 3U) == 0U) &&
            (lda_lo4 == (size_t(rows_A) >> 2)) &&
            ((lda_lo4 % general::threads_x_colwise_full_tiled) == 0U) &&
            vector_load_aligned;

        if (aligned_full) {
            const dim3 grid(
                unsigned(lda_lo4 / general::threads_x_colwise_full_tiled),
                (cols_A + general::threads_y_colwise_full_tiled - 1U) / general::threads_y_colwise_full_tiled);

            general::scaling_colwise_full_tiled_aligned_kernel<T, BACKEND, NUM_MODULI, CONJ>
                <<<grid, threads, 0, stream>>>(
                    rows_A, cols_A, A, lda, A_lo, lda_lo4, incA_lo4, sftA);

        } else {
            const dim3 grid(
                (unsigned(lda_lo4) + general::threads_x_colwise_full_tiled - 1U) / general::threads_x_colwise_full_tiled,
                (cols_A + general::threads_y_colwise_full_tiled - 1U) / general::threads_y_colwise_full_tiled);

            general::scaling_colwise_full_tiled_kernel<T, BACKEND, NUM_MODULI, CONJ>
                <<<grid, threads, 0, stream>>>(
                    rows_A, cols_A, A, lda, A_lo, lda_lo4, incA_lo4, sftA);
        }

    } else {

        const unsigned rows4 = (rows_A + 3U) >> 2;

        unsigned cols_grid = cols_A;
        if constexpr (UPLO == CUBLAS_FILL_MODE_LOWER) {
            cols_grid = (cols_A < rows_A) ? cols_A : rows_A;
        }

        if (rows4 == 0U || cols_grid == 0U) return;

        const dim3 grid(
            (rows4 + general::threads_x_colwise_full_tiled - 1U) / general::threads_x_colwise_full_tiled,
            (cols_grid + general::threads_y_colwise_full_tiled - 1U) / general::threads_y_colwise_full_tiled);

        general::scaling_colwise_tri_tiled_kernel<T, BACKEND, NUM_MODULI, UPLO, DIAG, CONJ>
            <<<grid, threads, 0, stream>>>(
                rows_A, cols_A, A, lda, A_lo, lda_lo4, incA_lo4, sftA);
    }
}

} // namespace gemmul8::scaling::accu
