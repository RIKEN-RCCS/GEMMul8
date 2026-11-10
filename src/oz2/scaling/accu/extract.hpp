#pragma once
#include "config.hpp"
#include "calc_sft.hpp"
#include "extract_rowwise.hpp"
#include "extract_colwise.hpp"

#include "../general/scaling_general_declaration.hpp"
#include "../general/helper_temp.hpp"

namespace gemmul8::scaling::accu {

template <typename T, Backend BACKEND, cublasFillMode_t UPLO, cublasDiagType_t DIAG>
void extract(
    const cudaStream_t stream,
    const cublasOperation_t op_A, const cublasSideMode_t side,
    const unsigned rows_A, const unsigned cols_A,
    const T *const A, const size_t lda,
    common::matptr_t<common::low_t<BACKEND>, common::isComplex<T>> &A_lo,
    const size_t lda_lo, const size_t incA_lo,
    int16_t *const sftA //
) {
    if (side == CUBLAS_SIDE_LEFT) {

        if (op_A == CUBLAS_OP_N) {

            // A: rows_A x cols_A -> A_lo: cols_A x rows_A
            using U               = common::underlying_t<T>;
            U *const partial_amax = general::temporary_memory<common::low_t<BACKEND>, U>(A_lo.ptr0);

            calc_sft_before_rowwise_launch<T, BACKEND, UPLO, DIAG>(
                stream, rows_A, cols_A, A, lda, sftA, partial_amax);

            if constexpr (UPLO == CUBLAS_FILL_MODE_FULL) {

                constexpr dim3 threads_extract_full(threads_x_extract_rowwise,
                                                    threads_y_extract_rowwise_full);
                dim3 grid((rows_A + common::TILE_DIM - 1) / common::TILE_DIM,
                          (lda_lo + common::TILE_DIM - 1) / common::TILE_DIM);

                common::launch_grid_y(grid, [&]<bool EXTENDED>(dim3 part) {
                    extract_rowwise_full_kernel<T, BACKEND, EXTENDED>
                        <<<part, threads_extract_full, 0, stream>>>(
                            rows_A, cols_A, A, lda, A_lo, lda_lo, sftA);
                });

            } else {

                general::memset_low_mats_async<T, BACKEND, 0U>(stream, A_lo, incA_lo);

                constexpr dim3 threads_extract_tri(threads_x_extract_rowwise,
                                                   threads_y_extract_rowwise_tri);
                dim3 grid((rows_A + common::TILE_DIM - 1) / common::TILE_DIM,
                          (cols_A + common::TILE_DIM - 1) / common::TILE_DIM);

                if constexpr (UPLO == CUBLAS_FILL_MODE_UPPER) {

                    common::launch_grid_y(grid, [&]<bool EXTENDED>(dim3 part) {
                        extract_rowwise_tri_kernel<true, T, BACKEND, DIAG, EXTENDED>
                            <<<part, threads_extract_tri, 0, stream>>>(
                                rows_A, cols_A, A, lda, A_lo, lda_lo, sftA);
                    });

                } else {

                    common::launch_grid_y(grid, [&]<bool EXTENDED>(dim3 part) {
                        extract_rowwise_tri_kernel<false, T, BACKEND, DIAG, EXTENDED>
                            <<<part, threads_extract_tri, 0, stream>>>(
                                rows_A, cols_A, A, lda, A_lo, lda_lo, sftA);
                    });
                }
            }

        } else {

            if constexpr (UPLO != CUBLAS_FILL_MODE_FULL) {
                general::memset_low_mats_async<T, BACKEND, 0U>(stream, A_lo, incA_lo);
            }

            // A: cols_A x rows_A -> A_lo: cols_A x rows_A
            extract_colwise<T, BACKEND, UPLO, DIAG>
                <<<rows_A, threads_accu, 0, stream>>>(
                    cols_A, A, lda, A_lo, lda_lo >> 2, sftA);
        }

    } else {

        if (op_A == CUBLAS_OP_N) {

            if constexpr (UPLO != CUBLAS_FILL_MODE_FULL) {
                general::memset_low_mats_async<T, BACKEND, 0U>(stream, A_lo, incA_lo);
            }

            // A: rows_A x cols_A -> A_lo: rows_A x cols_A
            extract_colwise<T, BACKEND, UPLO, DIAG>
                <<<cols_A, threads_accu, 0, stream>>>(
                    rows_A, A, lda, A_lo, lda_lo >> 2, sftA);

        } else {

            // A: cols_A x rows_A -> A_lo: rows_A x cols_A
            using U               = common::underlying_t<T>;
            U *const partial_amax = general::temporary_memory<common::low_t<BACKEND>, U>(A_lo.ptr0);

            calc_sft_before_rowwise_launch<T, BACKEND, UPLO, DIAG>(
                stream, cols_A, rows_A, A, lda, sftA, partial_amax);

            if constexpr (UPLO == CUBLAS_FILL_MODE_FULL) {

                constexpr dim3 threads_extract_full(threads_x_extract_rowwise,
                                                    threads_y_extract_rowwise_full);
                dim3 grid((cols_A + common::TILE_DIM - 1) / common::TILE_DIM,
                          (lda_lo + common::TILE_DIM - 1) / common::TILE_DIM);

                common::launch_grid_y(grid, [&]<bool EXTENDED>(dim3 part) {
                    extract_rowwise_full_kernel<T, BACKEND, EXTENDED>
                        <<<part, threads_extract_full, 0, stream>>>(
                            cols_A, rows_A, A, lda, A_lo, lda_lo, sftA);
                });

            } else {

                general::memset_low_mats_async<T, BACKEND, 0U>(stream, A_lo, incA_lo);

                constexpr dim3 threads_extract_tri(threads_x_extract_rowwise,
                                                   threads_y_extract_rowwise_tri);
                dim3 grid((cols_A + common::TILE_DIM - 1) / common::TILE_DIM,
                          (rows_A + common::TILE_DIM - 1) / common::TILE_DIM);

                if constexpr (UPLO == CUBLAS_FILL_MODE_UPPER) {

                    common::launch_grid_y(grid, [&]<bool EXTENDED>(dim3 part) {
                        extract_rowwise_tri_kernel<true, T, BACKEND, DIAG, EXTENDED>
                            <<<part, threads_extract_tri, 0, stream>>>(
                                cols_A, rows_A, A, lda, A_lo, lda_lo, sftA);
                    });

                } else {

                    common::launch_grid_y(grid, [&]<bool EXTENDED>(dim3 part) {
                        extract_rowwise_tri_kernel<false, T, BACKEND, DIAG, EXTENDED>
                            <<<part, threads_extract_tri, 0, stream>>>(
                                cols_A, rows_A, A, lda, A_lo, lda_lo, sftA);
                    });
                }
            }
        }
    }
}

} // namespace gemmul8::scaling::accu
