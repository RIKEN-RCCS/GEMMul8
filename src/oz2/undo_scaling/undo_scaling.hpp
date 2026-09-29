#pragma once

#include "config.hpp"
#include "rescale.hpp"
#include "undo_scaling_declaration.hpp"
#include "final_reduction.hpp"
#include "predicates.hpp"
#include "scalar.hpp"

namespace gemmul8::undo_scaling {

template <Backend BACKEND, unsigned NUM_MODULI, cublasFillMode_t UPLO, bool HERMITIAN>
__device__ __forceinline__ bool prepare_crt_tile(
    const unsigned m,
    const unsigned n,
    const size_t ldc_mid,
    const size_t incC_mid,
    const crt_tail<BACKEND> tail,
    int32_t (*opposite)[HERMITIAN ? 33 : 1] //
) {
    if constexpr (HERMITIAN) {
        if constexpr (UPLO == CUBLAS_FILL_MODE_UPPER) {
            if (blockIdx.x > blockIdx.y) return false;
        } else if constexpr (UPLO == CUBLAS_FILL_MODE_LOWER) {
            if (blockIdx.x < blockIdx.y) return false;
        }

        constexpr unsigned IDX = common::table::active_index<BACKEND, NUM_MODULI, NUM_MODULI - 1, true>;

        if (tail.ptr0 != nullptr && tail.ptr1 == nullptr) {
#pragma unroll
            for (unsigned j = 0; j < 32; j += 8) {
                const unsigned tr = blockIdx.y * 32 + threadIdx.x;
                const unsigned tc = blockIdx.x * 32 + threadIdx.y + j;

                opposite[threadIdx.y + j][threadIdx.x] =
                    (tr < n && tc < m)
                        ? raw_product_residue<BACKEND, IDX, true>(tail.ptr0, tc * ldc_mid + tr, incC_mid)
                        : 0;
            }
            __syncthreads();
        }
    }
    return true;
}

//------------------------------
// General kernel
//------------------------------
template <typename T, Backend BACKEND, unsigned NUM_MODULI,
          cublasFillMode_t UPLO, bool isTRTRMM,
          typename TAlpha, typename TBeta,
          bool HERMITIAN = false, bool GROUPED = false>
__global__ void undo_scaling_kernel(
    const TAlpha alpha, const TBeta beta,
    const unsigned m, const unsigned n,
    const common::mid_t<BACKEND, common::isComplex<T>> *const __restrict__ C_mid,
    const size_t ldc_mid, const size_t incC_mid,
    T *const __restrict__ C, const size_t ldc,
    const crt_tail<BACKEND> tail,
    const int16_t *const __restrict__ sftA,
    const int16_t *const __restrict__ sftB //
) {
    using TCrt = std::conditional_t<common::isComplex<T>, cuDoubleComplex, double>;

    constexpr bool NONGROUPED_HERM  = HERMITIAN && !GROUPED;
    constexpr unsigned tile_columns = NONGROUPED_HERM ? 32U : threads_y;
    constexpr unsigned step_columns = NONGROUPED_HERM ? 8U : threads_y;
    const unsigned row              = blockIdx.x * threads_x + threadIdx.x;
    __shared__ int32_t opposite[NONGROUPED_HERM ? 32 : 1][NONGROUPED_HERM ? 33 : 1];

    if (!prepare_crt_tile<BACKEND, NUM_MODULI, UPLO, NONGROUPED_HERM>(m, n, ldc_mid, incC_mid, tail, opposite)) return;

#pragma unroll 1
    for (unsigned j = 0; j < tile_columns; j += step_columns) {
        const unsigned col    = blockIdx.y * tile_columns + threadIdx.y + j;
        int32_t tail_opposite = 0;
        if constexpr (NONGROUPED_HERM) {
            if (tail.ptr0 != nullptr && tail.ptr1 == nullptr) {
                tail_opposite = opposite[threadIdx.x][threadIdx.y + j];
            }
        }
        if (row >= m || col >= n) continue;

        if constexpr (UPLO == CUBLAS_FILL_MODE_UPPER && !isTRTRMM) {
            if (row > col) continue;
        } else if constexpr (UPLO == CUBLAS_FILL_MODE_LOWER && !isTRTRMM) {
            if (row < col) continue;
        }

        if constexpr (UPLO == CUBLAS_FILL_MODE_UPPER && isTRTRMM) {
            if (row <= col) {
                const size_t idx_C_mid = col * ldc_mid + row;

                const TCrt AB_crt = reconstruct_from_crt < TCrt, BACKEND, NUM_MODULI, GROUPED, GROUPED && !HERMITIAN > (C_mid + idx_C_mid, incC_mid, tail, idx_C_mid, tail_opposite);

                const int sft        = int(sftA[row]) + int(sftB[col]);
                const TCrt AB_scaled = rescale_crt<TCrt>(AB_crt, sft);
                const T AB           = common::Tcast<TCrt, T>(AB_scaled);

                const size_t idx_C = col * ldc + row;
                const auto alpha_v = alpha.get();
                if constexpr (std::is_same_v<TBeta, void *>) {
                    C[idx_C] = common::Tmul<scalar_t<TAlpha>, T>(alpha_v, AB);
                } else {
                    const auto beta_v = beta.get();
                    C[idx_C]          = axpby<T>(alpha_v, AB, beta_v, C + idx_C);
                }
            } else {
                const size_t idx_C = col * ldc + row;
                if constexpr (std::is_same_v<TBeta, void *>) {
                    C[idx_C] = common::Tconst<T>::zero();
                } else {
                    const auto beta_v = beta.get();
                    C[idx_C]          = scale_or_zero<T>(beta_v, C + idx_C);
                }
            }
        } else if constexpr (UPLO == CUBLAS_FILL_MODE_LOWER && isTRTRMM) {
            if (row >= col) {
                const size_t idx_C_mid = col * ldc_mid + row;

                const TCrt AB_crt = reconstruct_from_crt < TCrt, BACKEND, NUM_MODULI, GROUPED, GROUPED && !HERMITIAN > (C_mid + idx_C_mid, incC_mid, tail, idx_C_mid, tail_opposite);

                const int sft        = int(sftA[row]) + int(sftB[col]);
                const TCrt AB_scaled = rescale_crt<TCrt>(AB_crt, sft);
                const T AB           = common::Tcast<TCrt, T>(AB_scaled);

                const size_t idx_C = col * ldc + row;
                const auto alpha_v = alpha.get();
                if constexpr (std::is_same_v<TBeta, void *>) {
                    C[idx_C] = common::Tmul<scalar_t<TAlpha>, T>(alpha_v, AB);
                } else {
                    const auto beta_v = beta.get();
                    C[idx_C]          = axpby<T>(alpha_v, AB, beta_v, C + idx_C);
                }
            } else {
                const size_t idx_C = col * ldc + row;
                if constexpr (std::is_same_v<TBeta, void *>) {
                    C[idx_C] = common::Tconst<T>::zero();
                } else {
                    const auto beta_v = beta.get();
                    C[idx_C]          = scale_or_zero<T>(beta_v, C + idx_C);
                }
            }
        } else {
            const size_t idx_C_mid = col * ldc_mid + row;

            const TCrt AB_crt = reconstruct_from_crt < TCrt, BACKEND, NUM_MODULI, GROUPED, GROUPED && !HERMITIAN > (C_mid + idx_C_mid, incC_mid, tail, idx_C_mid, tail_opposite);

            const int sft        = int(sftA[row]) + int(sftB[col]);
            const TCrt AB_scaled = rescale_crt<TCrt>(AB_crt, sft);
            const T AB           = common::Tcast<TCrt, T>(AB_scaled);

            const size_t idx_C = col * ldc + row;
            const auto alpha_v = alpha.get();
            if constexpr (HERMITIAN) {
                static_assert(common::isComplex<T> && !isTRTRMM);
                if (row == col) {
                    T value;
                    if constexpr (std::is_same_v<TBeta, void *>) {
                        value = common::Tmul<scalar_t<TAlpha>, T>(alpha_v, AB);
                    } else {
                        const auto beta_v = beta.get();
                        if (is_zero(beta_v)) {
                            value = common::Tmul<scalar_t<TAlpha>, T>(alpha_v, AB);
                        } else {
                            const T old{C[idx_C].x, common::underlying_t<T>(0)};
                            value = common::Taxpby<T, scalar_t<TAlpha>, scalar_t<TBeta>>(alpha_v, AB, beta_v, old);
                        }
                    }
                    value.y  = common::underlying_t<T>(0);
                    C[idx_C] = value;
                    continue;
                }
            }
            if constexpr (std::is_same_v<TBeta, void *>) {
                C[idx_C] = common::Tmul<scalar_t<TAlpha>, T>(alpha_v, AB);
            } else {
                const auto beta_v = beta.get();
                C[idx_C]          = axpby<T>(alpha_v, AB, beta_v, C + idx_C);
            }
        }
    }
}

//------------------------------
// Special kernel for alpha in {1, -1}, beta in {-1, 0, 1}
//------------------------------
template <typename T, Backend BACKEND, unsigned NUM_MODULI,
          cublasFillMode_t UPLO, bool isTRTRMM,
          int ALPHA, int BETA,
          bool HERMITIAN = false, bool GROUPED = false>
__global__ void undo_scaling_kernel_special(
    const unsigned m, const unsigned n,
    const common::mid_t<BACKEND, common::isComplex<T>> *const __restrict__ C_mid,
    const size_t ldc_mid, const size_t incC_mid,
    T *const __restrict__ C, const size_t ldc,
    const crt_tail<BACKEND> tail,
    const int16_t *const __restrict__ sftA,
    const int16_t *const __restrict__ sftB //
) {
    using TCrt = std::conditional_t<common::isComplex<T>, cuDoubleComplex, double>;

    constexpr bool NONGROUPED_HERM  = HERMITIAN && !GROUPED;
    constexpr unsigned tile_columns = NONGROUPED_HERM ? 32U : threads_y;
    constexpr unsigned step_columns = NONGROUPED_HERM ? 8U : threads_y;
    const unsigned row              = blockIdx.x * threads_x + threadIdx.x;
    __shared__ int32_t opposite[NONGROUPED_HERM ? 32 : 1][NONGROUPED_HERM ? 33 : 1];

    if (!prepare_crt_tile<BACKEND, NUM_MODULI, UPLO, NONGROUPED_HERM>(m, n, ldc_mid, incC_mid, tail, opposite)) return;

#pragma unroll 1
    for (unsigned j = 0; j < tile_columns; j += step_columns) {
        const unsigned col    = blockIdx.y * tile_columns + threadIdx.y + j;
        int32_t tail_opposite = 0;
        if constexpr (NONGROUPED_HERM) {
            if (tail.ptr0 != nullptr && tail.ptr1 == nullptr) {
                tail_opposite = opposite[threadIdx.x][threadIdx.y + j];
            }
        }
        if (row >= m || col >= n) continue;

        if constexpr (UPLO == CUBLAS_FILL_MODE_UPPER && !isTRTRMM) {
            if (row > col) continue;
        } else if constexpr (UPLO == CUBLAS_FILL_MODE_LOWER && !isTRTRMM) {
            if (row < col) continue;
        }

        if constexpr (UPLO == CUBLAS_FILL_MODE_UPPER && isTRTRMM) {
            if (row <= col) {
                const size_t idx_C_mid = col * ldc_mid + row;

                const TCrt AB_crt = reconstruct_from_crt < TCrt, BACKEND, NUM_MODULI, GROUPED, GROUPED && !HERMITIAN > (C_mid + idx_C_mid, incC_mid, tail, idx_C_mid, tail_opposite);

                const int sft        = int(sftA[row]) + int(sftB[col]);
                const TCrt AB_scaled = rescale_crt<TCrt>(AB_crt, sft);
                const T AB           = common::Tcast<TCrt, T>(AB_scaled);

                const size_t idx_C = col * ldc + row;
                C[idx_C]           = Taxpby_special<T, ALPHA, BETA>(AB, C + idx_C);
            } else {
                const size_t idx_C = col * ldc + row;
                C[idx_C]           = Tmul_special<T, BETA>(C + idx_C);
            }
        } else if constexpr (UPLO == CUBLAS_FILL_MODE_LOWER && isTRTRMM) {
            if (row >= col) {
                const size_t idx_C_mid = col * ldc_mid + row;

                const TCrt AB_crt = reconstruct_from_crt < TCrt, BACKEND, NUM_MODULI, GROUPED, GROUPED && !HERMITIAN > (C_mid + idx_C_mid, incC_mid, tail, idx_C_mid, tail_opposite);

                const int sft        = int(sftA[row]) + int(sftB[col]);
                const TCrt AB_scaled = rescale_crt<TCrt>(AB_crt, sft);
                const T AB           = common::Tcast<TCrt, T>(AB_scaled);

                const size_t idx_C = col * ldc + row;
                C[idx_C]           = Taxpby_special<T, ALPHA, BETA>(AB, C + idx_C);
            } else {
                const size_t idx_C = col * ldc + row;
                C[idx_C]           = Tmul_special<T, BETA>(C + idx_C);
            }
        } else {
            const size_t idx_C_mid = col * ldc_mid + row;

            const TCrt AB_crt = reconstruct_from_crt < TCrt, BACKEND, NUM_MODULI, GROUPED, GROUPED && !HERMITIAN > (C_mid + idx_C_mid, incC_mid, tail, idx_C_mid, tail_opposite);

            const int sft        = int(sftA[row]) + int(sftB[col]);
            const TCrt AB_scaled = rescale_crt<TCrt>(AB_crt, sft);
            const T AB           = common::Tcast<TCrt, T>(AB_scaled);

            const size_t idx_C = col * ldc + row;
            if constexpr (HERMITIAN) {
                static_assert(common::isComplex<T> && !isTRTRMM);
                if (row == col) {
                    using U       = common::underlying_t<T>;
                    const U value = Taxpby_special<U, ALPHA, BETA>(AB.x, &C[idx_C].x);
                    C[idx_C]      = T{value, U(0)};
                    continue;
                }
            }
            C[idx_C] = Taxpby_special<T, ALPHA, BETA>(AB, C + idx_C);
        }
    }
}

//------------------------------
// Launcher
//------------------------------
template <typename T, typename TAlpha, typename TBeta,
          Backend BACKEND, unsigned NUM_MODULI,
          cublasFillMode_t UPLO, bool isTRTRMM, bool HERMITIAN, bool GROUPED>
void undo_scaling_impl(
    const cudaStream_t stream,
    const unsigned m, const unsigned n,
    common::mid_t<BACKEND, common::isComplex<T>> *C_mid,
    const size_t ldc_mid, const size_t incC_mid,
    T *const C, const size_t ldc,
    const int16_t *const sftA, const int16_t *const sftB,
    const TAlpha *const alpha, const TBeta *const beta,
    const crt_tail<BACKEND> tail //
) {
    constexpr bool NONGROUPED_HERM  = HERMITIAN && !GROUPED;
    constexpr unsigned tile_columns = NONGROUPED_HERM ? 32U : threads_y;
    constexpr dim3 threads(threads_x, NONGROUPED_HERM ? 8U : threads_y);
    const dim3 grid((m + threads_x - 1) / threads_x,
                    (n + tile_columns - 1) / tile_columns);

    if (beta == nullptr) {
        const bool alpha_dev = is_device_pointer(alpha);

        if (alpha_dev) {
            using alpha_t   = DeviceScalar<TAlpha>;
            alpha_t alpha_d = alpha_t(alpha);

            undo_scaling_kernel<T, BACKEND, NUM_MODULI, UPLO, isTRTRMM, alpha_t, void *, HERMITIAN, GROUPED>
                <<<grid, threads, 0, stream>>>(
                    alpha_d, nullptr, m, n, C_mid, ldc_mid, incC_mid, C, ldc, tail, sftA, sftB);
            return;
        }

        const TAlpha alpha_v = *alpha;

        if (is_one_h(alpha_v)) {
            undo_scaling_kernel_special<T, BACKEND, NUM_MODULI, UPLO, isTRTRMM, 1, 0, HERMITIAN, GROUPED>
                <<<grid, threads, 0, stream>>>(
                    m, n, C_mid, ldc_mid, incC_mid, C, ldc, tail, sftA, sftB);
            return;
        }
        if (is_mone_h(alpha_v)) {
            undo_scaling_kernel_special<T, BACKEND, NUM_MODULI, UPLO, isTRTRMM, -1, 0, HERMITIAN, GROUPED>
                <<<grid, threads, 0, stream>>>(
                    m, n, C_mid, ldc_mid, incC_mid, C, ldc, tail, sftA, sftB);
            return;
        }

        using alpha_t   = HostScalar<TAlpha>;
        alpha_t alpha_h = alpha_t(alpha_v);

        undo_scaling_kernel<T, BACKEND, NUM_MODULI, UPLO, isTRTRMM, alpha_t, void *, HERMITIAN, GROUPED>
            <<<grid, threads, 0, stream>>>(
                alpha_h, nullptr, m, n, C_mid, ldc_mid, incC_mid, C, ldc, tail, sftA, sftB);

    } else {

        const bool alpha_dev = is_device_pointer(alpha);
        const bool beta_dev  = is_device_pointer(beta);

        if (alpha_dev || beta_dev) {
            if (!(alpha_dev && beta_dev)) {
                assert(false && "alpha and beta must both be host pointers or both be device-accessible pointers");
                return;
            }

            using alpha_t = DeviceScalar<TAlpha>;
            using beta_t  = DeviceScalar<TBeta>;

            alpha_t alpha_d = alpha_t(alpha);
            beta_t beta_d   = beta_t(beta);

            undo_scaling_kernel<T, BACKEND, NUM_MODULI, UPLO, isTRTRMM, alpha_t, beta_t, HERMITIAN, GROUPED>
                <<<grid, threads, 0, stream>>>(
                    alpha_d, beta_d, m, n, C_mid, ldc_mid, incC_mid, C, ldc, tail, sftA, sftB);
            return;
        }

        const TAlpha alpha_v = *alpha;
        const TBeta beta_v   = *beta;

        if (is_one_h(alpha_v)) {
            if (is_zero_h(beta_v)) {
                undo_scaling_kernel_special<T, BACKEND, NUM_MODULI, UPLO, isTRTRMM, 1, 0, HERMITIAN, GROUPED>
                    <<<grid, threads, 0, stream>>>(
                        m, n, C_mid, ldc_mid, incC_mid, C, ldc, tail, sftA, sftB);
                return;
            }

            if (is_one_h(beta_v)) {
                undo_scaling_kernel_special<T, BACKEND, NUM_MODULI, UPLO, isTRTRMM, 1, 1, HERMITIAN, GROUPED>
                    <<<grid, threads, 0, stream>>>(
                        m, n, C_mid, ldc_mid, incC_mid, C, ldc, tail, sftA, sftB);
                return;
            }

            if (is_mone_h(beta_v)) {
                undo_scaling_kernel_special<T, BACKEND, NUM_MODULI, UPLO, isTRTRMM, 1, -1, HERMITIAN, GROUPED>
                    <<<grid, threads, 0, stream>>>(
                        m, n, C_mid, ldc_mid, incC_mid, C, ldc, tail, sftA, sftB);
                return;
            }
        }

        if (is_mone_h(alpha_v)) {
            if (is_zero_h(beta_v)) {
                undo_scaling_kernel_special<T, BACKEND, NUM_MODULI, UPLO, isTRTRMM, -1, 0, HERMITIAN, GROUPED>
                    <<<grid, threads, 0, stream>>>(
                        m, n, C_mid, ldc_mid, incC_mid, C, ldc, tail, sftA, sftB);
                return;
            }

            if (is_one_h(beta_v)) {
                undo_scaling_kernel_special<T, BACKEND, NUM_MODULI, UPLO, isTRTRMM, -1, 1, HERMITIAN, GROUPED>
                    <<<grid, threads, 0, stream>>>(
                        m, n, C_mid, ldc_mid, incC_mid, C, ldc, tail, sftA, sftB);
                return;
            }

            if (is_mone_h(beta_v)) {
                undo_scaling_kernel_special<T, BACKEND, NUM_MODULI, UPLO, isTRTRMM, -1, -1, HERMITIAN, GROUPED>
                    <<<grid, threads, 0, stream>>>(
                        m, n, C_mid, ldc_mid, incC_mid, C, ldc, tail, sftA, sftB);
                return;
            }
        }

        if (is_zero_h(beta_v)) {
            using alpha_t = HostScalar<TAlpha>;
            undo_scaling_kernel<T, BACKEND, NUM_MODULI, UPLO, isTRTRMM, alpha_t, void *, HERMITIAN, GROUPED>
                <<<grid, threads, 0, stream>>>(
                    alpha_t(alpha_v), nullptr, m, n, C_mid, ldc_mid, incC_mid,
                    C, ldc, tail, sftA, sftB);
            return;
        }

        using alpha_t = HostScalar<TAlpha>;
        using beta_t  = HostScalar<TBeta>;

        alpha_t alpha_h = alpha_t(alpha_v);
        beta_t beta_h   = beta_t(beta_v);

        undo_scaling_kernel<T, BACKEND, NUM_MODULI, UPLO, isTRTRMM, alpha_t, beta_t, HERMITIAN, GROUPED>
            <<<grid, threads, 0, stream>>>(
                alpha_h, beta_h, m, n, C_mid, ldc_mid, incC_mid, C, ldc, tail, sftA, sftB);
    }
}

template <typename T, typename TAlpha, typename TBeta,
          Backend BACKEND, unsigned NUM_MODULI,
          cublasFillMode_t UPLO, bool isTRTRMM, bool HERMITIAN>
void undo_scaling(
    cudaStream_t stream,
    unsigned m, unsigned n,
    common::mid_t<BACKEND, common::isComplex<T>> *C_mid,
    size_t ldc_mid, size_t incC_mid, T *C, size_t ldc,
    const int16_t *sftA, const int16_t *sftB,
    const TAlpha *alpha, const TBeta *beta,
    crt_tail<BACKEND> tail //
) {
    if (tail.grouped) {
        undo_scaling_impl<T, TAlpha, TBeta, BACKEND, NUM_MODULI, UPLO, isTRTRMM, HERMITIAN, true>(
            stream, m, n, C_mid, ldc_mid, incC_mid, C, ldc, sftA, sftB, alpha, beta, tail);
    } else {
        undo_scaling_impl<T, TAlpha, TBeta, BACKEND, NUM_MODULI, UPLO, isTRTRMM, HERMITIAN, false>(
            stream, m, n, C_mid, ldc_mid, incC_mid, C, ldc, sftA, sftB, alpha, beta, tail);
    }
}

} // namespace gemmul8::undo_scaling
