#pragma once

#include "config.hpp"
#include "rescale.hpp"
#include "undo_scaling_declaration.hpp"
#include "final_reduction.hpp"
#include "predicates.hpp"
#include "scalar.hpp"

namespace gemmul8::undo_scaling {

//------------------------------
// General kernel
//------------------------------
template <typename T, Backend BACKEND, unsigned NUM_MODULI, typename TP,
          cublasFillMode_t UPLO, bool isTRTRMM, typename TAlpha, typename TBeta, bool HERMITIAN = false>
__global__ void undo_scaling_kernel(
    const TAlpha alpha, const TBeta beta,
    const unsigned m, const unsigned n,
    const common::mid_t<BACKEND, common::isComplex<T>> *const __restrict__ C_mid,
    const size_t ldc_mid, const size_t incC_mid,
    T *const __restrict__ C, const size_t ldc,
    const TP P, const double invP,
    const int16_t *const __restrict__ sftA,
    const int16_t *const __restrict__ sftB //
) {
    using TCrt = std::conditional_t<common::isComplex<T>, cuDoubleComplex, double>;

    const unsigned row = blockIdx.x * blockDim.x + threadIdx.x;
    const unsigned col = blockIdx.y * blockDim.y + threadIdx.y;
    if (row >= m || col >= n) return;

    if constexpr (UPLO == CUBLAS_FILL_MODE_UPPER && !isTRTRMM) {
        if (row > col) return;
    } else if constexpr (UPLO == CUBLAS_FILL_MODE_LOWER && !isTRTRMM) {
        if (row < col) return;
    }

    if constexpr (UPLO == CUBLAS_FILL_MODE_UPPER && isTRTRMM) {
        if (row <= col) {
            const size_t idx_C_mid = col * ldc_mid + row;

            const TCrt AB_crt = reconstruct_from_crt<TCrt, BACKEND, NUM_MODULI, TP>(
                C_mid + idx_C_mid, incC_mid, P, invP);

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

            const TCrt AB_crt = reconstruct_from_crt<TCrt, BACKEND, NUM_MODULI, TP>(
                C_mid + idx_C_mid, incC_mid, P, invP);

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

        const TCrt AB_crt = reconstruct_from_crt<TCrt, BACKEND, NUM_MODULI, TP>(
            C_mid + idx_C_mid, incC_mid, P, invP);

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
                return;
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

//------------------------------
// Special kernel for alpha in {1, -1}, beta in {-1, 0, 1}
//------------------------------
template <typename T, Backend BACKEND, unsigned NUM_MODULI, typename TP,
          cublasFillMode_t UPLO, bool isTRTRMM, int ALPHA, int BETA, bool HERMITIAN = false>
__global__ void undo_scaling_kernel_special(
    const unsigned m, const unsigned n,
    const common::mid_t<BACKEND, common::isComplex<T>> *const __restrict__ C_mid,
    const size_t ldc_mid, const size_t incC_mid,
    T *const __restrict__ C, const size_t ldc,
    const TP P, const double invP,
    const int16_t *const __restrict__ sftA,
    const int16_t *const __restrict__ sftB //
) {
    using TCrt = std::conditional_t<common::isComplex<T>, cuDoubleComplex, double>;

    const unsigned row = blockIdx.x * blockDim.x + threadIdx.x;
    const unsigned col = blockIdx.y * blockDim.y + threadIdx.y;
    if (row >= m || col >= n) return;

    if constexpr (UPLO == CUBLAS_FILL_MODE_UPPER && !isTRTRMM) {
        if (row > col) return;
    } else if constexpr (UPLO == CUBLAS_FILL_MODE_LOWER && !isTRTRMM) {
        if (row < col) return;
    }

    if constexpr (UPLO == CUBLAS_FILL_MODE_UPPER && isTRTRMM) {
        if (row <= col) {
            const size_t idx_C_mid = col * ldc_mid + row;

            const TCrt AB_crt = reconstruct_from_crt<TCrt, BACKEND, NUM_MODULI, TP>(
                C_mid + idx_C_mid, incC_mid, P, invP);

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

            const TCrt AB_crt = reconstruct_from_crt<TCrt, BACKEND, NUM_MODULI, TP>(
                C_mid + idx_C_mid, incC_mid, P, invP);

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

        const TCrt AB_crt = reconstruct_from_crt<TCrt, BACKEND, NUM_MODULI, TP>(
            C_mid + idx_C_mid, incC_mid, P, invP);

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
                return;
            }
        }
        C[idx_C] = Taxpby_special<T, ALPHA, BETA>(AB, C + idx_C);
    }
}

//------------------------------
// Launcher
//------------------------------
template <typename T, typename TAlpha, typename TBeta,
          Backend BACKEND, unsigned NUM_MODULI, cublasFillMode_t UPLO, bool isTRTRMM, bool HERMITIAN>
void undo_scaling(
    const cudaStream_t stream,
    const unsigned m, const unsigned n,
    common::mid_t<BACKEND, common::isComplex<T>> *C_mid,
    const size_t ldc_mid, const size_t incC_mid,
    T *const C, const size_t ldc,
    const int16_t *const sftA, const int16_t *const sftB,
    const TAlpha *const alpha, const TBeta *const beta //
) {
    constexpr dim3 threads(threads_x, threads_y);
    const dim3 grid((m + threads_x - 1) / threads_x,
                    (n + threads_y - 1) / threads_y);

    constexpr bool is_float = std::is_same_v<common::underlying_t<T>, float>;
    constexpr bool small_NM = NUM_MODULI <= common::threshold<BACKEND, common::isComplex<T>>::P_is_double;
    using TP                = std::conditional_t<is_float || small_NM, double, double2>;
    const TP P              = common::table::get_P<BACKEND, TP, common::isComplex<T>>(NUM_MODULI);
    const double invP       = common::table::get_invP<BACKEND, common::isComplex<T>>(NUM_MODULI);

    if (beta == nullptr) {

        const bool alpha_dev = is_device_pointer(alpha);

        if (alpha_dev) {
            using alpha_t   = DeviceScalar<TAlpha>;
            alpha_t alpha_d = alpha_t(alpha);

            undo_scaling_kernel<T, BACKEND, NUM_MODULI, TP, UPLO, isTRTRMM, alpha_t, void *, HERMITIAN>
                <<<grid, threads, 0, stream>>>(
                    alpha_d, nullptr, m, n, C_mid, ldc_mid, incC_mid, C, ldc, P, invP, sftA, sftB);
            return;
        }

        const TAlpha alpha_v = *alpha;

        if (is_one_h(alpha_v)) {
            undo_scaling_kernel_special<T, BACKEND, NUM_MODULI, TP, UPLO, isTRTRMM, 1, 0, HERMITIAN>
                <<<grid, threads, 0, stream>>>(
                    m, n, C_mid, ldc_mid, incC_mid, C, ldc, P, invP, sftA, sftB);
            return;
        }
        if (is_mone_h(alpha_v)) {
            undo_scaling_kernel_special<T, BACKEND, NUM_MODULI, TP, UPLO, isTRTRMM, -1, 0, HERMITIAN>
                <<<grid, threads, 0, stream>>>(
                    m, n, C_mid, ldc_mid, incC_mid, C, ldc, P, invP, sftA, sftB);
            return;
        }

        using alpha_t   = HostScalar<TAlpha>;
        alpha_t alpha_h = alpha_t(alpha_v);

        undo_scaling_kernel<T, BACKEND, NUM_MODULI, TP, UPLO, isTRTRMM, alpha_t, void *, HERMITIAN>
            <<<grid, threads, 0, stream>>>(
                alpha_h, nullptr, m, n, C_mid, ldc_mid, incC_mid, C, ldc, P, invP, sftA, sftB);

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

            undo_scaling_kernel<T, BACKEND, NUM_MODULI, TP, UPLO, isTRTRMM, alpha_t, beta_t, HERMITIAN>
                <<<grid, threads, 0, stream>>>(
                    alpha_d, beta_d, m, n, C_mid, ldc_mid, incC_mid, C, ldc, P, invP, sftA, sftB);
            return;
        }

        const TAlpha alpha_v = *alpha;
        const TBeta beta_v   = *beta;

        if (is_one_h(alpha_v)) {
            if (is_zero_h(beta_v)) {
                undo_scaling_kernel_special<T, BACKEND, NUM_MODULI, TP, UPLO, isTRTRMM, 1, 0, HERMITIAN>
                    <<<grid, threads, 0, stream>>>(
                        m, n, C_mid, ldc_mid, incC_mid, C, ldc, P, invP, sftA, sftB);
                return;
            }

            if (is_one_h(beta_v)) {
                undo_scaling_kernel_special<T, BACKEND, NUM_MODULI, TP, UPLO, isTRTRMM, 1, 1, HERMITIAN>
                    <<<grid, threads, 0, stream>>>(
                        m, n, C_mid, ldc_mid, incC_mid, C, ldc, P, invP, sftA, sftB);
                return;
            }

            if (is_mone_h(beta_v)) {
                undo_scaling_kernel_special<T, BACKEND, NUM_MODULI, TP, UPLO, isTRTRMM, 1, -1, HERMITIAN>
                    <<<grid, threads, 0, stream>>>(
                        m, n, C_mid, ldc_mid, incC_mid, C, ldc, P, invP, sftA, sftB);
                return;
            }
        }

        if (is_mone_h(alpha_v)) {
            if (is_zero_h(beta_v)) {
                undo_scaling_kernel_special<T, BACKEND, NUM_MODULI, TP, UPLO, isTRTRMM, -1, 0, HERMITIAN>
                    <<<grid, threads, 0, stream>>>(
                        m, n, C_mid, ldc_mid, incC_mid, C, ldc, P, invP, sftA, sftB);
                return;
            }

            if (is_one_h(beta_v)) {
                undo_scaling_kernel_special<T, BACKEND, NUM_MODULI, TP, UPLO, isTRTRMM, -1, 1, HERMITIAN>
                    <<<grid, threads, 0, stream>>>(
                        m, n, C_mid, ldc_mid, incC_mid, C, ldc, P, invP, sftA, sftB);
                return;
            }

            if (is_mone_h(beta_v)) {
                undo_scaling_kernel_special<T, BACKEND, NUM_MODULI, TP, UPLO, isTRTRMM, -1, -1, HERMITIAN>
                    <<<grid, threads, 0, stream>>>(
                        m, n, C_mid, ldc_mid, incC_mid, C, ldc, P, invP, sftA, sftB);
                return;
            }
        }

        if (is_zero_h(beta_v)) {
            using alpha_t = HostScalar<TAlpha>;
            undo_scaling_kernel<T, BACKEND, NUM_MODULI, TP, UPLO, isTRTRMM, alpha_t, void *, HERMITIAN>
                <<<grid, threads, 0, stream>>>(
                    alpha_t(alpha_v), nullptr, m, n, C_mid, ldc_mid, incC_mid,
                    C, ldc, P, invP, sftA, sftB);
            return;
        }

        using alpha_t = HostScalar<TAlpha>;
        using beta_t  = HostScalar<TBeta>;

        alpha_t alpha_h = alpha_t(alpha_v);
        beta_t beta_h   = beta_t(beta_v);

        undo_scaling_kernel<T, BACKEND, NUM_MODULI, TP, UPLO, isTRTRMM, alpha_t, beta_t, HERMITIAN>
            <<<grid, threads, 0, stream>>>(
                alpha_h, beta_h, m, n, C_mid, ldc_mid, incC_mid, C, ldc, P, invP, sftA, sftB);
    }
}

} // namespace gemmul8::undo_scaling
