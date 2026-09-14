#pragma once

#include "config.hpp"
#include "rescale.hpp"
#include "final_reduction.hpp"
#include "predicates.hpp"
#include "scalar.hpp"

namespace gemmul8::undo_scaling {

// alpha * (Dij + Dji)
template <typename T, typename TAlpha>
__device__ __forceinline__ T syr2k_update(const TAlpha alpha, const T Dij, const T Dji) {
    const T sum = common::Tadd<T>(Dij, Dji);
    return common::Tmul<TAlpha, T>(alpha, sum);
}

// alpha * Dij + conj(alpha) * conj(Dji)
template <typename T, typename TAlpha>
__device__ __forceinline__ T her2k_update(const TAlpha alpha, const T Dij, const T Dji) {
    if constexpr (common::isComplex<TAlpha>) {
        return common::Taxpby<T, TAlpha, TAlpha>(
            alpha, Dij, common::conj<TAlpha, true>(alpha), common::conj<T, true>(Dji));
    } else {
        return syr2k_update<T, TAlpha>(alpha, Dij, common::conj<T, true>(Dji));
    }
}

template <bool HER2K, typename T, typename TAlpha>
__device__ __forceinline__ T r2k_update(const TAlpha alpha, const T Dij, const T Dji) {
    if constexpr (HER2K) {
        return her2k_update<T, TAlpha>(alpha, Dij, Dji);
    } else {
        return syr2k_update<T, TAlpha>(alpha, Dij, Dji);
    }
}

// C_new = update + beta*C_old
template <bool HER2K, typename T, typename TBeta>
__device__ __forceinline__ T r2k_add_beta(const T update, const TBeta beta, const T *C_old, const bool diag) {
    T C_new = update;
    if (!is_zero(beta)) C_new = common::Taxpy<T, TBeta>(beta, *C_old, update);
    if constexpr (HER2K) {
        using U = common::underlying_t<T>;
        C_new.y = (diag) ? common::Tconst<U>::zero() : C_new.y;
    }
    return C_new;
}

// c_new = alpha*D + beta*c_old
template <bool HER2K, typename T, int ALPHA, int BETA>
__device__ __forceinline__ T r2k_update_special_offdiag(const T Dij, const T Dji, const T *C_old) {
    const T D = common::Tadd<T>(Dij, common::conj<T, HER2K>(Dji));
    return Taxpby_special<T, ALPHA, BETA>(D, C_old);
}

template <bool HER2K, typename T, int ALPHA, int BETA>
__device__ __forceinline__ T r2k_update_special(const T Dij, const T Dji, const T *C_old, const bool diag) {
    T C_new = r2k_update_special_offdiag<HER2K, T, ALPHA, BETA>(Dij, Dji, C_old);
    if constexpr (HER2K) {
        using U = common::underlying_t<T>;
        C_new.y = (diag) ? common::Tconst<U>::zero() : C_new.y;
    }
    return C_new;
}

template <bool HER2K, typename T, Backend BACKEND, unsigned NUM_MODULI, typename TP,
          cublasFillMode_t UPLO, typename TAlpha, typename TBeta>
__global__ void undo_scaling_r2k_offdiag_kernel(
    const TAlpha alpha, const TBeta beta,
    const unsigned n,
    const common::mid_t<BACKEND, common::isComplex<T>> *const __restrict__ C_mid,
    const size_t ldc_mid, const size_t incC_mid,
    T *const __restrict__ C, const size_t ldc,
    const TP P, const double invP,
    const int16_t *const __restrict__ sftA,
    const int16_t *const __restrict__ sftB //
) {
    static_assert(!HER2K || common::isComplex<T>, "HER2K requires complex T.");

    using TCrt = std::conditional_t<common::isComplex<T>, cuDoubleComplex, double>;

    __shared__ T tile_Dji[threads_x_r2k][threads_x_r2k + 1];

    const unsigned tile_i = blockIdx.x;
    const unsigned tile_j = blockIdx.y;
    if constexpr (UPLO == CUBLAS_FILL_MODE_UPPER) {
        if (tile_i >= tile_j) return;
    } else {
        if (tile_i <= tile_j) return;
    }

    const unsigned rowBase = tile_i * threads_x_r2k;
    const unsigned colBase = tile_j * threads_x_r2k;

#pragma unroll 1
    for (unsigned j = 0; j < threads_x_r2k; j += threads_y_r2k) {
        const unsigned row = colBase + threadIdx.x;
        const unsigned col = rowBase + threadIdx.y + j;

        T Dji = common::Tconst<T>::zero();

        if (row < n && col < n) {
            const size_t idx_ji = col * ldc_mid + row;

            const TCrt Dji_crt = reconstruct_from_crt<TCrt, BACKEND, NUM_MODULI, TP>(
                C_mid + idx_ji, incC_mid, P, invP);

            const int32_t sft_ji  = int32_t(sftA[row]) + int32_t(sftB[col]);
            const TCrt Dji_scaled = rescale_crt<TCrt>(Dji_crt, sft_ji);
            Dji                   = common::Tcast<TCrt, T>(Dji_scaled);
        }

        tile_Dji[threadIdx.x][threadIdx.y + j] = Dji;
    }
    __syncthreads();

#pragma unroll 1
    for (unsigned j = 0; j < threads_x_r2k; j += threads_y_r2k) {
        const unsigned row = rowBase + threadIdx.x;
        const unsigned col = colBase + threadIdx.y + j;

        if (row < n && col < n) {
            const size_t idx_ij = col * ldc_mid + row;

            const TCrt Dij_crt = reconstruct_from_crt<TCrt, BACKEND, NUM_MODULI, TP>(
                C_mid + idx_ij, incC_mid, P, invP);

            const int32_t sft_ij  = int32_t(sftA[row]) + int32_t(sftB[col]);
            const TCrt Dij_scaled = rescale_crt<TCrt>(Dij_crt, sft_ij);
            const T Dij           = common::Tcast<TCrt, T>(Dij_scaled);
            const T Dji           = tile_Dji[threadIdx.y + j][threadIdx.x];

            const auto alpha_v = alpha.get();
            const T update     = r2k_update<HER2K, T, scalar_t<TAlpha>>(alpha_v, Dij, Dji);

            const size_t idx_C = col * ldc + row;
            const auto beta_v  = beta.get();
            C[idx_C]           = r2k_add_beta<false, T, scalar_t<TBeta>>(update, beta_v, C + idx_C, false);
        }
    }
}

template <bool HER2K, typename T, Backend BACKEND, unsigned NUM_MODULI, typename TP,
          cublasFillMode_t UPLO, int ALPHA, int BETA>
__global__ void undo_scaling_r2k_offdiag_kernel_special(
    const unsigned n,
    const common::mid_t<BACKEND, common::isComplex<T>> *const __restrict__ C_mid,
    const size_t ldc_mid, const size_t incC_mid,
    T *const __restrict__ C, const size_t ldc,
    const TP P, const double invP,
    const int16_t *const __restrict__ sftA,
    const int16_t *const __restrict__ sftB //
) {
    static_assert(!HER2K || common::isComplex<T>, "HER2K requires complex T.");

    using TCrt = std::conditional_t<common::isComplex<T>, cuDoubleComplex, double>;

    __shared__ T tile_Dji[threads_x_r2k][threads_x_r2k + 1];

    const unsigned tile_i = blockIdx.x;
    const unsigned tile_j = blockIdx.y;
    if constexpr (UPLO == CUBLAS_FILL_MODE_UPPER) {
        if (tile_i >= tile_j) return;
    } else {
        if (tile_i <= tile_j) return;
    }

    const unsigned rowBase = tile_i * threads_x_r2k;
    const unsigned colBase = tile_j * threads_x_r2k;

#pragma unroll 1
    for (unsigned j = 0; j < threads_x_r2k; j += threads_y_r2k) {
        const unsigned row = colBase + threadIdx.x;
        const unsigned col = rowBase + threadIdx.y + j;

        T Dji = common::Tconst<T>::zero();

        if (row < n && col < n) {
            const size_t idx_ji = col * ldc_mid + row;

            const TCrt Dji_crt = reconstruct_from_crt<TCrt, BACKEND, NUM_MODULI, TP>(
                C_mid + idx_ji, incC_mid, P, invP);

            const int32_t sft_ji  = int32_t(sftA[row]) + int32_t(sftB[col]);
            const TCrt Dji_scaled = rescale_crt<TCrt>(Dji_crt, sft_ji);
            Dji                   = common::Tcast<TCrt, T>(Dji_scaled);
        }

        tile_Dji[threadIdx.x][threadIdx.y + j] = Dji;
    }
    __syncthreads();

#pragma unroll 1
    for (unsigned j = 0; j < threads_x_r2k; j += threads_y_r2k) {
        const unsigned row = rowBase + threadIdx.x;
        const unsigned col = colBase + threadIdx.y + j;

        if (row < n && col < n) {
            const size_t idx_ij = col * ldc_mid + row;

            const TCrt Dij_crt = reconstruct_from_crt<TCrt, BACKEND, NUM_MODULI, TP>(
                C_mid + idx_ij, incC_mid, P, invP);

            const int32_t sft_ij  = int32_t(sftA[row]) + int32_t(sftB[col]);
            const TCrt Dij_scaled = rescale_crt<TCrt>(Dij_crt, sft_ij);
            const T Dij           = common::Tcast<TCrt, T>(Dij_scaled);
            const T Dji           = tile_Dji[threadIdx.y + j][threadIdx.x];

            const size_t idx_C = col * ldc + row;
            C[idx_C]           = r2k_update_special_offdiag<HER2K, T, ALPHA, BETA>(Dij, Dji, C + idx_C);
        }
    }
}

template <bool HER2K, typename T, Backend BACKEND, unsigned NUM_MODULI, typename TP,
          cublasFillMode_t UPLO, typename TAlpha, typename TBeta>
__global__ void undo_scaling_r2k_diag_kernel(
    const TAlpha alpha, const TBeta beta,
    const unsigned n,
    const common::mid_t<BACKEND, common::isComplex<T>> *const __restrict__ C_mid,
    const size_t ldc_mid, const size_t incC_mid,
    T *const __restrict__ C, const size_t ldc,
    const TP P, const double invP,
    const int16_t *const __restrict__ sftA,
    const int16_t *const __restrict__ sftB //
) {
    static_assert(!HER2K || common::isComplex<T>, "HER2K requires complex T.");

    using TCrt = std::conditional_t<common::isComplex<T>, cuDoubleComplex, double>;

    __shared__ T tile_D[threads_x_r2k][threads_x_r2k + 1];

    const unsigned base = blockIdx.x * threads_x_r2k;

#pragma unroll 1
    for (unsigned j = 0; j < threads_x_r2k; j += threads_y_r2k) {
        const unsigned row = base + threadIdx.x;
        const unsigned col = base + threadIdx.y + j;

        T Dij = common::Tconst<T>::zero();

        if (row < n && col < n) {
            const size_t idx = col * ldc_mid + row;

            const TCrt Dij_crt = reconstruct_from_crt<TCrt, BACKEND, NUM_MODULI, TP>(
                C_mid + idx, incC_mid, P, invP);

            const int32_t sft_ij  = int32_t(sftA[row]) + int32_t(sftB[col]);
            const TCrt Dij_scaled = rescale_crt<TCrt>(Dij_crt, sft_ij);
            Dij                   = common::Tcast<TCrt, T>(Dij_scaled);
        }

        tile_D[threadIdx.x][threadIdx.y + j] = Dij;
    }
    __syncthreads();

#pragma unroll 1
    for (unsigned j = 0; j < threads_x_r2k; j += threads_y_r2k) {
        const unsigned row = base + threadIdx.x;
        const unsigned col = base + threadIdx.y + j;

        if (row >= n || col >= n) continue;
        if constexpr (UPLO == CUBLAS_FILL_MODE_UPPER) {
            if (row > col) continue;
        } else {
            if (row < col) continue;
        }

        const T Dij        = tile_D[threadIdx.x][threadIdx.y + j];
        const T Dji        = tile_D[threadIdx.y + j][threadIdx.x];
        const auto alpha_v = alpha.get();
        const T update     = r2k_update<HER2K, T, scalar_t<TAlpha>>(alpha_v, Dij, Dji);

        const size_t idx_C = col * ldc + row;
        const auto beta_v  = beta.get();
        C[idx_C]           = r2k_add_beta<HER2K, T, scalar_t<TBeta>>(update, beta_v, C + idx_C, row == col);
    }
}

template <bool HER2K, typename T, Backend BACKEND, unsigned NUM_MODULI, typename TP,
          cublasFillMode_t UPLO, int ALPHA, int BETA>
__global__ void undo_scaling_r2k_diag_kernel_special(
    const unsigned n,
    const common::mid_t<BACKEND, common::isComplex<T>> *const __restrict__ C_mid,
    const size_t ldc_mid, const size_t incC_mid,
    T *const __restrict__ C, const size_t ldc,
    const TP P, const double invP,
    const int16_t *const __restrict__ sftA,
    const int16_t *const __restrict__ sftB //
) {
    static_assert(!HER2K || common::isComplex<T>, "HER2K requires complex T.");

    using TCrt = std::conditional_t<common::isComplex<T>, cuDoubleComplex, double>;

    __shared__ T tile_D[threads_x_r2k][threads_x_r2k + 1];

    const unsigned base = blockIdx.x * threads_x_r2k;

#pragma unroll 1
    for (unsigned j = 0; j < threads_x_r2k; j += threads_y_r2k) {
        const unsigned row = base + threadIdx.x;
        const unsigned col = base + threadIdx.y + j;

        T Dij = common::Tconst<T>::zero();

        if (row < n && col < n) {
            const size_t idx = col * ldc_mid + row;

            const TCrt Dij_crt = reconstruct_from_crt<TCrt, BACKEND, NUM_MODULI, TP>(
                C_mid + idx, incC_mid, P, invP);

            const int32_t sft_ij  = int32_t(sftA[row]) + int32_t(sftB[col]);
            const TCrt Dij_scaled = rescale_crt<TCrt>(Dij_crt, sft_ij);
            Dij                   = common::Tcast<TCrt, T>(Dij_scaled);
        }

        tile_D[threadIdx.x][threadIdx.y + j] = Dij;
    }
    __syncthreads();

#pragma unroll 1
    for (unsigned j = 0; j < threads_x_r2k; j += threads_y_r2k) {
        const unsigned row = base + threadIdx.x;
        const unsigned col = base + threadIdx.y + j;

        if (row >= n || col >= n) continue;
        if constexpr (UPLO == CUBLAS_FILL_MODE_UPPER) {
            if (row > col) continue;
        } else {
            if (row < col) continue;
        }

        const T Dij = tile_D[threadIdx.x][threadIdx.y + j];
        const T Dji = tile_D[threadIdx.y + j][threadIdx.x];

        const size_t idx_C = col * ldc + row;
        C[idx_C]           = r2k_update_special<HER2K, T, ALPHA, BETA>(Dij, Dji, C + idx_C, row == col);
    }
}

template <bool HER2K, typename T, typename TAlpha, typename TBeta,
          Backend BACKEND, unsigned NUM_MODULI, cublasFillMode_t UPLO>
inline void undo_scaling_r2k_launch(
    const cudaStream_t stream,
    const unsigned n,
    common::mid_t<BACKEND, common::isComplex<T>> *C_mid,
    const size_t ldc_mid, const size_t incC_mid,
    T *const C, const size_t ldc,
    const int16_t *const sftA, const int16_t *const sftB,
    const TAlpha *const alpha, const TBeta *const beta //
) {
    constexpr dim3 threads(threads_x_r2k, threads_y_r2k);
    const unsigned nt = (n + threads_x_r2k - 1) / threads_x_r2k;
    const dim3 grid_offdiag(nt, nt);
    const dim3 grid_diag(nt);

    constexpr bool is_float = std::is_same_v<common::underlying_t<T>, float>;
    constexpr bool small_NM = NUM_MODULI <= common::threshold<BACKEND, common::isComplex<T>>::P_is_double;
    using TP                = std::conditional_t<is_float || small_NM, double, double2>;
    const TP P              = common::table::get_P<BACKEND, TP, common::isComplex<T>>(NUM_MODULI);
    const double invP       = common::table::get_invP<BACKEND, common::isComplex<T>>(NUM_MODULI);

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

        undo_scaling_r2k_offdiag_kernel<HER2K, T, BACKEND, NUM_MODULI, TP, UPLO, alpha_t, beta_t>
            <<<grid_offdiag, threads, 0, stream>>>(
                alpha_d, beta_d, n, C_mid, ldc_mid, incC_mid, C, ldc, P, invP, sftA, sftB);

        undo_scaling_r2k_diag_kernel<HER2K, T, BACKEND, NUM_MODULI, TP, UPLO, alpha_t, beta_t>
            <<<grid_diag, threads, 0, stream>>>(
                alpha_d, beta_d, n, C_mid, ldc_mid, incC_mid, C, ldc, P, invP, sftA, sftB);
        return;
    }

    const TAlpha alpha_v = *alpha;
    const TBeta beta_v   = *beta;

    if (is_one_h(alpha_v)) {
        if (is_zero_h(beta_v)) {
            undo_scaling_r2k_offdiag_kernel_special<HER2K, T, BACKEND, NUM_MODULI, TP, UPLO, 1, 0>
                <<<grid_offdiag, threads, 0, stream>>>(
                    n, C_mid, ldc_mid, incC_mid, C, ldc, P, invP, sftA, sftB);

            undo_scaling_r2k_diag_kernel_special<HER2K, T, BACKEND, NUM_MODULI, TP, UPLO, 1, 0>
                <<<grid_diag, threads, 0, stream>>>(
                    n, C_mid, ldc_mid, incC_mid, C, ldc, P, invP, sftA, sftB);
            return;
        }

        if (is_one_h(beta_v)) {
            undo_scaling_r2k_offdiag_kernel_special<HER2K, T, BACKEND, NUM_MODULI, TP, UPLO, 1, 1>
                <<<grid_offdiag, threads, 0, stream>>>(
                    n, C_mid, ldc_mid, incC_mid, C, ldc, P, invP, sftA, sftB);

            undo_scaling_r2k_diag_kernel_special<HER2K, T, BACKEND, NUM_MODULI, TP, UPLO, 1, 1>
                <<<grid_diag, threads, 0, stream>>>(
                    n, C_mid, ldc_mid, incC_mid, C, ldc, P, invP, sftA, sftB);
            return;
        }

        if (is_mone_h(beta_v)) {
            undo_scaling_r2k_offdiag_kernel_special<HER2K, T, BACKEND, NUM_MODULI, TP, UPLO, 1, -1>
                <<<grid_offdiag, threads, 0, stream>>>(
                    n, C_mid, ldc_mid, incC_mid, C, ldc, P, invP, sftA, sftB);

            undo_scaling_r2k_diag_kernel_special<HER2K, T, BACKEND, NUM_MODULI, TP, UPLO, 1, -1>
                <<<grid_diag, threads, 0, stream>>>(
                    n, C_mid, ldc_mid, incC_mid, C, ldc, P, invP, sftA, sftB);
            return;
        }
    }

    if (is_mone_h(alpha_v)) {
        if (is_zero_h(beta_v)) {
            undo_scaling_r2k_offdiag_kernel_special<HER2K, T, BACKEND, NUM_MODULI, TP, UPLO, -1, 0>
                <<<grid_offdiag, threads, 0, stream>>>(
                    n, C_mid, ldc_mid, incC_mid, C, ldc, P, invP, sftA, sftB);

            undo_scaling_r2k_diag_kernel_special<HER2K, T, BACKEND, NUM_MODULI, TP, UPLO, -1, 0>
                <<<grid_diag, threads, 0, stream>>>(
                    n, C_mid, ldc_mid, incC_mid, C, ldc, P, invP, sftA, sftB);
            return;
        }

        if (is_one_h(beta_v)) {
            undo_scaling_r2k_offdiag_kernel_special<HER2K, T, BACKEND, NUM_MODULI, TP, UPLO, -1, 1>
                <<<grid_offdiag, threads, 0, stream>>>(
                    n, C_mid, ldc_mid, incC_mid, C, ldc, P, invP, sftA, sftB);

            undo_scaling_r2k_diag_kernel_special<HER2K, T, BACKEND, NUM_MODULI, TP, UPLO, -1, 1>
                <<<grid_diag, threads, 0, stream>>>(
                    n, C_mid, ldc_mid, incC_mid, C, ldc, P, invP, sftA, sftB);
            return;
        }

        if (is_mone_h(beta_v)) {
            undo_scaling_r2k_offdiag_kernel_special<HER2K, T, BACKEND, NUM_MODULI, TP, UPLO, -1, -1>
                <<<grid_offdiag, threads, 0, stream>>>(
                    n, C_mid, ldc_mid, incC_mid, C, ldc, P, invP, sftA, sftB);

            undo_scaling_r2k_diag_kernel_special<HER2K, T, BACKEND, NUM_MODULI, TP, UPLO, -1, -1>
                <<<grid_diag, threads, 0, stream>>>(
                    n, C_mid, ldc_mid, incC_mid, C, ldc, P, invP, sftA, sftB);
            return;
        }
    }

    using alpha_t = HostScalar<TAlpha>;
    using beta_t  = HostScalar<TBeta>;

    alpha_t alpha_h = alpha_t(alpha_v);
    beta_t beta_h   = beta_t(beta_v);

    undo_scaling_r2k_offdiag_kernel<HER2K, T, BACKEND, NUM_MODULI, TP, UPLO, alpha_t, beta_t>
        <<<grid_offdiag, threads, 0, stream>>>(
            alpha_h, beta_h, n, C_mid, ldc_mid, incC_mid, C, ldc, P, invP, sftA, sftB);

    undo_scaling_r2k_diag_kernel<HER2K, T, BACKEND, NUM_MODULI, TP, UPLO, alpha_t, beta_t>
        <<<grid_diag, threads, 0, stream>>>(
            alpha_h, beta_h, n, C_mid, ldc_mid, incC_mid, C, ldc, P, invP, sftA, sftB);
}

template <typename T, typename TAlpha, typename TBeta,
          Backend BACKEND, unsigned NUM_MODULI, cublasFillMode_t UPLO>
void undo_scaling_syr2k(
    const cudaStream_t stream,
    const unsigned n,
    common::mid_t<BACKEND, common::isComplex<T>> *C_mid,
    const size_t ldc_mid, const size_t incC_mid,
    T *const C, const size_t ldc,
    const int16_t *const sftA, const int16_t *const sftB,
    const TAlpha *alpha, const TBeta *beta //
) {
    undo_scaling_r2k_launch<false, T, TAlpha, TBeta, BACKEND, NUM_MODULI, UPLO>(
        stream, n, C_mid, ldc_mid, incC_mid, C, ldc, sftA, sftB, alpha, beta);
}

template <typename T, typename TAlpha, typename TBeta,
          Backend BACKEND, unsigned NUM_MODULI, cublasFillMode_t UPLO>
void undo_scaling_her2k(
    const cudaStream_t stream,
    const unsigned n,
    common::mid_t<BACKEND, common::isComplex<T>> *C_mid,
    const size_t ldc_mid, const size_t incC_mid,
    T *const C, const size_t ldc,
    const int16_t *const sftA, const int16_t *const sftB,
    const TAlpha *alpha, const TBeta *beta //
) {
    static_assert(common::isComplex<T>, "undo_scaling_her2k requires complex T.");

    undo_scaling_r2k_launch<true, T, TAlpha, TBeta, BACKEND, NUM_MODULI, UPLO>(
        stream, n, C_mid, ldc_mid, incC_mid, C, ldc, sftA, sftB, alpha, beta);
}

} // namespace gemmul8::undo_scaling
