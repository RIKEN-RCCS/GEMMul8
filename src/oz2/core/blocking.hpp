#pragma once
#include "../common/common.hpp"
#include "../undo_scaling/predicates.hpp"
#include "../undo_scaling/scalar.hpp"

namespace gemmul8::oz2::core::blocking {

inline void add_timer(std::vector<double> &dst, const std::vector<double> &src) {
    const size_t n = std::min<size_t>(dst.size(), src.size());
    for (size_t i = 0; i < n; ++i) dst[i] += src[i];
}

inline constexpr cublasFillMode_t flip_uplo(cublasFillMode_t uplo) noexcept {
    return (uplo == CUBLAS_FILL_MODE_UPPER) ? CUBLAS_FILL_MODE_LOWER : CUBLAS_FILL_MODE_UPPER;
}

inline constexpr cublasFillMode_t effective_uplo(cublasFillMode_t uplo, cublasOperation_t op) noexcept {
    return (op == CUBLAS_OP_N) ? uplo : flip_uplo(uplo);
}

inline constexpr bool triangular_block_nonzero(cublasFillMode_t effective, size_t row, size_t col) noexcept {
    return (effective == CUBLAS_FILL_MODE_UPPER) ? (row <= col) : (row >= col);
}

struct BlockRange {
    size_t begin = 0;
    size_t end   = 0;
};

// Range of block offsets p for which op(A)(i,p) * op(B)(p,j) can be nonzero.
inline constexpr BlockRange triangular_product_p_range(
    cublasFillMode_t uplo_A, cublasFillMode_t uplo_B, size_t i, size_t j, size_t n, size_t block_size) noexcept {

    if (uplo_A == CUBLAS_FILL_MODE_UPPER) {
        if (uplo_B == CUBLAS_FILL_MODE_UPPER) {
            if (i > j) return {};
            return {i, std::min<size_t>(n, j + block_size)};
        }
        return {std::max<size_t>(i, j), n};
    }

    if (uplo_B == CUBLAS_FILL_MODE_UPPER) {
        return {0, std::min<size_t>(n, std::min<size_t>(i, j) + block_size)};
    }

    if (j > i) return {};
    return {j, std::min<size_t>(n, i + block_size)};
}

// Pointer to the logical block op(A)(row:row+m, col:col+n).
template <typename T>
inline const T *matrix_block_ptr(const T *A, size_t lda, cublasOperation_t op, size_t row, size_t col) noexcept {
    return (op == CUBLAS_OP_N) ? (A + row + col * lda) : (A + col + row * lda);
}

template <typename T>
struct StoredBlock {
    const T *ptr         = nullptr;
    cublasOperation_t op = CUBLAS_OP_N;
};

template <bool HERMITIAN, typename T>
inline StoredBlock<T> stored_structured_block(const T *S, size_t lds, cublasFillMode_t uplo, size_t row, size_t col) noexcept {
    const bool direct = (uplo == CUBLAS_FILL_MODE_UPPER) ? (row < col) : (row > col);
    if (direct) return {S + row + col * lds, CUBLAS_OP_N};
    return {S + col + row * lds, HERMITIAN ? CUBLAS_OP_C : CUBLAS_OP_T};
}

template <typename T, bool NEGATIVE>
__host__ __device__ __forceinline__ T unit_scalar_value() {
    if constexpr (NEGATIVE) {
        return common::Tconst<T>::mone();
    } else {
        return common::Tconst<T>::one();
    }
}

template <typename T, bool NEGATIVE>
__global__ void set_unit_scalar_kernel(T *dst) {
    if (threadIdx.x == 0 && blockIdx.x == 0) *dst = unit_scalar_value<T, NEGATIVE>();
}

template <typename T>
__global__ void conjugate_scalar_kernel(const T *src, T *dst) {
    if (threadIdx.x == 0 && blockIdx.x == 0) *dst = common::conj<T, true>(*src);
}

template <typename TC>
__host__ __device__ __forceinline__ TC complex_from_real(common::underlying_t<TC> x) {
    static_assert(common::isComplex<TC>);
    return TC{x, common::underlying_t<TC>(0)};
}

template <typename TC>
__global__ void real_to_complex_scalar_kernel(const common::underlying_t<TC> *src, TC *dst) {
    if (threadIdx.x == 0 && blockIdx.x == 0) *dst = complex_from_real<TC>(*src);
}

template <typename T>
inline T conjugate_host(T x) noexcept {
    static_assert(common::isComplex<T>);
    x.y = -x.y;
    return x;
}


template <typename T, bool NEGATIVE>
struct UnitScalar {
    T host    = unit_scalar_value<T, NEGATIVE>();
    T *device = nullptr;

    const T *get(const void *reference, cudaStream_t stream) {
        if (!reference || !undo_scaling::is_device_pointer(reference)) return &host;
        if (device) return device;

        if (cudaMallocAsync(reinterpret_cast<void **>(&device), sizeof(T), stream) != cudaSuccess) {
            device = nullptr;
            return nullptr;
        }
        set_unit_scalar_kernel<T, NEGATIVE><<<1, 1, 0, stream>>>(device);
        return device;
    }

    void release(cudaStream_t stream) noexcept {
        if (device) cudaFreeAsync(device, stream);
        device = nullptr;
    }
};

template <typename T>
using OneScalar = UnitScalar<T, false>;

template <typename T>
using MinusOneScalar = UnitScalar<T, true>;

template <typename T>
struct ConjugateScalar {
    static_assert(common::isComplex<T>);
    T host{};
    T *device = nullptr;

    const T *get(const T *src, cudaStream_t stream) {
        if (!undo_scaling::is_device_pointer(src)) {
            host = conjugate_host(*src);
            return &host;
        }
        if (!device && cudaMallocAsync(reinterpret_cast<void **>(&device), sizeof(T), stream) != cudaSuccess) {
            device = nullptr;
            return nullptr;
        }
        conjugate_scalar_kernel<T><<<1, 1, 0, stream>>>(src, device);
        return device;
    }

    void release(cudaStream_t stream) noexcept {
        if (device) cudaFreeAsync(device, stream);
        device = nullptr;
    }
};

template <typename TC>
struct RealToComplexScalar {
    static_assert(common::isComplex<TC>);
    using U = common::underlying_t<TC>;
    TC host{};
    TC *device = nullptr;

    const TC *get(const U *src, cudaStream_t stream) {
        if (src == nullptr) return nullptr;
        if (!undo_scaling::is_device_pointer(src)) {
            host = complex_from_real<TC>(*src);
            return &host;
        }
        if (!device && cudaMallocAsync(reinterpret_cast<void **>(&device), sizeof(TC), stream) != cudaSuccess) {
            device = nullptr;
            return nullptr;
        }
        real_to_complex_scalar_kernel<TC><<<1, 1, 0, stream>>>(src, device);
        return device;
    }

    void release(cudaStream_t stream) noexcept {
        if (device) cudaFreeAsync(device, stream);
        device = nullptr;
    }
};

template <typename T, typename BetaScalar>
__global__ void scale_block_kernel(unsigned m, unsigned n, BetaScalar beta, T *C, size_t ldc) {
    const unsigned row = blockIdx.x * blockDim.x + threadIdx.x;
    const unsigned col = blockIdx.y * blockDim.y + threadIdx.y;
    if (row >= m || col >= n) return;
    C[col * ldc + row] = common::Tmul<undo_scaling::scalar_t<BetaScalar>, T>(beta.get(), C[col * ldc + row]);
}

template <typename T>
__global__ void zero_block_kernel(unsigned m, unsigned n, T *C, size_t ldc) {
    const unsigned row = blockIdx.x * blockDim.x + threadIdx.x;
    const unsigned col = blockIdx.y * blockDim.y + threadIdx.y;
    if (row < m && col < n) C[col * ldc + row] = common::Tconst<T>::zero();
}

template <typename T>
inline void scale_block(cudaStream_t stream, size_t m, size_t n, const T *beta, T *C, size_t ldc) {
    constexpr dim3 threads(32, 8);
    const dim3 grid((m + threads.x - 1) / threads.x, (n + threads.y - 1) / threads.y);

    if (!beta) {
        zero_block_kernel<T><<<grid, threads, 0, stream>>>(unsigned(m), unsigned(n), C, ldc);
    } else if (undo_scaling::is_device_pointer(beta)) {
        undo_scaling::DeviceScalar<T> b(beta);
        scale_block_kernel<T><<<grid, threads, 0, stream>>>(unsigned(m), unsigned(n), b, C, ldc);
    } else {
        undo_scaling::HostScalar<T> b(*beta);
        scale_block_kernel<T><<<grid, threads, 0, stream>>>(unsigned(m), unsigned(n), b, C, ldc);
    }
}

template <class Diagonal, class OffDiagonal>
inline void for_each_rankk_tile(
    cublasFillMode_t uplo,
    size_t n, size_t k, size_t nB, size_t kB,
    Diagonal &&diagonal, OffDiagonal &&offdiag //
) {
    for (size_t i = 0; i < n; i += nB) {
        const size_t ni      = std::min<size_t>(nB, n - i);
        const size_t j_begin = (uplo == CUBLAS_FILL_MODE_UPPER) ? i : 0;
        const size_t j_end   = (uplo == CUBLAS_FILL_MODE_UPPER) ? n : std::min(n, i + nB);

        for (size_t j = j_begin; j < j_end; j += nB) {
            const size_t nj = std::min<size_t>(nB, n - j);
            bool first      = true;
            for (size_t p = 0; p < k; p += kB) {
                const size_t kp = std::min<size_t>(kB, k - p);
                if (i == j) {
                    diagonal(i, ni, p, kp, first);
                } else {
                    offdiag(i, ni, j, nj, p, kp, first);
                }
                first = false;
            }
        }
    }
}

template <bool TRIANGULAR, class Diagonal, class OffDiagonal>
inline void for_each_structured_tile(
    cublasSideMode_t side, cublasFillMode_t effective,
    size_t s, size_t f, size_t sB, size_t fB,
    Diagonal &&diagonal, OffDiagonal &&offdiag //
) {
    if (side == CUBLAS_SIDE_LEFT) {
        for (size_t i = 0; i < s; i += sB) {
            const size_t si = std::min(sB, s - i);
            for (size_t j = 0; j < f; j += fB) {
                const size_t fj = std::min<size_t>(fB, f - j);
                diagonal(i, si, j, fj);

                if constexpr (TRIANGULAR) {
                    const size_t p_begin = (effective == CUBLAS_FILL_MODE_UPPER) ? (i + sB) : 0;
                    const size_t p_end   = (effective == CUBLAS_FILL_MODE_UPPER) ? s : i;
                    for (size_t p = p_begin; p < p_end; p += sB) {
                        offdiag(i, si, p, std::min<size_t>(sB, s - p), j, fj);
                    }
                } else {
                    for (size_t p = 0; p < s; p += sB) {
                        if (p != i) offdiag(i, si, p, std::min<size_t>(sB, s - p), j, fj);
                    }
                }
            }
        }
    } else {
        for (size_t j = 0; j < s; j += sB) {
            const size_t sj = std::min(sB, s - j);
            for (size_t i = 0; i < f; i += fB) {
                const size_t fi = std::min<size_t>(fB, f - i);
                diagonal(j, sj, i, fi);

                if constexpr (TRIANGULAR) {
                    const size_t p_begin = (effective == CUBLAS_FILL_MODE_UPPER) ? 0 : (j + sB);
                    const size_t p_end   = (effective == CUBLAS_FILL_MODE_UPPER) ? j : s;
                    for (size_t p = p_begin; p < p_end; p += sB) {
                        offdiag(j, sj, p, std::min<size_t>(sB, s - p), i, fi);
                    }
                } else {
                    for (size_t p = 0; p < s; p += sB) {
                        if (p != j) offdiag(j, sj, p, std::min<size_t>(sB, s - p), i, fi);
                    }
                }
            }
        }
    }
}

} // namespace gemmul8::oz2::core::blocking
