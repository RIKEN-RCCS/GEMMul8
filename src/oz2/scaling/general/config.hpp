#pragma once
#include "../../common/common.hpp"

#define GEMMUL8_FP8_UNPACK(...) __VA_ARGS__

#if defined(__CUDA_ARCH__) && __CUDA_ARCH__ >= 1000 && __CUDA_ARCH__ < 1200
    #define GEMMUL8_FP8_USE_LIMIT (!common::isComplex<T> && NUM_MODULI < 10U && (NUM_MODULI != 8U || std::is_same_v<T, float>))
#elif defined(__CUDA_ARCH__)
    #define GEMMUL8_FP8_USE_LIMIT false
#else
    #define GEMMUL8_FP8_USE_LIMIT true
#endif

#if defined(__CUDACC__)
    #define GEMMUL8_FP8_KERNEL_EX(NAME, PARAMS, ARGS, LIMITED, LIMIT, EXTRA)                                                         \
        template <GEMMUL8_FP8_UNPACK PARAMS, class... KernelArgs> requires(!(((LIMITED) && GEMMUL8_FP8_USE_LIMIT) || (EXTRA) != 0U)) \
        __global__ void NAME(KernelArgs... args) { NAME##_body<GEMMUL8_FP8_UNPACK ARGS>(args...); }                                  \
        template <GEMMUL8_FP8_UNPACK PARAMS, class... KernelArgs> requires(((LIMITED) && GEMMUL8_FP8_USE_LIMIT) || (EXTRA) != 0U)    \
        __global__ __maxnreg__((EXTRA) != 0U ? (EXTRA) : (LIMIT)) void NAME(KernelArgs... args) { NAME##_body<GEMMUL8_FP8_UNPACK ARGS>(args...); }
#else
    #define GEMMUL8_FP8_KERNEL_EX(NAME, PARAMS, ARGS, LIMITED, LIMIT, EXTRA) \
        template <GEMMUL8_FP8_UNPACK PARAMS, class... KernelArgs>            \
        __global__ void NAME(KernelArgs... args) { NAME##_body<GEMMUL8_FP8_UNPACK ARGS>(args...); }
#endif

#define GEMMUL8_FP8_KERNEL(NAME, PARAMS, ARGS, LIMITED, LIMIT) \
    GEMMUL8_FP8_KERNEL_EX(NAME, PARAMS, ARGS, LIMITED, LIMIT, 0U)

#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ == 890 || __CUDA_ARCH__ >= 1000)
    #define GEMMUL8_FP8_COLUMN_LIMIT 40
#else
    #define GEMMUL8_FP8_COLUMN_LIMIT 32
#endif

namespace gemmul8::scaling::general {

enum class fp8_kernel_kind {
    row,
    vector,
    accurate,
    fast
};

#if defined(__CUDA_ARCH__) && __CUDA_ARCH__ >= 1000 && __CUDA_ARCH__ < 1200

template <typename T, unsigned N, fp8_kernel_kind K>
inline constexpr unsigned fp8_refined_registers = [] {
    if constexpr (std::is_same_v<T, float> && N == 20U &&
                  (K == fp8_kernel_kind::accurate ||
                   K == fp8_kernel_kind::fast ||
                   K == fp8_kernel_kind::vector)) {
        return 40U;
    } else if constexpr (std::is_same_v<T, cuFloatComplex> &&
                         (N == 14U || N == 15U) &&
                         K == fp8_kernel_kind::vector) {
        return N == 14U ? 51U : 49U;
    } else {
        return 0U;
    }
}();

template <typename T, unsigned N, bool Aligned = false>
inline constexpr unsigned fp8_vector_registers =
    Aligned && N == 20U ? 40U : 32U;

template <typename T, unsigned N>
inline constexpr unsigned fp8_scalar_registers =
    std::is_same_v<T, float> || N == 20U ? 32U : 0U;

template <typename T, unsigned N>
inline constexpr unsigned fp8_column_registers =
    std::is_same_v<T, double> && (N == 12U || N == 16U || N == 17U)                                ? 38U
    : (std::is_same_v<T, double> && N >= 12U && N < 20U) || (std::is_same_v<T, float> && N >= 19U) ? 40U
                                                                                                   : 32U;

template <typename T, unsigned N>
inline constexpr unsigned fp8_complex_column_registers =
    std::is_same_v<T, cuFloatComplex> && N == 4U ? 40U : std::is_same_v<T, cuFloatComplex> && N == 17U ? 56U
                                                     : std::is_same_v<T, cuDoubleComplex> && N == 10U  ? 48U
                                                                                                       : 0U;

#elif defined(__CUDA_ARCH__)

template <typename T, unsigned N, fp8_kernel_kind K>
inline constexpr unsigned fp8_refined_registers = [] {
    if constexpr (std::is_same_v<T, double> && N == 20U &&
                  (K == fp8_kernel_kind::accurate ||
                   K == fp8_kernel_kind::row)) {
        return 40U;
    } else {
        return 0U;
    }
}();

template <typename T, unsigned N, bool Aligned = false>
inline constexpr unsigned fp8_vector_registers =
    std::is_same_v<T, double> && N == 20U && Aligned                                        ? 34U
    : std::is_same_v<T, float> && (N == 7U || N == 10U)                                     ? 40U
    : N == 18U && std::is_same_v<T, float>                                                  ? 32U
    : std::is_same_v<T, double> || N == 12U || (N == 15U && GEMMUL8_FP8_COLUMN_LIMIT == 40) ? 32U
    : N == 19U || (GEMMUL8_FP8_COLUMN_LIMIT == 40 && N >= 12U)                              ? GEMMUL8_FP8_COLUMN_LIMIT
                                                                                            : 0U;

template <typename T, unsigned N>
inline constexpr unsigned fp8_scalar_registers =
    std::is_same_v<T, float> && (N == 15U || (GEMMUL8_FP8_COLUMN_LIMIT == 40 && (N == 17U || N == 18U || N == 20U)))
        ? GEMMUL8_FP8_COLUMN_LIMIT
        : 0U;

template <typename T, unsigned N>
inline constexpr unsigned fp8_column_registers = GEMMUL8_FP8_COLUMN_LIMIT;

template <typename T, unsigned N>
inline constexpr unsigned fp8_complex_column_registers =
    std::is_same_v<T, cuDoubleComplex> && N == 3U ? 40U : 0U;

#else

template <typename T, unsigned N, fp8_kernel_kind K>
inline constexpr unsigned fp8_refined_registers = 0U;

template <typename T, unsigned N, bool Aligned = false>
inline constexpr unsigned fp8_vector_registers =
    std::is_same_v<T, double> || N == 12U || (N == 15U && GEMMUL8_FP8_COLUMN_LIMIT == 40) ? 32U
    : N == 19U || (GEMMUL8_FP8_COLUMN_LIMIT == 40 && N >= 12U)                            ? GEMMUL8_FP8_COLUMN_LIMIT
                                                                                          : 0U;

template <typename T, unsigned N>
inline constexpr unsigned fp8_scalar_registers =
    std::is_same_v<T, float> && (N == 15U || (GEMMUL8_FP8_COLUMN_LIMIT == 40 && (N == 17U || N == 18U || N == 20U)))
        ? GEMMUL8_FP8_COLUMN_LIMIT
        : 0U;

template <typename T, unsigned N>
inline constexpr unsigned fp8_column_registers = GEMMUL8_FP8_COLUMN_LIMIT;

template <typename T, unsigned N>
inline constexpr unsigned fp8_complex_column_registers = 0U;

#endif

inline constexpr unsigned threads_x_general                            = common::TILE_DIM;
template <Backend BACKEND> inline constexpr unsigned threads_y_rowwise = 4U;
template <> inline constexpr unsigned threads_y_rowwise<Backend::FP8>  = 16U;
template <typename T> inline constexpr unsigned threads_y_symm_hemm    = (common::isComplex<T>) ? 16U : 8U;

inline constexpr unsigned threads_x_colwise_full_tiled = 32U;
inline constexpr unsigned threads_y_colwise_full_tiled = 4U;

} // namespace gemmul8::scaling::general
