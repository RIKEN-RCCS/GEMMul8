#include "../../include/worksize.hpp"

#include "../oz2/gemm/worksize.hpp"
#include "../oz2/symm/worksize.hpp"
#include "../oz2/syrk/worksize.hpp"
#include "../oz2/syr2k/worksize.hpp"
#include "../oz2/syrkx/worksize.hpp"
#include "../oz2/hemm/worksize.hpp"
#include "../oz2/herk/worksize.hpp"
#include "../oz2/her2k/worksize.hpp"
#include "../oz2/herkx/worksize.hpp"
#include "../oz2/trmm/worksize.hpp"
#include "../oz2/trtrmm/worksize.hpp"
#include "../oz2/trsm/worksize.hpp"
#include "../config/config.hpp"

namespace gemmul8 {

#define WORKSIZE_OZ1(F)              \
    if constexpr (FUNC == Func::F) { \
        return 0;                    \
    }

#define WORKSIZE_OZ2(F)                               \
    if constexpr (FUNC == Func::F) {                  \
        return oz2::F::workSize<is_Complex, BACKEND>( \
            m, n, k, NUM_MODULI,                      \
            enable_skip_scalA, enable_skip_scalB,     \
            workSizeA, workSizeB, fastmode);          \
    }

template <bool is_Complex, Backend BACKEND, Func FUNC>
size_t workSize(
    size_t m, size_t n, size_t k, int NUM_MODULI,
    bool enable_skip_scalA, bool enable_skip_scalB,
    size_t *workSizeA, size_t *workSizeB, bool fastmode //
) {
    if (NUM_MODULI > 0) {
        WORKSIZE_OZ2(gemm)
        WORKSIZE_OZ2(symm)
        WORKSIZE_OZ2(syrk)
        WORKSIZE_OZ2(syr2k)
        WORKSIZE_OZ2(syrkx)
        WORKSIZE_OZ2(hemm)
        WORKSIZE_OZ2(herk)
        WORKSIZE_OZ2(her2k)
        WORKSIZE_OZ2(herkx)
        WORKSIZE_OZ2(trmm)
        WORKSIZE_OZ2(trtrmm)
    } else {
        WORKSIZE_OZ1(gemm)
        WORKSIZE_OZ1(symm)
        WORKSIZE_OZ1(syrk)
        WORKSIZE_OZ1(syr2k)
        WORKSIZE_OZ1(syrkx)
        WORKSIZE_OZ1(hemm)
        WORKSIZE_OZ1(herk)
        WORKSIZE_OZ1(her2k)
        WORKSIZE_OZ1(herkx)
        WORKSIZE_OZ1(trmm)
        WORKSIZE_OZ1(trtrmm)
    }

    return 0;
}

template <typename T, Backend BACKEND>
inline size_t workSizeTrsm_impl(
    cublasSideMode_t side,
    size_t m, size_t n,
    int NUM_MODULI,
    int nB_override //
) {
    if (side == CUBLAS_SIDE_LEFT) {
        return oz2::trsm::workSize_left<T, BACKEND>(m, n, NUM_MODULI, nB_override);
    }
    if (side == CUBLAS_SIDE_RIGHT) {
        return oz2::trsm::workSize_right<T, BACKEND>(m, n, NUM_MODULI, nB_override);
    }
    return 0;
}

template <typename T, Backend BACKEND>
size_t workSizeTrsm(
    cublasHandle_t handle,
    cublasSideMode_t side,
    size_t m, size_t n,
    int NUM_MODULI //
) {
    return workSizeTrsm_impl<T, BACKEND>(
        side, m, n, NUM_MODULI,
        config::get_config(handle).block_size_trsm);
}

template <typename T, Backend BACKEND>
size_t workSizeTrsmLt(
    cublasLtHandle_t handle,
    cublasSideMode_t side,
    size_t m, size_t n,
    int NUM_MODULI //
) {
    return workSizeTrsm_impl<T, BACKEND>(
        side, m, n, NUM_MODULI,
        config::get_configLt(handle).block_size_trsm);
}

#undef WORKSIZE_OZ2
#undef WORKSIZE_OZ1

} // namespace gemmul8
