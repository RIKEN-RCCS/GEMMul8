#include "../undo_scaling_syr2k_her2k.hpp"

namespace gemmul8::undo_scaling {

namespace {

#if !defined(GEMMUL8_INST_TYPE)
    #error "GEMMUL8_INST_TYPE is not defined"
#endif
using T = GEMMUL8_INST_TYPE;

#if !defined(GEMMUL8_INST_TYPE_ALPHA)
    #error "GEMMUL8_INST_TYPE_ALPHA is not defined"
#endif
using TALPHA = GEMMUL8_INST_TYPE_ALPHA;

#if !defined(GEMMUL8_INST_TYPE_BETA)
    #error "GEMMUL8_INST_TYPE_BETA is not defined"
#endif
using TBETA = GEMMUL8_INST_TYPE_BETA;

#if !defined(GEMMUL8_INST_BACKEND)
    #error "GEMMUL8_INST_BACKEND is not defined"
#endif
inline constexpr Backend BE = Backend::GEMMUL8_INST_BACKEND;

#if !defined(GEMMUL8_INST_FILLMODE)
    #error "GEMMUL8_INST_FILLMODE is not defined"
#endif
inline constexpr cublasFillMode_t UPLO = GEMMUL8_INST_FILLMODE;

} // namespace

static_assert(common::isComplex<T>, "This instantiation requires complex T.");
static_assert(UPLO == CUBLAS_FILL_MODE_UPPER || UPLO == CUBLAS_FILL_MODE_LOWER,
              "UPLO must be CUBLAS_FILL_MODE_UPPER or CUBLAS_FILL_MODE_LOWER.");

#define INSTANTIATE_THIS(NM)                                          \
    template void undo_scaling_her2k<T, TALPHA, TBETA, BE, NM, UPLO>( \
        const cudaStream_t,                                           \
        const unsigned,                                               \
        common::mid_t<BE, common::isComplex<T>> *,                    \
        const size_t, const size_t,                                   \
        T *const, const size_t,                                       \
        const int16_t *const, const int16_t *const,                   \
        const TALPHA *const, const TBETA *const, const crt_tail<BE>);

INSTANTIATE_THIS(2U)
INSTANTIATE_THIS(3U)
INSTANTIATE_THIS(4U)
INSTANTIATE_THIS(5U)
INSTANTIATE_THIS(6U)
INSTANTIATE_THIS(7U)
INSTANTIATE_THIS(8U)
INSTANTIATE_THIS(9U)
INSTANTIATE_THIS(10U)
INSTANTIATE_THIS(11U)
INSTANTIATE_THIS(12U)
INSTANTIATE_THIS(13U)
INSTANTIATE_THIS(14U)
INSTANTIATE_THIS(15U)
INSTANTIATE_THIS(16U)
INSTANTIATE_THIS(17U)
INSTANTIATE_THIS(18U)
INSTANTIATE_THIS(19U)
INSTANTIATE_THIS(20U)

#undef INSTANTIATE_THIS

} // namespace gemmul8::undo_scaling
