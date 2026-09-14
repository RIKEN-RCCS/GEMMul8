#include "../undo_scaling.hpp"

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

#if defined(GEMMUL8_INST_HERMITIAN) && GEMMUL8_INST_HERMITIAN
    #define INSTANTIATE_HERMITIAN(NM)                                            \
        template void undo_scaling<T, TALPHA, TBETA, BE, NM, UPLO, false, true>( \
            const cudaStream_t, const unsigned, const unsigned,                  \
            common::mid_t<BE, common::isComplex<T>> *,                           \
            const size_t, const size_t, T *const, const size_t,                  \
            const int16_t *const, const int16_t *const,                          \
            const TALPHA *const, const TBETA *const);
#else
    #define INSTANTIATE_HERMITIAN(NM)
#endif

#define INSTANTIATE_THIS(NM)                                           \
    INSTANTIATE_HERMITIAN(NM)                                          \
    template void undo_scaling<T, TALPHA, TBETA, BE, NM, UPLO, true>(  \
        const cudaStream_t,                                            \
        const unsigned, const unsigned,                                \
        common::mid_t<BE, common::isComplex<T>> *,                     \
        const size_t, const size_t,                                    \
        T *const, const size_t,                                        \
        const int16_t *const, const int16_t *const,                    \
        const TALPHA *const, const TBETA *const);                      \
    template void undo_scaling<T, TALPHA, TBETA, BE, NM, UPLO, false>( \
        const cudaStream_t,                                            \
        const unsigned, const unsigned,                                \
        common::mid_t<BE, common::isComplex<T>> *,                     \
        const size_t, const size_t,                                    \
        T *const, const size_t,                                        \
        const int16_t *const, const int16_t *const,                    \
        const TALPHA *const, const TBETA *const);

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
#undef INSTANTIATE_HERMITIAN

} // namespace gemmul8::undo_scaling
