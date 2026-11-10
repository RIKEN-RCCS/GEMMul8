#include "../scaling_rowwise.hpp"

namespace gemmul8::scaling::general {

namespace {

#if !defined(GEMMUL8_INST_TYPE)
    #error "GEMMUL8_INST_TYPE is not defined"
#endif
using T = GEMMUL8_INST_TYPE;

#if !defined(GEMMUL8_INST_BACKEND)
    #error "GEMMUL8_INST_BACKEND is not defined"
#endif
inline constexpr Backend BE = Backend::GEMMUL8_INST_BACKEND;

#if !defined(GEMMUL8_INST_FILLMODE)
    #error "GEMMUL8_INST_FILLMODE is not defined"
#endif
inline constexpr cublasFillMode_t UPLO = GEMMUL8_INST_FILLMODE;

#if !defined(GEMMUL8_INST_DIAG)
    #error "GEMMUL8_INST_DIAG is not defined"
#endif
inline constexpr cublasDiagType_t DIAG = GEMMUL8_INST_DIAG;

#if !defined(GEMMUL8_INST_CONJ)
    #error "GEMMUL8_INST_CONJ is not defined"
#endif
inline constexpr bool CONJ = GEMMUL8_INST_CONJ;

#ifndef GEMMUL8_INST_HERK
    #define GEMMUL8_INST_HERK false
#endif
inline constexpr bool HERK = GEMMUL8_INST_HERK;

} // namespace

static_assert(UPLO != CUBLAS_FILL_MODE_FULL || DIAG == CUBLAS_DIAG_NON_UNIT,
              "When UPLO is CUBLAS_FILL_MODE_FULL, DIAG must be CUBLAS_DIAG_NON_UNIT.");
static_assert(common::isComplex<T> || !CONJ,
              "CONJ=true requires complex T.");

#define INSTANTIATE_THIS(NM)                                          \
    template void scaling_rowwise<T, BE, NM, UPLO, DIAG, CONJ, HERK>( \
        const cudaStream_t,                                           \
        const unsigned, const unsigned,                               \
        const T *const, const size_t,                                 \
        common::matptr_t<common::low_t<BE>, common::isComplex<T>>,    \
        const size_t, const size_t,                                   \
        int16_t *const);

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

} // namespace gemmul8::scaling::general
