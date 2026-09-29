#include "../mod_hi2mid.hpp"

namespace gemmul8::mod {

namespace {

#if !defined(GEMMUL8_INST_BACKEND)
    #error "GEMMUL8_INST_BACKEND is not defined"
#endif
inline constexpr Backend BE = Backend::GEMMUL8_INST_BACKEND;

#if !defined(GEMMUL8_INST_FILLMODE)
    #error "GEMMUL8_INST_FILLMODE is not defined"
#endif
inline constexpr cublasFillMode_t UPLO = GEMMUL8_INST_FILLMODE;

} // namespace

static_assert(UPLO == CUBLAS_FILL_MODE_UPPER || UPLO == CUBLAS_FILL_MODE_LOWER,
              "UPLO must be CUBLAS_FILL_MODE_UPPER or CUBLAS_FILL_MODE_LOWER.");

template void mod_hi2mid_AHA<BE, UPLO, true>(
    const cudaStream_t,
    const unsigned,
    const size_t, const unsigned,
    common::matptr_t<common::hi_t<BE>, true> &,
    common::mid_t<BE, true> *,
    const unsigned);

template void mod_hi2mid_AHA<BE, UPLO, false>(
    const cudaStream_t,
    const unsigned,
    const size_t, const unsigned,
    common::matptr_t<common::hi_t<BE>, true> &,
    common::mid_t<BE, true> *,
    const unsigned);

} // namespace gemmul8::mod
