#include "../mod_hi2mid.hpp"

namespace gemmul8::mod {

namespace {

#if !defined(GEMMUL8_INST_BACKEND)
    #error "GEMMUL8_INST_BACKEND is not defined"
#endif
inline constexpr Backend BE = Backend::GEMMUL8_INST_BACKEND;

#if !defined(GEMMUL8_INST_COMPLEX)
    #error "GEMMUL8_INST_COMPLEX is not defined"
#endif
inline constexpr bool COMPLEX = GEMMUL8_INST_COMPLEX;

#if !defined(GEMMUL8_INST_FILLMODE)
    #error "GEMMUL8_INST_FILLMODE is not defined"
#endif
inline constexpr cublasFillMode_t UPLO = GEMMUL8_INST_FILLMODE;

} // namespace

template void mod_hi2mid<BE, COMPLEX, UPLO>(
    const cudaStream_t,
    const unsigned,
    const size_t, const unsigned,
    common::matptr_t<common::hi_t<BE>, COMPLEX> &,
    common::mid_t<BE, COMPLEX> *);

} // namespace gemmul8::mod
