#include "../mod_reduce_matprod.hpp"

namespace gemmul8::mod {

namespace {

#if !defined(GEMMUL8_INST_BACKEND)
    #error "GEMMUL8_INST_BACKEND is not defined"
#endif
inline constexpr Backend BE = Backend::GEMMUL8_INST_BACKEND;

} // namespace

template void mod_reduce_matprod<BE>(
    const cudaStream_t,
    common::hi_t<BE> *,
    const int,
    const int,
    const size_t,
    const unsigned //
);

template void mod_reduce_matprod_strided<BE>(
    const cudaStream_t,
    common::hi_t<BE> *,
    const int,
    const int,
    const size_t,
    const int64_t,
    const int,
    const unsigned //
);

template void mod_reduce_matprod_pointer_and_advance<BE>(
    const cudaStream_t,
    void **,
    void **,
    void **,
    const int,
    const int,
    const size_t,
    const int,
    const int,
    const unsigned //
);

} // namespace gemmul8::mod
