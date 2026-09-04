#include "../worksize_impl.hpp"

namespace gemmul8 {

namespace {

#if !defined(GEMMUL8_INST_BACKEND)
    #error "GEMMUL8_INST_BACKEND is not defined"
#endif
inline constexpr Backend BE = Backend::GEMMUL8_INST_BACKEND;

#if !defined(GEMMUL8_INST_COMPLEX)
    #error "GEMMUL8_INST_COMPLEX is not defined"
#endif
#if GEMMUL8_INST_COMPLEX != 0 && GEMMUL8_INST_COMPLEX != 1
    #error "GEMMUL8_INST_COMPLEX must be 0 or 1"
#endif
inline constexpr bool COMPLEX = (GEMMUL8_INST_COMPLEX != 0);

} // namespace

template size_t workSize<COMPLEX, BE, Func::gemm>(
    size_t, size_t, size_t, int, bool, bool, size_t *, size_t *, bool);

template size_t workSize<COMPLEX, BE, Func::symm>(
    size_t, size_t, size_t, int, bool, bool, size_t *, size_t *, bool);

template size_t workSize<COMPLEX, BE, Func::syrk>(
    size_t, size_t, size_t, int, bool, bool, size_t *, size_t *, bool);

template size_t workSize<COMPLEX, BE, Func::syr2k>(
    size_t, size_t, size_t, int, bool, bool, size_t *, size_t *, bool);

template size_t workSize<COMPLEX, BE, Func::syrkx>(
    size_t, size_t, size_t, int, bool, bool, size_t *, size_t *, bool);

template size_t workSize<COMPLEX, BE, Func::trmm>(
    size_t, size_t, size_t, int, bool, bool, size_t *, size_t *, bool);

template size_t workSize<COMPLEX, BE, Func::trtrmm>(
    size_t, size_t, size_t, int, bool, bool, size_t *, size_t *, bool);

#if GEMMUL8_INST_COMPLEX
template size_t workSize<COMPLEX, BE, Func::hemm>(
    size_t, size_t, size_t, int, bool, bool, size_t *, size_t *, bool);

template size_t workSize<COMPLEX, BE, Func::herk>(
    size_t, size_t, size_t, int, bool, bool, size_t *, size_t *, bool);

template size_t workSize<COMPLEX, BE, Func::her2k>(
    size_t, size_t, size_t, int, bool, bool, size_t *, size_t *, bool);

template size_t workSize<COMPLEX, BE, Func::herkx>(
    size_t, size_t, size_t, int, bool, bool, size_t *, size_t *, bool);

template size_t workSizeTrsm<cuFloatComplex, BE>(
    cublasHandle_t, cublasSideMode_t, size_t, size_t, int);

template size_t workSizeTrsmLt<cuFloatComplex, BE>(
    cublasLtHandle_t, cublasSideMode_t, size_t, size_t, int);

template size_t workSizeTrsm<cuDoubleComplex, BE>(
    cublasHandle_t, cublasSideMode_t, size_t, size_t, int);

template size_t workSizeTrsmLt<cuDoubleComplex, BE>(
    cublasLtHandle_t, cublasSideMode_t, size_t, size_t, int);

#else

template size_t workSizeTrsm<float, BE>(
    cublasHandle_t, cublasSideMode_t, size_t, size_t, int);

template size_t workSizeTrsmLt<float, BE>(
    cublasLtHandle_t, cublasSideMode_t, size_t, size_t, int);

template size_t workSizeTrsm<double, BE>(
    cublasHandle_t, cublasSideMode_t, size_t, size_t, int);

template size_t workSizeTrsmLt<double, BE>(
    cublasLtHandle_t, cublasSideMode_t, size_t, size_t, int);

#endif

} // namespace gemmul8
