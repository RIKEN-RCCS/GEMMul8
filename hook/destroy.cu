/**
 * Destroy handle and workspace
 */
#include "common.hpp"

// =======================
// Hook: cublasDestroy
// =======================
extern "C" cublasStatus_t cublasDestroy_v2(cublasHandle_t handle) {
#ifdef __CUDA_ARCH__
    return CUBLAS_STATUS_NOT_SUPPORTED;
#else
    using Fn       = cublasStatus_t (*)(cublasHandle_t);
    static Fn real = gemmul8::hook::load_real<Fn>(STR(cublasDestroy_v2));
    if (!real) return CUBLAS_STATUS_NOT_INITIALIZED;

    if (gemmul8::hook::inside_hook()) {
        return real(handle);
    }

    gemmul8::hook::HookGuard guard;
    gemmul8::hook::cleanup_work(handle);

    const cublasStatus_t status = real(handle);
    if (status == CUBLAS_STATUS_SUCCESS) {
        gemmul8::clear_config(handle);
    }
    return status;
#endif
}
