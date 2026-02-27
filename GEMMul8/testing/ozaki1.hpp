#pragma once

namespace ozaki1 {

template <typename T1, typename T2> __forceinline__ T1 ceildiv(T1 n, T2 d) { return (n + d - 1) / d; }

size_t workSize(
    int m, int n, int k, int batchCount, bool isComplex,
    cudaEmulationMantissaControl mantissaControl, int maxMantissaBitCount //
) {
    constexpr double MULTIPLIER = 1.25;

    int mult      = isComplex ? 2 : 1;
    int numSlices = ceildiv(maxMantissaBitCount + 1, 8);

    int padded_m     = ceildiv(m, 1024) * 1024;
    int padded_n     = ceildiv(n, 1024) * 1024;
    int padded_k     = ceildiv(k, 128) * 128;
    int num_blocks_k = ceildiv(k, 64);

    size_t gemm_workspace = sizeof(int8_t) *
                            ((size_t)padded_m * padded_k + (size_t)padded_n * padded_k) * mult * numSlices;

    gemm_workspace += sizeof(int32_t) * ((size_t)padded_m + padded_n) * mult;
    if (isComplex) {
        gemm_workspace += sizeof(double) * (size_t)m * n * mult * mult;
    }

    size_t adp_workspace = 0;
    if (mantissaControl == CUDA_EMULATION_MANTISSA_CONTROL_DYNAMIC) {
        adp_workspace = sizeof(int32_t) * ((size_t)m * num_blocks_k + (size_t)n * num_blocks_k + (size_t)m * n) * mult;
    }

    constexpr size_t CONSTANT_SIZE = 128 * 1024 * 1024;
    return (size_t)(std::max(gemm_workspace, adp_workspace) * batchCount * MULTIPLIER) + CONSTANT_SIZE;
}

void setting(cublasHandle_t cublasH, int m, int n, int k, int num_slice) {
    cudaEmulationMantissaControl_t mControl = CUDA_EMULATION_MANTISSA_CONTROL_FIXED;
    int mantissaBitCount                    = num_slice * 8 - 1;
    size_t lwork_oz1                        = ozaki1::workSize(m, n, k, 1, false, mControl, mantissaBitCount);
    void *dwork_oz1;
    cudaMalloc(reinterpret_cast<void **>(&dwork_oz1), lwork_oz1);
    cublasSetWorkspace(cublasH, dwork_oz1, lwork_oz1);
    cublasSetMathMode(cublasH, CUBLAS_FP64_EMULATED_FIXEDPOINT_MATH);
    cublasSetEmulationStrategy(cublasH, CUBLAS_EMULATION_STRATEGY_EAGER);
    cublasSetFixedPointEmulationMantissaControl(cublasH, mControl);
    cublasSetFixedPointEmulationMaxMantissaBitCount(cublasH, mantissaBitCount);
}

} // namespace ozaki1
