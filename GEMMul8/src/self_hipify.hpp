#pragma once

#if defined(__HIPCC__)
    #include <hip/hip_runtime.h>
    #include <hipblas/hipblas.h>
    #include <hipblaslt/hipblaslt.h>
    #include <hip/hip_complex.h>
    #include <hip/hip_fp8.h>

    #if defined(__HIP_DEVICE_COMPILE__) && !defined(__CUDA_ARCH__)
        #define __CUDA_ARCH__ 0
    #endif
    #define CUBLAS_VER_MAJOR 0

    // cuBLAS
    #if defined(HIPBLAS_V2)
        #define cublasGemmEx   hipblasGemmEx
        #define cublasCgemm_v2 hipblasCgemm
        #define cublasZgemm_v2 hipblasZgemm
    #else
        #define cublasGemmEx   hipblasGemmEx_v2
        #define cublasCgemm_v2 hipblasCgemm_v2
        #define cublasZgemm_v2 hipblasZgemm_v2
    #endif
    #define cublasSgemm_v2                 hipblasSgemm
    #define cublasDgemm_v2                 hipblasDgemm
    #define cublasGetStream                hipblasGetStream
    #define cublasGetStream_v2             hipblasGetStream
    #define CUBLAS_OP_N                    HIPBLAS_OP_N
    #define CUBLAS_OP_T                    HIPBLAS_OP_T
    #define CUBLAS_OP_C                    HIPBLAS_OP_C
    #define CUDA_R_8F_E4M3                 HIP_R_8F_E4M3
    #define CUDA_R_8I                      HIP_R_8I
    #define CUDA_R_32I                     HIP_R_32I
    #define CUDA_R_32F                     HIP_R_32F
    #define CUDA_R_64F                     HIP_R_64F
    #define CUDA_C_32F                     HIP_C_32F
    #define CUDA_C_64F                     HIP_C_64F
    #define CUBLAS_COMPUTE_32I             HIPBLAS_COMPUTE_32I
    #define CUBLAS_GEMM_DEFAULT            HIPBLAS_GEMM_DEFAULT
    #define CUBLAS_COMPUTE_32F             HIPBLAS_COMPUTE_32F
    #define CUBLAS_COMPUTE_64F             HIPBLAS_COMPUTE_64F
    #define CUBLAS_COMPUTE_32F_FAST_TF32   HIPBLAS_COMPUTE_32F_FAST_TF32
    #define CUBLAS_STATUS_EXECUTION_FAILED HIPBLAS_STATUS_EXECUTION_FAILED
    #define CUBLAS_STATUS_INVALID_VALUE    HIPBLAS_STATUS_INVALID_VALUE
    #define CUBLAS_STATUS_INTERNAL_ERROR   HIPBLAS_STATUS_INTERNAL_ERROR
    #define cublasCreate                   hipblasCreate
    #define cublasDestroy                  hipblasDestroy
    #define cublasDestroy_v2               hipblasDestroy
    #define cublasHandle_t                 hipblasHandle_t
    #define cublasOperation_t              hipblasOperation_t
    #define cublasStatus_t                 hipblasStatus_t
    #define cublasComputeType_t            hipblasComputeType_t
    #define cublasGemmAlgo_t               hipblasGemmAlgo_t
    #define CUBLAS_STATUS_SUCCESS          HIPBLAS_STATUS_SUCCESS
    #define CUBLAS_STATUS_NOT_SUPPORTED    HIPBLAS_STATUS_NOT_SUPPORTED
    #define CUBLAS_STATUS_ALLOC_FAILED     HIPBLAS_STATUS_ALLOC_FAILED
    #define CUBLAS_STATUS_NOT_INITIALIZED  HIPBLAS_STATUS_NOT_INITIALIZED

    // cuBLASLt
    #define cublasLtCreate                           hipblasLtCreate
    #define cublasLtDestroy                          hipblasLtDestroy
    #define CUBLASLT_MATMUL_DESC_TRANSA              HIPBLASLT_MATMUL_DESC_TRANSA
    #define CUBLASLT_MATMUL_DESC_TRANSB              HIPBLASLT_MATMUL_DESC_TRANSB
    #define CUBLASLT_MATMUL_PREF_MAX_WORKSPACE_BYTES HIPBLASLT_MATMUL_PREF_MAX_WORKSPACE_BYTES
    #define cublasLtHandle_t                         hipblasLtHandle_t
    #define cublasLtMatmulDesc_t                     hipblasLtMatmulDesc_t
    #define cublasLtMatrixLayout_t                   hipblasLtMatrixLayout_t
    #define cublasLtMatmulPreference_t               hipblasLtMatmulPreference_t
    #define cublasLtMatmulHeuristicResult_t          hipblasLtMatmulHeuristicResult_t
    #define cublasLtMatmulDescCreate                 hipblasLtMatmulDescCreate
    #define cublasLtMatmulDescSetAttribute           hipblasLtMatmulDescSetAttribute
    #define cublasLtMatmulPreferenceCreate           hipblasLtMatmulPreferenceCreate
    #define cublasLtMatmulPreferenceSetAttribute     hipblasLtMatmulPreferenceSetAttribute
    #define cublasLtMatrixLayoutCreate               hipblasLtMatrixLayoutCreate
    #define cublasLtMatmulAlgoGetHeuristic           hipblasLtMatmulAlgoGetHeuristic
    #define cublasLtMatmul                           hipblasLtMatmul
    #define cublasLtMatmulPreferenceDestroy          hipblasLtMatmulPreferenceDestroy
    #define cublasLtMatrixLayoutDestroy              hipblasLtMatrixLayoutDestroy
    #define cublasLtMatmulDescDestroy                hipblasLtMatmulDescDestroy

    // CUDA
    #define __nv_fp8_e4m3              __hip_fp8_e4m3_fnuz
    #define __nv_fp8x2_e4m3            __hip_fp8x2_e4m3_fnuz
    #define __nv_fp8x4_e4m3            __hip_fp8x4_e4m3_fnuz
    #define __nv_fp8_storage_t         __hip_fp8_storage_t
    #define __nv_cvt_float_to_fp8      __hip_cvt_float_to_fp8
    #define __NV_SATFINITE             __HIP_SATFINITE
    #define __NV_E4M3                  __HIP_E4M3_FNUZ
    #define cudaDeviceSynchronize      hipDeviceSynchronize
    #define cudaStreamSynchronize      hipStreamSynchronize
    #define cudaDeviceProp             hipDeviceProp_t
    #define cudaGetDeviceProperties    hipGetDeviceProperties
    #define cudaMalloc                 hipMalloc
    #define cudaMallocAsync            hipMallocAsync
    #define cudaMemcpy                 hipMemcpy
    #define cudaMemcpyAsync            hipMemcpyAsync
    #define cudaMemcpyDeviceToHost     hipMemcpyDeviceToHost
    #define cudaFree                   hipFree
    #define cudaFreeAsync              hipFreeAsync
    #define cudaError_t                hipError_t
    #define cudaSuccess                hipSuccess
    #define cudaDataType               hipDataType
    #define cudaGetErrorString         hipGetErrorString
    #define cudaStream_t               hipStream_t
    #define cudaPointerAttributes      hipPointerAttribute_t
    #define cudaPointerGetAttributes   hipPointerGetAttributes
    #define cudaMemoryTypeHost         hipMemoryTypeHost
    #define cudaMemoryTypeUnregistered hipMemoryTypeUnregistered
    #define cudaEvent_t                hipEvent_t
    #define cudaEventCreate            hipEventCreate
    #define cudaEventCreateWithFlags   hipEventCreateWithFlags
    #define cudaEventRecord            hipEventRecord
    #define cudaStreamWaitEvent        hipStreamWaitEvent
    #define cudaEventDestroy           hipEventDestroy
    #define cudaEventDisableTiming     hipEventDisableTiming
    #define cudaSetDevice              hipSetDevice
    #define cudaMemcpyKind             hipMemcpyKind
    #define cudaMemcpyHostToDevice     hipMemcpyHostToDevice
    #define cudaMemcpyDeviceToDevice   hipMemcpyDeviceToDevice
    #define cudaMemcpyHostToHost       hipMemcpyHostToHost
    #define cudaMemcpyDefault          hipMemcpyDefault
    #define cudaGetLastError           hipGetLastError
    #define cudaEventSynchronize       hipEventSynchronize
    #define cudaEventElapsedTime       hipEventElapsedTime

    #define cudaMemcpyToSymbol(symbol, src, count)                            hipMemcpyToSymbol(HIP_SYMBOL(symbol), src, count)
    #define cudaMemcpyToSymbolAsync(symbol, src, count, offset, kind, stream) hipMemcpyToSymbolAsync(HIP_SYMBOL(symbol), src, count, offset, kind, stream)
    #define __shfl_down_sync(mask, val, offset, width)                        __shfl_down(val, offset, width)

    // cuComplex
    #define cuComplex              hipFloatComplex
    #define cuFloatComplex         hipFloatComplex
    #define cuDoubleComplex        hipDoubleComplex
    #define cuCreal                hipCreal
    #define cuCrealf               hipCrealf
    #define cuCimag                hipCimag
    #define cuCimagf               hipCimagf
    #define make_cuComplex         make_hipFloatComplex
    #define make_cuFloatComplex    make_hipFloatComplex
    #define make_cuDoubleComplex   make_hipDoubleComplex
    #define cuCmul                 hipCmul
    #define cuCmulf                hipCmulf
    #define cuCadd                 hipCadd
    #define cuCaddf                hipCaddf
    #define cuCfma                 hipCfma
    #define cuCfmaf                hipCfmaf
    #define cuCabs                 hipCabs
    #define cuCsub                 hipCsub
    #define cuCsubf                hipCsubf
    #define cuComplexFloatToDouble hipComplexFloatToDouble
    #define cuComplexDoubleToFloat hipComplexDoubleToFloat

    #define STR_MACRO(x) #x
    #define STR(x)       STR_MACRO(x)

#else
    // No change for cuda
    #define STR(x) #x
#endif
