#pragma once

#if defined(__HIPCC__)
    #include <hip/hip_fp8.h>

    #if defined(__HIP_DEVICE_COMPILE__) && !defined(__CUDA_ARCH__)
        #define __CUDA_ARCH__ 0
    #endif
    #define CUBLAS_VER_MAJOR 0

    // cuBLAS
    #if defined(HIPBLAS_V2)

        #define cublasSgemm                   hipblasSgemm
        #define cublasSgemm_v2                hipblasSgemm
        #define cublasSgemm_64                hipblasSgemm_64
        #define cublasSgemm_v2_64             hipblasSgemm_64
        #define cublasDgemm                   hipblasDgemm
        #define cublasDgemm_v2                hipblasDgemm
        #define cublasDgemm_64                hipblasDgemm_64
        #define cublasDgemm_v2_64             hipblasDgemm_64
        #define cublasCgemm                   hipblasCgemm
        #define cublasCgemm_v2                hipblasCgemm
        #define cublasCgemm_64                hipblasCgemm_64
        #define cublasCgemm_v2_64             hipblasCgemm_64
        #define cublasZgemm                   hipblasZgemm
        #define cublasZgemm_v2                hipblasZgemm
        #define cublasZgemm_64                hipblasZgemm_64
        #define cublasZgemm_v2_64             hipblasZgemm_64
        #define cublasGemmEx                  hipblasGemmEx
        #define cublasGemmEx_64               hipblasGemmEx_64
        #define cublasGemmStridedBatchedEx    hipblasGemmStridedBatchedEx
        #define cublasGemmStridedBatchedEx_64 hipblasGemmStridedBatchedEx_64

        #define cublasSsymm       hipblasSsymm
        #define cublasSsymm_v2    hipblasSsymm
        #define cublasSsymm_64    hipblasSsymm_64
        #define cublasSsymm_v2_64 hipblasSsymm_64
        #define cublasDsymm       hipblasDsymm
        #define cublasDsymm_v2    hipblasDsymm
        #define cublasDsymm_64    hipblasDsymm_64
        #define cublasDsymm_v2_64 hipblasDsymm_64
        #define cublasCsymm       hipblasCsymm
        #define cublasCsymm_v2    hipblasCsymm
        #define cublasCsymm_64    hipblasCsymm_64
        #define cublasCsymm_v2_64 hipblasCsymm_64
        #define cublasZsymm       hipblasZsymm
        #define cublasZsymm_v2    hipblasZsymm
        #define cublasZsymm_64    hipblasZsymm_64
        #define cublasZsymm_v2_64 hipblasZsymm_64

        #define cublasSsyrk       hipblasSsyrk
        #define cublasSsyrk_v2    hipblasSsyrk
        #define cublasSsyrk_64    hipblasSsyrk_64
        #define cublasSsyrk_v2_64 hipblasSsyrk_64
        #define cublasDsyrk       hipblasDsyrk
        #define cublasDsyrk_v2    hipblasDsyrk
        #define cublasDsyrk_64    hipblasDsyrk_64
        #define cublasDsyrk_v2_64 hipblasDsyrk_64
        #define cublasCsyrk       hipblasCsyrk
        #define cublasCsyrk_v2    hipblasCsyrk
        #define cublasCsyrk_64    hipblasCsyrk_64
        #define cublasCsyrk_v2_64 hipblasCsyrk_64
        #define cublasZsyrk       hipblasZsyrk
        #define cublasZsyrk_v2    hipblasZsyrk
        #define cublasZsyrk_64    hipblasZsyrk_64
        #define cublasZsyrk_v2_64 hipblasZsyrk_64

        #define cublasSsyr2k       hipblasSsyr2k
        #define cublasSsyr2k_v2    hipblasSsyr2k
        #define cublasSsyr2k_64    hipblasSsyr2k_64
        #define cublasSsyr2k_v2_64 hipblasSsyr2k_64
        #define cublasDsyr2k       hipblasDsyr2k
        #define cublasDsyr2k_v2    hipblasDsyr2k
        #define cublasDsyr2k_64    hipblasDsyr2k_64
        #define cublasDsyr2k_v2_64 hipblasDsyr2k_64
        #define cublasCsyr2k       hipblasCsyr2k
        #define cublasCsyr2k_v2    hipblasCsyr2k
        #define cublasCsyr2k_64    hipblasCsyr2k_64
        #define cublasCsyr2k_v2_64 hipblasCsyr2k_64
        #define cublasZsyr2k       hipblasZsyr2k
        #define cublasZsyr2k_v2    hipblasZsyr2k
        #define cublasZsyr2k_64    hipblasZsyr2k_64
        #define cublasZsyr2k_v2_64 hipblasZsyr2k_64

        #define cublasSsyrkx    hipblasSsyrkx
        #define cublasSsyrkx_64 hipblasSsyrkx_64
        #define cublasDsyrkx    hipblasDsyrkx
        #define cublasDsyrkx_64 hipblasDsyrkx_64
        #define cublasCsyrkx    hipblasCsyrkx
        #define cublasCsyrkx_64 hipblasCsyrkx_64
        #define cublasZsyrkx    hipblasZsyrkx
        #define cublasZsyrkx_64 hipblasZsyrkx_64

        #define cublasChemm       hipblasChemm
        #define cublasChemm_v2    hipblasChemm
        #define cublasChemm_64    hipblasChemm_64
        #define cublasChemm_v2_64 hipblasChemm_64
        #define cublasZhemm       hipblasZhemm
        #define cublasZhemm_v2    hipblasZhemm
        #define cublasZhemm_64    hipblasZhemm_64
        #define cublasZhemm_v2_64 hipblasZhemm_64

        #define cublasCherk       hipblasCherk
        #define cublasCherk_v2    hipblasCherk
        #define cublasCherk_64    hipblasCherk_64
        #define cublasCherk_v2_64 hipblasCherk_64
        #define cublasZherk       hipblasZherk
        #define cublasZherk_v2    hipblasZherk
        #define cublasZherk_64    hipblasZherk_64
        #define cublasZherk_v2_64 hipblasZherk_64

        #define cublasCher2k       hipblasCher2k
        #define cublasCher2k_v2    hipblasCher2k
        #define cublasCher2k_64    hipblasCher2k_64
        #define cublasCher2k_v2_64 hipblasCher2k_64
        #define cublasZher2k       hipblasZher2k
        #define cublasZher2k_v2    hipblasZher2k
        #define cublasZher2k_64    hipblasZher2k_64
        #define cublasZher2k_v2_64 hipblasZher2k_64

        #define cublasCherkx    hipblasCherkx
        #define cublasCherkx_64 hipblasCherkx_64
        #define cublasZherkx    hipblasZherkx
        #define cublasZherkx_64 hipblasZherkx_64

        #define cublasStrmm       hipblasStrmm
        #define cublasStrmm_v2    hipblasStrmm
        #define cublasStrmm_64    hipblasStrmm_64
        #define cublasStrmm_v2_64 hipblasStrmm_64
        #define cublasDtrmm       hipblasDtrmm
        #define cublasDtrmm_v2    hipblasDtrmm
        #define cublasDtrmm_64    hipblasDtrmm_64
        #define cublasDtrmm_v2_64 hipblasDtrmm_64
        #define cublasCtrmm       hipblasCtrmm
        #define cublasCtrmm_v2    hipblasCtrmm
        #define cublasCtrmm_64    hipblasCtrmm_64
        #define cublasCtrmm_v2_64 hipblasCtrmm_64
        #define cublasZtrmm       hipblasZtrmm
        #define cublasZtrmm_v2    hipblasZtrmm
        #define cublasZtrmm_64    hipblasZtrmm_64
        #define cublasZtrmm_v2_64 hipblasZtrmm_64

        #define cublasStrsm       hipblasStrsm
        #define cublasStrsm_v2    hipblasStrsm
        #define cublasStrsm_64    hipblasStrsm_64
        #define cublasStrsm_v2_64 hipblasStrsm_64
        #define cublasDtrsm       hipblasDtrsm
        #define cublasDtrsm_v2    hipblasDtrsm
        #define cublasDtrsm_64    hipblasDtrsm_64
        #define cublasDtrsm_v2_64 hipblasDtrsm_64
        #define cublasCtrsm       hipblasCtrsm
        #define cublasCtrsm_v2    hipblasCtrsm
        #define cublasCtrsm_64    hipblasCtrsm_64
        #define cublasCtrsm_v2_64 hipblasCtrsm_64
        #define cublasZtrsm       hipblasZtrsm
        #define cublasZtrsm_v2    hipblasZtrsm
        #define cublasZtrsm_64    hipblasZtrsm_64
        #define cublasZtrsm_v2_64 hipblasZtrsm_64

    #else

        #define cublasSgemm                   hipblasSgemm
        #define cublasSgemm_v2                hipblasSgemm
        #define cublasSgemm_64                hipblasSgemm_64
        #define cublasSgemm_v2_64             hipblasSgemm_64
        #define cublasDgemm                   hipblasDgemm
        #define cublasDgemm_v2                hipblasDgemm
        #define cublasDgemm_64                hipblasDgemm_64
        #define cublasDgemm_v2_64             hipblasDgemm_64
        #define cublasCgemm                   hipblasCgemm_v2
        #define cublasCgemm_v2                hipblasCgemm_v2
        #define cublasCgemm_64                hipblasCgemm_v2_64
        #define cublasCgemm_v2_64             hipblasCgemm_v2_64
        #define cublasZgemm                   hipblasZgemm_v2
        #define cublasZgemm_v2                hipblasZgemm_v2
        #define cublasZgemm_64                hipblasZgemm_v2_64
        #define cublasZgemm_v2_64             hipblasZgemm_v2_64
        #define cublasGemmEx                  hipblasGemmEx_v2
        #define cublasGemmEx_64               hipblasGemmEx_v2_64
        #define cublasGemmStridedBatchedEx    hipblasGemmStridedBatchedEx_v2
        #define cublasGemmStridedBatchedEx_64 hipblasGemmStridedBatchedEx_v2_64

        #define cublasSsymm       hipblasSsymm
        #define cublasSsymm_v2    hipblasSsymm
        #define cublasSsymm_64    hipblasSsymm_64
        #define cublasSsymm_v2_64 hipblasSsymm_64
        #define cublasDsymm       hipblasDsymm
        #define cublasDsymm_v2    hipblasDsymm
        #define cublasDsymm_64    hipblasDsymm_64
        #define cublasDsymm_v2_64 hipblasDsymm_64
        #define cublasCsymm       hipblasCsymm_v2
        #define cublasCsymm_v2    hipblasCsymm_v2
        #define cublasCsymm_64    hipblasCsymm_v2_64
        #define cublasCsymm_v2_64 hipblasCsymm_v2_64
        #define cublasZsymm       hipblasZsymm_v2
        #define cublasZsymm_v2    hipblasZsymm_v2
        #define cublasZsymm_64    hipblasZsymm_v2_64
        #define cublasZsymm_v2_64 hipblasZsymm_v2_64

        #define cublasSsyrk       hipblasSsyrk
        #define cublasSsyrk_v2    hipblasSsyrk
        #define cublasSsyrk_64    hipblasSsyrk_64
        #define cublasSsyrk_v2_64 hipblasSsyrk_64
        #define cublasDsyrk       hipblasDsyrk
        #define cublasDsyrk_v2    hipblasDsyrk
        #define cublasDsyrk_64    hipblasDsyrk_64
        #define cublasDsyrk_v2_64 hipblasDsyrk_64
        #define cublasCsyrk       hipblasCsyrk_v2
        #define cublasCsyrk_v2    hipblasCsyrk_v2
        #define cublasCsyrk_64    hipblasCsyrk_v2_64
        #define cublasCsyrk_v2_64 hipblasCsyrk_v2_64
        #define cublasZsyrk       hipblasZsyrk_v2
        #define cublasZsyrk_v2    hipblasZsyrk_v2
        #define cublasZsyrk_64    hipblasZsyrk_v2_64
        #define cublasZsyrk_v2_64 hipblasZsyrk_v2_64

        #define cublasSsyr2k       hipblasSsyr2k
        #define cublasSsyr2k_v2    hipblasSsyr2k
        #define cublasSsyr2k_64    hipblasSsyr2k_64
        #define cublasSsyr2k_v2_64 hipblasSsyr2k_64
        #define cublasDsyr2k       hipblasDsyr2k
        #define cublasDsyr2k_v2    hipblasDsyr2k
        #define cublasDsyr2k_64    hipblasDsyr2k_64
        #define cublasDsyr2k_v2_64 hipblasDsyr2k_64
        #define cublasCsyr2k       hipblasCsyr2k_v2
        #define cublasCsyr2k_v2    hipblasCsyr2k_v2
        #define cublasCsyr2k_64    hipblasCsyr2k_v2_64
        #define cublasCsyr2k_v2_64 hipblasCsyr2k_v2_64
        #define cublasZsyr2k       hipblasZsyr2k_v2
        #define cublasZsyr2k_v2    hipblasZsyr2k_v2
        #define cublasZsyr2k_64    hipblasZsyr2k_v2_64
        #define cublasZsyr2k_v2_64 hipblasZsyr2k_v2_64

        #define cublasSsyrkx    hipblasSsyrkx
        #define cublasSsyrkx_64 hipblasSsyrkx_64
        #define cublasDsyrkx    hipblasDsyrkx
        #define cublasDsyrkx_64 hipblasDsyrkx_64
        #define cublasCsyrkx    hipblasCsyrkx_v2
        #define cublasCsyrkx_64 hipblasCsyrkx_v2_64
        #define cublasZsyrkx    hipblasZsyrkx_v2
        #define cublasZsyrkx_64 hipblasZsyrkx_v2_64

        #define cublasChemm       hipblasChemm_v2
        #define cublasChemm_v2    hipblasChemm_v2
        #define cublasChemm_64    hipblasChemm_v2_64
        #define cublasChemm_v2_64 hipblasChemm_v2_64
        #define cublasZhemm       hipblasZhemm_v2
        #define cublasZhemm_v2    hipblasZhemm_v2
        #define cublasZhemm_64    hipblasZhemm_v2_64
        #define cublasZhemm_v2_64 hipblasZhemm_v2_64

        #define cublasCherk       hipblasCherk_v2
        #define cublasCherk_v2    hipblasCherk_v2
        #define cublasCherk_64    hipblasCherk_v2_64
        #define cublasCherk_v2_64 hipblasCherk_v2_64
        #define cublasZherk       hipblasZherk_v2
        #define cublasZherk_v2    hipblasZherk_v2
        #define cublasZherk_64    hipblasZherk_v2_64
        #define cublasZherk_v2_64 hipblasZherk_v2_64

        #define cublasCher2k       hipblasCher2k_v2
        #define cublasCher2k_v2    hipblasCher2k_v2
        #define cublasCher2k_64    hipblasCher2k_v2_64
        #define cublasCher2k_v2_64 hipblasCher2k_v2_64
        #define cublasZher2k       hipblasZher2k_v2
        #define cublasZher2k_v2    hipblasZher2k_v2
        #define cublasZher2k_64    hipblasZher2k_v2_64
        #define cublasZher2k_v2_64 hipblasZher2k_v2_64

        #define cublasCherkx    hipblasCherkx_v2
        #define cublasCherkx_64 hipblasCherkx_v2_64
        #define cublasZherkx    hipblasZherkx_v2
        #define cublasZherkx_64 hipblasZherkx_v2_64

        #define cublasStrmm       hipblasStrmm
        #define cublasStrmm_v2    hipblasStrmm
        #define cublasStrmm_64    hipblasStrmm_64
        #define cublasStrmm_v2_64 hipblasStrmm_64
        #define cublasDtrmm       hipblasDtrmm
        #define cublasDtrmm_v2    hipblasDtrmm
        #define cublasDtrmm_64    hipblasDtrmm_64
        #define cublasDtrmm_v2_64 hipblasDtrmm_64
        #define cublasCtrmm       hipblasCtrmm_v2
        #define cublasCtrmm_v2    hipblasCtrmm_v2
        #define cublasCtrmm_64    hipblasCtrmm_v2_64
        #define cublasCtrmm_v2_64 hipblasCtrmm_v2_64
        #define cublasZtrmm       hipblasZtrmm_v2
        #define cublasZtrmm_v2    hipblasZtrmm_v2
        #define cublasZtrmm_64    hipblasZtrmm_v2_64
        #define cublasZtrmm_v2_64 hipblasZtrmm_v2_64

        #define cublasStrsm       hipblasStrsm
        #define cublasStrsm_v2    hipblasStrsm
        #define cublasStrsm_64    hipblasStrsm_64
        #define cublasStrsm_v2_64 hipblasStrsm_64
        #define cublasDtrsm       hipblasDtrsm
        #define cublasDtrsm_v2    hipblasDtrsm
        #define cublasDtrsm_64    hipblasDtrsm_64
        #define cublasDtrsm_v2_64 hipblasDtrsm_64
        #define cublasCtrsm       hipblasCtrsm_v2
        #define cublasCtrsm_v2    hipblasCtrsm_v2
        #define cublasCtrsm_64    hipblasCtrsm_v2_64
        #define cublasCtrsm_v2_64 hipblasCtrsm_v2_64
        #define cublasZtrsm       hipblasZtrsm_v2
        #define cublasZtrsm_v2    hipblasZtrsm_v2
        #define cublasZtrsm_64    hipblasZtrsm_v2_64
        #define cublasZtrsm_v2_64 hipblasZtrsm_v2_64
    #endif

    #define cublasSetStream         hipblasSetStream
    #define cublasSetStream_v2      hipblasSetStream
    #define cublasGetStream         hipblasGetStream
    #define cublasGetStream_v2      hipblasGetStream
    #define cublasCreate            hipblasCreate
    #define cublasCreate_v2         hipblasCreate
    #define cublasDestroy           hipblasDestroy
    #define cublasDestroy_v2        hipblasDestroy
    #define cublasSetWorkspace      hipblasSetWorkspace
    #define cublasSetWorkspace_v2   hipblasSetWorkspace
    #define cublasSetPointerMode    hipblasSetPointerMode
    #define cublasSetPointerMode_v2 hipblasSetPointerMode
    #define cublasGetPointerMode    hipblasGetPointerMode
    #define cublasGetPointerMode_v2 hipblasGetPointerMode

    #define cublasHandle_t      hipblasHandle_t
    #define cublasOperation_t   hipblasOperation_t
    #define cublasStatus_t      hipblasStatus_t
    #define cublasComputeType_t hipblasComputeType_t
    #define cublasGemmAlgo_t    hipblasGemmAlgo_t
    #define cublasPointerMode_t hipblasPointerMode_t

    #define CUBLAS_SIDE_LEFT               HIPBLAS_SIDE_LEFT
    #define CUBLAS_SIDE_RIGHT              HIPBLAS_SIDE_RIGHT
    #define CUBLAS_FILL_MODE_LOWER         HIPBLAS_FILL_MODE_LOWER
    #define CUBLAS_FILL_MODE_UPPER         HIPBLAS_FILL_MODE_UPPER
    #define CUBLAS_FILL_MODE_FULL          HIPBLAS_FILL_MODE_FULL
    #define CUBLAS_DIAG_UNIT               HIPBLAS_DIAG_UNIT
    #define CUBLAS_DIAG_NON_UNIT           HIPBLAS_DIAG_NON_UNIT
    #define CUBLAS_OP_N                    HIPBLAS_OP_N
    #define CUBLAS_OP_T                    HIPBLAS_OP_T
    #define CUBLAS_OP_C                    HIPBLAS_OP_C
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
    #define CUBLAS_STATUS_SUCCESS          HIPBLAS_STATUS_SUCCESS
    #define CUBLAS_STATUS_NOT_SUPPORTED    HIPBLAS_STATUS_NOT_SUPPORTED
    #define CUBLAS_STATUS_ALLOC_FAILED     HIPBLAS_STATUS_ALLOC_FAILED
    #define CUBLAS_STATUS_NOT_INITIALIZED  HIPBLAS_STATUS_NOT_INITIALIZED
    #define CUBLAS_POINTER_MODE_HOST       HIPBLAS_POINTER_MODE_HOST
    #define CUBLAS_POINTER_MODE_DEVICE     HIPBLAS_POINTER_MODE_DEVICE

    #define cublasSideMode_t hipblasSideMode_t
    #define cublasFillMode_t hipblasFillMode_t
    #define cublasDiagType_t hipblasDiagType_t

    #define GPU_ARCH_ID_gfx940  940
    #define GPU_ARCH_ID_gfx941  941
    #define GPU_ARCH_ID_gfx942  942
    #define GPU_ARCH_ID_IMPL(x) GPU_ARCH_ID_##x
    #define GPU_ARCH_ID(x)      GPU_ARCH_ID_IMPL(x)
    #if defined(GPU_ARCH) && (GPU_ARCH_ID(GPU_ARCH) == 940 || GPU_ARCH_ID(GPU_ARCH) == 941 || GPU_ARCH_ID(GPU_ARCH) == 942)
        #define FP8_FNUZ        1
        #define CUDA_R_8F_E4M3  HIP_R_8F_E4M3_FNUZ
        #define __nv_fp8_e4m3   __hip_fp8_e4m3_fnuz
        #define __nv_fp8x2_e4m3 __hip_fp8x2_e4m3_fnuz
        #define __nv_fp8x4_e4m3 __hip_fp8x4_e4m3_fnuz
        #define __NV_E4M3       __HIP_E4M3_FNUZ
    #else
        #define FP8_FNUZ        0
        #define CUDA_R_8F_E4M3  HIP_R_8F_E4M3
        #define __nv_fp8_e4m3   __hip_fp8_e4m3
        #define __nv_fp8x2_e4m3 __hip_fp8x2_e4m3
        #define __nv_fp8x4_e4m3 __hip_fp8x4_e4m3
        #define __NV_E4M3       __HIP_E4M3
    #endif

    // cuBLASLt
    #define cublasLtHandle_t                hipblasLtHandle_t
    #define cublasLtMatmulDesc_t            hipblasLtMatmulDesc_t
    #define cublasLtMatrixLayout_t          hipblasLtMatrixLayout_t
    #define cublasLtMatmulPreference_t      hipblasLtMatmulPreference_t
    #define cublasLtMatmulHeuristicResult_t hipblasLtMatmulHeuristicResult_t

    #define cublasLtCreate                       hipblasLtCreate
    #define cublasLtDestroy                      hipblasLtDestroy
    #define cublasLtMatmulDescCreate             hipblasLtMatmulDescCreate
    #define cublasLtMatmulDescSetAttribute       hipblasLtMatmulDescSetAttribute
    #define cublasLtMatmulPreferenceCreate       hipblasLtMatmulPreferenceCreate
    #define cublasLtMatmulPreferenceSetAttribute hipblasLtMatmulPreferenceSetAttribute
    #define cublasLtMatrixLayoutCreate           hipblasLtMatrixLayoutCreate
    #define cublasLtMatmulAlgoGetHeuristic       hipblasLtMatmulAlgoGetHeuristic
    #define cublasLtMatmul                       hipblasLtMatmul
    #define cublasLtMatmulPreferenceDestroy      hipblasLtMatmulPreferenceDestroy
    #define cublasLtMatrixLayoutDestroy          hipblasLtMatrixLayoutDestroy
    #define cublasLtMatmulDescDestroy            hipblasLtMatmulDescDestroy
    #define cublasLtMatrixLayoutSetAttribute     hipblasLtMatrixLayoutSetAttribute

    #define CUBLASLT_MATMUL_DESC_TRANSA                 HIPBLASLT_MATMUL_DESC_TRANSA
    #define CUBLASLT_MATMUL_DESC_TRANSB                 HIPBLASLT_MATMUL_DESC_TRANSB
    #define CUBLASLT_MATMUL_PREF_MAX_WORKSPACE_BYTES    HIPBLASLT_MATMUL_PREF_MAX_WORKSPACE_BYTES
    #define CUBLASLT_MATRIX_LAYOUT_ROWS                 HIPBLASLT_MATRIX_LAYOUT_ROWS
    #define CUBLASLT_MATRIX_LAYOUT_COLS                 HIPBLASLT_MATRIX_LAYOUT_COLS
    #define CUBLASLT_MATRIX_LAYOUT_LD                   HIPBLASLT_MATRIX_LAYOUT_LD
    #define CUBLASLT_MATRIX_LAYOUT_BATCH_COUNT          HIPBLASLT_MATRIX_LAYOUT_BATCH_COUNT
    #define CUBLASLT_MATRIX_LAYOUT_STRIDED_BATCH_OFFSET HIPBLASLT_MATRIX_LAYOUT_STRIDED_BATCH_OFFSET
    #define CUBLASLT_MATMUL_MATRIX_SCALE_SCALAR_32F     HIPBLASLT_MATMUL_MATRIX_SCALE_SCALAR_32F
    #define CUBLASLT_MATMUL_DESC_A_SCALE_MODE           HIPBLASLT_MATMUL_DESC_A_SCALE_MODE
    #define CUBLASLT_MATMUL_DESC_B_SCALE_MODE           HIPBLASLT_MATMUL_DESC_B_SCALE_MODE

    // CUDA
    #define __nv_fp8_storage_t    __hip_fp8_storage_t
    #define __nv_fp8x4_storage_t  __hip_fp8x4_storage_t
    #define cudaError_t           hipError_t
    #define cudaDataType_t        hipDataType
    #define cudaStream_t          hipStream_t
    #define cudaPointerAttributes hipPointerAttribute_t
    #define cudaEvent_t           hipEvent_t
    #define cudaDeviceProp        hipDeviceProp_t

    #define __nv_cvt_float_to_fp8 __hip_cvt_float_to_fp8
    #define __NV_SATFINITE        __HIP_SATFINITE

    #define cudaDeviceSynchronize             hipDeviceSynchronize
    #define cudaStreamCreate                  hipStreamCreate
    #define cudaStreamDestroy                 hipStreamDestroy
    #define cudaStreamSynchronize             hipStreamSynchronize
    #define cudaMemGetInfo                    hipMemGetInfo
    #define cudaGetDeviceProperties           hipGetDeviceProperties
    #define cudaGetDevice                     hipGetDevice
    #define cudaDeviceGetAttribute            hipDeviceGetAttribute
    #define cudaDriverGetVersion              hipDriverGetVersion
    #define cudaRuntimeGetVersion             hipRuntimeGetVersion
    #define cudaDevAttrComputeCapabilityMajor hipDeviceAttributeComputeCapabilityMajor
    #define cudaDevAttrComputeCapabilityMinor hipDeviceAttributeComputeCapabilityMinor
    #define cudaMalloc                        hipMalloc
    #define cudaMallocAsync                   hipMallocAsync
    #define cudaMemcpy                        hipMemcpy
    #define cudaMemcpyAsync                   hipMemcpyAsync
    #define cudaMemcpy2DAsync                 hipMemcpy2DAsync
    #define cudaMemsetAsync                   hipMemsetAsync
    #define cudaMemset2DAsync                 hipMemset2DAsync
    #define cudaMemcpyDeviceToHost            hipMemcpyDeviceToHost
    #define cudaFree                          hipFree
    #define cudaFreeAsync                     hipFreeAsync
    #define cudaSuccess                       hipSuccess
    #define cudaDataType                      hipDataType
    #define cudaGetErrorString                hipGetErrorString
    #define cudaPointerGetAttributes          hipPointerGetAttributes
    #define cudaMemoryTypeHost                hipMemoryTypeHost
    #define cudaMemoryTypeUnregistered        hipMemoryTypeUnregistered
    #define cudaEventCreate                   hipEventCreate
    #define cudaEventCreateWithFlags          hipEventCreateWithFlags
    #define cudaEventRecord                   hipEventRecord
    #define cudaStreamWaitEvent               hipStreamWaitEvent
    #define cudaEventDestroy                  hipEventDestroy
    #define cudaEventDisableTiming            hipEventDisableTiming
    #define cudaSetDevice                     hipSetDevice
    #define cudaMemcpyKind                    hipMemcpyKind
    #define cudaMemcpyHostToDevice            hipMemcpyHostToDevice
    #define cudaMemcpyDeviceToDevice          hipMemcpyDeviceToDevice
    #define cudaMemcpyHostToHost              hipMemcpyHostToHost
    #define cudaMemcpyDefault                 hipMemcpyDefault
    #define cudaGetLastError                  hipGetLastError
    #define cudaEventSynchronize              hipEventSynchronize
    #define cudaEventElapsedTime              hipEventElapsedTime
    #define cudaDeviceReset                   hipDeviceReset

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
    #define FP8_FNUZ 0
    #define STR(x)   #x // No change for cuda
#endif
