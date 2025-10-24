# GEMMul8<!-- omit in toc -->

GEMMul8 (GEMMulate): GEMM emulation using int8 matrix engines based on Ozaki Scheme II

GEMMul8 is a library for emulating high-precision matrix multiplication (SGEMM/DGEMM) using low-precision INT8 matrix engines.
It is based on the Ozaki Scheme II, enabling bit-wise reproducible results with superior performance and power efficiency compared to native floating-point implementations.

- [Technical Overview](#technical-overview)
- [Build](#build)
  - [make options](#make-options)
  - [Example](#example)
    - [CUDA build](#cuda-build)
    - [HIP build](#hip-build)
- [Running Sample Codes](#running-sample-codes)
  - [How to Run](#how-to-run)
  - [make options](#make-options-1)
  - [Example](#example-1)
- [Usage](#usage)
  - [1. Direct Usage (Normal mode)](#1-direct-usage-normal-mode)
    - [Example: run emulation for the CUDA backend](#example-run-emulation-for-the-cuda-backend)
    - [Arguments of `gemmul8::gemm`](#arguments-of-gemmul8gemm)
    - [Behavior of `skip_scalA` / `skip_scalB`](#behavior-of-skip_scala--skip_scalb)
    - [Example: How to skip scaling step](#example-how-to-skip-scaling-step)
  - [2. Hijack cuBLAS/hipBLAS GEMM (Hook Mode)](#2-hijack-cublashipblas-gemm-hook-mode)
    - [Interception target](#interception-target)
    - [How to enable the hook](#how-to-enable-the-hook)
    - [Configure emulation parameters via environment variables](#configure-emulation-parameters-via-environment-variables)
    - [How to change environment variables programmatically](#how-to-change-environment-variables-programmatically)
- [Numerical results](#numerical-results)
  - [Environments](#environments)
  - [Accuracy](#accuracy)
  - [Throughput performance](#throughput-performance)
  - [Power efficiency](#power-efficiency)
  - [Supplementary Experiment on a B200](#supplementary-experiment-on-a-b200)
- [Acknowledgment](#acknowledgment)
  - [Assistance with debugging](#assistance-with-debugging)
  - [Assistance with experiments on a B200 GPU](#assistance-with-experiments-on-a-b200-gpu)
- [Contact (Responsible Developer)](#contact-responsible-developer)
- [References](#references)
- [Citations](#citations)
- [License](#license)

## Technical Overview

GEMMul8 implements GEMM emulation based on Ozaki scheme II, which utilizes the Chinese Remainder Theorem (CRT).
A larger number of moduli (`num_moduli`) for the CRT results in higher precision at the cost of increased computation time.

This project supports both **CUDA** and **HIP** backends for GPU computation.

The GEMM emulation ensures bit-wise numerical reproducibility.
This is achieved through a combination of exact, error-free computations and a fixed order of operations.

However, please note that this reproducibility guarantee may not extend to environments using different toolkits or compilers.
Discrepancies can arise from varying propagation rules for Inf/NaN values or differences in the precision of underlying math library functions.

## Build

Run `make` in the `GEMMUl8/GEMMul8` directory to compile all source files.

```bash
make -j8
```

To rebuild from scratch, run `make clean && make -j8`.

### make options

| Option          | Default           | Description                                                                                           |
| :-------------- | :---------------- | :---------------------------------------------------------------------------------------------------- |
| `CUDA_PATH`     | `/usr/local/cuda` | Path to your CUDA toolkit installation. Used for CUDA backends.                                       |
| `HIP_PATH`      | `/opt/rocm`       | Path to your HIP (ROCm) toolkit installation. Used for HIP backends.                                  |
| `BACKEND`       | `auto`            | Select GPU backend: `cuda`, `hip`, or `auto` (auto-detect).                                           |
| `GPU_ARCH`      | `auto`            | Target GPU architecture.<br>Examples: `80` (A100), `90` (H100), `gfx90a` (MI250X), `gfx942` (MI300X). |
| `ozIMMU_EF`     | `no`              | Set to `yes` to use ozIMMU in the sample code.                                                        |
| `ozIMMU_EF_DIR` | `../../ozIMMU_EF` | Path to ozIMMU implementation. Used if `ozIMMU_EF=yes`.                                               |
| `cuMpSGEMM`     | `no`              | Set to `yes` to use cuMpSGEMM in the sample code.                                                     |
| `cuMpSGEMM_DIR` | `../../cuMpSGEMM` | Path to cuMpSGEMM. Used if `cuMpSGEMM=yes`.                                                           |
| `HIJACK`        | `no`              | Set to `yes` to hijack cuBLAS GEMM calls with emulation in the sample code.                           |

> [!NOTE]
>
> - If you enable optional modules (ozIMMU_EF=yes, cuMpSGEMM=yes), please clone and build their repositories first.
>   - [cuMpSGEMM - CUDA Mutable-precision SGEMM](https://github.com/enp1s0/cuMpSGEMM)
>   - [ozIMMU - DGEMM on Int8 Tensor Core](https://github.com/enp1s0/ozIMMU)
>   - See also [Accelerator for ozIMMU](https://github.com/RIKEN-RCCS/accelerator_for_ozIMMU)
>   - If you use these, please ensure compliance with their respective license terms.
> - `BACKEND=auto` will attempt to detect your GPU vendor automatically.
> - `GPU_ARCH=auto` will automatically detect and use the appropriate compute capability or architecture for your GPU.
> - Target GPU architecture can be found from e.g., [NVIDIA GPU CC](https://developer.nvidia.com/cuda-gpus) or [AMD hardware specs](https://rocm.docs.amd.com/en/latest/reference/gpu-arch-specs.html).

### Example

#### CUDA build

Build for an NVIDIA H100/H200 GPU (Compute Capability 9.0)

```bash
make -j8 BACKEND=cuda CUDA_PATH=/usr/local/cuda GPU_ARCH=90
```

#### HIP build

Build for an AMD MI300X GPU (gfx942 architecture)

```bash
make -j8 BACKEND=hip HIP_PATH=/opt/rocm GPU_ARCH=gfx942
```

## Running Sample Codes

After building the project, you can run sample codes from the testing directory.

### How to Run

Navigate to the `testing` directory and use `make` to run tests for different precisions.

### make options

| Mode     | Value        | Description                                                     |
| :------- | :----------- | :-------------------------------------------------------------- |
| `test_f` | (no)         | Tests for SGEMM                                                 |
| `test_d` | (no)         | Tests for DGEMM                                                 |
| `MODE`   | `accuracy`   | Tests numerical accuracy (maximum element-wise relative error). |
|          | `flops`      | Measures TFLOPS with square matrices.                           |
|          | `watt`       | Measures watt and GFLOPS/watt with square matrices.             |
|          | `flops_rect` | Measures TFLOPS with rectangular matrices.                      |
|          | `watt_rect`  | Measures watt and GFLOPS/watt with rectangular matrices.        |
|          | `all`        | Runs all available test modes.                                  |

### Example

```bash
# Run accuracy test (DGEMM)
make test_d MODE="accuracy"
```

```bash
# Run all tests (SGEMM)
make test_f MODE="all"
```

```bash
# Run accuracy, flops, & watt tests (SGEMM & DGEMM)
make test_f test_d MODE="accuracy flops watt"
```

## Usage

This library provides two ways to use the GEMM emulation.

### 1. Direct Usage (Normal mode)

Call GEMMul8 functions explicitly from your source code.
This gives you fine-grained control over the emulation parameters.

#### Example: run emulation for the CUDA backend

```cpp
#include "gemmul8.hpp"

// 1. Create a handle to the cuBLAS library context
cublasHandle_t cublas_handle;
cublasCreate(&cublas_handle);

// 2. Settings
const unsigned num_moduli = 14u;  // Accuracy knob: 2 <= num_moduli <= 20
const bool fastmode = true;       // true (fast mode) or false (accurate mode)

// 3. Allocate workspace
const size_t worksize = gemmul8::workSize(m, n, k, num_moduli); // calculate required memory (Byte)
void *work;
cudaMalloc(&work, worksize);

// 4. (Optional) Create a vector to store timing breakdown
std::vector<double> time_breakdown(4, 0.0);

// 5. Run emulation
// The function returns a vector with execution times for each phase.
time_breakdown = gemmul8::gemm(cublas_handle,
                               CUBLAS_OP_N, CUBLAS_OP_N,
                               m, n, k,
                               &alpha, devA, lda,
                               devB, ldb,
                               &beta, devC, ldc,
                               num_moduli, fastmode, work);

// 6. Free workspace
cudaFree(work);

// 7. Destroy a handle
cublasDestroy(cublas_andle);
```

#### Arguments of `gemmul8::gemm`

The arguments for `gemmul8::gemm` closely match the standard [cublas&lt;t&gt;gemm](https://docs.nvidia.com/cuda/cublas/#cublas-t-gemm) and [hipblasXgemm](https://rocm.docs.amd.com/projects/hipBLAS/en/develop/reference/hipblas-api-functions.html#hipblasxgemm-batched-stridedbatched) interfaces, with additional parameters that control internal preprocessing and workspace reuse.

GEMMul8 provides both a workspace query function (`workSize`) and an extended GEMM computation interface (`gemm`), as shown below:

```cpp
namespace gemmul8 {

// workSize returns the required workspace size in bytes.
size_t workSize(
    const size_t m,                         // Number of rows of C
    const size_t n,                         // Number of columns of C
    const size_t k,                         // Inner dimension <= 2^17
    const unsigned num_moduli,              // #moduli, 2 <= num_moduli <= 20
    const bool enable_skip_scalA = false,   // [option] Reserve extra space for A to allow skip_scalA
    const bool enable_skip_scalB = false,   // [option] Reserve extra space for B to allow skip_scalB
    size_t *workSizeA            = nullptr, // [option] Output: workspace size used for A8i and sftA
    size_t *workSizeB            = nullptr  // [option] Output: workspace size used for B8i and sftB
);

// gemm returns computation time in second of each computational phase
// gemm returns computation time in second of each computational phase
#if defined(__NVCC__)
template <typename T> std::vector<double> gemm(
    cublasHandle_t handle,                  // Handle to the cuBLAS library context
    const cublasOperation_t op_A,           // CUBLAS_OP_N or CUBLAS_OP_T
    const cublasOperation_t op_B,           // CUBLAS_OP_N or CUBLAS_OP_T
    const size_t m,                         // Number of rows of C
    const size_t n,                         // Number of columns of C
    const size_t k,                         // Inner dimension <= 2^17
    const T *alpha,                         // Scaling factor for op(A)*op(B)
    const T *const A,                       // 1-D device array of dimensions lda*k (CUBLAS_OP_N) or lda*m (CUBLAS_OP_T)
    const size_t lda,                       // Leading dimension of A
    const T *const B,                       // 1-D device array of dimensions ldb*n (CUBLAS_OP_N) or ldb*k (CUBLAS_OP_T)
    const size_t ldb,                       // Leading dimension of B
    const T *beta,                          // Scaling factor for C
    T *const C,                             // 1-D device array of dimensions ldc*n
    const size_t ldc,                       // Leading dimension of C
    const unsigned num_moduli,              // #moduli, 2 <= num_moduli <= 20
    const bool fastmode,                    // false (accurate mode) or true (fast mode)
    void *const work,                       // Preallocated workspace
    void *const workA            = nullptr, // [optional] Separate workspace for A (if nullptr, uses work)
    void *const workB            = nullptr, // [optional] Separate workspace for B (if nullptr, uses work)
    const bool enable_skip_scalA = false,   // [optional] Enables scaling-skip mechanism for A
    const bool enable_skip_scalB = false,   // [optional] Enables scaling-skip mechanism for B
    bool skip_scalA              = false,   // [optional] If true, skip preprocessing for A
    bool skip_scalB              = false    // [optional] If true, skip preprocessing for B
);
#endif

#if defined(__HIPCC__)
template <typename T> std::vector<double> gemm(
    hipblasHandle_t handle,                 // Handle to the cuBLAS library context
    const hipblasOperation_t op_A,          // CUBLAS_OP_N or CUBLAS_OP_T
    const hipblasOperation_t op_B,          // CUBLAS_OP_N or CUBLAS_OP_T
    const size_t m,                         // Number of rows of C
    const size_t n,                         // Number of columns of C
    const size_t k,                         // Inner dimension <= 2^17
    const T *alpha,                         // Scaling factor for op(A)*op(B)
    const T *const A,                       // 1-D device array of dimensions lda*k (CUBLAS_OP_N) or lda*m (CUBLAS_OP_T)
    const size_t lda,                       // Leading dimension of A
    const T *const B,                       // 1-D device array of dimensions ldb*n (CUBLAS_OP_N) or ldb*k (CUBLAS_OP_T)
    const size_t ldb,                       // Leading dimension of B
    const T *beta,                          // Scaling factor for C
    T *const C,                             // 1-D device array of dimensions ldc*n
    const size_t ldc,                       // Leading dimension of C
    const unsigned num_moduli,              // #moduli, 2 <= num_moduli <= 20
    const bool fastmode,                    // false (accurate mode) or true (fast mode)
    void *const work,                       // Preallocated workspace
    void *const workA            = nullptr, // [optional] Separate workspace for A (if nullptr, uses work)
    void *const workB            = nullptr, // [optional] Separate workspace for B (if nullptr, uses work)
    const bool enable_skip_scalA = false,   // [optional] Enables scaling-skip mechanism for A
    const bool enable_skip_scalB = false,   // [optional] Enables scaling-skip mechanism for B
    bool skip_scalA              = false,   // [optional] If true, skip preprocessing for A
    bool skip_scalB              = false    // [optional] If true, skip preprocessing for B
);
#endif

} // namespace gemmul8
```

#### Behavior of `skip_scalA` / `skip_scalB`

- `gemmul8::gemm` internally converts the input matrices `A` and `B` into INT8 format and performs modular multiplications across multiple moduli.
- The conversion (scaling) step can be **skipped** in consecutive GEMM calls if the same matrices are reused, allowing for substantial performance gains.
- If `enable_skip_scal{A|B} = true`, additional workspace is preallocated so that the scaled INT8 representation of `A`/`B` can be retained between calls.
- If both `enable_skip_scal{A|B} && skip_scal{A|B} = true`, the scaling step for `A`/`B` is **actually skipped**, and previously prepared data are reused for faster execution.
  This mode assumes that the contents of `A`/`B` in device memory remain unchanged.
- This mechanism is particularly effective when repeatedly multiplying the same `A` (or `B`) with different `B` (or `A`) matrices.

> [!NOTE]
> When using `skip_scalA` / `skip_scalB`, the preprocessing step that converts `A`/`B` into its internal INT8 representation is skipped.
> For correctness, the following conditions **must all hold** between consecutive GEMM calls:
>
> 1. The dimensions (`M`/`K` for `A`, `K`/`N` for `B`) must be identical to those in the previous call.
> 2. The operation type (`CUBLAS_OP_N` / `CUBLAS_OP_T`) for `A`/`B` must be the same as before.
> 3. The value of `num_moduli` must remain unchanged.
> 4. The `fastmode` setting must be identical to that of the previous call.
> 5. The contents of `A`/`B` in device memory must not be modified between calls.

> [!CAUTION]
> If any of these conditions differ, the cached scaled data become invalid, and skipping must **not** be used.
> In such cases, set `skip_scalA=false` / `skip_scalB=false`.

> [!CAUTION]
> This skip mechanism is designed for repeated GEMM calls with identical A or B.
> Use it only when you are certain that the input matrices and configuration have not changed.
> When in doubt, disable skipping to ensure correctness.

#### Example: How to skip scaling step

```cpp
#include "gemmul8.hpp"

// 1. Create a handle to the cuBLAS library context
cublasHandle_t cublas_handle;
cublasCreate(&cublas_handle);

// 2. Settings
const unsigned num_moduli = 14u;
const bool fastmode = false;

// 3. Matrix shapes
const size_t m1 = 64, n1 = 10, k1 = 10; // 1st GEMM: 64×10 × 10×10
const size_t m2 = 20, n2 = 10, k2 = 10; // 2nd GEMM: 20×10 × 10×10

// 4. Allocate workspace
size_t worksizeA, worksizeB;
const size_t worksize = gemmul8::workSize(
    std::max(m1, m2), std::max(n1, n2), std::max(k1, k2),
    num_moduli, false, true, &worksizeA, &worksizeB);

void *work;
cudaMalloc(&work, worksize);

int8_t *workA = reinterpret_cast<int8_t *>(work); // fixed workspace for A
int8_t *workB = workA + offsetA;                  // fixed workspace for B
int8_t *work_rem = workB + offsetB;               // remaining workspace

// 5. Run GEMM (first call: scaling performed)
gemmul8::gemm(cublas_handle,
              CUBLAS_OP_N, CUBLAS_OP_N,
              m1, n1, k1,
              &alpha, devA1, lda,
              devB, ldb,
              &beta, devC, ldc,
              num_moduli, fastmode, (void*)work_rem, (void*)workA, (void*)workB,
              false, true, false, false);

// 6. Reuse scaled A (second call: skip_scalB = true)
gemmul8::gemm(cublas_handle,
              CUBLAS_OP_N, CUBLAS_OP_N,
              m1, n1, k1,
              &alpha, devA1, lda,
              devB, ldb,
              &beta, devC, ldc,
              num_moduli, fastmode, (void*)work_rem, (void*)workA, (void*)workB,
              false, true, false, true);

// 7. Free workspace
cudaFree(work);

// 8. Destroy a handle
cublasDestroy(cublas_andle);
```

### 2. Hijack cuBLAS/hipBLAS GEMM (Hook Mode)

Intercept standard GEMM calls automatically without modifying the application source code.

#### Interception target

- `cublasSgemm`, `cublasDgemm`, `cublasSgemm_v2`, `cublasDgemm_v2`, `cublasGemmEx`
- `hipblasSgemm`, `hipblasDgemm`, `hipblasGemmEx`, `hipblasGemmEx_v2`

#### How to enable the hook

1. Build the library.
2. Set the `LD_PRELOAD` environment variable.

```bash
export LD_PRELOAD=/path-to-GEMMul8/lib/libgemmul8.so
```

3. Run your application.

#### Configure emulation parameters via environment variables

```bash
export GEMMUL8_NUM_MOD_D=15
export GEMMUL8_NUM_MOD_S=7
export GEMMUL8_FASTMODE_D=1
export GEMMUL8_FASTMODE_S=0
export GEMMUL8_MAX_M=4096
export GEMMUL8_MAX_N=4096
export GEMMUL8_MAX_K=4096
export GEMMUL8_MAX_NUM_MOD=18
export GEMMUL8_SKIP_SCALE_A=1
export GEMMUL8_SKIP_SCALE_B=1
```

| Variable               | Default | Applies to | Description                                                                          |
| :--------------------- | :------ | :--------- | :----------------------------------------------------------------------------------- |
| `GEMMUL8_NUM_MOD_D`    | `2`     | DGEMM      | Number of moduli (`unsigned num_moduli`) used in DGEMM emulation. Controls accuracy. |
| `GEMMUL8_NUM_MOD_S`    | `2`     | SGEMM      | Number of moduli (`unsigned num_moduli`) used in SGEMM emulation. Controls accuracy. |
| `GEMMUL8_FASTMODE_D`   | `0`     | DGEMM      | Enables fast mode (`1` = fast mode, `0` = accurate mode).                            |
| `GEMMUL8_FASTMODE_S`   | `0`     | SGEMM      | Enables fast mode (`1` = fast mode, `0` = accurate mode).                            |
| `GEMMUL8_MAX_M`        | `0`     | both       | Maximum value of `M` used to preallocate workspace memory.                           |
| `GEMMUL8_MAX_N`        | `0`     | both       | Maximum value of `N` used to preallocate workspace memory.                           |
| `GEMMUL8_MAX_K`        | `0`     | both       | Maximum value of `K` used to preallocate workspace memory.                           |
| `GEMMUL8_MAX_NUM_MOD`  | `2`     | both       | Maximum number of moduli used when computing the size of the preallocated workspace. |
| `GEMMUL8_SKIP_SCALE_A` | `0`     | both       | Enables skipping redundant preprocessing for `A` (`1` = enable, `0` = disable).      |
| `GEMMUL8_SKIP_SCALE_B` | `0`     | both       | Enables skipping redundant preprocessing for `B` (`1` = enable, `0` = disable).      |

- This hook mode preallocates a single large GPU workspace on demand and resizes as needed.
- Each `cublasHandle_t` maintains an independent workspace.
- The workspace is automatically released when the corresponding handle is destroyed via `cublasDestroy()` or `cublasDestroy_v2()`.
- Workspace size is determined as `max(wsmax, ws, pre_ws)`, where
  - `wsmax  = workSize(GEMMUL8_MAX_M, GEMMUL8_MAX_N, GEMMUL8_MAX_K, GEMMUL8_MAX_NUM_MOD)`
  - `ws     = workSize(m, n, k, GEMMUL8_NUM_MOD_{D|S})`
  - `pre_ws =` the size of the previously allocated workspace for the same `cublasHandle_t`
- If `GEMMUL8_MAX_{M|N|K|NUM_MOD}` variables are set appropriately, the workspace size becomes fixed to `wsmax`, avoiding costly reallocations during subsequent GEMM calls.
- The workspace is **never shrunk** automatically; it only grows when a larger allocation is required.
- When `GEMMUL8_SKIP_SCALE_{A|B}=1`, redundant preprocessing for `A`/`B` is skipped (see below).

> [!IMPORTANT]
>
> - `GEMMUL8_SKIP_SCALE_{A|B}=1` allows consecutive GEMM calls using the same matrix pointers (`A`/`B`) to reuse already-scaled intermediate data.
> - Automatic skip of conversion from `A` or `B` to INT8 is enabled when:
>   1. `GEMMUL8_SKIP_SCALE_{A|B}=1`.
>   2. The same computation mode (fast mode or accurate mode) is used for both the previous and current calls.
>   3. The same `cublasHandle_t` is used for both the previous and current calls.
>   4. The same device pointer (`A` or `B`) is used for both the previous and current calls.
>   5. The same matrix shapes (`M`/`K`/`lda` for `A`, `K`/`N`/`ldb` for `B`) are used for both the previous and current calls.
>   6. The same operation flag (`transa` for `A`, `transb` for `B`) is used for both the previous and current calls.
>   7. The same number of moduli (`GEMMUL8_NUM_MOD_{D|S}`) is used for both the previous and current calls.
>   8. The same workspace size (`ws`) is used for both the previous and current calls.

> [!TIP]
> To ensure the last condition is met, it is recommended to set the `GEMMUL8_MAX_{M|N|K|NUM_MOD}` variables.
> This fixes the workspace size to wsmax, preventing it from changing between calls.

> [!CAUTION]
> ⚠️Note: Skip scaling assumes that the matrices `A` or `B` remain unchanged in GPU memory.  
> If their contents are modified between GEMM calls, set `GEMMUL8_SKIP_SCALE_{A|B}=0` to ensure correctness.

#### How to change environment variables programmatically

You can also set these environment variables programmatically from within your code using setenv.

```cpp
char num_moduli[12];

// Run emulation with num_moduli = 15 & fastmode = true
snprintf(num_moduli, sizeof(num_moduli), "%u", 15u);
setenv("GEMMUL8_NUM_MOD_D", num_moduli, 1);
setenv("GEMMUL8_FASTMODE_D", "1", 1);
cublasDgemm_v2(...);

// Run emulation with num_moduli = 18 & fastmode = false
snprintf(num_moduli, sizeof(num_moduli), "%u", 18u);
setenv("GEMMUL8_NUM_MOD_D", num_moduli, 1);
setenv("GEMMUL8_FASTMODE_D", "0", 1);
cublasDgemm_v2(...);
```

## Numerical results

The constant $\phi$ controls the difficulty of matrix multiplication (exponent distribution of input matrices).
The difficulty of $\phi = 0.5$ is comparable to that of matrix multiplication in HPL.

See all numerical results in the separate repository: [GEMMul8_numerical_results](https://github.com/UCHINO-Yuki/GEMMul8_numerical_results)

### Environments

- NVIDIA A100 SXM4 GPU with an AMD EPYC 7713 CPU , CUDA Toolkit 12.9.86, and gcc 11.5.0
- NVIDIA GH200 Grace Hopper Superchip with CUDA Toolkit 12.9.86, and gcc 11.5.0
- NVIDIA GeForce RTX 5080 GPU with an AMD Ryzen 9 7950X CPU, CUDA Toolkit 12.9.86, and gcc 11.5.0
- NVIDIA B200 GPU with an Intel Xeon Platinum 8570 CPU, CUDA Toolkit 12.8.93, and gcc 13.3.0

### Accuracy

![accuracy_dgemm](https://raw.githubusercontent.com/UCHINO-Yuki/GEMMul8_numerical_results/main/results_02/fig/oz2_results_df_accuracy.png)
_Accuracy of DGEMM (top) and SGEMM (bottom) emulation for $m=n=1024$ on GH200. Solid lines represent results for $k=1024$, and dashed lines for $k=16384$._

### Throughput performance

![throughput_dgemm](https://raw.githubusercontent.com/UCHINO-Yuki/GEMMul8_numerical_results/main/results_02/fig/oz2_results_d_time.png)
_Throughput performance of DGEMM emulation on A100 (top), GH200 (middle), and RTX 5080 (bottom)._

![throughput_sgemm](https://raw.githubusercontent.com/UCHINO-Yuki/GEMMul8_numerical_results/main/results_02/fig/oz2_results_f_time.png)
_Throughput performance of SGEMM emulation on A100 (top), GH200 (middle), and RTX 5080 (bottom)._

### Power efficiency

![power_dgemm](https://raw.githubusercontent.com/UCHINO-Yuki/GEMMul8_numerical_results/main/results_02/fig/oz2_results_d_watt.png)
_Power efficiency of DGEMM emulation on A100 (top), GH200 (middle), and RTX 5080 (bottom)._

![power_sgemm](https://raw.githubusercontent.com/UCHINO-Yuki/GEMMul8_numerical_results/main/results_02/fig/oz2_results_f_watt.png)
_Power efficiency of SGEMM emulation on A100 (top), GH200 (middle), and RTX 5080 (bottom)._

### Supplementary Experiment on a B200

- $m=n=k=16384$
- DGEMM:
  - Native DGEMM: 36.0 TFLOPS
  - Emulation (Fast Mode, 14-17 moduli): 98.9-123 TFLOPS
  - Emulation (Accurate Mode, 14-17 moduli): 92.5-110 TFLOPS
- SGEMM:
  - Native SGEMM: 67.2 TFLOPS
  - Emulation (Fast Mode, 7-9 moduli): 197--259 TFLOPS
  - Emulation (Accurate Mode, 6-8 moduli): 194--240 TFLOPS

## Acknowledgment

> [!CAUTION]
> Please do not contact the individuals listed below regarding this code.

### Assistance with debugging

- Patrick Gutsche (École Normale Supérieure de Lyon, France)
- Prajval Kumar (Indian Institute of Science and Education Research, India)
- Dr. William Dawson (RIKEN Center for Computational Science, Japan)
- Dr. Toshiyuki Imamura (RIKEN Center for Computational Science, Japan)

(Affiliations as of 2025)

### Assistance with experiments on a B200 GPU

- Dr. Qianxiang Ma (RIKEN Center for Computational Science, Japan)
- Prof. Rio Yokota (Institute of Science Tokyo, Japan)

(Affiliations as of 2025)

## Contact (Responsible Developer)

- Yuki Uchino (RIKEN Center for Computational Science, Japan)
- yuki.uchino.fe (at) riken.jp

## References

- Ootomo, H., & Yokota, R. (2022). Recovering single precision accuracy from Tensor Cores while surpassing the FP32 theoretical peak performance. The International Journal of High Performance Computing Applications, 36(4), 475-491, [doi.org/10.1177/10943420221090256](https://doi.org/10.1177/10943420221090256).
- Ootomo, H., Manabe, H., Harada, K., & Yokota, R. (2023). Quantum Circuit Simulation by SGEMM Emulation on Tensor Cores and Automatic Precision Selection. In High Performance Computing (pp. 259-276). Springer, [doi.org/10.1007/978-3-031-32041-5_14](https://doi.org/10.1007/978-3-031-32041-5_14).
- Ootomo, H., Ozaki, K., & Yokota, R. (2024). DGEMM on integer matrix multiplication unit. The International Journal of High Performance Computing Applications, 38(4), 297-313, [https://doi.org/10.1177/10943420241239588](https://doi.org/10.1177/10943420241239588).
- Uchino, Y., Ozaki, K., & Imamura, T. (2025). Performance enhancement of the Ozaki Scheme on integer matrix multiplication unit. The International Journal of High Performance Computing Applications, 39(3), 462-476, [doi.org/10.1177/10943420241313064](https://doi.org/10.1177/10943420241313064).

## Citations

```
@misc{ozaki2025ozakischemeiigemmoriented,
      title={Ozaki Scheme II: A GEMM-oriented emulation of floating-point matrix multiplication using an integer modular technique},
      author={Katsuhisa Ozaki and Yuki Uchino and Toshiyuki Imamura},
      year={2025},
      eprint={2504.08009},
      archivePrefix={arXiv},
      primaryClass={cs.MS},
      url={https://arxiv.org/abs/2504.08009},
}

@misc{uchino2025highperformancepowerefficientemulationmatrix,
      title={High-Performance and Power-Efficient Emulation of Matrix Multiplication using INT8 Matrix Engines},
      author={Yuki Uchino and Katsuhisa Ozaki and Toshiyuki Imamura},
      year={2025},
      eprint={2508.03984},
      archivePrefix={arXiv},
      primaryClass={cs.DC},
      url={https://arxiv.org/abs/2508.03984},
}

and

TBD
```

## License

This project is licensed under the MIT License. See the LICENSE file for details.
