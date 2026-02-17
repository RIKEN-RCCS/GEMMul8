# GEMMul8<!-- omit in toc -->

GEMMul8 (GEMMulate): GEMM emulation using INT8/FP8 matrix engines based on the Ozaki Scheme II

GEMMul8 is a library for emulating high-precision matrix multiplication (SGEMM, DGEMM, CGEMM, and ZGEMM) using low-precision matrix engines, including INT8 and FP8.
The library is based on the Ozaki Scheme II and supports selectable INT8- or FP8-based emulation backends within each GEMM routine.
This design enables bit-wise reproducible results while achieving superior performance and power efficiency compared to native floating-point implementations.

- [Technical Overview](#technical-overview)
- [Build](#build)
  - [make options](#make-options)
  - [Example](#example)
    - [CUDA build](#cuda-build)
    - [HIP build](#hip-build)
- [Running Test Codes](#running-test-codes)
  - [How to Run](#how-to-run)
  - [MODE options](#mode-options)
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

GEMMul8 supports two emulation backends:

- INT8 backend: uses standard BLAS handle (cuBLAS/hipBLAS handle) or Lt handle (cuBLASLt/hipBLASLt handle).
- FP8 backend: uses Lt handle (cuBLASLt/hipBLASLt handle).

## Build

Run `make` in the `GEMMUl8/GEMMul8` directory to compile all source files.

```bash
make -j8
```

To rebuild from scratch, run `make clean && make -j8`.

### make options

| Option      | Default           | Description                                                                                           |
| :---------- | :---------------- | :---------------------------------------------------------------------------------------------------- |
| `CUDA_PATH` | `/usr/local/cuda` | Path to your CUDA toolkit installation. Used for CUDA backends.                                       |
| `HIP_PATH`  | `/opt/rocm`       | Path to your HIP (ROCm) toolkit installation. Used for HIP backends.                                  |
| `BACKEND`   | `auto`            | Select GPU backend: `cuda`, `hip`, or `auto` (auto-detect).                                           |
| `GPU_ARCH`  | `auto`            | Target GPU architecture.<br>Examples: `80` (A100), `90` (H100), `gfx90a` (MI250X), `gfx942` (MI300X). |

> [!NOTE]
>
> - `BACKEND=auto` will attempt to detect your GPU vendor automatically.
> - `GPU_ARCH=auto` will automatically detect and use the appropriate compute capability or architecture for your GPU.
> - Target GPU architecture can be found from e.g., [CUDA GPU Compute Capability](https://developer.nvidia.com/cuda-gpus) or [AMD GPU hardware specifications](https://rocm.docs.amd.com/en/latest/reference/gpu-arch-specs.html).

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

## Running Test Codes

After building the project, you can run test codes from the testing directory.

### How to Run

Navigate to the `testing` directory and use `make run MODE="test_modes"` to run tests for different precisions.

### MODE options

| Value      | Description                                                     |
| :--------- | :-------------------------------------------------------------- |
| `SGEMM`    | Runs SGEMM tests.                                               |
| `DGEMM`    | Runs DGEMM tests.                                               |
| `CGEMM`    | Runs CGEMM tests.                                               |
| `ZGEMM`    | Runs ZGEMM tests.                                               |
| `accuracy` | Tests numerical accuracy (maximum element-wise relative error). |
| `flops`    | Measures TFLOPS with square matrices.                           |
| `watt`     | Measures watt and GFLOPS/watt with square matrices.             |
| `all`      | Runs `accuracy`, `flops`, and `watt`.                           |

### Example

```bash
# Run accuracy test (DGEMM)
make run MODE="DGEMM accuracy"
```

```bash
# Run accuracy, flops, & watt tests (SGEMM)
make run MODE="SGEMM all"
```

```bash
# Run accuracy & flops tests (SGEMM & DGEMM)
make run MODE="SGEMM DGEMM accuracy flops"
```

## Usage

This library provides two ways to use the GEMM emulation.

### 1. Direct Usage (Normal mode)

Call GEMMul8 functions explicitly from your source code.
This gives you fine-grained control over the emulation parameters.

#### Example: run emulation for the CUDA backend

See the sample code in `sample/`.

#### Arguments of `gemmul8::gemm`

The arguments for `gemmul8::gemm` closely match the standard [cublas&lt;t&gt;gemm](https://docs.nvidia.com/cuda/cublas/#cublas-t-gemm) and [hipblasXgemm](https://rocm.docs.amd.com/projects/hipBLAS/en/develop/reference/hipblas-api-functions.html#hipblasxgemm-batched-stridedbatched) interfaces, with additional parameters that control internal preprocessing and workspace reuse.

GEMMul8 provides both a workspace query function (`workSize`) and an extended GEMM computation interface (`gemm`).
See `include/gemmul8.hpp` for the full function signatures:

- `gemmul8::workSize<is_Complex, Backend>(...)`
- `gemmul8::gemm<T, Backend>(...)` (BLAS handle / Lt handle overloads)

#### Behavior of `skip_scalA` / `skip_scalB`

- `gemmul8::gemm` internally preprocesses the input matrices `A` and `B` into a backend-specific low-precision representation (INT8 or FP8) and performs modular multiplications across multiple moduli.
- This preprocessing step can be **skipped** in consecutive GEMM calls if the same matrices are reused, allowing for substantial performance gains.
- If `enable_skip_scal{A|B} = true`, additional workspace is reserved so that the preprocessed representation of `A`/`B` can be retained between calls.
- If both `enable_skip_scal{A|B} && skip_scal{A|B} = true`, the preprocessing step for `A`/`B` is **actually skipped**, and previously prepared data are reused for faster execution.
  This mode assumes that the contents of `A`/`B` in device memory remain unchanged.
- This mechanism is particularly effective when repeatedly multiplying the same `A` (or `B`) with different `B` (or `A`) matrices.

> [!NOTE]
> When using `skip_scalA` / `skip_scalB`, the preprocessing step that converts `A`/`B` into an internal backend-specific representation (INT8/FP8) is skipped.
> For correctness, the following conditions **must all hold** between consecutive GEMM calls:
>
> 1. The dimensions (`M`/`K` for `A`, `K`/`N` for `B`) must be identical to those in the previous call.
> 2. The operation type (`CUBLAS_OP_N` / `CUBLAS_OP_T`) for `A`/`B` must be the same as before.
> 3. The value of `num_moduli` must remain unchanged.
> 4. The `fastmode` setting must be identical to that of the previous call.
> 5. The contents of `A`/`B` in device memory must not be modified between calls.
> 6. The selected emulation backend must be identical in both calls (INT8 or FP8).

> [!NOTE]
>
> - GEMMul8 does not verify the contents of `A`/`B` when skipping; correctness requires the user to ensure immutability.
> - `alpha` / `beta` follow the BLAS handle pointer mode (host or device). GEMMul8 does not override the pointer mode.

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
    num_moduli, true, true, &worksizeA, &worksizeB);

void *work;
cudaMalloc(&work, worksize);

// NOTE: worksizeA/worksizeB are byte sizes for the dedicated A/B work areas.
const size_t offsetA = worksizeA;
const size_t offsetB = worksizeB;

int8_t *workA    = reinterpret_cast<int8_t *>(work); // dedicated workspace for A
int8_t *workB    = workA + offsetA;                  // dedicated workspace for B
int8_t *work_rem = workB + offsetB;                  // remaining workspace

// 5. Run GEMM (first call: preprocessing performed)
gemmul8::gemm(cublas_handle,
              CUBLAS_OP_N, CUBLAS_OP_N,
              m1, n1, k1,
              &alpha, devA1, lda,
              devB, ldb,
              &beta, devC, ldc,
              num_moduli, fastmode, (void*)work_rem, (void*)workA, (void*)workB,
              true, true, false, false);

// 6. Reuse preprocessed B (second call: skip_scalB = true)
gemmul8::gemm(cublas_handle,
              CUBLAS_OP_N, CUBLAS_OP_N,
              m1, n1, k1,
              &alpha, devA1, lda,
              devB, ldb,
              &beta, devC, ldc,
              num_moduli, fastmode, (void*)work_rem, (void*)workA, (void*)workB,
              true, true, false, true);

// 7. Free workspace
cudaFree(work);

// 8. Destroy a handle
cublasDestroy(cublas_handle);
```

### 2. Hijack cuBLAS/hipBLAS GEMM (Hook Mode)

Intercept standard GEMM calls automatically without modifying the application source code.

#### Interception target

- `cublasSgemm`, `cublasDgemm`, `cublasCgemm`, `cublasZgemm`
- `cublasSgemm_v2`, `cublasDgemm_v2`, `cublasCgemm_v2`, `cublasZgemm_v2`
- `cublasGemmEx`
- `hipblasSgemm`, `hipblasDgemm`, `hipblasCgemm`, `hipblasZgemm`
- `hipblasGemmEx`, `hipblasGemmEx_v2`

#### How to enable the hook

1. Build the library.
2. Set the `LD_PRELOAD` environment variable.

```bash
export LD_PRELOAD=/path-to-GEMMul8/lib/libgemmul8.so
```

3. Run your application.

#### Configure emulation parameters via environment variables

```bash
export GEMMUL8_BACKEND=INT8
export GEMMUL8_NUM_MOD_D=15
export GEMMUL8_NUM_MOD_S=7
export GEMMUL8_NUM_MOD_Z=15
export GEMMUL8_NUM_MOD_C=7
export GEMMUL8_FASTMODE_D=0
export GEMMUL8_FASTMODE_S=1
export GEMMUL8_FASTMODE_Z=0
export GEMMUL8_FASTMODE_C=1
export GEMMUL8_MAXWS_BACKEND=INT8
export GEMMUL8_MAX_M=4096
export GEMMUL8_MAX_N=4096
export GEMMUL8_MAX_K=4096
export GEMMUL8_MAX_NUM_MOD=18
export GEMMUL8_SKIP_SCALE_A=1
export GEMMUL8_SKIP_SCALE_B=1
```

| Variable                | Default | Applies to | Description                                                                                               |
| :---------------------- | :------ | :--------- | :-------------------------------------------------------------------------------------------------------- |
| `GEMMUL8_BACKEND`       | `0`     | all        | Selects the emulation backend (`0` or `INT8` = INT8-based emulation, `1` or `FP8` = FP8-based emulation). |
| `GEMMUL8_NUM_MOD_D`     | `0`     | DGEMM      | Number of moduli used in DGEMM emulation. When num_moduli < 2 or 20 < num_moduli, native DGEMM is used.   |
| `GEMMUL8_NUM_MOD_S`     | `0`     | SGEMM      | Number of moduli used in SGEMM emulation. When num_moduli < 2 or 13 < num_moduli, native SGEMM is used.   |
| `GEMMUL8_NUM_MOD_Z`     | `0`     | ZGEMM      | Number of moduli used in ZGEMM emulation. When num_moduli < 2 or 20 < num_moduli, native ZGEMM is used.   |
| `GEMMUL8_NUM_MOD_C`     | `0`     | CGEMM      | Number of moduli used in CGEMM emulation. When num_moduli < 2 or 13 < num_moduli, native CGEMM is used.   |
| `GEMMUL8_FASTMODE_D`    | `0`     | DGEMM      | Enables fast mode (`1` = fast mode, `0` = accurate mode).                                                 |
| `GEMMUL8_FASTMODE_S`    | `0`     | SGEMM      | Enables fast mode (`1` = fast mode, `0` = accurate mode).                                                 |
| `GEMMUL8_FASTMODE_Z`    | `0`     | ZGEMM      | Enables fast mode (`1` = fast mode, `0` = accurate mode).                                                 |
| `GEMMUL8_FASTMODE_C`    | `0`     | CGEMM      | Enables fast mode (`1` = fast mode, `0` = accurate mode).                                                 |
| `GEMMUL8_MAXWS_BACKEND` | `0`     | all        | Max workspace calc target (`0` or `INT8`, `1` or `FP8`, `2` or `BOTH`). Default is INT8.                  |
| `GEMMUL8_MAX_M`         | `0`     | all        | Maximum value of `M` used to preallocate workspace memory.                                                |
| `GEMMUL8_MAX_N`         | `0`     | all        | Maximum value of `N` used to preallocate workspace memory.                                                |
| `GEMMUL8_MAX_K`         | `0`     | all        | Maximum value of `K` used to preallocate workspace memory.                                                |
| `GEMMUL8_MAX_NUM_MOD`   | `2`     | all        | Maximum number of moduli used when computing the size of the preallocated workspace.                      |
| `GEMMUL8_SKIP_SCALE_A`  | `0`     | all        | Enables skipping redundant preprocessing for `A` (`1` = enable, `0` = disable).                           |
| `GEMMUL8_SKIP_SCALE_B`  | `0`     | all        | Enables skipping redundant preprocessing for `B` (`1` = enable, `0` = disable).                           |

- This hook mode maintains an independent workspace **per BLAS handle** (`cublasHandle_t` / `hipblasHandle_t`).
- For each handle, the hook allocates **three independent GPU work buffers**:
  - `workA`: cache area for preprocessed `A`
  - `workB`: cache area for preprocessed `B`
  - `workC`: remaining workspace used by the emulation routine
- Each buffer follows a **grow-only** policy:
  it is resized upward on demand and is never shrunk automatically.
- Allocation/free use stream-ordered APIs (`cudaMallocAsync/cudaFreeAsync` or HIP equivalents) on the **current stream**.
- When the same handle is used with different CUDA/HIP streams across calls, the hook enforces ordering by inserting
  an event dependency (`eventRecord` on the previous stream -> `streamWaitEvent` on the current stream).
  This prevents freeing/reallocating a workspace buffer before the previous stream has finished using it.
- Workspace sizing logic (per buffer):
  - Let `(needA, needB, needC)` be computed from `gemmul8::workSize(...)` for the current GEMM shape and parameters.
  - Let `(prevA, prevB, prevC)` be the previously allocated sizes for the same handle.
  - If `GEMMUL8_SKIP_SCALE_A=1` and/or `GEMMUL8_SKIP_SCALE_B=1`, the hook may also enforce a per-process maximum size
    `(maxA, maxB, maxC)` computed from `GEMMUL8_MAX_{M,N,K,NUM_MOD}` and `GEMMUL8_MAXWS_BACKEND`.
- The requested sizes become:
  - `reqA = max(prevA, needA, maxA)` if `GEMMUL8_SKIP_SCALE_A=1`, else `max(prevA, needA)`
  - `reqB = max(prevB, needB, maxB)` if `GEMMUL8_SKIP_SCALE_B=1`, else `max(prevB, needB)`
  - `reqC = max(prevC, needC, maxC)` if either skip flag is enabled, else `max(prevC, needC)`
- The workspaces are released when the corresponding handle is destroyed
  (intercepted via `cublasDestroy_v2` on CUDA or the mapped `hipblasDestroy` on HIP).
- When `GEMMUL8_SKIP_SCALE_{A|B}=1`, redundant preprocessing for `A`/`B` is skipped (see below).

> [!IMPORTANT]
>
> - `GEMMUL8_SKIP_SCALE_{A|B}=1` enables **automatic reuse** of already-preprocessed intermediate data for `A`/`B`
>   within the hook, when it is safe according to the cache conditions below.
> - The decision is based on **pointer identity and cached metadata only**.
>   The hook does **not** verify the contents of `A`/`B`.
>
> Automatic skipping for `A` or `B` is enabled only when **all** of the following hold between consecutive calls:
>
> 1. `GEMMUL8_SKIP_SCALE_A=1` (for skipping `A`) and/or `GEMMUL8_SKIP_SCALE_B=1` (for skipping `B`).
> 2. The emulation path is taken in both calls:
>    - SGEMM / CGEMM: `2 <= GEMMUL8_NUM_MOD_{S|C} <= 13`
>    - DGEMM / ZGEMM: `2 <= GEMMUL8_NUM_MOD_{D|Z} <= 20`
> 3. The same BLAS handle is used (cache is **per-handle**).
> 4. The same emulation backend is used (`GEMMUL8_BACKEND`: INT8 or FP8).
> 5. The same `fastmode` setting is used.
> 6. The same `num_moduli` is used.
> 7. The inner dimension `K` is identical.
> 8. (For skipping `A`) the following are identical:
>    - `A` device pointer
>    - `M`, `lda`, and `transa`
>    - the internal cached workspace pointer for `A` (`workA`) is unchanged
> 9. (For skipping `B`) the following are identical:
>    - `B` device pointer
>    - `N`, `ldb`, and `transb`
>    - the internal cached workspace pointer for `B` (`workB`) is unchanged
>
> If any condition differs, the hook performs preprocessing again for that operand.

> [!TIP]
> To keep the internal workspace pointers (`workA`/`workB`) stable across calls, avoid workspace reallocation.
> Setting `GEMMUL8_MAX_M`, `GEMMUL8_MAX_N`, `GEMMUL8_MAX_K`, `GEMMUL8_MAX_NUM_MOD` appropriately helps the hook allocate sufficiently large buffers early, reducing pointer changes.
> If you may switch `GEMMUL8_BACKEND` at runtime, set `GEMMUL8_MAXWS_BACKEND=BOTH`.

> [!CAUTION]
> Skip scaling assumes that the contents of `A` or `B` remain unchanged in GPU memory.
> If `A`/`B` data are modified between GEMM calls, do not rely on skipping (set `GEMMUL8_SKIP_SCALE_{A|B}=0`).

> [!NOTE]
> `GEMMUL8_MAX_*` and `GEMMUL8_MAXWS_BACKEND` are read only once (on first use) to compute `wsmax`.
> Other variables (e.g., `GEMMUL8_NUM_MOD_*`, `GEMMUL8_FASTMODE_*`, `GEMMUL8_BACKEND`, `GEMMUL8_SKIP_SCALE_*`) are read at each GEMM call.

#### How to change environment variables programmatically

You can also set these environment variables programmatically from within your code using setenv.

```cpp
// Run emulation with num_moduli = 15 & fastmode = true
char num_moduli[12];
snprintf(num_moduli, sizeof(num_moduli), "%u", 15u);
setenv("GEMMUL8_NUM_MOD_D", num_moduli, 1);
setenv("GEMMUL8_FASTMODE_D", "1", 1);
cublasDgemm_v2(...);
```

## Numerical results

See numerical results in the separate repository: [GEMMul8_numerical_results](https://github.com/UCHINO-Yuki/GEMMul8_numerical_results)

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

```bibtex
@inproceedings{10.1145/3731599.3767539,
  author = {Uchino, Yuki and Ozaki, Katsuhisa and Imamura, Toshiyuki},
  title = {High-Performance and Power-Efficient Emulation of Matrix Multiplication using INT8 Matrix Engines},
  year = {2025},
  isbn = {9798400718717},
  publisher = {Association for Computing Machinery},
  address = {New York, NY, USA},
  url = {https://doi.org/10.1145/3731599.3767539},
  doi = {10.1145/3731599.3767539},
  booktitle = {Proceedings of the SC '25 Workshops of the International Conference for High Performance Computing, Networking, Storage and Analysis},
  pages = {1824-1831},
  numpages = {8},
  series = {SC Workshops '25}
}
```

```bibtex
@misc{ozaki2025ozakischemeiigemmoriented,
    title={Ozaki Scheme II: A GEMM-oriented emulation of floating-point matrix multiplication using an integer modular technique},
    author={Katsuhisa Ozaki and Yuki Uchino and Toshiyuki Imamura},
    year={2025},
    eprint={2504.08009},
    archivePrefix={arXiv},
    primaryClass={cs.MS},
    url={https://arxiv.org/abs/2504.08009},
}
```

```bibtex
@misc{uchino2025emulationcomplexmatrixmultiplication,
    title={Emulation of Complex Matrix Multiplication based on the Chinese Remainder Theorem},
    author={Yuki Uchino and Qianxiang Ma and Toshiyuki Imamura and Katsuhisa Ozaki and Patrick Lars Gutsche},
    year={2025},
    eprint={2512.08321},
    archivePrefix={arXiv},
    primaryClass={cs.DC},
    url={https://arxiv.org/abs/2512.08321},
}
```

## License

This project is licensed under the MIT License. See the LICENSE file for details.
