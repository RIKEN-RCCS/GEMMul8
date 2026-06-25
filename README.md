# GEMMul8<!-- omit in toc -->

GEMMul8 (GEMMulate/ジェミュレート): GEMM emulation and its extension to BLAS-like matrix operations using INT8/FP8 matrix engines

GEMMul8 is a library for emulating GEMM using low-precision matrix engines, including INT8 and FP8.
The current version extends this GEMM emulation framework to several BLAS-like Level-3 matrix operations, including symmetric, Hermitian, triangular, and triangular-solve routines.

The library is based on the Ozaki Scheme II and supports selectable INT8- or FP8-based emulation backends within each supported routine.
This design enables bit-wise reproducible results while using low-precision matrix engines for high-throughput computation.

- [Technical Overview](#technical-overview)
- [Supported operations](#supported-operations)
- [Build](#build)
  - [make options](#make-options)
  - [Example](#example)
    - [CUDA build](#cuda-build)
    - [HIP build](#hip-build)
- [Running Test Codes](#running-test-codes)
  - [Test options](#test-options)
  - [Routine options](#routine-options)
  - [Precision options](#precision-options)
  - [Disable options](#disable-options)
  - [BLAS parameter options](#blas-parameter-options)
  - [Examples](#examples)
- [Usage](#usage)
  - [1. Direct Usage (Normal mode)](#1-direct-usage-normal-mode)
    - [Example: run emulation for the CUDA backend](#example-run-emulation-for-the-cuda-backend)
    - [Public API](#public-api)
    - [Return value](#return-value)
    - [Workspace query](#workspace-query)
    - [TRSM implementation and block-size control](#trsm-implementation-and-block-size-control)
    - [Behavior of `skip_scalA` / `skip_scalB`](#behavior-of-skip_scala--skip_scalb)
    - [Example: GEMM with skip scaling](#example-gemm-with-skip-scaling)
  - [2. Hijack cuBLAS/hipBLAS routines (Hook Mode)](#2-hijack-cublashipblas-routines-hook-mode)
    - [Interception targets](#interception-targets)
    - [Ex-routine dispatch policy](#ex-routine-dispatch-policy)
    - [How to enable the hook](#how-to-enable-the-hook)
    - [Configure emulation parameters via environment variables](#configure-emulation-parameters-via-environment-variables)
    - [Max-workspace preallocation](#max-workspace-preallocation)
    - [Hook workspace and skip-scaling behavior](#hook-workspace-and-skip-scaling-behavior)
    - [How to change environment variables programmatically](#how-to-change-environment-variables-programmatically)
- [Numerical results](#numerical-results)
- [Acknowledgment](#acknowledgment)
  - [Assistance with debugging](#assistance-with-debugging)
  - [Assistance with preliminary experiments](#assistance-with-preliminary-experiments)
- [Contact (Responsible Developer)](#contact-responsible-developer)
- [References](#references)
- [Citations](#citations)
- [License](#license)

## Technical Overview

GEMMul8 implements high-precision emulation of BLAS-like matrix operations based on Ozaki Scheme II, which utilizes the Chinese Remainder Theorem (CRT).
A larger number of moduli (`num_moduli`) for the CRT results in higher precision at the cost of increased computation time.

The current implementation supports both **CUDA** and **HIP** backends.

GEMMul8 supports two low-precision emulation backends:

- INT8 backend: uses standard BLAS handle (cuBLAS/hipBLAS handle) or Lt handle (cuBLASLt/hipBLASLt handle).
- FP8 backend: uses Lt handle (cuBLASLt/hipBLASLt handle).

> [!CAUTION]
>
> This library does not support FP8-based emulation on Hopper architectures.

## Supported operations

GEMMul8 currently provides the following BLAS-like operations.

| Routine              | Operation type                                  |
| :------------------- | :---------------------------------------------- |
| `gemm`, `gemmLt`     | general matrix-matrix multiplication            |
| `symm`, `symmLt`     | symmetric matrix-matrix multiplication          |
| `syrk`, `syrkLt`     | symmetric rank-k update                         |
| `syr2k`, `syr2kLt`   | symmetric rank-2k update                        |
| `syrkx`, `syrkxLt`   | symmetric rank-k update with two input matrices |
| `hemm`, `hemmLt`     | Hermitian matrix-matrix multiplication          |
| `herk`, `herkLt`     | Hermitian rank-k update                         |
| `her2k`, `her2kLt`   | Hermitian rank-2k update                        |
| `herkx`, `herkxLt`   | Hermitian rank-k update with two input matrices |
| `trmm`, `trmmLt`     | triangular matrix-matrix multiplication         |
| `trsm`, `trsmLt`     | triangular solve with multiple right-hand sides |
| `trtrmm`, `trtrmmLt` | triangular-by-triangular matrix multiplication  |

The Hermitian routines are intended for complex arithmetic.

## Build

Run `make` in the project root directory to build both the static and shared libraries.

```bash
make -j$(nproc)
```

This creates:

- `lib/libgemmul8.a`
- `lib/libgemmul8.so`

To rebuild from scratch:

```bash
make clean
make -j$(nproc)
```

### make options

| Option      | Default           | Description                                                                                            |
| :---------- | :---------------- | :----------------------------------------------------------------------------------------------------- |
| `CUDA_PATH` | `/usr/local/cuda` | Path to your CUDA toolkit installation. Used for CUDA backends.                                        |
| `HIP_PATH`  | `/opt/rocm`       | Path to your HIP (ROCm) toolkit installation. Used for HIP backends.                                   |
| `BACKEND`   | `auto`            | Select GPU backend: `cuda`, `hip`, or `auto` (auto-detect).                                            |
| `GPU_ARCH`  | `auto`            | Target GPU architecture.<br>Examples: `90` (H100), `100` (B200), `gfx90a` (MI250X), `gfx942` (MI300X). |
| `TEMPDIR`   | `build/tmp`       | Temporary directory used by the compiler.                                                              |

> [!NOTE]
>
> - `BACKEND=auto` will attempt to detect your GPU vendor automatically.
> - `GPU_ARCH=auto` will automatically detect and use the appropriate compute capability or architecture for your GPU.
> - Target GPU architecture can be found from e.g., [CUDA GPU Compute Capability](https://developer.nvidia.com/cuda-gpus) or [AMD GPU hardware specifications](https://rocm.docs.amd.com/en/latest/reference/gpu-arch-specs.html).

### Example

#### CUDA build

Build for an NVIDIA H100/H200 GPU (Compute Capability 9.0)

```bash
make -j$(nproc) BACKEND=cuda CUDA_PATH=/usr/local/cuda GPU_ARCH=90
```

#### HIP build

Build for an AMD MI300X GPU (gfx942 architecture)

```bash
make -j$(nproc) BACKEND=hip HIP_PATH=/opt/rocm GPU_ARCH=gfx942
```

## Running Test Codes

After building the library, the test program can be built and run from the `test/` directory.

```bash
cd test
make -j$(nproc)
```

The test executable accepts three groups of options:

```bash
make run MODE="<test-option>... <routine-option>... <precision-option>... [disable-option]..."
```

### Test options

| Option               | Description                                 |
| :------------------- | :------------------------------------------ |
| `accuracy_square`    | Run accuracy tests for square matrices      |
| `accuracy_rectangle` | Run accuracy tests for rectangular matrices |
| `time_square`        | Run timing tests for square matrices        |
| `time_rectangle`     | Run timing tests for rectangular matrices   |

### Routine options

| Option   | Description |
| :------- | :---------- |
| `GEMM`   | Run GEMM    |
| `SYMM`   | Run SYMM    |
| `SYRK`   | Run SYRK    |
| `SYR2K`  | Run SYR2K   |
| `SYRKX`  | Run SYRKX   |
| `HEMM`   | Run HEMM    |
| `HERK`   | Run HERK    |
| `HER2K`  | Run HER2K   |
| `HERKX`  | Run HERKX   |
| `TRMM`   | Run TRMM    |
| `TRSM`   | Run TRSM    |
| `TRTRMM` | Run TRTRMM  |

### Precision options

| Option | Description            |
| :----- | :--------------------- |
| `S`    | Run FP32 real tests    |
| `D`    | Run FP64 real tests    |
| `C`    | Run FP32 complex tests |
| `Z`    | Run FP64 complex tests |

### Disable options

| Option           | Description                 |
| :--------------- | :-------------------------- |
| `no_Ozaki2_INT8` | Disable Ozaki-II INT8 tests |
| `no_Ozaki2_FP8`  | Disable Ozaki-II FP8 tests  |
| `no_Ozaki1_INT8` | Disable Ozaki-I INT8 tests  |

### BLAS parameter options

By default, the test driver runs all supported combinations of BLAS parameters for each selected routine.
The following options can be used to restrict the tested parameter combinations.

| Option       | Values                   | Applies to                                                                         |
| :----------- | :----------------------- | :--------------------------------------------------------------------------------- |
| `trans=...`  | `all`, `N`, `T`, `C`     | `SYRK`, `SYR2K`, `SYRKX`, `HERK`, `HER2K`, `HERKX`, `TRMM`, `TRSM`                 |
| `transA=...` | `all`, `N`, `T`, `C`     | `GEMM`, `TRTRMM`                                                                   |
| `transB=...` | `all`, `N`, `T`, `C`     | `GEMM`, `TRTRMM`                                                                   |
| `uplo=...`   | `all`, `upper`, `lower`  | `SYRK`, `SYR2K`, `SYRKX`, `HERK`, `HER2K`, `HERKX`, `SYMM`, `HEMM`, `TRMM`, `TRSM` |
| `uploA=...`  | `all`, `upper`, `lower`  | `TRTRMM`                                                                           |
| `uploB=...`  | `all`, `upper`, `lower`  | `TRTRMM`                                                                           |
| `diag=...`   | `all`, `nonunit`, `unit` | `TRMM`, `TRSM`                                                                     |
| `diagA=...`  | `all`, `nonunit`, `unit` | `TRTRMM`                                                                           |
| `diagB=...`  | `all`, `nonunit`, `unit` | `TRTRMM`                                                                           |
| `side=...`   | `all`, `left`, `right`   | `SYMM`, `HEMM`, `TRMM`, `TRSM`                                                     |

Short aliases are also accepted:

| Parameter         | Aliases  |
| :---------------- | :------- |
| `all`             | `A`      |
| `upper`, `lower`  | `U`, `L` |
| `left`, `right`   | `L`, `R` |
| `nonunit`, `unit` | `N`, `U` |

For `SYRK`, `SYR2K`, and `SYRKX`, only `trans=N` and `trans=T` are used.
For `HERK`, `HER2K`, and `HERKX`, only `trans=N` and `trans=C` are used.

### Examples

```bash
# Run only non-transposed FP64 GEMM accuracy tests
make run MODE="accuracy_rectangle GEMM D transA=N transB=N"

# Run lower-triangular SYRK timing tests only
make run MODE="time_square SYRK D uplo=lower trans=N"

# Run left-side upper-triangular TRSM timing tests with non-unit diagonal
make run MODE="time_square TRSM D side=left uplo=upper trans=N diag=nonunit"

# Run one TRTRMM parameter subset
make run MODE="time_square TRTRMM Z uploA=upper uploB=lower transA=N transB=C diagA=nonunit diagB=unit"
```

## Usage

This library provides two ways to use GEMMul8:

1. Direct usage: explicitly call `gemmul8::gemm`, `gemmul8::syrk`, `gemmul8::trsm`, etc.
2. Hook mode: intercept existing BLAS routine calls through `LD_PRELOAD`.

### 1. Direct Usage (Normal mode)

Call GEMMul8 functions explicitly from your source code.
This gives you fine-grained control over the emulation parameters.

#### Example: run emulation for the CUDA backend

See the sample code in `sample/`.

#### Public API

Include the umbrella header:

```cpp
#include "gemmul8.hpp"
```

Each routine follows the corresponding cuBLAS/hipBLAS argument convention as closely as possible, with additional GEMMul8-specific arguments.
See `include/gemm.hpp`, `include/symm.hpp`, etc. for the full function signatures.

The TRSM block-size control API is declared in `include/trsm.hpp`:

```cpp
void gemmul8::set_block_size_trsm(const int nB) noexcept;
int gemmul8::get_block_size_trsm() noexcept;
```

> [!NOTE]
>
> `gemmul8::trsm` and `gemmul8::trsmLt` follow the BLAS `trsm` convention:
>
> ```text
> B := X
> ```
>
> That is, the input/output matrix `B` is overwritten in place by the solution matrix.
> For left-side solve:
>
> ```text
> op(A) * X = alpha * B
> ```
>
> For right-side solve:
>
> ```text
> X * op(A) = alpha * B
> ```

#### Return value

When `work != nullptr`, GEMMul8 routines execute the requested operation and return elapsed times in seconds.

For most routines, the returned vector contains four internal phase timings:

```text
t[0]: scaling and quantization
t[1]: low-precision matrix multiplication
t[2]: re-quantization of matrix products
t[3]: final CRT reduction and undo scaling
```

For `trsm`, the returned vector has a different meaning:

```text
t[0]: standard BLAS TRSM phase
t[1]: GEMMul8 GEMM phase
```

When `work == nullptr`, the routine is used as a workspace-query call.
The requested BLAS-like operation is not executed; instead, GEMMul8 computes and returns the required workspace sizes in bytes.

For most routines:

```text
t[0]: total workspace size
t[1]: workspace size associated with A
t[2]: workspace size associated with B
```

For one-input routines such as `syrk` and `herk`, the B-associated workspace may be absent or unused.
For `trsm`, only `t[0]`, the total workspace size, is returned.

#### Workspace query

GEMMul8 provides two ways to query the required workspace size.

1. Call the lightweight query functions:

- `gemmul8::workSize`
- `gemmul8::workSizeTrsm`

2. Call the corresponding GEMMul8 routine with `work == nullptr`.

- In this mode, the operation itself is not executed.
- The returned vector contains workspace sizes in bytes, using the same convention described in [Return value](#return-value).

See `include/worksize.hpp` for the full function signatures.
The compact size arguments are interpreted as follows.

| Function                            | Recommended workspace query |
| :---------------------------------- | :-------------------------- |
| `gemm`                              | `workSize(m, n, k, ...)`    |
| `symm`, `hemm` with `side == LEFT`  | `workSize(m, n, m, ...)`    |
| `symm`, `hemm` with `side == RIGHT` | `workSize(m, n, n, ...)`    |
| `syrk`, `herk`                      | `workSize(n, n, k, ...)`    |
| `syr2k`, `her2k`                    | `workSize(n, n, k, ...)`    |
| `syrkx`, `herkx`                    | `workSize(n, n, k, ...)`    |
| `trmm` with `side == LEFT`          | `workSize(m, n, m, ...)`    |
| `trmm` with `side == RIGHT`         | `workSize(m, n, n, ...)`    |
| `trtrmm`                            | `workSize(n, n, n, ...)`    |

For `trsm`, use the dedicated query function `gemmul8::workSizeTrsm`.
`workSizeTrsm()` depends on `side`, `m`, `n`, `num_moduli`, the selected backend, the element type, and the current TRSM block-size.
When a custom TRSM block size is used, call `gemmul8::set_block_size_trsm(nB)` before calling `workSizeTrsm()` and before allocating the workspace.

#### TRSM implementation and block-size control

`gemmul8::trsm` and `gemmul8::trsmLt` use a blocked triangular-solve algorithm internally.

The implementation combines:

- standard cuBLAS/hipBLAS TRSM for triangular solves on diagonal blocks, and
- `gemmul8::gemm` / `gemmul8::gemmLt` for updates to the remaining blocks.

The internal block size can be controlled with:

```cpp
gemmul8::set_block_size_trsm(const int nB);
```

If `set_block_size_trsm(nB)` has not been called, or if the value set by `set_block_size_trsm(nB)` is non-positive, GEMMul8 automatically selects the TRSM block size from the detected GPU architecture and backend.

> [!NOTE]
>
> The automatically selected TRSM block size is a heuristic default and is not guaranteed to be the fastest setting.
> For performance tuning, benchmark several block sizes and set a custom value with `gemmul8::set_block_size_trsm(nB)`.

A positive value passed to `set_block_size_trsm(nB)` is used as the block size for subsequent `trsm()` and `trsmLt()` calls.
The setting is process-global and also affects the workspace size returned by `workSizeTrsm()`.
Therefore, when using a custom block size, call `set_block_size_trsm(nB)` before calling `workSizeTrsm()` and before allocating the workspace.

> [!CAUTION]
>
> `get_block_size_trsm()` only returns the value explicitly set by `set_block_size_trsm()`.
> It does not report the block size automatically selected by GEMMul8.
> Therefore, if get_block_size_trsm() returns a non-positive value, it means that automatic block-size selection is enabled, not that the internally selected block size is non-positive.

#### Behavior of `skip_scalA` / `skip_scalB`

This skip mechanism is designed for repeated calls to GEMMul8 routines that reuse the same input matrix `A` and/or `B`.
It applies to routines that expose `workA`, `workB`, `enable_skip_scalA`, `enable_skip_scalB`, `skip_scalA`, and `skip_scalB`.

> [!NOTE]
>
> It does not apply to `trsm`, because the current `trsm` interface does not expose skip-scaling arguments.

- Most routines internally preprocess the input matrices `A` and/or `B` into a backend-specific low-precision representation (INT8/FP8) and perform modular multiplications across multiple moduli.
- This preprocessing step can be **skipped** in consecutive calls if the same matrices are reused, allowing for substantial performance gains.
- If `enable_skip_scal{A|B} = true`, additional workspace is reserved so that the preprocessed representation of `A`/`B` can be retained between calls.
- If `enable_skip_scal{A|B} = true && skip_scal{A|B} = true`, the preprocessing step for `A`/`B` is **actually skipped**, and previously prepared data are reused for faster execution.

> [!NOTE]
>
> When using `skip_scalA` / `skip_scalB`, the preprocessing step that converts `A`/`B` into an internal backend-specific representation (INT8/FP8) is skipped.
> For correctness, the following conditions **must all hold** between consecutive calls:
>
> 1. The effective dimensions of the reused input matrix must be identical to those in the previous call.
> 2. The operation type (`CUBLAS_OP_N` / `CUBLAS_OP_C` / `CUBLAS_OP_T`) for `A`/`B` must be the same as before.
> 3. The value of `num_moduli` must remain unchanged.
> 4. The `fastmode` setting must be identical to that of the previous call.
> 5. The contents of `A`/`B` in device memory must not be modified between calls.
> 6. The selected emulation backend must be identical in both calls (INT8 or FP8).
>
> GEMMul8 does not verify the contents of `A`/`B` when skipping; correctness requires the user to ensure immutability.
>
> If any of these conditions differ, the cached scaled data become invalid, and skipping must **not** be used.
> In such cases, set `skip_scalA=false` / `skip_scalB=false`.

> [!CAUTION]
>
> This skip mechanism is designed for repeated routine calls with identical A or B.
> Use it only when you are certain that the input matrices and configuration have not changed.
> When in doubt, disable skipping to ensure correctness.

#### Example: GEMM with skip scaling

The following example demonstrates skip scaling for `gemm`.
The same concept applies to other routines that expose `workA`, `workB`, and skip-scaling arguments, but the effective matrix dimensions must be chosen according to the routine.

```cpp
#include "gemmul8.hpp"

// 1. Create a handle to the cuBLAS library context
cublasHandle_t cublas_handle;
cublasCreate(&cublas_handle);

// 2. Settings
const unsigned num_moduli = 14u;
const bool fastmode = false;

bool enable_skip_scalA = false;
bool enable_skip_scalB = true;
bool skip_scalA = false;
bool skip_scalB = false;

// 3. Matrix shapes
const size_t m1 = 64, n1 = 10, k1 = 10; // 1st GEMM: 64×10 × 10×10
const size_t m2 = 20, n2 = 10, k2 = 10; // 2nd GEMM: 20×10 × 10×10

// 4. Allocate workspace
size_t worksizeA, worksizeB;
const size_t worksize = gemmul8::workSize(
    std::max(m1, m2), std::max(n1, n2), std::max(k1, k2), num_moduli,
    enable_skip_scalA, enable_skip_scalB, &worksizeA, &worksizeB);

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
              enable_skip_scalA, enable_skip_scalB, skip_scalA, skip_scalB);

// 6. Reuse preprocessed B (second call: skip_scalB = true)
skip_scalB = true;
gemmul8::gemm(cublas_handle,
              CUBLAS_OP_N, CUBLAS_OP_N,
              m1, n1, k1,
              &alpha, devA1, lda,
              devB, ldb,
              &beta, devC, ldc,
              num_moduli, fastmode, (void*)work_rem, (void*)workA, (void*)workB,
              enable_skip_scalA, enable_skip_scalB, skip_scalA, skip_scalB);

// 7. Free workspace
cudaFree(work);

// 8. Destroy a handle
cublasDestroy(cublas_handle);
```

### 2. Hijack cuBLAS/hipBLAS routines (Hook Mode)

Intercept standard cuBLAS/hipBLAS routine calls automatically without modifying the application source code.
The hook path is intended to support all GEMMul8-supported routines except `trtrmm`, because `trtrmm` is a GEMMul8-specific extension and has no standard cuBLAS/hipBLAS routine to intercept.

#### Interception targets

The hook mode intercepts selected standard cuBLAS/hipBLAS entry points and routes matching calls to GEMMul8 emulation.

The hook targets are exact-symbol based. If a `v2` symbol exists in cuBLAS, GEMMul8 hooks the `v2` / `v2_64` symbols rather than the legacy non-`v2` symbols.

Batched, strided-batched, and grouped-batched routines are not hook targets.

| Family | CUDA hook targets                                                                                              | HIP hook targets                                          |
| :----- | :------------------------------------------------------------------------------------------------------------- | :-------------------------------------------------------- |
| GEMM   | `cublas{S,D,C,Z}gemm_v2` / `_64`<br>`cublasGemmEx` / `_64`<br>                                                 | `hipblas{S,D,C,Z}gemm` / `_64`<br>`hipblasGemmEx` / `_64` |
| GEMM   | `cublas{C,Z}gemm3m` / `_64`<br>`cublasSgemmEx` / `_64`<br>`cublasCgemmEx` / `_64`<br>`cublasCgemm3mEx` / `_64` | not supported                                             |
| SYMM   | `cublas{S,D,C,Z}symm_v2` / `_64`                                                                               | `hipblas{S,D,C,Z}symm` / `_64`                            |
| SYRK   | `cublas{S,D,C,Z}syrk_v2` / `_64`<br>`cublasCsyrkEx` / `_64`<br>`cublasCsyrk3mEx` / `_64`                       | `hipblas{S,D,C,Z}syrk` / `_64`                            |
| SYR2K  | `cublas{S,D,C,Z}syr2k_v2` / `_64`                                                                              | `hipblas{S,D,C,Z}syr2k` / `_64`                           |
| SYRKX  | `cublas{S,D,C,Z}syrkx` / `_64`                                                                                 | `hipblas{S,D,C,Z}syrkx` / `_64`                           |
| HEMM   | `cublas{C,Z}hemm_v2` / `_64`                                                                                   | `hipblas{C,Z}hemm` / `_64`                                |
| HERK   | `cublas{C,Z}herk_v2` / `_64`<br>`cublasCherkEx` / `_64`<br>`cublasCherk3mEx` / `_64`                           | `hipblas{C,Z}herk` / `_64`                                |
| HER2K  | `cublas{C,Z}her2k_v2` / `_64`                                                                                  | `hipblas{C,Z}her2k` / `_64`                               |
| HERKX  | `cublas{C,Z}herkx` / `_64`                                                                                     | `hipblas{C,Z}herkx` / `_64`                               |
| TRMM   | `cublas{S,D,C,Z}trmm_v2` / `_64`                                                                               | `hipblas{S,D,C,Z}trmm` / `_64`                            |
| TRSM   | `cublas{S,D,C,Z}trsm_v2` / `_64`                                                                               | `hipblas{S,D,C,Z}trsm` / `_64`                            |

`trtrmm` / `trtrmmLt` are not hook targets because they are GEMMul8-specific routines rather than standard cuBLAS/hipBLAS routines.

> [!NOTE]
>
> For cuBLAS routines with `_v2` variants, the table lists the actual hook targets.  
> In user code, however, the `_v2` suffix is usually unnecessary because `cublas_v2.h` maps the non-`_v2` routine names to the corresponding `_v2` entry points.

#### Ex-routine dispatch policy

The hook intercepts several Ex routines, but GEMMul8 emulation is applied only to the same-type FP32/FP64 cases supported by the direct GEMMul8 API.

For `cublasGemmEx` / `cublasGemmEx_64`:

| Input/output types                       | compute type         | Hook behavior                  |
| :--------------------------------------- | :------------------- | :----------------------------- |
| `CUDA_R_32F`, `CUDA_R_32F`, `CUDA_R_32F` | `CUBLAS_COMPUTE_32F` | GEMMul8 FP32 real emulation    |
| `CUDA_R_64F`, `CUDA_R_64F`, `CUDA_R_64F` | `CUBLAS_COMPUTE_64F` | GEMMul8 FP64 real emulation    |
| `CUDA_C_32F`, `CUDA_C_32F`, `CUDA_C_32F` | `CUBLAS_COMPUTE_32F` | GEMMul8 FP32 complex emulation |
| `CUDA_C_64F`, `CUDA_C_64F`, `CUDA_C_64F` | `CUBLAS_COMPUTE_64F` | GEMMul8 FP64 complex emulation |
| otherwise                                | any                  | native cuBLAS/hipBLAS fallback |

For CUDA-only `cublasSgemmEx`, GEMMul8 emulation is used only when `A`, `B`, and `C` are all `CUDA_R_32F`.

For CUDA-only `cublasCgemmEx` and `cublasCgemm3mEx`, GEMMul8 emulation is used only when `A`, `B`, and `C` are all `CUDA_C_32F`.

Other mixed-precision Ex cases, such as FP16, BF16, TF32, INT8, and mixed input/output types, are forwarded to the native cuBLAS routine.

#### How to enable the hook

1. Build the library.
2. Set the `LD_PRELOAD` environment variable.

```bash
export LD_PRELOAD=/<path-to-GEMMul8>/lib/libgemmul8.so
```

3. Run your application.

#### Configure emulation parameters via environment variables

Hook mode uses operation-specific environment variables. The general form is:

```text
GEMMUL8_BACKEND_<OP>
GEMMUL8_NUM_MOD_<S|D|C|Z>_<OP>
GEMMUL8_FASTMODE_<S|D|C|Z>_<OP>
```

where `<OP>` is one of:

```text
GEMM
SYMM_LEFT, SYMM_RIGHT
SYRK
SYR2K
SYRKX
HEMM_LEFT, HEMM_RIGHT
HERK
HER2K
HERKX
TRMM_LEFT, TRMM_RIGHT
TRSM_LEFT, TRSM_RIGHT
```

For side-dependent routines, the suffix depends on the runtime `side` argument.

> [!CAUTION]
>
> `trtrmm` is not a hook target because it is a GEMMul8-specific routine and has no standard cuBLAS/hipBLAS routine to intercept.

```bash
# GEMM, operation-specific form
export GEMMUL8_BACKEND_GEMM=INT8
export GEMMUL8_NUM_MOD_D_GEMM=15
export GEMMUL8_FASTMODE_D_GEMM=1

# GEMM, backward-compatible form
export GEMMUL8_BACKEND=INT8
export GEMMUL8_NUM_MOD_D=15
export GEMMUL8_FASTMODE_D=1

# SYMM with side == LEFT
export GEMMUL8_BACKEND_SYMM_LEFT=INT8
export GEMMUL8_NUM_MOD_D_SYMM_LEFT=15
export GEMMUL8_FASTMODE_D_SYMM_LEFT=1

# SYRK
export GEMMUL8_BACKEND_SYRK=INT8
export GEMMUL8_NUM_MOD_D_SYRK=15
export GEMMUL8_FASTMODE_D_SYRK=1

# TRMM with side == RIGHT
export GEMMUL8_BACKEND_TRMM_RIGHT=FP8
export GEMMUL8_NUM_MOD_D_TRMM_RIGHT=12
export GEMMUL8_FASTMODE_D_TRMM_RIGHT=1

# TRSM with side == LEFT
export GEMMUL8_BACKEND_TRSM_LEFT=INT8
export GEMMUL8_NUM_MOD_D_TRSM_LEFT=10
export GEMMUL8_FASTMODE_D_TRSM_LEFT=0

# Global skip-scaling switches
export GEMMUL8_SKIP_SCALE_A=1
export GEMMUL8_SKIP_SCALE_B=1
```

| Variable pattern          | Default | Description                                                                                     |
| :------------------------ | :------ | :---------------------------------------------------------------------------------------------- |
| `GEMMUL8_BACKEND_<OP>`    | `INT8`  | Selects the emulation backend. `0` or `INT8` = INT8 backend; `1` or `FP8` = FP8 backend.        |
| `GEMMUL8_NUM_MOD_S_<OP>`  | `0`     | Number of moduli for FP32 real routines. Native BLAS is used if outside `[2, 13]`.              |
| `GEMMUL8_NUM_MOD_D_<OP>`  | `0`     | Number of moduli for FP64 real routines. Native BLAS is used if outside `[2, 20]`.              |
| `GEMMUL8_NUM_MOD_C_<OP>`  | `0`     | Number of moduli for FP32 complex routines. Native BLAS is used if outside `[2, 13]`.           |
| `GEMMUL8_NUM_MOD_Z_<OP>`  | `0`     | Number of moduli for FP64 complex routines. Native BLAS is used if outside `[2, 20]`.           |
| `GEMMUL8_FASTMODE_S_<OP>` | `0`     | Fast mode switch for FP32 real routines. `1` = fast mode; `0` = accurate mode.                  |
| `GEMMUL8_FASTMODE_D_<OP>` | `0`     | Fast mode switch for FP64 real routines. `1` = fast mode; `0` = accurate mode.                  |
| `GEMMUL8_FASTMODE_C_<OP>` | `0`     | Fast mode switch for FP32 complex routines. `1` = fast mode; `0` = accurate mode.               |
| `GEMMUL8_FASTMODE_Z_<OP>` | `0`     | Fast mode switch for FP64 complex routines. `1` = fast mode; `0` = accurate mode.               |
| `GEMMUL8_SKIP_SCALE_A`    | `0`     | Global switch that enables reuse of preprocessed/scaled `A` when the operand cache key matches. |
| `GEMMUL8_SKIP_SCALE_B`    | `0`     | Global switch that enables reuse of preprocessed/scaled `B` when the operand cache key matches. |

#### Max-workspace preallocation

GEMMul8 normally grows hook workspaces on demand. To stabilize workspace addresses, avoid reallocating workspace, and improve skip-scaling reuse, define the maximum BLAS size arguments for the operations that will be used.

| Operation suffix | Required size variables                                          | Internal workspace query         |
| :--------------- | :--------------------------------------------------------------- | :------------------------------- |
| `GEMM`           | `GEMMUL8_MAX_M_GEMM`, `GEMMUL8_MAX_N_GEMM`, `GEMMUL8_MAX_K_GEMM` | `workSize(m, n, k, ...)`         |
| `SYMM_LEFT`      | `GEMMUL8_MAX_M_SYMM_LEFT`, `GEMMUL8_MAX_N_SYMM_LEFT`             | `workSize(m, n, m, ...)`         |
| `SYMM_RIGHT`     | `GEMMUL8_MAX_M_SYMM_RIGHT`, `GEMMUL8_MAX_N_SYMM_RIGHT`           | `workSize(m, n, n, ...)`         |
| `SYRK`           | `GEMMUL8_MAX_N_SYRK`, `GEMMUL8_MAX_K_SYRK`                       | `workSize(n, n, k, ...)`         |
| `SYR2K`          | `GEMMUL8_MAX_N_SYR2K`, `GEMMUL8_MAX_K_SYR2K`                     | `workSize(n, n, k, ...)`         |
| `SYRKX`          | `GEMMUL8_MAX_N_SYRKX`, `GEMMUL8_MAX_K_SYRKX`                     | `workSize(n, n, k, ...)`         |
| `HEMM_LEFT`      | `GEMMUL8_MAX_M_HEMM_LEFT`, `GEMMUL8_MAX_N_HEMM_LEFT`             | `workSize(m, n, m, ...)`         |
| `HEMM_RIGHT`     | `GEMMUL8_MAX_M_HEMM_RIGHT`, `GEMMUL8_MAX_N_HEMM_RIGHT`           | `workSize(m, n, n, ...)`         |
| `HERK`           | `GEMMUL8_MAX_N_HERK`, `GEMMUL8_MAX_K_HERK`                       | `workSize(n, n, k, ...)`         |
| `HER2K`          | `GEMMUL8_MAX_N_HER2K`, `GEMMUL8_MAX_K_HER2K`                     | `workSize(n, n, k, ...)`         |
| `HERKX`          | `GEMMUL8_MAX_N_HERKX`, `GEMMUL8_MAX_K_HERKX`                     | `workSize(n, n, k, ...)`         |
| `TRMM_LEFT`      | `GEMMUL8_MAX_M_TRMM_LEFT`, `GEMMUL8_MAX_N_TRMM_LEFT`             | `workSize(m, n, m, ...)`         |
| `TRMM_RIGHT`     | `GEMMUL8_MAX_M_TRMM_RIGHT`, `GEMMUL8_MAX_N_TRMM_RIGHT`           | `workSize(m, n, n, ...)`         |
| `TRSM_LEFT`      | `GEMMUL8_MAX_M_TRSM_LEFT`, `GEMMUL8_MAX_N_TRSM_LEFT`             | `workSizeTrsm(LEFT, m, n, ...)`  |
| `TRSM_RIGHT`     | `GEMMUL8_MAX_M_TRSM_RIGHT`, `GEMMUL8_MAX_N_TRSM_RIGHT`           | `workSizeTrsm(RIGHT, m, n, ...)` |

Additional max-workspace variables are also operation-specific:

| Variable pattern             | Default | Description                                                                                                        |
| :--------------------------- | :------ | :----------------------------------------------------------------------------------------------------------------- |
| `GEMMUL8_MAXWS_BACKEND_<OP>` | `INT8`  | Backend used for max-workspace calculation. `0` or `INT8` = INT8, `1` or `FP8` = FP8, `2` or `BOTH` = max of both. |
| `GEMMUL8_MAX_NUM_MOD_<OP>`   | `2`     | Number of moduli used for max-workspace calculation.                                                               |

> [!NOTE]
>
> For GEMM only, the following old names are also accepted when the corresponding `_GEMM` variables are not defined:
>
> ```text
> GEMMUL8_BACKEND
> GEMMUL8_NUM_MOD_S, GEMMUL8_NUM_MOD_D, GEMMUL8_NUM_MOD_C, GEMMUL8_NUM_MOD_Z
> GEMMUL8_FASTMODE_S, GEMMUL8_FASTMODE_D, GEMMUL8_FASTMODE_C, GEMMUL8_FASTMODE_Z
> GEMMUL8_MAXWS_BACKEND
> GEMMUL8_MAX_M, GEMMUL8_MAX_N, GEMMUL8_MAX_K, GEMMUL8_MAX_NUM_MOD
> ```

Example:

```bash
# GEMM max-workspace, operation-specific form
export GEMMUL8_MAXWS_BACKEND_GEMM=BOTH
export GEMMUL8_MAX_M_GEMM=32768
export GEMMUL8_MAX_N_GEMM=32768
export GEMMUL8_MAX_K_GEMM=32768
export GEMMUL8_MAX_NUM_MOD_GEMM=15

# SYMM_LEFT max-workspace
export GEMMUL8_MAXWS_BACKEND_SYMM_LEFT=INT8
export GEMMUL8_MAX_M_SYMM_LEFT=32768
export GEMMUL8_MAX_N_SYMM_LEFT=32768
export GEMMUL8_MAX_NUM_MOD_SYMM_LEFT=15

# SYRK max-workspace
export GEMMUL8_MAXWS_BACKEND_SYRK=INT8
export GEMMUL8_MAX_N_SYRK=32768
export GEMMUL8_MAX_K_SYRK=32768
export GEMMUL8_MAX_NUM_MOD_SYRK=15

# TRMM_RIGHT max-workspace
export GEMMUL8_MAXWS_BACKEND_TRMM_RIGHT=BOTH
export GEMMUL8_MAX_M_TRMM_RIGHT=32768
export GEMMUL8_MAX_N_TRMM_RIGHT=32768
export GEMMUL8_MAX_NUM_MOD_TRMM_RIGHT=12

# TRSM_LEFT max-workspace
export GEMMUL8_MAXWS_BACKEND_TRSM_LEFT=INT8
export GEMMUL8_MAX_M_TRSM_LEFT=32768
export GEMMUL8_MAX_N_TRSM_LEFT=32768
export GEMMUL8_MAX_NUM_MOD_TRSM_LEFT=10
```

#### Hook workspace and skip-scaling behavior

This hook mode maintains an independent workspace per BLAS handle (`cublasHandle_t` / `hipblasHandle_t`).

For each handle, the hook allocates GPU work buffers used by the emulation routine.

- For routines that support skip scaling, the hook may keep separate `workA` and/or `workB` cache areas for preprocessed input matrices.
- The remaining workspace is used as the routine's internal work buffer.
- For `trsm`, the current direct interface uses a single workspace and does not expose `workA` or `workB`.

Each buffer follows a grow-only policy: it is resized upward on demand and is never shrunk automatically.

Allocation/free use stream-ordered APIs (`cudaMallocAsync/cudaFreeAsync` or HIP equivalents) on the current stream. When the same handle is used with different CUDA/HIP streams across calls, the hook enforces ordering by inserting an event dependency (`eventRecord` on the previous stream -> `streamWaitEvent` on the current stream).

The workspaces are released when the corresponding handle is destroyed.

> [!IMPORTANT]
>
> `GEMMUL8_SKIP_SCALE_A=1` and/or `GEMMUL8_SKIP_SCALE_B=1` enables automatic reuse of already-preprocessed intermediate data for `A` and/or `B` within the hook, when it is safe according to the cache conditions below.
>
> The decision is based on pointer identity and cached metadata only. The hook does not verify the contents of `A` or `B`.

Automatic skipping for `A` or `B` is enabled only when all of the following hold between consecutive calls:

1. `GEMMUL8_SKIP_SCALE_A=1` and/or `GEMMUL8_SKIP_SCALE_B=1`.
2. The emulation path is taken in both calls.
3. The same BLAS handle is used.
4. The same emulation backend is used.
5. The same `fastmode` setting is used.
6. The same `num_moduli` is used.
7. The effective dimensions and operation flags of the reused operand are identical.
8. The reused operand has the same device pointer and leading dimension as before.
9. The internal cached workspace pointer for the reused operand is unchanged.

If any condition differs, the hook performs preprocessing again for that operand.

> [!TIP]
>
> To keep internal workspace pointers stable across calls, define the operation-specific maximum-size variables for the operations that will be used.
> For example,
>
> - use `GEMMUL8_MAX_M_GEMM`, `GEMMUL8_MAX_N_GEMM`, `GEMMUL8_MAX_K_GEMM`, and `GEMMUL8_MAX_NUM_MOD_GEMM` for GEMM;
> - use `GEMMUL8_MAX_M_TRMM_RIGHT`, `GEMMUL8_MAX_N_TRMM_RIGHT`, and `GEMMUL8_MAX_NUM_MOD_TRMM_RIGHT` for right-side TRMM.
>
> If you may switch backend for an operation at runtime, set `GEMMUL8_MAXWS_BACKEND_<OP>=BOTH`.

> [!CAUTION]
>
> Skip scaling assumes that the contents of `A` or `B` remain unchanged in GPU memory. If `A` or `B` data are modified between routine calls, do not rely on skipping.

> [!NOTE]
>
> `GEMMUL8_MAX_*_<OP>`, `GEMMUL8_MAXWS_BACKEND_<OP>`, and `GEMMUL8_MAX_NUM_MOD_<OP>` are read only once on first hook use to compute the maximum workspace sizes.
>
> Runtime variables such as `GEMMUL8_NUM_MOD_<S|D|C|Z>_<OP>`, `GEMMUL8_FASTMODE_<S|D|C|Z>_<OP>`, `GEMMUL8_BACKEND_<OP>`, and `GEMMUL8_SKIP_SCALE_*` are read at each intercepted routine call.

#### How to change environment variables programmatically

You can also set these environment variables programmatically from within your code using setenv.

```cpp
// Run GEMM emulation with Backend = INT8, num_moduli = 15 & fastmode = true
char num_moduli[12];
snprintf(num_moduli, sizeof(num_moduli), "%u", 15u);

setenv("GEMMUL8_BACKEND_GEMM", "INT8", 1);
setenv("GEMMUL8_NUM_MOD_D_GEMM", num_moduli, 1);
setenv("GEMMUL8_FASTMODE_D_GEMM", "1", 1);

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

### Assistance with preliminary experiments

The following individuals helped conduct preliminary performance experiments on B200 systems at Yokota Lab:

- Dr. Qianxiang Ma (RIKEN Center for Computational Science, Japan)
- Prof. Rio Yokota (Institute of Science Tokyo, Japan)

The following individuals helped conduct preliminary experiments on the B200 environment of SAKURAONE, SAKURA internet Inc.'s managed HPC cluster service:

- Takeshi Yamashita (SAKURA internet Inc., Japan)
- Fumikazu Konishi (SAKURA internet Inc., Japan)

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
@inproceedings{10.23919/ISC.2026.11520500,
    author={Uchino, Yuki and Ma, Qianxiang and Imamura, Toshiyuki and Ozaki, Katsuhisa and Gutsche, Patrick Lars},
    booktitle={ISC High Performance 2026 Research Paper Proceedings (41st International Conference)}, 
    title={Emulation of Complex Matrix Multiplication based on the Chinese Remainder Theorem}, 
    year={2026},
    volume={},
    number={},
    pages={1-12},
  url = {https://doi.org/10.23919/ISC.2026.11520500},
    doi={10.23919/ISC.2026.11520500}
}
```

```bibtex
@misc{uchino2026doubleprecisionmatrixmultiplicationemulation,
      title={Double-Precision Matrix Multiplication Emulation via Ozaki-II Scheme with FP8 Quantization},
      author={Yuki Uchino and Katsuhisa Ozaki and Toshiyuki Imamura},
      year={2026},
      eprint={2603.10634},
      archivePrefix={arXiv},
      primaryClass={cs.DC},
      url={https://arxiv.org/abs/2603.10634},
}
```

## License

This project is licensed under the MIT License. See the LICENSE file for details.
