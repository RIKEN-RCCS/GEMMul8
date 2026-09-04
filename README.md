# GEMMul8<!-- omit in toc -->

GEMMul8 (GEMMulate/ジェミュレート): GEMM emulation and its extension to BLAS-like matrix operations using INT8/FP8 matrix engines

GEMMul8 is a library for emulating GEMM using low-precision matrix engines, including INT8 and FP8.
The current version extends this GEMM emulation framework to several BLAS-like Level-3 matrix operations, including symmetric, Hermitian, triangular, and triangular-solve routines.

The library is based on the Ozaki Scheme II and supports selectable INT8- or FP8-based emulation backends within each supported routine.
This design enables bit-wise reproducible results while using low-precision matrix engines for high-throughput computation.

- [Technical Overview](#technical-overview)
- [Requirements](#requirements)
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
  - [Memory saving options](#memory-saving-options)
  - [BLAS parameter options](#blas-parameter-options)
  - [Examples](#examples)
- [Usage](#usage)
  - [1. Direct Usage (Normal mode)](#1-direct-usage-normal-mode)
    - [Example: run emulation for the CUDA backend](#example-run-emulation-for-the-cuda-backend)
    - [Public API](#public-api)
    - [Handle-local execution settings](#handle-local-execution-settings)
    - [Return value](#return-value)
    - [Workspace query](#workspace-query)
    - [Memory-saving mode](#memory-saving-mode)
    - [TRSM implementation and block-size control](#trsm-implementation-and-block-size-control)
    - [Behavior of `skip_scalA` / `skip_scalB`](#behavior-of-skip_scala--skip_scalb)
    - [Example: GEMM with skip scaling](#example-gemm-with-skip-scaling)
  - [2. Hijack cuBLAS/hipBLAS routines (Hook Mode)](#2-hijack-cublashipblas-routines-hook-mode)
    - [Interception targets](#interception-targets)
    - [Ex-routine dispatch policy](#ex-routine-dispatch-policy)
    - [How to enable the hook](#how-to-enable-the-hook)
    - [Configure emulation parameters via environment variables](#configure-emulation-parameters-via-environment-variables)
    - [Handle-local settings and environment-variable overrides](#handle-local-settings-and-environment-variable-overrides)
    - [Max-workspace preallocation](#max-workspace-preallocation)
    - [Hook workspace, memory-saving, and skip-scaling behavior](#hook-workspace-memory-saving-and-skip-scaling-behavior)
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

As a practical rule of thumb, the following settings typically provide accuracy comparable to cuBLAS INT8-based fixed-point emulation with `mantissaBitCount = 55`, corresponding to INT8-based Ozaki Scheme I with 7 slices.

| Backend | `num_moduli` | `fastmode`         |
| :------ | :----------- | :----------------- |
| INT8    | 14 or 15     | `true` (fast mode) |
| FP8     | 10 or 11     | `true` (fast mode) |

> [!NOTE]
>
> These values are practical starting points, not accuracy guarantees.
> The required number of moduli depends on the input matrices and the target application.

## Requirements

- Linux
- GNU Make 3.81 or later
- C++20-capable C++ compiler

- NVIDIA CUDA backend:
  - CUDA Toolkit 12.9 or later
  - cuBLAS and cuBLASLt
  - NVML and cuRAND for building the test programs

- AMD HIP backend:
  - ROCm 7.0 or later
  - hipBLAS and hipBLASLt
  - AMD SMI and hipRAND for building the test programs

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

The Makefile automatically detects the appropriate build configuration.
If `make` fails, try setting the following [`make` options](#make-options).

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

### Memory saving options

| Option              | Value    | Description                   |
| :------------------ | :------- | :---------------------------- |
| `memory_saving=...` | `0`, `1` | `0` = disabled; `1` = enabled |
| `max_memory=...`    |          | workspace-size limit in bytes |

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

# Run only non-transposed FP64 GEMM timing tests with workspace-size limit: 8 GiB
make run MODE="time_square GEMM D transA=N transB=N memory_saving=1 max_memory=8589934592"

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

Handle-local execution settings are declared in `include/config.hpp`.
For CUDA, the relevant APIs are:

```cpp
// cuBLAS handle
void gemmul8::set_memory_saving(cublasHandle_t handle, bool enable) noexcept;
bool gemmul8::get_memory_saving(cublasHandle_t handle) noexcept;

void gemmul8::set_max_worksize(cublasHandle_t handle, size_t bytes) noexcept;
size_t gemmul8::get_max_worksize(cublasHandle_t handle) noexcept;

void gemmul8::set_block_size_trsm(cublasHandle_t handle, int nB) noexcept;
int gemmul8::get_block_size_trsm(cublasHandle_t handle) noexcept;

void gemmul8::clear_config(cublasHandle_t handle) noexcept;

// cuBLASLt handle
void gemmul8::set_memory_savingLt(cublasLtHandle_t handle, bool enable) noexcept;
bool gemmul8::get_memory_savingLt(cublasLtHandle_t handle) noexcept;

void gemmul8::set_max_worksizeLt(cublasLtHandle_t handle, size_t bytes) noexcept;
size_t gemmul8::get_max_worksizeLt(cublasLtHandle_t handle) noexcept;

void gemmul8::set_block_size_trsmLt(cublasLtHandle_t handle, int nB) noexcept;
int gemmul8::get_block_size_trsmLt(cublasLtHandle_t handle) noexcept;

void gemmul8::clear_configLt(cublasLtHandle_t handle) noexcept;
```

The corresponding HIP APIs use `hipblasHandle_t` / `hipblasLtHandle_t`.

> [!IMPORTANT]
>
> `set_block_size_trsm()` is handle-local and requires a BLAS handle.
> This is a breaking API change from GEMMul8 v3.2.0 and earlier releases, in which the TRSM block size was configured globally.

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

#### Handle-local execution settings

The memory-saving flag, workspace-size limit, and TRSM block-size override are associated with individual BLAS/Lt handles.
Different handles can therefore use different settings concurrently.

The default settings are:

- memory saving: disabled;
- maximum workspace size: 12 GiB;
- TRSM block-size override: non-positive, which selects automatic block-size selection.

A getter returns the setting currently associated with the specified handle.
For the TRSM block size, a non-positive value means that automatic architecture/backend-dependent selection is enabled; the getter does not report the automatically selected internal block size.

> [!CAUTION]
>
> In direct mode, if GEMMul8-specific handle-local settings have been used, call `clear_config(handle)` or `clear_configLt(handle)` before destroying the corresponding BLAS/Lt handle.
>
> In hook mode, GEMMul8 removes the configuration associated with an intercepted BLAS handle automatically when the handle is destroyed.

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
- `gemmul8::workSizeTrsmLt`

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

For `trsm`, use `workSizeTrsm(handle, ...)` for a standard BLAS handle or `workSizeTrsmLt(handle, ...)` for an Lt handle.
The returned workspace size depends on `side`, `m`, `n`, `num_moduli`, the selected backend, the element type, and the TRSM block-size setting associated with the specified handle.

When using a custom TRSM block size, configure the same handle before querying the workspace:

```cpp
gemmul8::set_block_size_trsm(handle, 2048);

const size_t worksize =
    gemmul8::workSizeTrsm<double, gemmul8::Backend::INT8>(
        handle, side, m, n, num_moduli);
```

> [!NOTE]
>
> The workspace-size limit configured by `set_max_worksize()` / `set_max_worksizeLt()` does **not** change the value returned by `workSize()`, `workSizeTrsm()`, or `workSizeTrsmLt()`.
> These query functions report the default/full workspace requirement.

#### Memory-saving mode

GEMMul8 can reduce the workspace required by BLAS-like operations by internally blocking the operation.
Memory saving is configured independently for each BLAS/Lt handle:

```cpp
gemmul8::set_memory_saving(handle, true);
gemmul8::set_max_worksize(handle, size_t(4) << 30); // 4 GiB
```

When memory saving is enabled:

- if the default/full workspace requirement does not exceed the configured limit, GEMMul8 uses the normal unblocked execution path;
- if the default/full workspace requirement exceeds the configured limit, GEMMul8 selects block sizes internally and executes the operation as a sequence of smaller BLAS-like operations;
- skip scaling/reuse is disabled, regardless of the `enable_skip_scalA`, `enable_skip_scalB`, `skip_scalA`, and `skip_scalB` arguments.

If internal blocking is required, `workA` and `workB` are not used by the blocked path.
The block workspace is taken entirely from `work`.

In general, a larger workspace-size limit allows larger blocks and can provide better performance.

A workspace-size limit of zero disables workspace limiting:

```cpp
gemmul8::set_max_worksize(handle, 0);
```

In this case, GEMMul8 uses the default/full workspace even if memory saving is enabled.
Skip scaling/reuse remains disabled while memory saving itself is enabled.

For a simple direct-mode memory-saving allocation, the caller can allocate a single `work` buffer:

```cpp
const size_t limit = size_t(4) << 30;
gemmul8::set_memory_saving(handle, true);
gemmul8::set_max_worksize(handle, limit);

const size_t full_worksize = gemmul8::workSize(m, n, k, num_moduli);

const size_t allocated_worksize = std::min(full_worksize, limit);

void *work = nullptr;
cudaMalloc(&work, allocated_worksize);

gemmul8::gemm(
    handle,
    transA, transB,
    m, n, k,
    &alpha, A, lda,
    B, ldb,
    &beta, C, ldc,
    num_moduli, fastmode,
    work);

cudaFree(work);
```

If `set_max_worksize(handle, 0)` is used instead, allocate the full workspace.

#### TRSM implementation and block-size control

`gemmul8::trsm` and `gemmul8::trsmLt` use a blocked triangular-solve algorithm internally.

The implementation combines:

- standard cuBLAS/hipBLAS TRSM for triangular solves on diagonal blocks, and
- `gemmul8::gemm` / `gemmul8::gemmLt` for updates to the remaining blocks.

The internal TRSM block size is configured independently for each handle:

```cpp
gemmul8::set_block_size_trsm(handle, nB);
gemmul8::set_block_size_trsmLt(ltHandle, nB);
```

A positive value selects the specified block size.
A non-positive value enables automatic architecture/backend-dependent selection.

The corresponding getter returns the explicitly configured value for that handle.
It does not return the automatically selected internal block size.

> [!NOTE]
>
> The automatically selected TRSM block size is a heuristic default and is not guaranteed to be the fastest setting.
> For performance tuning, benchmark several block sizes and set a custom value for the corresponding handle.

The block-size setting also affects the workspace size returned by `workSizeTrsm()` / `workSizeTrsmLt()`.
Therefore, when using a custom block size, set it before querying and allocating the TRSM workspace.

When memory saving is enabled, the existing TRSM solve blocking is retained.
The GEMM updates are executed through GEMMul8's GEMM path and are themselves memory-saving blocked when their workspace requirement exceeds the configured limit.

#### Behavior of `skip_scalA` / `skip_scalB`

This skip mechanism is designed for repeated calls to GEMMul8 routines that reuse the same input matrix `A` and/or `B`.
It applies to routines that expose `workA`, `workB`, `enable_skip_scalA`, `enable_skip_scalB`, `skip_scalA`, and `skip_scalB`.

> [!NOTE]
>
> It does not apply to `trsm`, because the current `trsm` interface does not expose skip-scaling arguments.

> [!IMPORTANT]
>
> Skip scaling/reuse is disabled whenever memory saving is enabled for the corresponding handle, even when the configured workspace-size limit is zero.

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

The following three variables are global hook-mode overrides of handle-local execution settings:

```text
GEMMUL8_MEMORY_SAVING
GEMMUL8_MAX_WORKSIZE
GEMMUL8_BLK_SIZE_TRSM
```

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

# Memory-saving override
export GEMMUL8_MEMORY_SAVING=1

# Workspace-size limit in bytes: 4 GiB
export GEMMUL8_MAX_WORKSIZE=4294967296

# Global TRSM block-size override
# 0 (default): automatic architecture/backend-dependent selection
# >0: use the specified block size
export GEMMUL8_BLK_SIZE_TRSM=2048

# Global skip-scaling switches
export GEMMUL8_SKIP_SCALE_A=1
export GEMMUL8_SKIP_SCALE_B=1
```

| Variable pattern          | Default | Description                                                                                                               |
| :------------------------ | :------ | :------------------------------------------------------------------------------------------------------------------------ |
| `GEMMUL8_BACKEND_<OP>`    | `INT8`  | Selects the emulation backend. `0` or `INT8` = INT8 backend; `1` or `FP8` = FP8 backend.                                  |
| `GEMMUL8_NUM_MOD_S_<OP>`  | `0`     | Number of moduli for FP32 real routines. Native BLAS is used if outside `[2, 13]`.                                        |
| `GEMMUL8_NUM_MOD_D_<OP>`  | `0`     | Number of moduli for FP64 real routines. Native BLAS is used if outside `[2, 20]`.                                        |
| `GEMMUL8_NUM_MOD_C_<OP>`  | `0`     | Number of moduli for FP32 complex routines. Native BLAS is used if outside `[2, 13]`.                                     |
| `GEMMUL8_NUM_MOD_Z_<OP>`  | `0`     | Number of moduli for FP64 complex routines. Native BLAS is used if outside `[2, 20]`.                                     |
| `GEMMUL8_FASTMODE_S_<OP>` | `1`     | Fast mode switch for FP32 real routines. `1` = fast mode; `0` = accurate mode.                                            |
| `GEMMUL8_FASTMODE_D_<OP>` | `1`     | Fast mode switch for FP64 real routines. `1` = fast mode; `0` = accurate mode.                                            |
| `GEMMUL8_FASTMODE_C_<OP>` | `1`     | Fast mode switch for FP32 complex routines. `1` = fast mode; `0` = accurate mode.                                         |
| `GEMMUL8_FASTMODE_Z_<OP>` | `1`     | Fast mode switch for FP64 complex routines. `1` = fast mode; `0` = accurate mode.                                         |
| `GEMMUL8_MEMORY_SAVING`   | unset   | Overrides the handle-local memory-saving setting when defined. `1` = enabled; `0` = disabled.                             |
| `GEMMUL8_MAX_WORKSIZE`    | unset   | Overrides the handle-local workspace-size limit when defined. Value is in bytes. `0` disables workspace limiting.         |
| `GEMMUL8_BLK_SIZE_TRSM`   | unset   | Overrides the handle-local TRSM block size when defined. `>0` uses the specified size; `<=0` enables automatic selection. |
| `GEMMUL8_SKIP_SCALE_A`    | `0`     | Global switch that enables reuse of preprocessed/scaled `A` when the operand cache key matches.                           |
| `GEMMUL8_SKIP_SCALE_B`    | `0`     | Global switch that enables reuse of preprocessed/scaled `B` when the operand cache key matches.                           |

#### Handle-local settings and environment-variable overrides

`GEMMUL8_MEMORY_SAVING`, `GEMMUL8_MAX_WORKSIZE`, and `GEMMUL8_BLK_SIZE_TRSM` override the corresponding handle-local setting only when the environment variable is explicitly defined.

If an environment variable is not defined, the current value configured through the public `gemmul8::set_*()` API is used.

When defined, the environment variable has precedence for matching hook calls.
These configuration variables are read repeatedly rather than being cached at process initialization:

- `GEMMUL8_MEMORY_SAVING` and `GEMMUL8_MAX_WORKSIZE` are read on intercepted GEMMul8 emulation calls;
- `GEMMUL8_BLK_SIZE_TRSM` is read on intercepted TRSM calls before the TRSM workspace is queried.

Therefore, changing one of these values during program execution takes effect from the next matching call.

> [!NOTE]
>
> The value of an explicitly defined configuration environment variable is copied into the handle-local configuration.
> Unsetting the environment variable does not restore the value that was present before the override.
> To change or reset the setting, explicitly set the desired value through the environment variable or the public setter.

Memory saving disables skip-scaling/reuse.
Therefore, while `GEMMUL8_MEMORY_SAVING=1` or the corresponding handle-local memory-saving setting is enabled, `GEMMUL8_SKIP_SCALE_A` and `GEMMUL8_SKIP_SCALE_B` do not cause preprocessing reuse.

#### Max-workspace preallocation

GEMMul8 normally grows hook workspaces on demand. To stabilize workspace addresses, avoid reallocating workspace, and improve skip-scaling reuse, define the maximum BLAS size arguments for the operations that will be used.

| Operation suffix | Required size variables                                          | Internal workspace query |
| :--------------- | :--------------------------------------------------------------- | :----------------------- |
| `GEMM`           | `GEMMUL8_MAX_M_GEMM`, `GEMMUL8_MAX_N_GEMM`, `GEMMUL8_MAX_K_GEMM` | `workSize(m, n, k, ...)` |
| `SYMM_LEFT`      | `GEMMUL8_MAX_M_SYMM_LEFT`, `GEMMUL8_MAX_N_SYMM_LEFT`             | `workSize(m, n, m, ...)` |
| `SYMM_RIGHT`     | `GEMMUL8_MAX_M_SYMM_RIGHT`, `GEMMUL8_MAX_N_SYMM_RIGHT`           | `workSize(m, n, n, ...)` |
| `SYRK`           | `GEMMUL8_MAX_N_SYRK`, `GEMMUL8_MAX_K_SYRK`                       | `workSize(n, n, k, ...)` |
| `SYR2K`          | `GEMMUL8_MAX_N_SYR2K`, `GEMMUL8_MAX_K_SYR2K`                     | `workSize(n, n, k, ...)` |
| `SYRKX`          | `GEMMUL8_MAX_N_SYRKX`, `GEMMUL8_MAX_K_SYRKX`                     | `workSize(n, n, k, ...)` |
| `HEMM_LEFT`      | `GEMMUL8_MAX_M_HEMM_LEFT`, `GEMMUL8_MAX_N_HEMM_LEFT`             | `workSize(m, n, m, ...)` |
| `HEMM_RIGHT`     | `GEMMUL8_MAX_M_HEMM_RIGHT`, `GEMMUL8_MAX_N_HEMM_RIGHT`           | `workSize(m, n, n, ...)` |
| `HERK`           | `GEMMUL8_MAX_N_HERK`, `GEMMUL8_MAX_K_HERK`                       | `workSize(n, n, k, ...)` |
| `HER2K`          | `GEMMUL8_MAX_N_HER2K`, `GEMMUL8_MAX_K_HER2K`                     | `workSize(n, n, k, ...)` |
| `HERKX`          | `GEMMUL8_MAX_N_HERKX`, `GEMMUL8_MAX_K_HERKX`                     | `workSize(n, n, k, ...)` |
| `TRMM_LEFT`      | `GEMMUL8_MAX_M_TRMM_LEFT`, `GEMMUL8_MAX_N_TRMM_LEFT`             | `workSize(m, n, m, ...)` |
| `TRMM_RIGHT`     | `GEMMUL8_MAX_M_TRMM_RIGHT`, `GEMMUL8_MAX_N_TRMM_RIGHT`           | `workSize(m, n, n, ...)` |
| `TRSM_LEFT`      | `GEMMUL8_MAX_M_TRSM_LEFT`, `GEMMUL8_MAX_N_TRSM_LEFT`             | TRSM workspace query     |
| `TRSM_RIGHT`     | `GEMMUL8_MAX_M_TRSM_RIGHT`, `GEMMUL8_MAX_N_TRSM_RIGHT`           | TRSM workspace query     |

Additional max-workspace variables are also operation-specific:

| Variable pattern             | Default | Description                                                                                                        |
| :--------------------------- | :------ | :----------------------------------------------------------------------------------------------------------------- |
| `GEMMUL8_MAXWS_BACKEND_<OP>` | `INT8`  | Backend used for max-workspace calculation. `0` or `INT8` = INT8, `1` or `FP8` = FP8, `2` or `BOTH` = max of both. |
| `GEMMUL8_MAX_NUM_MOD_<OP>`   | `2`     | Number of moduli used for max-workspace calculation.                                                               |

> [!NOTE]
>
> `GEMMUL8_MAX_WORKSIZE` and the max-workspace preallocation variables have different purposes.
>
> - `GEMMUL8_MAX_WORKSIZE` is the workspace-size limit used by the memory-saving mechanism.
> - `GEMMUL8_MAX_M_*`, `GEMMUL8_MAX_N_*`, `GEMMUL8_MAX_K_*`, `GEMMUL8_MAX_NUM_MOD_*`, and `GEMMUL8_MAXWS_BACKEND_*` are used to preallocate hook workspaces for selected maximum problem sizes.

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

#### Hook workspace, memory-saving, and skip-scaling behavior

Hook mode maintains an independent workspace per BLAS handle (`cublasHandle_t` / `hipblasHandle_t`).

For each handle, the hook allocates GPU work buffers used by the emulation routine.

- For routines that support skip scaling, the hook may keep separate `workA` and/or `workB` cache areas for preprocessed input matrices.
- The remaining workspace is used as the routine's internal work buffer.
- For `trsm`, the direct interface uses a single workspace and does not expose `workA` or `workB`.

Under normal execution, hook workspaces grow on demand and are reused across calls.

When memory saving is active, however, GEMMul8 may release or reallocate existing work buffers in order to satisfy the configured workspace-size limit.
If internal blocking is required, the dedicated `workA` and `workB` buffers are not used by the blocked path; the block workspace is allocated as a single internal work buffer.

Allocation/free use stream-ordered APIs (`cudaMallocAsync/cudaFreeAsync` or HIP equivalents) on the current stream.
When the same handle is used with different CUDA/HIP streams across calls, the hook enforces ordering by inserting an event dependency (`eventRecord` on the previous stream -> `streamWaitEvent` on the current stream).

The workspaces are released when the corresponding handle is destroyed.

> [!IMPORTANT]
>
> Skip scaling/reuse is disabled while memory saving is enabled.
> Therefore, `GEMMUL8_SKIP_SCALE_A=1` and `GEMMUL8_SKIP_SCALE_B=1` have no preprocessing-reuse effect while memory saving is enabled for the corresponding handle.
>
> When memory saving is disabled, `GEMMUL8_SKIP_SCALE_A=1` and/or `GEMMUL8_SKIP_SCALE_B=1` enables automatic reuse of already-preprocessed intermediate data for `A` and/or `B` within the hook, when it is safe according to the cache conditions below.
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
> To keep internal workspace pointers stable across calls when memory saving is disabled, define the operation-specific maximum-size variables for the operations that will be used.
> If you may switch backend for an operation at runtime, set `GEMMUL8_MAXWS_BACKEND_<OP>=BOTH`.

> [!CAUTION]
>
> Skip scaling assumes that the contents of `A` or `B` remain unchanged in GPU memory.
> If `A` or `B` data are modified between routine calls, do not rely on skipping.

> [!NOTE]
>
> `GEMMUL8_MAX_*_<OP>`, `GEMMUL8_MAXWS_BACKEND_<OP>`, and `GEMMUL8_MAX_NUM_MOD_<OP>` are read only once on first hook use to compute the maximum workspace sizes.
>
> Runtime variables such as `GEMMUL8_NUM_MOD_<S|D|C|Z>_<OP>`, `GEMMUL8_FASTMODE_<S|D|C|Z>_<OP>`, `GEMMUL8_BACKEND_<OP>`, and `GEMMUL8_SKIP_SCALE_*` are read at each intercepted routine call.
>
> `GEMMUL8_MEMORY_SAVING`, `GEMMUL8_MAX_WORKSIZE`, and `GEMMUL8_BLK_SIZE_TRSM` are also read dynamically as described above.

#### How to change environment variables programmatically

You can set these environment variables programmatically using `setenv()`.

```cpp
// Run GEMM emulation with Backend = INT8, num_moduli = 15 & fastmode = true
char num_moduli[12];
snprintf(num_moduli, sizeof(num_moduli), "%u", 15u);

setenv("GEMMUL8_BACKEND_GEMM", "INT8", 1);
setenv("GEMMUL8_NUM_MOD_D_GEMM", num_moduli, 1);
setenv("GEMMUL8_FASTMODE_D_GEMM", "1", 1);

cublasDgemm_v2(...);
```

Memory-saving settings can also be changed during execution:

```cpp
setenv("GEMMUL8_MEMORY_SAVING", "1", 1);
setenv("GEMMUL8_MAX_WORKSIZE", "4294967296", 1); // 4 GiB

cublasDgemm_v2(...); // 4-GiB workspace-size limit

setenv("GEMMUL8_MAX_WORKSIZE", "8589934592", 1); // 8 GiB

cublasDgemm_v2(...); // next call uses the 8-GiB limit
```

The TRSM block-size override is likewise reread:

```cpp
setenv("GEMMUL8_BLK_SIZE_TRSM", "1024", 1);
cublasDtrsm_v2(...); // block size 1024

setenv("GEMMUL8_BLK_SIZE_TRSM", "2048", 1);
cublasDtrsm_v2(...); // next TRSM call uses 2048

setenv("GEMMUL8_BLK_SIZE_TRSM", "0", 1);
cublasDtrsm_v2(...); // automatic block-size selection
```

The public handle-local setters can also be changed during execution when the corresponding environment variable is not defined:

```cpp
gemmul8::set_memory_saving(handle, true);
gemmul8::set_max_worksize(handle, size_t(4) << 30);
cublasDgemm_v2(...);

gemmul8::set_max_worksize(handle, size_t(8) << 30);
cublasDgemm_v2(...);
```

If a corresponding environment variable is explicitly defined, that environment-variable value overrides the setter value for matching hook calls.

## Numerical results

See numerical results in the separate repository: [GEMMul8_numerical_results](https://github.com/UCHINO-Yuki/GEMMul8_numerical_results)

## Acknowledgment

> [!CAUTION]
> Please do not contact the individuals listed below regarding this code.

### Assistance with debugging

- Patrick Gutsche (École Normale Supérieure de Lyon, France; affiliation as of 2025)
- Prajval Kumar (Indian Institute of Science and Education Research, India; affiliation as of 2025)
- Dr. William Dawson (RIKEN Center for Computational Science, Japan; affiliation as of 2025)
- Dr. Toshiyuki Imamura (RIKEN Center for Computational Science, Japan; affiliation as of 2025)

### Assistance with preliminary experiments

The following individuals helped conduct preliminary performance experiments on B200 systems at Yokota Lab:

- Dr. Qianxiang Ma (RIKEN Center for Computational Science, Japan; affiliation as of 2025)
- Prof. Rio Yokota (Institute of Science Tokyo, Japan; affiliation as of 2025)

The following individuals helped conduct preliminary experiments on the B200 environment of SAKURAONE, SAKURA internet Inc.'s managed HPC cluster service:

- Takeshi Yamashita (SAKURA internet Inc., Japan; affiliation as of 2025)
- Fumikazu Konishi (SAKURA internet Inc., Japan; affiliation as of 2025)

## Contact (Responsible Developer)

- Yuki Uchino (RIKEN Center for Computational Science, Japan)
- yuki.uchino.fe (at) riken.jp

## References

- Ootomo, H., & Yokota, R. (2022). Recovering single precision accuracy from Tensor Cores while surpassing the FP32 theoretical peak performance. The International Journal of High Performance Computing Applications, 36(4), 475-491, [doi.org/10.1177/10943420221090256](https://doi.org/10.1177/10943420221090256).
- Ootomo, H., Manabe, H., Harada, K., & Yokota, R. (2023). Quantum Circuit Simulation by SGEMM Emulation on Tensor Cores and Automatic Precision Selection. In High Performance Computing (pp. 259-276). Springer, [doi.org/10.1007/978-3-031-32041-5_14](https://doi.org/10.1007/978-3-031-32041-5_14).
- Ootomo, H., Ozaki, K., & Yokota, R. (2024). DGEMM on integer matrix multiplication unit. The International Journal of High Performance Computing Applications, 38(4), 297-313, [https://doi.org/10.1177/10943420241239588](https://doi.org/10.1177/10943420241239588).
- Uchino, Y., Ozaki, K., & Imamura, T. (2025). Performance enhancement of the Ozaki Scheme on integer matrix multiplication unit. The International Journal of High Performance Computing Applications, 39(3), 462-476, [doi.org/10.1177/10943420241313064](https://doi.org/10.1177/10943420241313064).
- Kawakami S. & Takahashi D. (2026). Improved Scaling for Fast Mode of Ozaki Scheme II, [doi.org/10.48550/arXiv.2606.29129](https://doi.org/10.48550/arXiv.2606.29129).
- Kawakami S. (2026). GEMMul8 (fork with improved fast mode scaling), GitHub, [https://github.com/kotatsumuri/GEMMul8](https://github.com/kotatsumuri/GEMMul8).
- Hayashi S., Mukunoki D., Hoshino T., Katagiri T. (2026). DGEMM with Ozaki Scheme I/II on FP4 Tensor Cores: A Base-13 E2M1 Limb Representation, [doi.org/10.48550/arXiv.2608.06812](https://doi.org/10.48550/arXiv.2608.06812).

## Citations

> [!NOTE]
>
> If you refer to the algorithm used in the fast mode, please also cite the following work:
> Kawakami S. & Takahashi D. (2026). Improved Scaling for Fast Mode of Ozaki Scheme II, [doi.org/10.48550/arXiv.2606.29129](https://doi.org/10.48550/arXiv.2606.29129).

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
@article{doi:10.1177/10943420261467787,
    author = {Katsuhisa Ozaki and Yuki Uchino and Toshiyuki Imamura},
    title ={Ozaki scheme II: A GEMM-oriented emulation of floating-point matrix multiplication using an integer modular technique},
    journal = {The International Journal of High Performance Computing Applications},
    year = {2026},
    doi = {10.1177/10943420261467787},
    URL = {https://doi.org/10.1177/10943420261467787},
    note = {OnlineFirst},
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
