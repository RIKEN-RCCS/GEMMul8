/**
 * Handle-local execution configuration
 * ------------------------------------
 * GEMMul8 associates execution settings with individual BLAS handles.
 * The settings are thread-safe and may be configured independently for
 * cuBLAS/cuBLASLt (or hipBLAS/hipBLASLt) handles.
 *
 * Direct mode lifetime note:
 *   If any handle-local setting is used, call clear_config(handle) for a
 *   cuBLAS/hipBLAS handle, or clear_configLt(handle) for a
 *   cuBLASLt/hipBLASLt handle, before destroying the corresponding handle.
 *   In hook mode, GEMMul8 clears the configuration automatically when the
 *   parent cublasDestroy()/hipblasDestroy() is intercepted.
 */
#pragma once
#include "types.hpp"

namespace gemmul8 {

//------------------------------
// CUDA
//------------------------------
#if defined(__CUDACC__)

/**
 * Enable or disable workspace-size reduction.
 *
 * If enable = false, GEMMul8 uses the default workspace size regardless
 * of the value set by set_max_worksize().
 *
 * If enable = true, GEMMul8 reduces the workspace size as necessary so
 * that it does not exceed the limit set by set_max_worksize().
 *
 * When memory saving is enabled, skip-scaling/reuse is disabled for the
 * corresponding handle, even if the skip-scaling arguments or hook-mode
 * environment variables request it.
 *
 * Memory saving is disabled by default.
 */
void set_memory_saving(cublasHandle_t handle, bool enable) noexcept;
void set_memory_savingLt(cublasLtHandle_t handle, bool enable) noexcept;

/**
 * Return whether workspace-size reduction is enabled.
 */
bool get_memory_saving(cublasHandle_t handle) noexcept;
bool get_memory_savingLt(cublasLtHandle_t handle) noexcept;

/**
 * Set the workspace-size limit used when memory saving is enabled.
 *
 * If bytes = 0, GEMMul8 uses the default workspace size regardless of whether
 * memory saving is enabled.
 *
 * Operations with a memory-saving blocked path internally use blocking to
 * satisfy the workspace-size limit.  In general, a larger value of bytes
 * allows more efficient execution by using larger blocks.
 *
 * The value does not affect the return values of workSize(), workSizeTrsm(),
 * or workSizeTrsmLt().
 * The default workspace-size limit is 12 GiB.
 */
void set_max_worksize(cublasHandle_t handle, size_t bytes) noexcept;
void set_max_worksizeLt(cublasLtHandle_t handle, size_t bytes) noexcept;

/**
 * Return the current workspace-size limit in bytes.
 */
size_t get_max_worksize(cublasHandle_t handle) noexcept;
size_t get_max_worksizeLt(cublasLtHandle_t handle) noexcept;

/**
 * Override the internal block size used by GEMMul8 TRSM for this handle.
 *
 * GEMMul8 TRSM uses a blocked algorithm internally.  By default, the block
 * size is selected automatically from the detected GPU architecture and
 * backend.  This function overrides that selection for subsequent trsm() and
 * trsmLt() calls.
 *
 * If nB > 0, the specified value is used as the TRSM block size.
 * If nB <= 0, the automatic architecture/backend-dependent block size is used.
 *
 * This setting also affects the workspace size returned by workSizeTrsm()
 * or workSizeTrsmLt(). Therefore, when overriding the block size, call
 * set_block_size_trsm(handle, nB) or set_block_size_trsmLt(handle, nB)
 * before querying and allocating the workspace.
 */
void set_block_size_trsm(cublasHandle_t handle, int nB) noexcept;
void set_block_size_trsmLt(cublasLtHandle_t handle, int nB) noexcept;

/**
 * Return the current TRSM block-size override.
 *
 * A positive value means that the returned value is used as the TRSM block
 * size for subsequent trsm() and trsmLt() calls.  A non-positive value means
 * that the automatic architecture/backend-dependent block size is used.
 */
int get_block_size_trsm(cublasHandle_t handle) noexcept;
int get_block_size_trsmLt(cublasLtHandle_t handle) noexcept;

/**
 * Remove all GEMMul8 configuration associated with handle.
 *
 * In direct mode, call this before destroying a handle for which GEMMul8
 * handle-local settings have been used.
 */
void clear_config(cublasHandle_t handle) noexcept;
void clear_configLt(cublasLtHandle_t handle) noexcept;

#endif

//------------------------------
// HIP
//------------------------------
#if defined(__HIPCC__)

void set_memory_saving(hipblasHandle_t handle, bool enable) noexcept;
void set_memory_savingLt(hipblasLtHandle_t handle, bool enable) noexcept;
bool get_memory_saving(hipblasHandle_t handle) noexcept;
bool get_memory_savingLt(hipblasLtHandle_t handle) noexcept;

void set_max_worksize(hipblasHandle_t handle, size_t bytes) noexcept;
void set_max_worksizeLt(hipblasLtHandle_t handle, size_t bytes) noexcept;
size_t get_max_worksize(hipblasHandle_t handle) noexcept;
size_t get_max_worksizeLt(hipblasLtHandle_t handle) noexcept;

void set_block_size_trsm(hipblasHandle_t handle, int nB) noexcept;
void set_block_size_trsmLt(hipblasLtHandle_t handle, int nB) noexcept;
int get_block_size_trsm(hipblasHandle_t handle) noexcept;
int get_block_size_trsmLt(hipblasLtHandle_t handle) noexcept;

void clear_config(hipblasHandle_t handle) noexcept;
void clear_configLt(hipblasLtHandle_t handle) noexcept;

#endif

} // namespace gemmul8
