#pragma once
#include "../oz2/common/include.hpp"
#include "../../include/config.hpp"

namespace gemmul8::config {

inline constexpr bool DEFAULT_MEMORY_SAVING  = false;
inline constexpr size_t DEFAULT_MAX_WORKSIZE = size_t(12) << 30; // 12 GiB
inline constexpr int DEFAULT_BLOCK_SIZE_TRSM = 0;

struct ConfigSnapshot {
    bool memory_saving  = DEFAULT_MEMORY_SAVING;
    size_t max_worksize = DEFAULT_MAX_WORKSIZE;
    int block_size_trsm = DEFAULT_BLOCK_SIZE_TRSM;
};

ConfigSnapshot get_config(cublasHandle_t handle) noexcept;
ConfigSnapshot get_configLt(cublasLtHandle_t handle) noexcept;

void set_memory_saving_impl(cublasHandle_t handle, bool enable) noexcept;
void set_memory_savingLt_impl(cublasLtHandle_t handle, bool enable) noexcept;

void set_max_worksize_impl(cublasHandle_t handle, size_t bytes) noexcept;
void set_max_worksizeLt_impl(cublasLtHandle_t handle, size_t bytes) noexcept;

void set_block_size_trsm_impl(cublasHandle_t handle, int nB) noexcept;
void set_block_size_trsmLt_impl(cublasLtHandle_t handle, int nB) noexcept;

void clear_config_impl(cublasHandle_t handle) noexcept;
void clear_configLt_impl(cublasLtHandle_t handle) noexcept;

void bind_config(cublasHandle_t parent, cublasLtHandle_t child) noexcept;
void unbind_configLt(cublasLtHandle_t child) noexcept;

} // namespace gemmul8::config
