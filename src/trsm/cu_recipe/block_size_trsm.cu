#include "../../config/config.hpp"

namespace gemmul8 {

void set_block_size_trsm(cublasHandle_t handle, int nB) noexcept {
    config::set_block_size_trsm_impl(handle, nB);
}

void set_block_size_trsmLt(cublasLtHandle_t handle, int nB) noexcept {
    config::set_block_size_trsmLt_impl(handle, nB);
}

int get_block_size_trsm(cublasHandle_t handle) noexcept {
    return config::get_config(handle).block_size_trsm;
}

int get_block_size_trsmLt(cublasLtHandle_t handle) noexcept {
    return config::get_configLt(handle).block_size_trsm;
}

} // namespace gemmul8
