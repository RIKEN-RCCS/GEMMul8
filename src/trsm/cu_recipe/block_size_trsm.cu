#include "../../../include/trsm.hpp"
#include <atomic>

namespace gemmul8 {

namespace {

std::atomic<int> &trsm_block_size_override_storage() noexcept {
    static std::atomic<int> nB = 0;
    return nB;
}

} // namespace

void set_block_size_trsm(const int nB) noexcept {
    trsm_block_size_override_storage().store(nB, std::memory_order_relaxed);
}

int get_block_size_trsm() noexcept {
    return trsm_block_size_override_storage().load(std::memory_order_relaxed);
}

} // namespace gemmul8
