#pragma once
#include "../common/common.hpp"

namespace bench::time::symm {

template <typename T>
void check_time(
    std::string &deviceName,
    std::string &dateTime,
    size_t memory_limit,
    cublasFillMode_t uplo,
    cublasSideMode_t side,
    const bool run_Ozaki2_I8,
    const bool run_Ozaki2_F8,
    const bool run_cuBLAS_FP64_emu,
    const bool is_square = false //
);

} // namespace bench::time::symm
