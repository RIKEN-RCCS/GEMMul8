#pragma once
#include "../common/common.hpp"

namespace bench::time::gemm {

template <typename T>
void check_time(
    std::string &deviceName,
    std::string &dateTime,
    size_t memory_limit,
    cublasOperation_t transa,
    cublasOperation_t transb,
    const bool run_Ozaki2_I8,
    const bool run_Ozaki2_F8,
    const bool run_Ozaki1_I8,
    const bool is_square = false //
);

} // namespace bench::time::gemm
