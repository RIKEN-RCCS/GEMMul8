#pragma once
#include "../common/common.hpp"

namespace bench::time::her2k {

template <typename T>
void check_time(
    std::string &deviceName,
    std::string &dateTime,
    size_t memory_limit,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    const bool run_Ozaki2_I8,
    const bool run_Ozaki2_F8,
    const bool run_Ozaki1_I8,
    const bool is_square = false //
);

} // namespace bench::time::her2k
