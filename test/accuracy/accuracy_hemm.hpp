#pragma once
#include "../common/common.hpp"

namespace bench::accuracy::hemm {

template <typename T>
void check_accuracy(
    std::string &deviceName,
    std::string &dateTime,
    size_t memory_limit,
    cublasFillMode_t uplo,
    cublasSideMode_t side,
    const bool run_Ozaki2_I8,
    const bool run_Ozaki2_F8,
    const bool run_Ozaki1_I8,
    const bool is_square = false);

} // namespace bench::accuracy::hemm
