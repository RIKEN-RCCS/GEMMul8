#pragma once
#include "../common/common.hpp"

namespace bench::accuracy::trmm {

template <typename T>
void check_accuracy(
    std::string &deviceName,
    std::string &dateTime,
    size_t memory_limit,
    cublasSideMode_t side,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    cublasDiagType_t diag,
    const bool run_Ozaki2_I8,
    const bool run_Ozaki2_F8,
    const bool run_Ozaki1_I8,
    const bool is_square = false);

} // namespace bench::accuracy::trmm
