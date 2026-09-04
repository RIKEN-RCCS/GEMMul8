#pragma once
#include "../common/common.hpp"

namespace bench::accuracy::trtrmm {

template <typename T>
void check_accuracy(
    std::string &deviceName,
    std::string &dateTime,
    size_t memory_limit,
    cublasFillMode_t uplo_A,
    cublasFillMode_t uplo_B,
    cublasOperation_t trans_A,
    cublasOperation_t trans_B,
    cublasDiagType_t diag_A,
    cublasDiagType_t diag_B,
    const bool run_Ozaki2_I8,
    const bool run_Ozaki2_F8,
    const bool run_Ozaki1_I8);

} // namespace bench::accuracy::trtrmm
