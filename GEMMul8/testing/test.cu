#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <chrono>
#include <type_traits>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdint>

#include "../include/gemmul8.hpp"
#include "self_hipify.hpp"
#include "eval.hpp"
#include "make_matrix.hpp"

#include "common.hpp"
#include "test_accuracy.hpp"
#include "test_flops.hpp"
#include "test_watt.hpp"

int main(int argc, char **argv) {
    std::chrono::system_clock::time_point start, stop;
    std::string deviceName = getDeviceName();
    std::string startTime  = getCurrentDateTime(start);

    bool run_accuracy = false;
    bool run_flops    = false;
    bool run_watt     = false;
    bool run_SGEMM    = false;
    bool run_DGEMM    = false;
    bool run_CGEMM    = false;
    bool run_ZGEMM    = false;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "accuracy") {
            run_accuracy = true;
        } else if (arg == "flops") {
            run_flops = true;
        } else if (arg == "watt") {
            run_watt = true;
        } else if (arg == "all") {
            run_accuracy = true;
            run_flops    = true;
            run_watt     = true;
        } else if (arg == "SGEMM") {
            run_SGEMM = true;
        } else if (arg == "DGEMM") {
            run_DGEMM = true;
        } else if (arg == "CGEMM") {
            run_CGEMM = true;
        } else if (arg == "ZGEMM") {
            run_ZGEMM = true;
        }
    }

    if (run_SGEMM) {
        if (run_accuracy) accuracy_check<float>(deviceName, startTime);
        if (run_flops) time_check<float>(deviceName, startTime);
        if (run_watt) watt_check<float>(deviceName, startTime);
    }

    if (run_DGEMM) {
        if (run_accuracy) accuracy_check<double>(deviceName, startTime);
        if (run_flops) time_check<double>(deviceName, startTime);
        if (run_watt) watt_check<double>(deviceName, startTime);
    }

    if (run_CGEMM) {
        if (run_accuracy) accuracy_check<cuFloatComplex>(deviceName, startTime);
        if (run_flops) time_check<cuFloatComplex>(deviceName, startTime);
        if (run_watt) watt_check<cuFloatComplex>(deviceName, startTime);
    }

    if (run_ZGEMM) {
        if (run_accuracy) accuracy_check<cuDoubleComplex>(deviceName, startTime);
        if (run_flops) time_check<cuDoubleComplex>(deviceName, startTime);
        if (run_watt) watt_check<cuDoubleComplex>(deviceName, startTime);
    }

    std::string endTime = getCurrentDateTime(stop);
    auto sec            = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count() * 1.e-9;
    std::cout << std::endl;
    std::cout << "=========================" << std::endl;
    std::cout << "start        : " << startTime << std::endl;
    std::cout << "end          : " << endTime << std::endl;
    std::cout << "elapsed time : " << sec << " [sec]" << std::endl;
    std::cout << "=========================" << std::endl;
    std::cout << std::endl;

    CHECK_CUDA(cudaDeviceReset());
    return 0;
}
