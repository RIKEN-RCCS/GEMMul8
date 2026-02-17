#pragma once
#include "self_hipify.hpp"
#include <chrono>
#include <ctime>
#include <fstream>
#include <functional>
#include <iostream>
#include <sstream>
#include <thread>
#include <vector>

namespace getWatt {

struct PowerProfile {
    double power;
    std::time_t timestamp;
};

inline void gpu_power_init() {
    nvmlReturn_t result = nvmlInit();
    if (result != NVML_SUCCESS) {
        std::cerr << "Failed to initialize Management Library: " << nvmlErrorString(result) << std::endl;
    }
}

inline double get_current_power(unsigned gpu_id) {
    nvmlDevice_t device;
    nvmlReturn_t result = nvmlDeviceGetHandleByIndex_v2(gpu_id, &device);
    if (result != NVML_SUCCESS) {
        std::cerr << "Failed to get GPU handle:" << nvmlErrorString(result) << std::endl;
        return 0.0;
    }
    unsigned int mw;
    result = nvmlDeviceGetPowerUsage(device, &mw);
    if (result != NVML_SUCCESS) {
        std::cerr << "Failed to get power usage: " << nvmlErrorString(result) << std::endl;
        return 0.0;
    }
    return mw / 1000.0;
}

std::vector<PowerProfile> getGpuPowerUsage(const std::function<void(void)> func, const std::time_t interval) //
{
    gpu_power_init();

    std::vector<PowerProfile> prof_res;
    unsigned gpu_id = 0;
    unsigned count  = 0;
    bool running    = true;

    std::thread worker([&]() {
        func();
        running = false;
    });

    const auto start = std::chrono::high_resolution_clock::now();
    do {
        const auto now_1          = std::chrono::high_resolution_clock::now();
        const auto elapsed_time_1 = std::chrono::duration_cast<std::chrono::microseconds>(now_1 - start).count();
        const auto power          = get_current_power(gpu_id);
        const auto now_2          = std::chrono::high_resolution_clock::now();
        const auto elapsed_time_2 = std::chrono::duration_cast<std::chrono::milliseconds>(now_2 - start).count();
        prof_res.push_back(PowerProfile{power, elapsed_time_1});

        using namespace std::chrono_literals;
        std::this_thread::sleep_for(std::chrono::milliseconds(std::max<std::time_t>(static_cast<int>(interval) * count, elapsed_time_2) - elapsed_time_2));
        count++;
    } while (running);

    worker.join();
    nvmlShutdown();
    return prof_res;
}

double get_integrated_power_consumption(const std::vector<PowerProfile> &list) {
    if (list.empty()) return 0.0;

    double power_consumption = 0.;
    for (unsigned i = 1; i < list.size(); i++) {
        const auto elapsed_time = (list[i].timestamp - list[i - 1].timestamp) * 1e-6;
        power_consumption += (list[i].power + list[i - 1].power) / 2 * elapsed_time;
    }
    return power_consumption;
}

double get_elapsed_time(const std::vector<PowerProfile> &list) {
    if (list.empty()) return 0.0;
    return (list.back().timestamp - list.front().timestamp) * 1.e-6;
}

//=================================================================
// Function returns power consumption Watt
//================================================================
std::vector<double> getWatt(const std::function<void(void)> func, const size_t m, const size_t n, const size_t k) {
    constexpr size_t duration_time = 10;
    size_t cnt                     = 0;
    std::vector<PowerProfile> powerUsages;
    powerUsages = getGpuPowerUsage(
        [&]() {
            cudaDeviceSynchronize();
            const auto start = std::chrono::system_clock::now();
            while (true) {
                func();
                if (((++cnt) % 10) == 0) {
                    cudaDeviceSynchronize();
                    const auto now          = std::chrono::system_clock::now();
                    const auto elapsed_time = std::chrono::duration_cast<std::chrono::microseconds>(now - start).count() * 1e-6;
                    if (elapsed_time > duration_time) {
                        break;
                    }
                }
            }
        },
        100);
    const double power          = get_integrated_power_consumption(powerUsages);
    const double elapsed_time   = get_elapsed_time(powerUsages);
    const double watt           = power / elapsed_time;
    const double flops_per_watt = 2.0 * m * n * k * cnt / power;
    std::vector<double> results{watt, flops_per_watt};
    return results;
}

} // namespace getWatt
