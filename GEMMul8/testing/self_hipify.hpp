#pragma once

#include "../src/self_hipify.hpp"

#if defined(__HIPCC__)
    #include <amd_smi/amdsmi.h>
    #include <hiprand/hiprand_kernel.h>

    #define curandState           hiprandState
    #define curand_init           hiprand_init
    #define curand_uniform_double hiprand_uniform_double
    #define curand_normal_double  hiprand_normal_double

    #define nvmlShutdown amdsmi_shut_down
    #define NVML_SUCCESS AMDSMI_STATUS_SUCCESS
    #define nvmlReturn_t amdsmi_status_t
    #define nvmlInit()   amdsmi_init(AMDSMI_INIT_AMD_GPUS)
    #define nvmlDevice_t amdsmi_processor_handle

namespace getWatt {
const char *nvmlErrorString(amdsmi_status_t result) {
    const char *error_string = nullptr;
    amdsmi_status_code_to_string(result, &error_string);
    return error_string;
}

amdsmi_status_t nvmlDeviceGetHandleByIndex_v2(unsigned gpu_id, amdsmi_processor_handle *device) {
    unsigned socket_count  = 0;
    amdsmi_status_t result = amdsmi_get_socket_handles(&socket_count, nullptr);
    if (result != AMDSMI_STATUS_SUCCESS) return result;

    std::vector<amdsmi_socket_handle> sockets(socket_count);
    result = amdsmi_get_socket_handles(&socket_count, &sockets[0]);
    if (result != AMDSMI_STATUS_SUCCESS) return result;

    unsigned device_count = 0;
    result                = amdsmi_get_processor_handles(sockets[0], &device_count, nullptr);

    std::vector<amdsmi_processor_handle> proc_handles(device_count);
    result = amdsmi_get_processor_handles(sockets[0], &device_count, &proc_handles[0]);
    if (result != AMDSMI_STATUS_SUCCESS) return result;

    *device = proc_handles[gpu_id];
    return AMDSMI_STATUS_SUCCESS;
}

amdsmi_status_t nvmlDeviceGetPowerUsage(amdsmi_processor_handle device, unsigned *mw) {
    amdsmi_power_info_t info{};
    amdsmi_status_t result = amdsmi_get_power_info(device, &info);
    if (result != AMDSMI_STATUS_SUCCESS) return result;
    *mw = (info.average_socket_power >= 10000)
              ? static_cast<unsigned>(info.current_socket_power * 1000.0)
              : static_cast<unsigned>(info.average_socket_power * 1000.0);
    return AMDSMI_STATUS_SUCCESS;
}
} // namespace getWatt

#elif defined(__NVCC__)
    #include <cuda_runtime.h>
    #include <cublas_v2.h>
    #include <cublasLt.h>
    #include <cuComplex.h>
    #include <curand_kernel.h>
    #include <nvml.h>
#endif
