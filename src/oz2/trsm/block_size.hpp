#pragma once
#include "../common/common.hpp"

namespace gemmul8::oz2::trsm {

template <typename T, Backend BACKEND>
inline int block_size_setting_trsm(size_t n, int arch, int nB_override) noexcept {
    if (nB_override > 0) return nB_override;

    constexpr bool isDouble = std::is_same_v<T, double> || std::is_same_v<T, cuDoubleComplex>;
    constexpr bool isFloat  = !isDouble;
    constexpr bool isINT8   = BACKEND == Backend::INT8;
    constexpr bool isFP8    = !isINT8;

    switch (arch) {
    case 89: {
        if constexpr (isFloat && isINT8) return n;
        if constexpr (isFloat && isFP8) return n;
        if constexpr (isDouble && isINT8) return 1024;
        if constexpr (isDouble && isFP8) return (n < 8192) ? 1024 : 2048;
        return 4096;
    }
    case 90: {
        if constexpr (isFloat && isINT8) return (n <= 4096) ? 2048 : 3072;
        if constexpr (isFloat && isFP8) return std::min<int>(8192, n);
        if constexpr (isDouble && isINT8) return (n <= 8192) ? n : int(common::padding(n / 2));
        if constexpr (isDouble && isFP8) return n;
        return 4096;
    }
    case 100: {
        if constexpr (isFloat && isINT8) return (n <= 4096) ? 2048 : 3072;
        if constexpr (isFloat && isFP8) return 4096;
        if constexpr (isDouble && isINT8) return (n <= 4096) ? 2048 : 3072;
        if constexpr (isDouble && isFP8) return (n <= 8192) ? 4096 : int(common::padding(n / 2));
        return 3072;
    }
    case 103: {
        if constexpr (isFloat && isINT8) return n;
        if constexpr (isFloat && isFP8) return 4096;
        if constexpr (isDouble && isINT8) return n;
        if constexpr (isDouble && isFP8) return (n <= 8192) ? 4096 : int(common::padding(n / 2));
        return 3072;
    }
    case 120: {
        if constexpr (isFloat && isINT8) return n;
        if constexpr (isFloat && isFP8) return n;
        if constexpr (isDouble && isINT8) return 1024;
        if constexpr (isDouble && isFP8) return 1024;
        return 1024;
    }
    case 121: {
        if constexpr (isFloat && isINT8) return n;
        if constexpr (isFloat && isFP8) return n;
        if constexpr (isDouble && isINT8) return 1024;
        if constexpr (isDouble && isFP8) return (n <= 4096) ? 1024 : 2048;
        return 1024;
    }
    default:
        return 4096;
    }
}

template <typename T, Backend BACKEND>
inline int block_size_trsm(size_t n, int &arch, int nB_override = 0) noexcept {
    if (arch == 0) {
#if defined(__CUDACC__)
    #if defined(GPU_ARCH)
        arch = GPU_ARCH;
    #else
        int dev = 0;
        if (cudaGetDevice(&dev) == cudaSuccess) {
            int major = 0, minor = 0;
            if (cudaDeviceGetAttribute(&major, cudaDevAttrComputeCapabilityMajor, dev) == cudaSuccess &&
                cudaDeviceGetAttribute(&minor, cudaDevAttrComputeCapabilityMinor, dev) == cudaSuccess) {
                arch = major * 10 + minor;
            }
        }
    #endif
#endif
    }
    return block_size_setting_trsm<T, BACKEND>(n, arch, nB_override);
}

} // namespace gemmul8::oz2::trsm
