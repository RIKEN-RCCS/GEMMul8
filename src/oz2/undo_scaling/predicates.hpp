#pragma once
#include "../common/common.hpp"

namespace gemmul8::undo_scaling {

template <typename T>
__host__ __device__ __forceinline__ bool is_zero(const T v) {
    if constexpr (common::isComplex<T>) {
        using U = common::underlying_t<T>;
        return v.x == U(0) && v.y == U(0);
    } else {
        return v == T(0);
    }
}

template <typename T>
__host__ __forceinline__ bool is_zero_h(const T v) {
    return is_zero(v);
}

template <typename T>
__host__ __forceinline__ bool is_one_h(const T v) {
    if constexpr (common::isComplex<T>) {
        using U = common::underlying_t<T>;
        return v.x == U(1) && v.y == U(0);
    } else {
        return v == T(1);
    }
}

template <typename T>
__host__ __forceinline__ bool is_mone_h(const T v) {
    if constexpr (common::isComplex<T>) {
        using U = common::underlying_t<T>;
        return v.x == U(-1) && v.y == U(0);
    } else {
        return v == T(-1);
    }
}

inline bool is_device_pointer(const void *ptr) {
    cudaPointerAttributes attr{};
    const cudaError_t err = cudaPointerGetAttributes(&attr, ptr);

    if (err != cudaSuccess) {
        cudaGetLastError();
        return false;
    }

    return attr.type != cudaMemoryTypeUnregistered &&
           attr.type != cudaMemoryTypeHost;
}

} // namespace gemmul8::undo_scaling
