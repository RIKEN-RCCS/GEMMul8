#pragma once
#include "../common/common.hpp"

namespace gemmul8::undo_scaling {

template <typename T>
struct HostScalar {
    T v;
    __device__ __forceinline__ T get() const { return v; }
    __host__ __device__ __forceinline__ constexpr HostScalar() : v(T(0)) {}
    __host__ __device__ __forceinline__ constexpr HostScalar(T v_) : v(v_) {}
};

template <typename T>
struct DeviceScalar {
    const T *p;
    __device__ __forceinline__ T get() const { return *p; }
    __host__ __device__ __forceinline__ constexpr DeviceScalar() : p(nullptr) {}
    __host__ __device__ __forceinline__ constexpr DeviceScalar(const T *p_) : p(p_) {}
};

template <typename T> struct scalar_type {
    using type = T;
};
template <> struct scalar_type<HostScalar<float>> {
    using type = float;
};
template <> struct scalar_type<HostScalar<double>> {
    using type = double;
};
template <> struct scalar_type<HostScalar<cuFloatComplex>> {
    using type = cuFloatComplex;
};
template <> struct scalar_type<HostScalar<cuDoubleComplex>> {
    using type = cuDoubleComplex;
};
template <> struct scalar_type<DeviceScalar<float>> {
    using type = float;
};
template <> struct scalar_type<DeviceScalar<double>> {
    using type = double;
};
template <> struct scalar_type<DeviceScalar<cuFloatComplex>> {
    using type = cuFloatComplex;
};
template <> struct scalar_type<DeviceScalar<cuDoubleComplex>> {
    using type = cuDoubleComplex;
};
template <typename T> using scalar_t = typename scalar_type<T>::type;

template <typename T, typename TAlpha, typename TBeta>
__device__ __forceinline__ T axpby(const TAlpha a, const T x, const TBeta b, const T *C) {
    if (is_zero(b)) return common::Tmul<TAlpha, T>(a, x);
    return common::Taxpby<T, TAlpha, TBeta>(a, x, b, *C);
}

template <typename T, typename TBeta>
__device__ __forceinline__ T scale_or_zero(const TBeta b, const T *C) {
    if (is_zero(b)) return common::Tconst<T>::zero();
    return common::Tmul<TBeta, T>(b, *C);
}

// C_new = alpha*D + beta*C_old
template <typename T, int ALPHA, int BETA>
__device__ __forceinline__ T Taxpby_special(const T D, const T *C_old) {
    if constexpr (ALPHA == 1) {
        if constexpr (BETA == 0) {
            return D;
        } else if constexpr (BETA == 1) {
            return common::Tadd<T>(*C_old, D);
        } else if constexpr (BETA == -1) {
            return common::Tsub<T>(D, *C_old);
        }
    } else if constexpr (ALPHA == -1) {
        if constexpr (BETA == 0) {
            return common::Tneg<T>(D);
        } else if constexpr (BETA == 1) {
            return common::Tsub<T>(*C_old, D);
        } else if constexpr (BETA == -1) {
            return common::Tsub<T>(common::Tneg<T>(D), *C_old);
        }
    }
    return common::Tconst<T>::zero();
}

template <typename T, int BETA>
__device__ __forceinline__ T Tmul_special(const T *C_old) {
    if constexpr (BETA == 1) {
        return *C_old;
    } else if constexpr (BETA == -1) {
        return common::Tneg<T>(*C_old);
    }
    return common::Tconst<T>::zero();
}

} // namespace gemmul8::undo_scaling
