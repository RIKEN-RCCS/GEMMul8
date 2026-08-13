#pragma once
#include "mod_core.hpp"
#include "config.hpp"

#include "../common/common.hpp"

namespace gemmul8::mod {

namespace {

template <unsigned IDX, bool nowrap = false>
__device__ __forceinline__ int32_t mod_f32x3_2_i32(const float C0, const float C1, const float C2) {
    //------------------------------
    // Transform 3 float integers
    //     isSqr<p> => sqrt(p_i)*(C0 + C1) + C2    with C0=Ahi*Blo, C1=Alo*Bhi, C2=Alo*Blo
    //    !isSqr<p> => 2401*C0 + 49*(C2-C0-C1) + C1 with C0=Ahi*Bhi, C1=Alo*Blo, C2=(Ahi+Alo)*(Bhi+Blo)
    // into 1 INT32
    //------------------------------
    const int32_t c0 = __float2int_rn(C0);
    const int32_t c1 = __float2int_rn(C1);
    const int32_t c2 = __float2int_rn(C2);
    const int32_t r0 = mod_small_nowrap<Backend::FP8, IDX>(c0);
    const int32_t r1 = mod_small_nowrap<Backend::FP8, IDX>(c1);
    const int32_t r2 = mod_small_nowrap<Backend::FP8, IDX>(c2);

    if constexpr (common::table::sqrt_moduli<IDX> != 0) {

        constexpr int32_t sqrt_p = common::table::sqrt_moduli<IDX>;
        const int32_t t          = sqrt_p * (r0 + r1) + r2;
        if constexpr (nowrap) return t;
        return mod_small<Backend::FP8, IDX>(t);

    } else {

        const int32_t t = (r0 * 2401) + ((r2 - r0 - r1) * 49) + r1;
        if constexpr (nowrap) return t;
        return mod_small<Backend::FP8, IDX>(t);
    }
}

template <Backend BACKEND, unsigned IDX>
__device__ __forceinline__ common::mid_t<BACKEND> mod_hi2mid_core(const int in) {
    return common::mid_t<BACKEND>(mod_small<BACKEND, IDX>(in));
}
template <>
__device__ __forceinline__ int8_t mod_hi2mid_core<Backend::INT8, 0U>(const int in) {
    return int8_t(in);
}

template <unsigned IDX>
__device__ __forceinline__ char4 mod_hi2mid_core_x4(const int4 in) {
    char4 out;
    out.x = mod_hi2mid_core<Backend::INT8, IDX>(in.x);
    out.y = mod_hi2mid_core<Backend::INT8, IDX>(in.y);
    out.z = mod_hi2mid_core<Backend::INT8, IDX>(in.z);
    out.w = mod_hi2mid_core<Backend::INT8, IDX>(in.w);
    return out;
}

template <unsigned IDX>
__device__ __forceinline__ short4 mod_hi2mid_core_x4(const float4 in0, const float4 in1, const float4 in2) {
    short4 out;
    out.x = int16_t(mod_f32x3_2_i32<IDX>(in0.x, in1.x, in2.x));
    out.y = int16_t(mod_f32x3_2_i32<IDX>(in0.y, in1.y, in2.y));
    out.z = int16_t(mod_f32x3_2_i32<IDX>(in0.z, in1.z, in2.z));
    out.w = int16_t(mod_f32x3_2_i32<IDX>(in0.w, in1.w, in2.w));
    return out;
}

template <unsigned IDX>
__device__ __forceinline__ int4 mod_hi2mid_core_nowrap_x4(const int4 in) {
    int4 out;
    out.x = mod_small_nowrap<Backend::INT8, IDX>(in.x);
    out.y = mod_small_nowrap<Backend::INT8, IDX>(in.y);
    out.z = mod_small_nowrap<Backend::INT8, IDX>(in.z);
    out.w = mod_small_nowrap<Backend::INT8, IDX>(in.w);
    return out;
}
template <>
__device__ __forceinline__ int4 mod_hi2mid_core_nowrap_x4<0U>(const int4 in) {
    return in;
}

template <unsigned IDX>
__device__ __forceinline__ int4 mod_hi2mid_core_nowrap_x4(const float4 in0, const float4 in1, const float4 in2) {
    int4 out;
    out.x = mod_f32x3_2_i32<IDX, true>(in0.x, in1.x, in2.x);
    out.y = mod_f32x3_2_i32<IDX, true>(in0.y, in1.y, in2.y);
    out.z = mod_f32x3_2_i32<IDX, true>(in0.z, in1.z, in2.z);
    out.w = mod_f32x3_2_i32<IDX, true>(in0.w, in1.w, in2.w);
    return out;
}

// complex Backend::INT8
template <unsigned IDX>
__device__ __forceinline__ short4 mod_hi2mid_device(
    const size_t idx,
    const int4 *__restrict__ C_hix4_1,
    const int4 *__restrict__ C_hix4_2,
    const int4 *__restrict__ C_hix4_3 //
) {
    const int4 ArBr   = mod_hi2mid_core_nowrap_x4<IDX>(C_hix4_1[idx]); // Ar x Br
    const int4 AiBi   = mod_hi2mid_core_nowrap_x4<IDX>(C_hix4_2[idx]); // Ai x Bi
    const int4 AriBri = mod_hi2mid_core_nowrap_x4<IDX>(C_hix4_3[idx]); // (Ar+Ai) x (Br+Bi)

    char2 C_mid_x;
    C_mid_x.x = mod_hi2mid_core<Backend::INT8, IDX>(ArBr.x - AiBi.x);            // real
    C_mid_x.y = mod_hi2mid_core<Backend::INT8, IDX>(AriBri.x - ArBr.x - AiBi.x); // imag

    char2 C_mid_y;
    C_mid_y.x = mod_hi2mid_core<Backend::INT8, IDX>(ArBr.y - AiBi.y);            // real
    C_mid_y.y = mod_hi2mid_core<Backend::INT8, IDX>(AriBri.y - ArBr.y - AiBi.y); // imag

    char2 C_mid_z;
    C_mid_z.x = mod_hi2mid_core<Backend::INT8, IDX>(ArBr.z - AiBi.z);            // real
    C_mid_z.y = mod_hi2mid_core<Backend::INT8, IDX>(AriBri.z - ArBr.z - AiBi.z); // imag

    char2 C_mid_w;
    C_mid_w.x = mod_hi2mid_core<Backend::INT8, IDX>(ArBr.w - AiBi.w);            // real
    C_mid_w.y = mod_hi2mid_core<Backend::INT8, IDX>(AriBri.w - ArBr.w - AiBi.w); // imag

    return short4{*reinterpret_cast<short *>(&C_mid_x),
                  *reinterpret_cast<short *>(&C_mid_y),
                  *reinterpret_cast<short *>(&C_mid_z),
                  *reinterpret_cast<short *>(&C_mid_w)};
}

// complex Backend::FP8
template <unsigned IDX>
__device__ __forceinline__ int4 mod_hi2mid_device(
    const size_t idx,
    const size_t sizeC4,
    const float4 *__restrict__ C_hix4_1,
    const float4 *__restrict__ C_hix4_2,
    const float4 *__restrict__ C_hix4_3 //
) {
    const int4 ArBr = mod_hi2mid_core_nowrap_x4<IDX>(
        C_hix4_1[idx], C_hix4_1[idx + sizeC4], C_hix4_1[idx + 2 * sizeC4]); // Ar x Br
    const int4 AiBi = mod_hi2mid_core_nowrap_x4<IDX>(
        C_hix4_2[idx], C_hix4_2[idx + sizeC4], C_hix4_2[idx + 2 * sizeC4]); // Ai x Bi
    const int4 AriBri = mod_hi2mid_core_nowrap_x4<IDX>(
        C_hix4_3[idx], C_hix4_3[idx + sizeC4], C_hix4_3[idx + 2 * sizeC4]); // (Ar+Ai) x (Br+Bi)

    short2 C_mid_x;
    C_mid_x.x = mod_hi2mid_core<Backend::FP8, IDX>(ArBr.x - AiBi.x);            // real
    C_mid_x.y = mod_hi2mid_core<Backend::FP8, IDX>(AriBri.x - ArBr.x - AiBi.x); // imag

    short2 C_mid_y;
    C_mid_y.x = mod_hi2mid_core<Backend::FP8, IDX>(ArBr.y - AiBi.y);            // real
    C_mid_y.y = mod_hi2mid_core<Backend::FP8, IDX>(AriBri.y - ArBr.y - AiBi.y); // imag

    short2 C_mid_z;
    C_mid_z.x = mod_hi2mid_core<Backend::FP8, IDX>(ArBr.z - AiBi.z);            // real
    C_mid_z.y = mod_hi2mid_core<Backend::FP8, IDX>(AriBri.z - ArBr.z - AiBi.z); // imag

    short2 C_mid_w;
    C_mid_w.x = mod_hi2mid_core<Backend::FP8, IDX>(ArBr.w - AiBi.w);            // real
    C_mid_w.y = mod_hi2mid_core<Backend::FP8, IDX>(AriBri.w - ArBr.w - AiBi.w); // imag

    return int4{*reinterpret_cast<int32_t *>(&C_mid_x),
                *reinterpret_cast<int32_t *>(&C_mid_y),
                *reinterpret_cast<int32_t *>(&C_mid_z),
                *reinterpret_cast<int32_t *>(&C_mid_w)};
}

// complex Backend::INT8
template <unsigned IDX, bool FLIP_IMAG = false>
__device__ __forceinline__ short4 mod_hi2mid_AHA_device(
    const size_t idx,
    const int4 *__restrict__ C_hix4_1,
    const int4 *__restrict__ C_hix4_2,
    const int4 *__restrict__ C_hix4_3 //
) {
    const int4 ArBi   = mod_hi2mid_core_nowrap_x4<IDX>(C_hix4_1[idx]); // Ar x Bi
    const int4 AiBr   = mod_hi2mid_core_nowrap_x4<IDX>(C_hix4_2[idx]); // Ai x Br
    const int4 AriBri = mod_hi2mid_core_nowrap_x4<IDX>(C_hix4_3[idx]); // (Ar+Ai) x (Br+Bi)

    char2 C_mid_x;
    C_mid_x.x = mod_hi2mid_core<Backend::INT8, IDX>(AriBri.x - ArBi.x - AiBr.x); // real
    if constexpr (FLIP_IMAG) {
        C_mid_x.y = mod_hi2mid_core<Backend::INT8, IDX>(AiBr.x - ArBi.x); // imag
    } else {
        C_mid_x.y = mod_hi2mid_core<Backend::INT8, IDX>(ArBi.x - AiBr.x); // imag
    }

    char2 C_mid_y;
    C_mid_y.x = mod_hi2mid_core<Backend::INT8, IDX>(AriBri.y - ArBi.y - AiBr.y); // real
    if constexpr (FLIP_IMAG) {
        C_mid_y.y = mod_hi2mid_core<Backend::INT8, IDX>(AiBr.y - ArBi.y); // imag
    } else {
        C_mid_y.y = mod_hi2mid_core<Backend::INT8, IDX>(ArBi.y - AiBr.y); // imag
    }

    char2 C_mid_z;
    C_mid_z.x = mod_hi2mid_core<Backend::INT8, IDX>(AriBri.z - ArBi.z - AiBr.z); // real
    if constexpr (FLIP_IMAG) {
        C_mid_z.y = mod_hi2mid_core<Backend::INT8, IDX>(AiBr.z - ArBi.z); // imag
    } else {
        C_mid_z.y = mod_hi2mid_core<Backend::INT8, IDX>(ArBi.z - AiBr.z); // imag
    }

    char2 C_mid_w;
    C_mid_w.x = mod_hi2mid_core<Backend::INT8, IDX>(AriBri.w - ArBi.w - AiBr.w); // real
    if constexpr (FLIP_IMAG) {
        C_mid_w.y = mod_hi2mid_core<Backend::INT8, IDX>(AiBr.w - ArBi.w); // imag
    } else {
        C_mid_w.y = mod_hi2mid_core<Backend::INT8, IDX>(ArBi.w - AiBr.w); // imag
    }

    return short4{*reinterpret_cast<short *>(&C_mid_x),
                  *reinterpret_cast<short *>(&C_mid_y),
                  *reinterpret_cast<short *>(&C_mid_z),
                  *reinterpret_cast<short *>(&C_mid_w)};
}

// complex Backend::FP8
template <unsigned IDX, bool FLIP_IMAG = false>
__device__ __forceinline__ int4 mod_hi2mid_AHA_device(
    const size_t idx,
    const size_t sizeC4,
    const float4 *__restrict__ C_hix4_1,
    const float4 *__restrict__ C_hix4_2,
    const float4 *__restrict__ C_hix4_3 //
) {
    const int4 ArBi = mod_hi2mid_core_nowrap_x4<IDX>(
        C_hix4_1[idx], C_hix4_1[idx + sizeC4], C_hix4_1[idx + 2 * sizeC4]); // Ar x Bi
    const int4 AiBr = mod_hi2mid_core_nowrap_x4<IDX>(
        C_hix4_2[idx], C_hix4_2[idx + sizeC4], C_hix4_2[idx + 2 * sizeC4]); // Ai x Br
    const int4 AriBri = mod_hi2mid_core_nowrap_x4<IDX>(
        C_hix4_3[idx], C_hix4_3[idx + sizeC4], C_hix4_3[idx + 2 * sizeC4]); // (Ar+Ai) x (Br+Bi)

    short2 C_mid_x;
    C_mid_x.x = mod_hi2mid_core<Backend::FP8, IDX>(AriBri.x - ArBi.x - AiBr.x); // real
    if constexpr (FLIP_IMAG) {
        C_mid_x.y = mod_hi2mid_core<Backend::FP8, IDX>(AiBr.x - ArBi.x); // imag
    } else {
        C_mid_x.y = mod_hi2mid_core<Backend::FP8, IDX>(ArBi.x - AiBr.x); // imag
    }

    short2 C_mid_y;
    C_mid_y.x = mod_hi2mid_core<Backend::FP8, IDX>(AriBri.y - ArBi.y - AiBr.y); // real
    if constexpr (FLIP_IMAG) {
        C_mid_y.y = mod_hi2mid_core<Backend::FP8, IDX>(AiBr.y - ArBi.y); // imag
    } else {
        C_mid_y.y = mod_hi2mid_core<Backend::FP8, IDX>(ArBi.y - AiBr.y); // imag
    }

    short2 C_mid_z;
    C_mid_z.x = mod_hi2mid_core<Backend::FP8, IDX>(AriBri.z - ArBi.z - AiBr.z); // real
    if constexpr (FLIP_IMAG) {
        C_mid_z.y = mod_hi2mid_core<Backend::FP8, IDX>(AiBr.z - ArBi.z); // imag
    } else {
        C_mid_z.y = mod_hi2mid_core<Backend::FP8, IDX>(ArBi.z - AiBr.z); // imag
    }

    short2 C_mid_w;
    C_mid_w.x = mod_hi2mid_core<Backend::FP8, IDX>(AriBri.w - ArBi.w - AiBr.w); // real
    if constexpr (FLIP_IMAG) {
        C_mid_w.y = mod_hi2mid_core<Backend::FP8, IDX>(AiBr.w - ArBi.w); // imag
    } else {
        C_mid_w.y = mod_hi2mid_core<Backend::FP8, IDX>(ArBi.w - AiBr.w); // imag
    }

    return int4{*reinterpret_cast<int32_t *>(&C_mid_x),
                *reinterpret_cast<int32_t *>(&C_mid_y),
                *reinterpret_cast<int32_t *>(&C_mid_z),
                *reinterpret_cast<int32_t *>(&C_mid_w)};
}

// real
template <Backend BACKEND, unsigned IDX>
__global__ void mod_hi2mid_ge_kernel(
    const size_t sizeC4,
    const common::hix4_t<BACKEND> *__restrict__ C_hix4,
    common::midx4_t<BACKEND, false> *__restrict__ C_midx4 //
) {
    const size_t idx = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (idx >= sizeC4) return;

    if constexpr (BACKEND == Backend::INT8) {
        C_midx4[idx] = mod_hi2mid_core_x4<IDX>(C_hix4[idx]);
    } else {
        C_midx4[idx] = mod_hi2mid_core_x4<IDX>(C_hix4[idx], C_hix4[idx + sizeC4], C_hix4[idx + 2U * sizeC4]);
    }
}

// real
template <Backend BACKEND, unsigned IDX, cublasFillMode_t UPLO>
__global__ void mod_hi2mid_tri_kernel(
    const unsigned n, const size_t sizeC4,
    const common::hix4_t<BACKEND> *__restrict__ C_hix4,
    common::midx4_t<BACKEND, false> *__restrict__ C_midx4,
    const size_t ldc4 //
) {
    const unsigned row4 = blockIdx.x * blockDim.x + threadIdx.x;
    const unsigned col  = blockIdx.y * blockDim.y + threadIdx.y;
    if (row4 >= ldc4 || col >= n) return;

    const unsigned row = row4 << 2;

    if constexpr (UPLO == CUBLAS_FILL_MODE_UPPER) {
        if (row > col) return;
    } else if constexpr (UPLO == CUBLAS_FILL_MODE_LOWER) {
        if (row + 3 < col) return;
    }

    const size_t idx = col * ldc4 + row4;

    if constexpr (BACKEND == Backend::INT8) {
        C_midx4[idx] = mod_hi2mid_core_x4<IDX>(C_hix4[idx]);
    } else {
        C_midx4[idx] = mod_hi2mid_core_x4<IDX>(C_hix4[idx], C_hix4[idx + sizeC4], C_hix4[idx + 2U * sizeC4]);
    }
}

// complex
template <Backend BACKEND, unsigned IDX>
__global__ void mod_hi2mid_ge_kernel(
    const size_t sizeC4,
    const common::hix4_t<BACKEND> *__restrict__ C_hix4_1,
    const common::hix4_t<BACKEND> *__restrict__ C_hix4_2,
    const common::hix4_t<BACKEND> *__restrict__ C_hix4_3,
    common::midx4_t<BACKEND, true> *__restrict__ C_midx4 //
) {
    const size_t idx = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (idx >= sizeC4) return;

    if constexpr (BACKEND == Backend::INT8) {
        C_midx4[idx] = mod_hi2mid_device<IDX>(idx, C_hix4_1, C_hix4_2, C_hix4_3);
    } else {
        C_midx4[idx] = mod_hi2mid_device<IDX>(idx, sizeC4, C_hix4_1, C_hix4_2, C_hix4_3);
    }
}

// complex
template <Backend BACKEND, unsigned IDX, cublasFillMode_t UPLO>
__global__ void mod_hi2mid_tri_kernel(
    const unsigned n, const size_t sizeC4,
    const common::hix4_t<BACKEND> *__restrict__ C_hix4_1,
    const common::hix4_t<BACKEND> *__restrict__ C_hix4_2,
    const common::hix4_t<BACKEND> *__restrict__ C_hix4_3,
    common::midx4_t<BACKEND, true> *__restrict__ C_midx4,
    const size_t ldc4 //
) {
    const unsigned row4 = blockIdx.x * blockDim.x + threadIdx.x;
    const unsigned col  = blockIdx.y * blockDim.y + threadIdx.y;
    if (row4 >= ldc4 || col >= n) return;

    const unsigned row = row4 << 2;

    if constexpr (UPLO == CUBLAS_FILL_MODE_UPPER) {
        if (row > col) return;
    } else if constexpr (UPLO == CUBLAS_FILL_MODE_LOWER) {
        if (row + 3 < col) return;
    }

    const size_t idx = col * ldc4 + row4;

    if constexpr (BACKEND == Backend::INT8) {
        C_midx4[idx] = mod_hi2mid_device<IDX>(idx, C_hix4_1, C_hix4_2, C_hix4_3);
    } else {
        C_midx4[idx] = mod_hi2mid_device<IDX>(idx, sizeC4, C_hix4_1, C_hix4_2, C_hix4_3);
    }
}

// complex
template <Backend BACKEND, unsigned IDX, cublasFillMode_t UPLO, bool FLIP_IMAG = false>
__global__ void mod_hi2mid_AHA_kernel(
    const unsigned n, const size_t sizeC4,
    const common::hix4_t<BACKEND> *__restrict__ C_hix4_1,
    const common::hix4_t<BACKEND> *__restrict__ C_hix4_2,
    const common::hix4_t<BACKEND> *__restrict__ C_hix4_3,
    common::midx4_t<BACKEND, true> *__restrict__ C_midx4,
    const size_t ldc4 //
) {
    const unsigned row4 = blockIdx.x * blockDim.x + threadIdx.x;
    const unsigned col  = blockIdx.y * blockDim.y + threadIdx.y;
    if (row4 >= ldc4 || col >= n) return;

    const unsigned row = row4 << 2;

    if constexpr (UPLO == CUBLAS_FILL_MODE_UPPER) {
        if (row > col) return;
    } else if constexpr (UPLO == CUBLAS_FILL_MODE_LOWER) {
        if (row + 3 < col) return;
    }

    const size_t idx = col * ldc4 + row4;

    if constexpr (BACKEND == Backend::INT8) {
        C_midx4[idx] = mod_hi2mid_AHA_device<IDX, FLIP_IMAG>(idx, C_hix4_1, C_hix4_2, C_hix4_3);
    } else {
        C_midx4[idx] = mod_hi2mid_AHA_device<IDX, FLIP_IMAG>(idx, sizeC4, C_hix4_1, C_hix4_2, C_hix4_3);
    }
}

template <Backend BACKEND, bool COMPLEX, unsigned IDX, cublasFillMode_t UPLO>
inline void mod_hi2mid_launch(
    const cudaStream_t stream,
    const size_t ldc, const unsigned n,
    common::matptr_t<common::hi_t<BACKEND>, COMPLEX> &C_hi,
    common::mid_t<BACKEND, COMPLEX> *C_mid //
) {
    using HI4  = common::hix4_t<BACKEND>;
    using MID4 = common::midx4_t<BACKEND, COMPLEX>;

    const size_t sizeC  = ldc * n;
    const size_t sizeC4 = sizeC >> 2;

    MID4 *C_midx4 = reinterpret_cast<MID4 *>(C_mid + IDX * sizeC);

    if constexpr (COMPLEX) {

        const HI4 *C_hix4_1 = reinterpret_cast<const HI4 *>(C_hi.ptr0);
        const HI4 *C_hix4_2 = reinterpret_cast<const HI4 *>(C_hi.ptr1);
        const HI4 *C_hix4_3 = reinterpret_cast<const HI4 *>(C_hi.ptr2);

        if constexpr (UPLO == CUBLAS_FILL_MODE_FULL) {

            const dim3 grid((sizeC4 + threads_1d - 1) / threads_1d);

            mod_hi2mid_ge_kernel<BACKEND, IDX>
                <<<grid, threads_1d, 0, stream>>>(
                    sizeC4, C_hix4_1, C_hix4_2, C_hix4_3, C_midx4);

        } else {

            const size_t ldc4        = ldc >> 2;
            const unsigned threads_y = select_threads_y<BACKEND>(sizeC);
            const dim3 threads(threads_x, threads_y);
            const dim3 grid((ldc4 + threads_x - 1) / threads_x,
                            (n + threads_y - 1) / threads_y);

            mod_hi2mid_tri_kernel<BACKEND, IDX, UPLO>
                <<<grid, threads, 0, stream>>>(
                    n, sizeC4, C_hix4_1, C_hix4_2, C_hix4_3, C_midx4, ldc4);
        }

    } else {

        const HI4 *C_hix4 = reinterpret_cast<const HI4 *>(C_hi.ptr0);

        if constexpr (UPLO == CUBLAS_FILL_MODE_FULL) {

            const dim3 grid((sizeC4 + threads_1d - 1) / threads_1d);

            mod_hi2mid_ge_kernel<BACKEND, IDX>
                <<<grid, threads_1d, 0, stream>>>(
                    sizeC4, C_hix4, C_midx4);

        } else {

            const size_t ldc4        = ldc >> 2;
            const unsigned threads_y = select_threads_y<BACKEND>(sizeC);
            const dim3 threads(threads_x, threads_y);
            const dim3 grid((ldc4 + threads_x - 1) / threads_x,
                            (n + threads_y - 1) / threads_y);

            mod_hi2mid_tri_kernel<BACKEND, IDX, UPLO>
                <<<grid, threads, 0, stream>>>(
                    n, sizeC4, C_hix4, C_midx4, ldc4);
        }
    }
}

template <Backend BACKEND, unsigned IDX, cublasFillMode_t UPLO, bool FLIP_IMAG = false>
inline void mod_hi2mid_AHA_launch(
    const cudaStream_t stream,
    const size_t ldc, const unsigned n,
    common::matptr_t<common::hi_t<BACKEND>, true> &C_hi,
    common::mid_t<BACKEND, true> *C_mid //
) {
    using HI4  = common::hix4_t<BACKEND>;
    using MID4 = common::midx4_t<BACKEND, true>;

    const size_t sizeC  = ldc * n;
    const size_t sizeC4 = sizeC >> 2;

    MID4 *C_midx4       = reinterpret_cast<MID4 *>(C_mid + IDX * sizeC);
    const HI4 *C_hix4_1 = reinterpret_cast<const HI4 *>(C_hi.ptr0);
    const HI4 *C_hix4_2 = reinterpret_cast<const HI4 *>(C_hi.ptr1);
    const HI4 *C_hix4_3 = reinterpret_cast<const HI4 *>(C_hi.ptr2);

    const size_t ldc4        = ldc >> 2;
    const unsigned threads_y = select_threads_y<BACKEND>(sizeC);
    const dim3 threads(threads_x, threads_y);
    const dim3 grid((ldc4 + threads_x - 1) / threads_x,
                    (n + threads_y - 1) / threads_y);

    mod_hi2mid_AHA_kernel<BACKEND, IDX, UPLO, FLIP_IMAG>
        <<<grid, threads, 0, stream>>>(
            n, sizeC4, C_hix4_1, C_hix4_2, C_hix4_3, C_midx4, ldc4);
}

} // namespace

template <Backend BACKEND, bool COMPLEX, cublasFillMode_t UPLO>
void mod_hi2mid(
    const cudaStream_t stream,
    const unsigned idx,
    const size_t ldc, const unsigned n,
    common::matptr_t<common::hi_t<BACKEND>, COMPLEX> &C_hi,
    common::mid_t<BACKEND, COMPLEX> *C_mid //
) {
    switch (idx) {
    case 0U: mod_hi2mid_launch<BACKEND, COMPLEX, 0U, UPLO>(stream, ldc, n, C_hi, C_mid); break;
    case 1U: mod_hi2mid_launch<BACKEND, COMPLEX, 1U, UPLO>(stream, ldc, n, C_hi, C_mid); break;
    case 2U: mod_hi2mid_launch<BACKEND, COMPLEX, 2U, UPLO>(stream, ldc, n, C_hi, C_mid); break;
    case 3U: mod_hi2mid_launch<BACKEND, COMPLEX, 3U, UPLO>(stream, ldc, n, C_hi, C_mid); break;
    case 4U: mod_hi2mid_launch<BACKEND, COMPLEX, 4U, UPLO>(stream, ldc, n, C_hi, C_mid); break;
    case 5U: mod_hi2mid_launch<BACKEND, COMPLEX, 5U, UPLO>(stream, ldc, n, C_hi, C_mid); break;
    case 6U: mod_hi2mid_launch<BACKEND, COMPLEX, 6U, UPLO>(stream, ldc, n, C_hi, C_mid); break;
    case 7U: mod_hi2mid_launch<BACKEND, COMPLEX, 7U, UPLO>(stream, ldc, n, C_hi, C_mid); break;
    case 8U: mod_hi2mid_launch<BACKEND, COMPLEX, 8U, UPLO>(stream, ldc, n, C_hi, C_mid); break;
    case 9U: mod_hi2mid_launch<BACKEND, COMPLEX, 9U, UPLO>(stream, ldc, n, C_hi, C_mid); break;
    case 10U: mod_hi2mid_launch<BACKEND, COMPLEX, 10U, UPLO>(stream, ldc, n, C_hi, C_mid); break;
    case 11U: mod_hi2mid_launch<BACKEND, COMPLEX, 11U, UPLO>(stream, ldc, n, C_hi, C_mid); break;
    case 12U: mod_hi2mid_launch<BACKEND, COMPLEX, 12U, UPLO>(stream, ldc, n, C_hi, C_mid); break;
    case 13U: mod_hi2mid_launch<BACKEND, COMPLEX, 13U, UPLO>(stream, ldc, n, C_hi, C_mid); break;
    case 14U: mod_hi2mid_launch<BACKEND, COMPLEX, 14U, UPLO>(stream, ldc, n, C_hi, C_mid); break;
    case 15U: mod_hi2mid_launch<BACKEND, COMPLEX, 15U, UPLO>(stream, ldc, n, C_hi, C_mid); break;
    case 16U: mod_hi2mid_launch<BACKEND, COMPLEX, 16U, UPLO>(stream, ldc, n, C_hi, C_mid); break;
    case 17U: mod_hi2mid_launch<BACKEND, COMPLEX, 17U, UPLO>(stream, ldc, n, C_hi, C_mid); break;
    case 18U: mod_hi2mid_launch<BACKEND, COMPLEX, 18U, UPLO>(stream, ldc, n, C_hi, C_mid); break;
    case 19U: mod_hi2mid_launch<BACKEND, COMPLEX, 19U, UPLO>(stream, ldc, n, C_hi, C_mid); break;
    default: break;
    }
}

template <Backend BACKEND, cublasFillMode_t UPLO, bool FLIP_IMAG>
void mod_hi2mid_AHA(
    const cudaStream_t stream,
    const unsigned idx,
    const size_t ldc, const unsigned n,
    common::matptr_t<common::hi_t<BACKEND>, true> &C_hi,
    common::mid_t<BACKEND, true> *C_mid //
) {
    switch (idx) {
    case 0U: mod_hi2mid_AHA_launch<BACKEND, 0U, UPLO, FLIP_IMAG>(stream, ldc, n, C_hi, C_mid); break;
    case 1U: mod_hi2mid_AHA_launch<BACKEND, 1U, UPLO, FLIP_IMAG>(stream, ldc, n, C_hi, C_mid); break;
    case 2U: mod_hi2mid_AHA_launch<BACKEND, 2U, UPLO, FLIP_IMAG>(stream, ldc, n, C_hi, C_mid); break;
    case 3U: mod_hi2mid_AHA_launch<BACKEND, 3U, UPLO, FLIP_IMAG>(stream, ldc, n, C_hi, C_mid); break;
    case 4U: mod_hi2mid_AHA_launch<BACKEND, 4U, UPLO, FLIP_IMAG>(stream, ldc, n, C_hi, C_mid); break;
    case 5U: mod_hi2mid_AHA_launch<BACKEND, 5U, UPLO, FLIP_IMAG>(stream, ldc, n, C_hi, C_mid); break;
    case 6U: mod_hi2mid_AHA_launch<BACKEND, 6U, UPLO, FLIP_IMAG>(stream, ldc, n, C_hi, C_mid); break;
    case 7U: mod_hi2mid_AHA_launch<BACKEND, 7U, UPLO, FLIP_IMAG>(stream, ldc, n, C_hi, C_mid); break;
    case 8U: mod_hi2mid_AHA_launch<BACKEND, 8U, UPLO, FLIP_IMAG>(stream, ldc, n, C_hi, C_mid); break;
    case 9U: mod_hi2mid_AHA_launch<BACKEND, 9U, UPLO, FLIP_IMAG>(stream, ldc, n, C_hi, C_mid); break;
    case 10U: mod_hi2mid_AHA_launch<BACKEND, 10U, UPLO, FLIP_IMAG>(stream, ldc, n, C_hi, C_mid); break;
    case 11U: mod_hi2mid_AHA_launch<BACKEND, 11U, UPLO, FLIP_IMAG>(stream, ldc, n, C_hi, C_mid); break;
    case 12U: mod_hi2mid_AHA_launch<BACKEND, 12U, UPLO, FLIP_IMAG>(stream, ldc, n, C_hi, C_mid); break;
    case 13U: mod_hi2mid_AHA_launch<BACKEND, 13U, UPLO, FLIP_IMAG>(stream, ldc, n, C_hi, C_mid); break;
    case 14U: mod_hi2mid_AHA_launch<BACKEND, 14U, UPLO, FLIP_IMAG>(stream, ldc, n, C_hi, C_mid); break;
    case 15U: mod_hi2mid_AHA_launch<BACKEND, 15U, UPLO, FLIP_IMAG>(stream, ldc, n, C_hi, C_mid); break;
    case 16U: mod_hi2mid_AHA_launch<BACKEND, 16U, UPLO, FLIP_IMAG>(stream, ldc, n, C_hi, C_mid); break;
    case 17U: mod_hi2mid_AHA_launch<BACKEND, 17U, UPLO, FLIP_IMAG>(stream, ldc, n, C_hi, C_mid); break;
    case 18U: mod_hi2mid_AHA_launch<BACKEND, 18U, UPLO, FLIP_IMAG>(stream, ldc, n, C_hi, C_mid); break;
    case 19U: mod_hi2mid_AHA_launch<BACKEND, 19U, UPLO, FLIP_IMAG>(stream, ldc, n, C_hi, C_mid); break;
    default: break;
    }
}

} // namespace gemmul8::mod
