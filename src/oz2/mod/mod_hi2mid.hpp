#pragma once
#include "reconstruct_fp8.hpp"
#include "complex_2m.hpp"
#include "config.hpp"
#include "../common/common.hpp"

namespace gemmul8::mod {

namespace {

template <Backend BACKEND, unsigned IDX>
__device__ __forceinline__ common::mid_t<BACKEND> mod_hi2mid_core(const int in) {
    if constexpr (common::table::moduli<BACKEND, IDX, false> == 256) {
        return int8_t(in);
    } else {
        return common::mid_t<BACKEND>(mod_small<BACKEND, IDX>(in));
    }
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

template <Backend BACKEND, unsigned IDX, bool FLIP_IMAG = false>
__device__ __forceinline__ auto pack_complex_2m(int32_t plus, int32_t minus) {
    const int2 v = reconstruct_complex_2m<BACKEND, IDX, FLIP_IMAG>(plus, minus);
    if constexpr (BACKEND == Backend::INT8) {
        const uint16_t bits = uint16_t(uint8_t(v.x)) | (uint16_t(uint8_t(v.y)) << 8);
        return static_cast<int16_t>(bits);
    } else {
        const uint32_t bits = uint32_t(uint16_t(v.x)) | (uint32_t(uint16_t(v.y)) << 16);
        return static_cast<int32_t>(bits);
    }
}

template <unsigned IDX, bool FLIP_IMAG = false>
__device__ __forceinline__ short4 mod_hi2mid_device(
    const size_t idx,
    const int4 *__restrict__ C_plus,
    const int4 *__restrict__ C_minus //
) {
    const int4 a = C_plus[idx];
    const int4 b = C_minus[idx];
    short4 out;
    out.x = pack_complex_2m<Backend::INT8, IDX, FLIP_IMAG>(
        mod_small<Backend::INT8, IDX, true>(a.x), mod_small<Backend::INT8, IDX, true>(b.x));
    out.y = pack_complex_2m<Backend::INT8, IDX, FLIP_IMAG>(
        mod_small<Backend::INT8, IDX, true>(a.y), mod_small<Backend::INT8, IDX, true>(b.y));
    out.z = pack_complex_2m<Backend::INT8, IDX, FLIP_IMAG>(
        mod_small<Backend::INT8, IDX, true>(a.z), mod_small<Backend::INT8, IDX, true>(b.z));
    out.w = pack_complex_2m<Backend::INT8, IDX, FLIP_IMAG>(
        mod_small<Backend::INT8, IDX, true>(a.w), mod_small<Backend::INT8, IDX, true>(b.w));
    return out;
}

template <unsigned IDX, bool FLIP_IMAG = false>
__device__ __forceinline__ int4 mod_hi2mid_device(
    const size_t idx,
    const size_t sizeC4,
    const float4 *__restrict__ C_plus,
    const float4 *__restrict__ C_minus //
) {
    const float4 a0 = C_plus[idx], a1 = C_plus[idx + sizeC4], a2 = C_plus[idx + 2 * sizeC4];
    const float4 b0 = C_minus[idx], b1 = C_minus[idx + sizeC4], b2 = C_minus[idx + 2 * sizeC4];
    int4 out;
    out.x = pack_complex_2m<Backend::FP8, IDX, FLIP_IMAG>(
        mod_f32x3_2_i32<IDX, false, true>(a0.x, a1.x, a2.x),
        mod_f32x3_2_i32<IDX, false, true>(b0.x, b1.x, b2.x));
    out.y = pack_complex_2m<Backend::FP8, IDX, FLIP_IMAG>(
        mod_f32x3_2_i32<IDX, false, true>(a0.y, a1.y, a2.y),
        mod_f32x3_2_i32<IDX, false, true>(b0.y, b1.y, b2.y));
    out.z = pack_complex_2m<Backend::FP8, IDX, FLIP_IMAG>(
        mod_f32x3_2_i32<IDX, false, true>(a0.z, a1.z, a2.z),
        mod_f32x3_2_i32<IDX, false, true>(b0.z, b1.z, b2.z));
    out.w = pack_complex_2m<Backend::FP8, IDX, FLIP_IMAG>(
        mod_f32x3_2_i32<IDX, false, true>(a0.w, a1.w, a2.w),
        mod_f32x3_2_i32<IDX, false, true>(b0.w, b1.w, b2.w));
    return out;
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
    common::midx4_t<BACKEND, true> *__restrict__ C_midx4 //
) {
    const size_t idx = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (idx >= sizeC4) return;

    if constexpr (BACKEND == Backend::INT8) {
        C_midx4[idx] = mod_hi2mid_device<IDX>(idx, C_hix4_1, C_hix4_2);
    } else {
        C_midx4[idx] = mod_hi2mid_device<IDX>(idx, sizeC4, C_hix4_1, C_hix4_2);
    }
}

// complex
template <Backend BACKEND, unsigned IDX, cublasFillMode_t UPLO>
__global__ void mod_hi2mid_tri_kernel(
    const unsigned n, const size_t sizeC4,
    const common::hix4_t<BACKEND> *__restrict__ C_hix4_1,
    const common::hix4_t<BACKEND> *__restrict__ C_hix4_2,
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
        C_midx4[idx] = mod_hi2mid_device<IDX>(idx, C_hix4_1, C_hix4_2);
    } else {
        C_midx4[idx] = mod_hi2mid_device<IDX>(idx, sizeC4, C_hix4_1, C_hix4_2);
    }
}

template <Backend BACKEND, unsigned IDX>
__device__ __forceinline__ int32_t herk_product_residue(
    const common::hi_t<BACKEND> *__restrict__ C_plus,
    const size_t idx,
    const size_t sizeC //
) {
    if constexpr (BACKEND == Backend::INT8) {
        return mod_small<BACKEND, IDX, true>(C_plus[idx]);
    } else {
        return mod_f32x3_2_i32<IDX, false, true>(
            C_plus[idx], C_plus[idx + sizeC], C_plus[idx + 2 * sizeC]);
    }
}

template <Backend BACKEND, unsigned IDX, cublasFillMode_t UPLO, bool FLIP_IMAG = false>
__global__ void mod_hi2mid_AHA_kernel(
    const unsigned n, const size_t ldc, const size_t sizeC,
    const common::hi_t<BACKEND> *__restrict__ C_plus,
    common::mid_t<BACKEND, true> *__restrict__ C_mid //
) {
    constexpr unsigned tile_size = 32;
    constexpr unsigned tile_rows = 8;
    if constexpr (UPLO == CUBLAS_FILL_MODE_UPPER) {
        if (blockIdx.x > blockIdx.y) return;
    } else if constexpr (UPLO == CUBLAS_FILL_MODE_LOWER) {
        if (blockIdx.x < blockIdx.y) return;
    }
    __shared__ int32_t opposite[tile_size][tile_size + 1];
    int32_t direct[tile_size / tile_rows];
    const unsigned row_base = blockIdx.x * tile_size;
    const unsigned col_base = blockIdx.y * tile_size;
    const unsigned tx       = threadIdx.x;
    const unsigned ty       = threadIdx.y;

#pragma unroll
    for (unsigned j = 0; j < tile_size; j += tile_rows) {
        const unsigned yy     = ty + j;
        const unsigned row    = row_base + tx;
        const unsigned col    = col_base + yy;
        const int32_t value   = (row < n && col < n)
                                    ? herk_product_residue<BACKEND, IDX>(C_plus, size_t(col) * ldc + row, sizeC)
                                    : 0;
        direct[j / tile_rows] = value;
        if (blockIdx.x == blockIdx.y) {
            opposite[yy][tx] = value;
        } else {
            const unsigned tr = col_base + tx;
            const unsigned tc = row_base + yy;
            opposite[yy][tx]  = (tr < n && tc < n)
                                    ? herk_product_residue<BACKEND, IDX>(C_plus, size_t(tc) * ldc + tr, sizeC)
                                    : 0;
        }
    }
    __syncthreads();

#pragma unroll
    for (unsigned j = 0; j < tile_size; j += tile_rows) {
        const unsigned yy  = ty + j;
        const unsigned row = row_base + tx;
        const unsigned col = col_base + yy;
        if (row >= n || col >= n) continue;
        if constexpr (UPLO == CUBLAS_FILL_MODE_UPPER) {
            if (row > col) continue;
        } else if constexpr (UPLO == CUBLAS_FILL_MODE_LOWER) {
            if (row < col) continue;
        }
        const int2 value               = reconstruct_complex_2m<BACKEND, IDX, FLIP_IMAG>(direct[j / tile_rows], opposite[tx][yy]);
        C_mid[size_t(col) * ldc + row] = {static_cast<common::mid_t<BACKEND>>(value.x),
                                          static_cast<common::mid_t<BACKEND>>(value.y)};
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

        if constexpr (UPLO == CUBLAS_FILL_MODE_FULL) {

            const dim3 grid((sizeC4 + threads_1d - 1) / threads_1d);

            mod_hi2mid_ge_kernel<BACKEND, IDX>
                <<<grid, threads_1d, 0, stream>>>(
                    sizeC4, C_hix4_1, C_hix4_2, C_midx4);

        } else {

            const size_t ldc4        = ldc >> 2;
            const unsigned threads_y = select_threads_y<BACKEND>(sizeC);
            const dim3 threads(threads_x, threads_y);
            const dim3 grid((ldc4 + threads_x - 1) / threads_x,
                            (n + threads_y - 1) / threads_y);

            mod_hi2mid_tri_kernel<BACKEND, IDX, UPLO>
                <<<grid, threads, 0, stream>>>(
                    n, sizeC4, C_hix4_1, C_hix4_2, C_midx4, ldc4);
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
    if (n == 0) return;
    const size_t sizeC = ldc * n;
    constexpr dim3 threads(32, 8);
    const dim3 grid((n + 31u) / 32u, (n + 31u) / 32u);

    mod_hi2mid_AHA_kernel<BACKEND, IDX, UPLO, FLIP_IMAG>
        <<<grid, threads, 0, stream>>>(
            n, ldc, sizeC, C_hi.ptr0, C_mid + IDX * sizeC);
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
