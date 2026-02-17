#pragma once

namespace complex {

//------------------------------
// C_mid = mod(C_hi, p), scalar
//------------------------------
// INT8
template <int IDX> __forceinline__ __device__ double2 conv_hi2mid_scal(
    const int32_t ArBr, const int32_t AiBi, const int32_t AriBri //
) {
    const int32_t ArBr_rem   = mod_small_nowrap<gemmul8::Backend::INT8, IDX>(ArBr);
    const int32_t AiBi_rem   = mod_small_nowrap<gemmul8::Backend::INT8, IDX>(AiBi);
    const int32_t AriBri_rem = mod_small_nowrap<gemmul8::Backend::INT8, IDX>(AriBri);
    const int32_t Ci         = mod_small<gemmul8::Backend::INT8, IDX>(AriBri_rem - ArBr_rem - AiBi_rem);
    const int32_t Cr         = mod_small<gemmul8::Backend::INT8, IDX>(ArBr_rem - AiBi_rem);
    return {__int2double_rn(Cr), __int2double_rn(Ci)};
}

template <> __forceinline__ __device__ double2 conv_hi2mid_scal<0>(
    const int32_t ArBr, const int32_t AiBi, const int32_t AriBri //
) {
    const int32_t Ci = mod_small<gemmul8::Backend::INT8, 0>(AriBri - ArBr - AiBi);
    const int32_t Cr = mod_small<gemmul8::Backend::INT8, 0>(ArBr - AiBi);
    return {__int2double_rn(Cr), __int2double_rn(Ci)};
}

// FP8
template <int IDX> __forceinline__ __device__ double2 conv_hi2mid_scal(
    const float ArBr0, const float ArBr1, const float ArBr2,
    const float AiBi0, const float AiBi1, const float AiBi2,
    const float AriBri0, const float AriBri1, const float AriBri2 //
) {
    const int32_t ArBr_rem   = mod_f32x3_2_i32_nowrap<IDX>(ArBr0, ArBr1, ArBr2);
    const int32_t AiBi_rem   = mod_f32x3_2_i32_nowrap<IDX>(AiBi0, AiBi1, AiBi2);
    const int32_t AriBri_rem = mod_f32x3_2_i32_nowrap<IDX>(AriBri0, AriBri1, AriBri2);
    const int32_t Ci         = mod_small<gemmul8::Backend::FP8, IDX>(AriBri_rem - ArBr_rem - AiBi_rem);
    const int32_t Cr         = mod_small<gemmul8::Backend::FP8, IDX>(ArBr_rem - AiBi_rem);
    return {__int2double_rn(Cr), __int2double_rn(Ci)};
}

//------------------------------
// C_mid = mod(C_hi, p), vector
//------------------------------
// INT8
template <int IDX> __global__ void conv_hi2mid_kernel(
    const size_t sizeC,                // padding(m*n)/4
    const int4 *__restrict__ C_hix4_1, // input
    const int4 *__restrict__ C_hix4_2, // input
    const int4 *__restrict__ C_hix4_3, // input
    short4 *__restrict__ C_midx4       // output
) {
    const size_t idx = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (idx >= sizeC) return;

    int4 ArBr   = C_hix4_1[idx]; // ArBr
    int4 AiBi   = C_hix4_2[idx]; // AiBi
    int4 AriBri = C_hix4_3[idx]; // (Ar+Ai)(Br+Bi)

    ArBr.x = mod_small_nowrap<gemmul8::Backend::INT8, IDX>(ArBr.x);
    ArBr.y = mod_small_nowrap<gemmul8::Backend::INT8, IDX>(ArBr.y);
    ArBr.z = mod_small_nowrap<gemmul8::Backend::INT8, IDX>(ArBr.z);
    ArBr.w = mod_small_nowrap<gemmul8::Backend::INT8, IDX>(ArBr.w);

    AiBi.x = mod_small_nowrap<gemmul8::Backend::INT8, IDX>(AiBi.x);
    AiBi.y = mod_small_nowrap<gemmul8::Backend::INT8, IDX>(AiBi.y);
    AiBi.z = mod_small_nowrap<gemmul8::Backend::INT8, IDX>(AiBi.z);
    AiBi.w = mod_small_nowrap<gemmul8::Backend::INT8, IDX>(AiBi.w);

    AriBri.x = mod_small_nowrap<gemmul8::Backend::INT8, IDX>(AriBri.x);
    AriBri.y = mod_small_nowrap<gemmul8::Backend::INT8, IDX>(AriBri.y);
    AriBri.z = mod_small_nowrap<gemmul8::Backend::INT8, IDX>(AriBri.z);
    AriBri.w = mod_small_nowrap<gemmul8::Backend::INT8, IDX>(AriBri.w);

    char2 C_mid_x;
    C_mid_x.y = mod_small<gemmul8::Backend::INT8, IDX>(AriBri.x - ArBr.x - AiBi.x); // imag
    C_mid_x.x = mod_small<gemmul8::Backend::INT8, IDX>(ArBr.x - AiBi.x);            // real

    char2 C_mid_y;
    C_mid_y.y = mod_small<gemmul8::Backend::INT8, IDX>(AriBri.y - ArBr.y - AiBi.y); // imag
    C_mid_y.x = mod_small<gemmul8::Backend::INT8, IDX>(ArBr.y - AiBi.y);            // real

    char2 C_mid_z;
    C_mid_z.y = mod_small<gemmul8::Backend::INT8, IDX>(AriBri.z - ArBr.z - AiBi.z); // imag
    C_mid_z.x = mod_small<gemmul8::Backend::INT8, IDX>(ArBr.z - AiBi.z);            // real

    char2 C_mid_w;
    C_mid_w.y = mod_small<gemmul8::Backend::INT8, IDX>(AriBri.w - ArBr.w - AiBi.w); // imag
    C_mid_w.x = mod_small<gemmul8::Backend::INT8, IDX>(ArBr.w - AiBi.w);            // real

    C_midx4[idx] = {*reinterpret_cast<short *>(&C_mid_x), *reinterpret_cast<short *>(&C_mid_y),
                    *reinterpret_cast<short *>(&C_mid_z), *reinterpret_cast<short *>(&C_mid_w)};
}

template <> __global__ void conv_hi2mid_kernel<0>(
    const size_t sizeC,                // padding(m*n)/4
    const int4 *__restrict__ C_hix4_1, // input
    const int4 *__restrict__ C_hix4_2, // input
    const int4 *__restrict__ C_hix4_3, // input
    short4 *__restrict__ C_midx4       // output
) {
    const size_t idx = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (idx >= sizeC) return;

    int4 ArBr   = C_hix4_1[idx]; // ArBr
    int4 AiBi   = C_hix4_2[idx]; // AiBi
    int4 AriBri = C_hix4_3[idx]; // (Ar+Ai)(Br+Bi)

    char2 C_mid_x;
    C_mid_x.y = int8_t(AriBri.x - ArBr.x - AiBi.x); // imag
    C_mid_x.x = int8_t(ArBr.x - AiBi.x);            // real

    char2 C_mid_y;
    C_mid_y.y = int8_t(AriBri.y - ArBr.y - AiBi.y); // imag
    C_mid_y.x = int8_t(ArBr.y - AiBi.y);            // real

    char2 C_mid_z;
    C_mid_z.y = int8_t(AriBri.z - ArBr.z - AiBi.z); // imag
    C_mid_z.x = int8_t(ArBr.z - AiBi.z);            // real

    char2 C_mid_w;
    C_mid_w.y = int8_t(AriBri.w - ArBr.w - AiBi.w); // imag
    C_mid_w.x = int8_t(ArBr.w - AiBi.w);            // real

    C_midx4[idx] = {*reinterpret_cast<short *>(&C_mid_x), *reinterpret_cast<short *>(&C_mid_y),
                    *reinterpret_cast<short *>(&C_mid_z), *reinterpret_cast<short *>(&C_mid_w)};
}

// FP8
template <int IDX> __global__ void conv_hi2mid_kernel(
    const size_t sizeC,                  // padding(m*n)/4
    const float4 *__restrict__ C_hix4_1, // input
    const float4 *__restrict__ C_hix4_2, // input
    const float4 *__restrict__ C_hix4_3, // input
    int4 *__restrict__ C_midx4           // output
) {
    const size_t idx = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (idx >= sizeC) return;

    // ArBr
    const float4 ArBr0 = C_hix4_1[idx];
    const float4 ArBr1 = C_hix4_1[idx + sizeC];
    const float4 ArBr2 = C_hix4_1[idx + 2 * sizeC];
    int4 ArBr;
    ArBr.x = mod_f32x3_2_i32_nowrap<IDX>(ArBr0.x, ArBr1.x, ArBr2.x);
    ArBr.y = mod_f32x3_2_i32_nowrap<IDX>(ArBr0.y, ArBr1.y, ArBr2.y);
    ArBr.z = mod_f32x3_2_i32_nowrap<IDX>(ArBr0.z, ArBr1.z, ArBr2.z);
    ArBr.w = mod_f32x3_2_i32_nowrap<IDX>(ArBr0.w, ArBr1.w, ArBr2.w);

    // AiBi
    const float4 AiBi0 = C_hix4_2[idx];
    const float4 AiBi1 = C_hix4_2[idx + sizeC];
    const float4 AiBi2 = C_hix4_2[idx + 2 * sizeC];
    int4 AiBi;
    AiBi.x = mod_f32x3_2_i32_nowrap<IDX>(AiBi0.x, AiBi1.x, AiBi2.x);
    AiBi.y = mod_f32x3_2_i32_nowrap<IDX>(AiBi0.y, AiBi1.y, AiBi2.y);
    AiBi.z = mod_f32x3_2_i32_nowrap<IDX>(AiBi0.z, AiBi1.z, AiBi2.z);
    AiBi.w = mod_f32x3_2_i32_nowrap<IDX>(AiBi0.w, AiBi1.w, AiBi2.w);

    // (Ar+Ai)(Br+Bi)
    const float4 AriBri0 = C_hix4_3[idx];
    const float4 AriBri1 = C_hix4_3[idx + sizeC];
    const float4 AriBri2 = C_hix4_3[idx + 2 * sizeC];
    int4 AriBri;
    AriBri.x = mod_f32x3_2_i32_nowrap<IDX>(AriBri0.x, AriBri1.x, AriBri2.x);
    AriBri.y = mod_f32x3_2_i32_nowrap<IDX>(AriBri0.y, AriBri1.y, AriBri2.y);
    AriBri.z = mod_f32x3_2_i32_nowrap<IDX>(AriBri0.z, AriBri1.z, AriBri2.z);
    AriBri.w = mod_f32x3_2_i32_nowrap<IDX>(AriBri0.w, AriBri1.w, AriBri2.w);

    short2 C_mid_x;
    C_mid_x.y = mod_small<gemmul8::Backend::FP8, IDX>(AriBri.x - ArBr.x - AiBi.x); // imag
    C_mid_x.x = mod_small<gemmul8::Backend::FP8, IDX>(ArBr.x - AiBi.x);            // real

    short2 C_mid_y;
    C_mid_y.y = mod_small<gemmul8::Backend::FP8, IDX>(AriBri.y - ArBr.y - AiBi.y); // imag
    C_mid_y.x = mod_small<gemmul8::Backend::FP8, IDX>(ArBr.y - AiBi.y);            // real

    short2 C_mid_z;
    C_mid_z.y = mod_small<gemmul8::Backend::FP8, IDX>(AriBri.z - ArBr.z - AiBi.z); // imag
    C_mid_z.x = mod_small<gemmul8::Backend::FP8, IDX>(ArBr.z - AiBi.z);            // real

    short2 C_mid_w;
    C_mid_w.y = mod_small<gemmul8::Backend::FP8, IDX>(AriBri.w - ArBr.w - AiBi.w); // imag
    C_mid_w.x = mod_small<gemmul8::Backend::FP8, IDX>(ArBr.w - AiBi.w);            // real

    C_midx4[idx] = {*reinterpret_cast<int32_t *>(&C_mid_x), *reinterpret_cast<int32_t *>(&C_mid_y),
                    *reinterpret_cast<int32_t *>(&C_mid_z), *reinterpret_cast<int32_t *>(&C_mid_w)};
}

//------------------------------
// Interface!!
//------------------------------
template <gemmul8::Backend backend>
__inline__ void conv_hi2mid(
    const cudaStream_t &stream,       //
    const unsigned i,                 //
    const size_t sizeC,               // padding(m*n) / 4
    hi_t<backend> *const *const C_hi, // input
    midx2_t<backend> *C_mid           // output
) {
    const hix4_t<backend> *C_hix4_1 = reinterpret_cast<hix4_t<backend> *>(C_hi[0]);
    const hix4_t<backend> *C_hix4_2 = reinterpret_cast<hix4_t<backend> *>(C_hi[1]);
    const hix4_t<backend> *C_hix4_3 = reinterpret_cast<hix4_t<backend> *>(C_hi[2]);
    midx8_t<backend> *C_midx4       = reinterpret_cast<midx8_t<backend> *>(C_mid);
    const size_t grid_conv_hi2mid   = (sizeC + threads_conv_hi2mid - 1) / threads_conv_hi2mid;

    switch (i) {
    case 0: conv_hi2mid_kernel<0><<<grid_conv_hi2mid, threads_conv_hi2mid, 0, stream>>>(sizeC, C_hix4_1, C_hix4_2, C_hix4_3, C_midx4); break;
    case 1: conv_hi2mid_kernel<1><<<grid_conv_hi2mid, threads_conv_hi2mid, 0, stream>>>(sizeC, C_hix4_1, C_hix4_2, C_hix4_3, C_midx4); break;
    case 2: conv_hi2mid_kernel<2><<<grid_conv_hi2mid, threads_conv_hi2mid, 0, stream>>>(sizeC, C_hix4_1, C_hix4_2, C_hix4_3, C_midx4); break;
    case 3: conv_hi2mid_kernel<3><<<grid_conv_hi2mid, threads_conv_hi2mid, 0, stream>>>(sizeC, C_hix4_1, C_hix4_2, C_hix4_3, C_midx4); break;
    case 4: conv_hi2mid_kernel<4><<<grid_conv_hi2mid, threads_conv_hi2mid, 0, stream>>>(sizeC, C_hix4_1, C_hix4_2, C_hix4_3, C_midx4); break;
    case 5: conv_hi2mid_kernel<5><<<grid_conv_hi2mid, threads_conv_hi2mid, 0, stream>>>(sizeC, C_hix4_1, C_hix4_2, C_hix4_3, C_midx4); break;
    case 6: conv_hi2mid_kernel<6><<<grid_conv_hi2mid, threads_conv_hi2mid, 0, stream>>>(sizeC, C_hix4_1, C_hix4_2, C_hix4_3, C_midx4); break;
    case 7: conv_hi2mid_kernel<7><<<grid_conv_hi2mid, threads_conv_hi2mid, 0, stream>>>(sizeC, C_hix4_1, C_hix4_2, C_hix4_3, C_midx4); break;
    case 8: conv_hi2mid_kernel<8><<<grid_conv_hi2mid, threads_conv_hi2mid, 0, stream>>>(sizeC, C_hix4_1, C_hix4_2, C_hix4_3, C_midx4); break;
    case 9: conv_hi2mid_kernel<9><<<grid_conv_hi2mid, threads_conv_hi2mid, 0, stream>>>(sizeC, C_hix4_1, C_hix4_2, C_hix4_3, C_midx4); break;
    case 10: conv_hi2mid_kernel<10><<<grid_conv_hi2mid, threads_conv_hi2mid, 0, stream>>>(sizeC, C_hix4_1, C_hix4_2, C_hix4_3, C_midx4); break;
    case 11: conv_hi2mid_kernel<11><<<grid_conv_hi2mid, threads_conv_hi2mid, 0, stream>>>(sizeC, C_hix4_1, C_hix4_2, C_hix4_3, C_midx4); break;
    case 12: conv_hi2mid_kernel<12><<<grid_conv_hi2mid, threads_conv_hi2mid, 0, stream>>>(sizeC, C_hix4_1, C_hix4_2, C_hix4_3, C_midx4); break;
    case 13: conv_hi2mid_kernel<13><<<grid_conv_hi2mid, threads_conv_hi2mid, 0, stream>>>(sizeC, C_hix4_1, C_hix4_2, C_hix4_3, C_midx4); break;
    case 14: conv_hi2mid_kernel<14><<<grid_conv_hi2mid, threads_conv_hi2mid, 0, stream>>>(sizeC, C_hix4_1, C_hix4_2, C_hix4_3, C_midx4); break;
    case 15: conv_hi2mid_kernel<15><<<grid_conv_hi2mid, threads_conv_hi2mid, 0, stream>>>(sizeC, C_hix4_1, C_hix4_2, C_hix4_3, C_midx4); break;
    case 16: conv_hi2mid_kernel<16><<<grid_conv_hi2mid, threads_conv_hi2mid, 0, stream>>>(sizeC, C_hix4_1, C_hix4_2, C_hix4_3, C_midx4); break;
    case 17: conv_hi2mid_kernel<17><<<grid_conv_hi2mid, threads_conv_hi2mid, 0, stream>>>(sizeC, C_hix4_1, C_hix4_2, C_hix4_3, C_midx4); break;
    case 18: conv_hi2mid_kernel<18><<<grid_conv_hi2mid, threads_conv_hi2mid, 0, stream>>>(sizeC, C_hix4_1, C_hix4_2, C_hix4_3, C_midx4); break;
    case 19: conv_hi2mid_kernel<19><<<grid_conv_hi2mid, threads_conv_hi2mid, 0, stream>>>(sizeC, C_hix4_1, C_hix4_2, C_hix4_3, C_midx4); break;
    default: break;
    }
}

} // namespace complex
