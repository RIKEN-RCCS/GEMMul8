#pragma once

namespace real {

//------------------------------
// C_mid = mod(C_hi, p)
//------------------------------
// INT8
template <int IDX> __global__ void conv_hi2mid_kernel(
    const size_t sizeC,              // padding(m*n)/4
    const int4 *__restrict__ C_hix4, // input
    char4 *__restrict__ C_midx4      // output
) {
    const size_t idx = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (idx >= sizeC) return;

    const int4 in = C_hix4[idx];
    char4 out;
    out.x = int8_t(mod_small<gemmul8::Backend::INT8, IDX>(in.x));
    out.y = int8_t(mod_small<gemmul8::Backend::INT8, IDX>(in.y));
    out.z = int8_t(mod_small<gemmul8::Backend::INT8, IDX>(in.z));
    out.w = int8_t(mod_small<gemmul8::Backend::INT8, IDX>(in.w));

    C_midx4[idx] = out;
}

// FP8
template <int IDX> __global__ void conv_hi2mid_kernel(
    const size_t sizeC,                // padding(m*n)/4
    const float4 *__restrict__ C_hix4, // input
    short4 *__restrict__ C_midx4       // output
) {
    const size_t idx = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (idx >= sizeC) return;

    const float4 in0 = C_hix4[idx];
    const float4 in1 = C_hix4[idx + sizeC];
    const float4 in2 = C_hix4[idx + 2 * sizeC];
    short4 out;
    out.x = int16_t(mod_f32x3_2_i32<IDX>(in0.x, in1.x, in2.x));
    out.y = int16_t(mod_f32x3_2_i32<IDX>(in0.y, in1.y, in2.y));
    out.z = int16_t(mod_f32x3_2_i32<IDX>(in0.z, in1.z, in2.z));
    out.w = int16_t(mod_f32x3_2_i32<IDX>(in0.w, in1.w, in2.w));

    C_midx4[idx] = out;
}

//------------------------------
// Interface!!
//------------------------------
template <gemmul8::Backend backend>
__inline__ void conv_hi2mid(
    const cudaStream_t &stream,      //
    const unsigned i,                //
    const size_t sizeC,              // padding(m*n) / 4
    const hi_t<backend> *const C_hi, // input
    mid_t<backend> *C_mid            // output
) {
    const hix4_t<backend> *C_hix4 = reinterpret_cast<const hix4_t<backend> *>(C_hi);
    midx4_t<backend> *C_midx4     = reinterpret_cast<midx4_t<backend> *>(C_mid);
    const size_t grid_conv_hi2mid = (sizeC + threads_conv_hi2mid - 1) / threads_conv_hi2mid;

    switch (i) {
    case 0: conv_hi2mid_kernel<0><<<grid_conv_hi2mid, threads_conv_hi2mid, 0, stream>>>(sizeC, C_hix4, C_midx4); break;
    case 1: conv_hi2mid_kernel<1><<<grid_conv_hi2mid, threads_conv_hi2mid, 0, stream>>>(sizeC, C_hix4, C_midx4); break;
    case 2: conv_hi2mid_kernel<2><<<grid_conv_hi2mid, threads_conv_hi2mid, 0, stream>>>(sizeC, C_hix4, C_midx4); break;
    case 3: conv_hi2mid_kernel<3><<<grid_conv_hi2mid, threads_conv_hi2mid, 0, stream>>>(sizeC, C_hix4, C_midx4); break;
    case 4: conv_hi2mid_kernel<4><<<grid_conv_hi2mid, threads_conv_hi2mid, 0, stream>>>(sizeC, C_hix4, C_midx4); break;
    case 5: conv_hi2mid_kernel<5><<<grid_conv_hi2mid, threads_conv_hi2mid, 0, stream>>>(sizeC, C_hix4, C_midx4); break;
    case 6: conv_hi2mid_kernel<6><<<grid_conv_hi2mid, threads_conv_hi2mid, 0, stream>>>(sizeC, C_hix4, C_midx4); break;
    case 7: conv_hi2mid_kernel<7><<<grid_conv_hi2mid, threads_conv_hi2mid, 0, stream>>>(sizeC, C_hix4, C_midx4); break;
    case 8: conv_hi2mid_kernel<8><<<grid_conv_hi2mid, threads_conv_hi2mid, 0, stream>>>(sizeC, C_hix4, C_midx4); break;
    case 9: conv_hi2mid_kernel<9><<<grid_conv_hi2mid, threads_conv_hi2mid, 0, stream>>>(sizeC, C_hix4, C_midx4); break;
    case 10: conv_hi2mid_kernel<10><<<grid_conv_hi2mid, threads_conv_hi2mid, 0, stream>>>(sizeC, C_hix4, C_midx4); break;
    case 11: conv_hi2mid_kernel<11><<<grid_conv_hi2mid, threads_conv_hi2mid, 0, stream>>>(sizeC, C_hix4, C_midx4); break;
    case 12: conv_hi2mid_kernel<12><<<grid_conv_hi2mid, threads_conv_hi2mid, 0, stream>>>(sizeC, C_hix4, C_midx4); break;
    case 13: conv_hi2mid_kernel<13><<<grid_conv_hi2mid, threads_conv_hi2mid, 0, stream>>>(sizeC, C_hix4, C_midx4); break;
    case 14: conv_hi2mid_kernel<14><<<grid_conv_hi2mid, threads_conv_hi2mid, 0, stream>>>(sizeC, C_hix4, C_midx4); break;
    case 15: conv_hi2mid_kernel<15><<<grid_conv_hi2mid, threads_conv_hi2mid, 0, stream>>>(sizeC, C_hix4, C_midx4); break;
    case 16: conv_hi2mid_kernel<16><<<grid_conv_hi2mid, threads_conv_hi2mid, 0, stream>>>(sizeC, C_hix4, C_midx4); break;
    case 17: conv_hi2mid_kernel<17><<<grid_conv_hi2mid, threads_conv_hi2mid, 0, stream>>>(sizeC, C_hix4, C_midx4); break;
    case 18: conv_hi2mid_kernel<18><<<grid_conv_hi2mid, threads_conv_hi2mid, 0, stream>>>(sizeC, C_hix4, C_midx4); break;
    case 19: conv_hi2mid_kernel<19><<<grid_conv_hi2mid, threads_conv_hi2mid, 0, stream>>>(sizeC, C_hix4, C_midx4); break;
    default: break;
    }
}

} // namespace real
