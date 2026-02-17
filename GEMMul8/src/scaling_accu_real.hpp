#pragma once

namespace real {
namespace accu {

template <int num_moduli>
__forceinline__ __device__ int32_t compute_sft(int32_t amax) {
    constexpr float log2P = table::log2P<gemmul8::Backend::INT8, num_moduli>; // fld(log2(P-1)/2 - 0.5)
    const float log2amax  = __log2f(__int2float_rn(amax));
    return __float2int_rd(__fmaf_rd(-0x1.0000060000000p-1F, log2amax, log2P));
}

template <int num_moduli>
__forceinline__ __device__ int32_t compute_sft(float amax) {
    constexpr float log2P = table::log2P<gemmul8::Backend::FP8, num_moduli>; // fld(log2(P-1)/2 - 0.5)
    const float log2amax  = __log2f(amax);
    return __float2int_rd(__fmaf_rd(-0x1.0000060000000p-1F, log2amax, log2P));
}

//------------------------------
// Determine row-wise shift values before A*B
//------------------------------
template <gemmul8::Backend backend, typename T>
__global__ void compute_sft_extract_kernel(
    const unsigned m,                // size(A,1)
    const unsigned k,                // size(A,2)
    const T *const __restrict__ A,   // input (lda * k)
    const size_t lda,                // leading dimension
    int16_t *const __restrict__ sftA // exponent of shift values
) {
    using U = underlying_t<T>;
    __shared__ U samax[TILE_DIM][TILE_DIM + 1];

    unsigned row_idx = blockIdx.x * TILE_DIM + threadIdx.x;
    const U amax     = find_amax_tile<T>(m, k, row_idx, A, lda, samax);

    row_idx = blockIdx.x * TILE_DIM + threadIdx.y;
    if (row_idx < m && threadIdx.x == 0) {
        sftA[row_idx] = maxUFP<backend> - Tilogb<U>(amax); // sftA(j)*|A(ij)| < 2^6 (INT8) or 2^8 (FP8)
    }
}

//------------------------------
// Extract first several-bit of abs(A^T)
//------------------------------
template <gemmul8::Backend backend, typename T>
__global__ void extract_A_lo_kernel(
    const unsigned m,                        // size(A,1)
    const unsigned k,                        // size(A,2)
    const T *const __restrict__ A,           // input (lda * k)
    const size_t lda,                        // leading dimension
    low_t<backend> *const __restrict__ A_lo, // output (lda_lo * (m+pad))
    const size_t lda_lo,                     // leading dimension
    int16_t *const __restrict__ sftA         // exponent of shift values
) {
    __shared__ low_t<backend> tile[TILE_DIM][TILE_DIM + 1];

    const auto rowBase = blockIdx.x * TILE_DIM;
    const auto colBase = blockIdx.y * TILE_DIM;

    const auto in_row = rowBase + threadIdx.x;
    const auto in_col = colBase + threadIdx.y;

    const int sft                  = (in_row < m) ? sftA[in_row] : 0;
    const T Atmp                   = (in_row < m && in_col < k) ? A[in_col * lda + in_row] : Tconst<T>::zero();
    tile[threadIdx.y][threadIdx.x] = upperBound_lo<backend, T>(Atmp, sft); // <= 2^6 (INT8) or 2^8 (FP8)
    __syncthreads();

    const auto out_row = colBase + threadIdx.x;
    const auto out_col = rowBase + threadIdx.y;
    if (out_col >= m || out_row >= lda_lo) return;

    A_lo[out_col * lda_lo + out_row] = tile[threadIdx.x][threadIdx.y];
}

//------------------------------
// Extract first several-bit of abs(B)
//------------------------------
template <gemmul8::Backend backend, typename T>
__global__ void extract_B_lo_kernel(
    const unsigned k,                          // size(B,1)
    const T *const __restrict__ B,             // input (ldb * n)
    const size_t ldb,                          // leading dimension
    lowx4_t<backend> *const __restrict__ B_lo, // output (ldb_lo / 4 * n)
    const size_t ldb_lo,                       // leading dimension ldb_lo / 4
    int16_t *const __restrict__ sftB           // exponent of shift values
) {
    __shared__ T shm[32];
    const auto col_idx             = blockIdx.x;
    const T *const __restrict__ in = B + col_idx * ldb;
    const T amax                   = find_amax<T>(in, k, shm);
    const int sft                  = maxUFP<backend> - Tilogb<T>(amax); // sftA(j)*|A(ij)| < 2^6 (INT8) or 2^8 (FP8)
    if (threadIdx.x == 0) {
        sftB[col_idx] = int16_t(sft);
    }

    lowx4_t<backend> *const __restrict__ out = B_lo + col_idx * ldb_lo;
    const unsigned kmax                      = k >> 2;
    unsigned i                               = threadIdx.x;

    for (; i < kmax; i += blockDim.x) {
        unsigned idx = i << 2;

        const T in0                = in[idx];
        const low_t<backend> out4x = upperBound_lo<backend, T>(in0, sft); // <= 2^6 (INT8), 2^8 (FP8)

        const T in1                = in[idx + 1];
        const low_t<backend> out4y = upperBound_lo<backend, T>(in1, sft); // <= 2^6 (INT8), 2^8 (FP8)

        const T in2                = in[idx + 2];
        const low_t<backend> out4z = upperBound_lo<backend, T>(in2, sft); // <= 2^6 (INT8), 2^8 (FP8)

        const T in3                = in[idx + 3];
        const low_t<backend> out4w = upperBound_lo<backend, T>(in3, sft); // <= 2^6 (INT8), 2^8 (FP8)

        out[i] = concat(out4x, out4y, out4z, out4w);
    }
    for (; i < ldb_lo; i += blockDim.x) {
        unsigned idx = i << 2;
        lowx4_t<backend> out4;

        const T in0                = (idx < k) ? in[idx] : Tconst<T>::zero();
        const low_t<backend> out4x = upperBound_lo<backend, T>(in0, sft); // <= 2^6 (INT8), 2^8 (FP8)

        const T in1                = (idx + 1 < k) ? in[idx + 1] : Tconst<T>::zero();
        const low_t<backend> out4y = upperBound_lo<backend, T>(in1, sft); // <= 2^6 (INT8), 2^8 (FP8)

        const T in2                = (idx + 2 < k) ? in[idx + 2] : Tconst<T>::zero();
        const low_t<backend> out4z = upperBound_lo<backend, T>(in2, sft); // <= 2^6 (INT8), 2^8 (FP8)

        const T in3                = (idx + 3 < k) ? in[idx + 3] : Tconst<T>::zero();
        const low_t<backend> out4w = upperBound_lo<backend, T>(in3, sft); // <= 2^6 (INT8), 2^8 (FP8)

        out[i] = concat(out4x, out4y, out4z, out4w);
    }
}

//------------------------------
// Determine row-wise shift values after A*B
//------------------------------
// INT8
template <int num_moduli>
__global__ void compute_sft_rowwise_kernel(
    const unsigned m,                       // size(C_hi,1)
    const unsigned n,                       // size(C_hi,2)
    const int32_t *const __restrict__ C_hi, // input (ldc_hi * n)
    const size_t ldc_hi,                    // leading dimension
    int16_t *const __restrict__ sft         // exponent of shift values
) {
    __shared__ int32_t samax[TILE_DIM][TILE_DIM + 1];

    unsigned row_idx   = blockIdx.x * TILE_DIM + threadIdx.x;
    const int32_t amax = find_max_tile(m, n, row_idx, C_hi, ldc_hi, samax);

    row_idx = blockIdx.x * TILE_DIM + threadIdx.y;
    if (row_idx < m && threadIdx.x == 0) {
        int sft_tmp = sft[row_idx];
        sft_tmp += compute_sft<num_moduli>(amax);
        sft[row_idx] = int16_t(-sft_tmp);
    }
}

// FP8
template <int num_moduli>
__global__ void compute_sft_rowwise_kernel(
    const unsigned m,                     // size(C_hi,1)
    const unsigned n,                     // size(C_hi,2)
    const unsigned k,                     // inner dimension
    const float *const __restrict__ C_hi, // input (ldc_hi * n)
    const size_t ldc_hi,                  // leading dimension
    int16_t *const __restrict__ sft       // exponent of shift values
) {
    __shared__ float samax[TILE_DIM][TILE_DIM + 1];

    unsigned row_idx = blockIdx.x * TILE_DIM + threadIdx.x;
    const float amax = find_max_tile(m, n, k, row_idx, C_hi, ldc_hi, samax);

    row_idx = blockIdx.x * TILE_DIM + threadIdx.y;
    if (row_idx < m && threadIdx.x == 0) {
        int sft_tmp = sft[row_idx];
        sft_tmp += compute_sft<num_moduli>(amax);
        sft[row_idx] = int16_t(-sft_tmp);
    }
}

//------------------------------
// Determine column-wise shift values after A*B
//------------------------------
// INT8
template <int num_moduli>
__global__ void compute_sft_colwise_kernel(
    const unsigned m,                       // size(C_hi,1)
    const int32_t *const __restrict__ C_hi, // input (ldc_hi * n)
    const size_t ldc_hi,                    // leading dimension
    int16_t *const __restrict__ sft         // exponent of shift values
) {
    __shared__ int32_t shm[32];
    const auto col_idx = blockIdx.x;
    const int32_t amax = find_max(C_hi + col_idx * ldc_hi, m, shm);

    if (threadIdx.x == 0) {
        int32_t sft_tmp = sft[col_idx];
        sft_tmp += compute_sft<num_moduli>(amax);
        sft[col_idx] = int16_t(-sft_tmp);
    }
}

// FP8
template <int num_moduli>
__global__ void compute_sft_colwise_kernel(
    const unsigned m,                     // size(C_hi,1)
    const unsigned k,                     // inner dimension
    const float *const __restrict__ C_hi, // input (ldc_hi * n)
    const size_t ldc_hi,                  // leading dimension
    int16_t *const __restrict__ sft       // exponent of shift values
) {
    __shared__ float shm[32];
    const auto col_idx = blockIdx.x;
    const float amax   = find_max(k, C_hi + col_idx * ldc_hi, m, shm);

    if (threadIdx.x == 0) {
        int32_t sft_tmp = sft[col_idx];
        sft_tmp += compute_sft<num_moduli>(amax);
        sft[col_idx] = int16_t(-sft_tmp);
    }
}

//------------------------------
// Convert trunc(B*diag(2^sftB)) to B_lo
//------------------------------
template <gemmul8::Backend backend, typename T, int num_moduli>
__global__ void scalingB_kernel(
    const unsigned k,                          // size(B,1)
    const size_t incB_lo,                      // ldb_lo / 4 * n
    const T *const __restrict__ B,             // input (ldb * n)
    const size_t ldb,                          // leading dimension
    lowx4_t<backend> *const __restrict__ B_lo, // output (ldb_lo / 4 * n)
    const size_t ldb_lo,                       // leading dimension ldb_lo / 4
    int16_t *const __restrict__ sftB           // exponent of shift values
) {
    const auto col_idx = blockIdx.x;
    int sft            = -sftB[col_idx];

    const T *const __restrict__ in     = B + col_idx * ldb;
    lowx4_t<backend> *__restrict__ out = B_lo + col_idx * ldb_lo;
    real::fast::scalingB_device<backend, T, num_moduli>(k, incB_lo, in, out, ldb_lo, sft);
}

//------------------------------
// Launcher!!
//------------------------------
template <gemmul8::Backend backend, typename T>
__forceinline__ void extract_lo_launch(
    const cudaStream_t stream,    //
    const cublasOperation_t op_A, // CUBLAS_OP_N or CUBLAS_OP_T
    const cublasOperation_t op_B, // CUBLAS_OP_N or CUBLAS_OP_T
    const size_t m,               // Number of rows of C
    const size_t n,               // Number of columns of C
    const size_t k,               // Inner dimension
    const T *const A,             // input
    const size_t lda,             // leading dimension
    low_t<backend> *const A_lo,   // output (lda_lo * (m+pad))
    const size_t lda_lo,          // leading dimension
    int16_t *const sftA,          // exponent of shift values for rows of A
    const T *const B,             // input
    const size_t ldb,             // leading dimension
    low_t<backend> *const B_lo,   // output (ldb_lo * n)
    const size_t ldb_lo,          // leading dimension
    int16_t *const sftB,          // exponent of shift values for cols of B
    const bool skip_scalA,        // false (unskip scaling_A) or true (skip scaling_A)
    const bool skip_scalB         // false (unskip scaling_B) or true (skip scaling_B)
) {
    if (!skip_scalA) {
        if (op_A == CUBLAS_OP_N) {
            // m*k -> k*m
            constexpr dim3 threads1(TILE_DIM, TILE_DIM);
            dim3 grid((m + (TILE_DIM - 1)) / TILE_DIM, (lda_lo + (TILE_DIM - 1)) / TILE_DIM);
            compute_sft_extract_kernel<backend, T><<<grid.x, threads1, 0, stream>>>(m, k, A, lda, sftA);
            extract_A_lo_kernel<backend, T><<<grid, threads1, 0, stream>>>(m, k, A, lda, A_lo, lda_lo, sftA);
        } else {
            // k*m -> k*m
            extract_B_lo_kernel<backend, T><<<m, threads_scaling, 0, stream>>>(k, A, lda,
                                                                               reinterpret_cast<lowx4_t<backend> *>(A_lo),
                                                                               lda_lo >> 2, sftA);
        }
    }

    if (!skip_scalB) {
        if (op_B == CUBLAS_OP_N) {
            // k*n -> k*n
            extract_B_lo_kernel<backend, T><<<n, threads_scaling, 0, stream>>>(k, B, ldb,
                                                                               reinterpret_cast<lowx4_t<backend> *>(B_lo),
                                                                               ldb_lo >> 2, sftB);
        } else {
            // n*k -> k*n
            constexpr dim3 threads1(TILE_DIM, TILE_DIM);
            dim3 grid((n + (TILE_DIM - 1)) / TILE_DIM, (ldb_lo + (TILE_DIM - 1)) / TILE_DIM);
            compute_sft_extract_kernel<backend, T><<<grid.x, threads1, 0, stream>>>(n, k, B, ldb, sftB);
            extract_A_lo_kernel<backend, T><<<grid, threads1, 0, stream>>>(n, k, B, ldb, B_lo, ldb_lo, sftB);
        }
    }
}

//------------------------------
// Launcher!!
//------------------------------
template <gemmul8::Backend backend, typename T, int num_moduli>
__forceinline__ void scaling_launch(
    const cudaStream_t stream,    //
    const cublasOperation_t op_A, // CUBLAS_OP_N or CUBLAS_OP_T
    const cublasOperation_t op_B, // CUBLAS_OP_N or CUBLAS_OP_T
    const size_t m,               // Number of rows of C
    const size_t n,               // Number of columns of C
    const size_t k,               // Inner dimension
    const T *const A,             // input
    const size_t lda,             // leading dimension
    low_t<backend> *const A_lo,   // output (lda_lo * (m+pad))
    const size_t lda_lo,          // leading dimension
    const size_t incA_lo,         // increment between the A_lo
    int16_t *const sftA,          // exponent of shift values for rows of A
    const T *const B,             // input
    const size_t ldb,             // leading dimension
    low_t<backend> *const B_lo,   // output (ldb_lo * n)
    const size_t ldb_lo,          // leading dimension
    const size_t incB_lo,         // increment between the B_lo
    int16_t *const sftB,          // exponent of shift values for cols of B
    hi_t<backend> *const C_hi,    // tmp (ldc_hi * n)
    const size_t ldc_hi,          // ((m + 15) >> 4) << 4
    const bool skip_scalA,        // false (unskip scaling_A) or true (skip scaling_A)
    const bool skip_scalB         // false (unskip scaling_B) or true (skip scaling_B)
) {
    if (!skip_scalA) {
        if (op_A == CUBLAS_OP_N) {
            // m*k -> k*m
            dim3 grid((m + (TILE_DIM - 1)) / TILE_DIM, (lda_lo + (TILE_DIM - 1)) / TILE_DIM);
            constexpr dim3 threads1(TILE_DIM, TILE_DIM);
            if constexpr (backend == gemmul8::Backend::INT8) {
                compute_sft_rowwise_kernel<num_moduli><<<grid.x, threads1, 0, stream>>>(m, n, C_hi, ldc_hi, sftA);
            } else {
                compute_sft_rowwise_kernel<num_moduli><<<grid.x, threads1, 0, stream>>>(m, n, k, C_hi, ldc_hi, sftA);
            }
            real::fast::scalingA_kernel<backend, T, num_moduli><<<grid, threads1, 0, stream>>>(m, k, incA_lo, A, lda, A_lo, lda_lo, sftA);
        } else {
            // k*m -> k*m
            constexpr dim3 threads1(TILE_DIM, TILE_DIM);
            if constexpr (backend == gemmul8::Backend::INT8) {
                compute_sft_rowwise_kernel<num_moduli><<<(m + (TILE_DIM - 1)) / TILE_DIM, threads1, 0, stream>>>(m, n, C_hi, ldc_hi, sftA);
            } else {
                compute_sft_rowwise_kernel<num_moduli><<<(m + (TILE_DIM - 1)) / TILE_DIM, threads1, 0, stream>>>(m, n, k, C_hi, ldc_hi, sftA);
            }
            scalingB_kernel<backend, T, num_moduli><<<m, threads_scaling, 0, stream>>>(k, incA_lo >> 2, A, lda,
                                                                                       reinterpret_cast<lowx4_t<backend> *>(A_lo),
                                                                                       lda_lo >> 2, sftA);
        }
    }

    if (!skip_scalB) {
        if constexpr (backend == gemmul8::Backend::INT8) {
            compute_sft_colwise_kernel<num_moduli><<<n, threads_scaling, 0, stream>>>(m, C_hi, ldc_hi, sftB);
        } else {
            compute_sft_colwise_kernel<num_moduli><<<n, threads_scaling, 0, stream>>>(m, k, C_hi, ldc_hi, sftB);
        }
        if (op_B == CUBLAS_OP_N) {
            // k*n -> k*n
            scalingB_kernel<backend, T, num_moduli><<<n, threads_scaling, 0, stream>>>(k, incB_lo >> 2, B, ldb,
                                                                                       reinterpret_cast<lowx4_t<backend> *>(B_lo),
                                                                                       ldb_lo >> 2, sftB);
        } else {
            // n*k -> k*n
            dim3 grid((n + (TILE_DIM - 1)) / TILE_DIM, (ldb_lo + (TILE_DIM - 1)) / TILE_DIM);
            constexpr dim3 threads2(TILE_DIM, TILE_DIM);
            real::fast::scalingA_kernel<backend, T, num_moduli><<<grid, threads2, 0, stream>>>(n, k, incB_lo, B, ldb, B_lo, ldb_lo, sftB);
        }
    }
}

//------------------------------
// Interface!!
//------------------------------
template <gemmul8::Backend backend, typename T>
__inline__ void scaling(
    const cudaStream_t stream,       //
    Handle_t &handle,                // Handle to the cuBLAS/cuBLASLt library context
    const cublasOperation_t op_A,    // CUBLAS_OP_N or CUBLAS_OP_T
    const cublasOperation_t op_B,    // CUBLAS_OP_N or CUBLAS_OP_T
    const size_t m,                  // Number of rows of C
    const size_t n,                  // Number of columns of C
    const size_t k,                  // Inner dimension
    const unsigned num_moduli,       // #moduli
    const T *const A,                // input
    const size_t lda,                // leading dimension
    low_t<backend> *const A_lo,      // output (lda_lo * m)
    low_t<backend> *const A_lo_high, // work/output (lda_lo * (m+pad))
    const size_t lda_lo,             // leading dimension
    const size_t incA_lo,            // increment between the A_lo
    int16_t *const sftA,             // exponent of shift values for rows of A
    const T *const B,                // input
    const size_t ldb,                // leading dimension
    low_t<backend> *const B_lo,      // output (ldb_lo * n)
    low_t<backend> *const B_lo_high, // work/output (lda_lo * n)
    const size_t ldb_lo,             // leading dimension
    const size_t incB_lo,            // increment between the B_lo
    int16_t *const sftB,             // exponent of shift values for cols of B
    hi_t<backend> *const C_hi,       // tmp (ldc_hi * n)
    const size_t ldc_hi,             // ((m + 15) >> 4) << 4
    const bool skip_scalA,           // false (unskip scaling_A) or true (skip scaling_A)
    const bool skip_scalB            // false (unskip scaling_B) or true (skip scaling_B)
) {
    // Extract first several-bit from A and B
    extract_lo_launch<backend, T>(stream, op_A, op_B, m, n, k,
                                  A, lda, A_lo_high, lda_lo, sftA,
                                  B, ldb, B_lo_high, ldb_lo, sftB,
                                  skip_scalA, skip_scalB);

    // C_hi := A_lo^T*B_lo
    constexpr hi_t<backend> one  = 1;
    constexpr hi_t<backend> zero = 0;
    if constexpr (backend == gemmul8::Backend::INT8) {
        gemm_low_prec_i8x1(stream, handle, ldc_hi, n, lda_lo,
                           &one,
                           A_lo_high, lda_lo,
                           B_lo_high, ldb_lo,
                           &zero,
                           C_hi, ldc_hi);
    } else {
        gemm_low_prec_f8x1(stream, handle, ldc_hi, n, lda_lo,
                           &one,
                           A_lo_high, lda_lo,
                           B_lo_high, ldb_lo,
                           &zero,
                           C_hi, ldc_hi);
    }

    // Convert A and B to INT8 matrices
    switch (num_moduli) {
    case 2: scaling_launch<backend, T, 2>(stream, op_A, op_B, m, n, k, A, lda, A_lo, lda_lo, incA_lo, sftA, B, ldb, B_lo, ldb_lo, incB_lo, sftB, C_hi, ldc_hi, skip_scalA, skip_scalB); break;
    case 3: scaling_launch<backend, T, 3>(stream, op_A, op_B, m, n, k, A, lda, A_lo, lda_lo, incA_lo, sftA, B, ldb, B_lo, ldb_lo, incB_lo, sftB, C_hi, ldc_hi, skip_scalA, skip_scalB); break;
    case 4: scaling_launch<backend, T, 4>(stream, op_A, op_B, m, n, k, A, lda, A_lo, lda_lo, incA_lo, sftA, B, ldb, B_lo, ldb_lo, incB_lo, sftB, C_hi, ldc_hi, skip_scalA, skip_scalB); break;
    case 5: scaling_launch<backend, T, 5>(stream, op_A, op_B, m, n, k, A, lda, A_lo, lda_lo, incA_lo, sftA, B, ldb, B_lo, ldb_lo, incB_lo, sftB, C_hi, ldc_hi, skip_scalA, skip_scalB); break;
    case 6: scaling_launch<backend, T, 6>(stream, op_A, op_B, m, n, k, A, lda, A_lo, lda_lo, incA_lo, sftA, B, ldb, B_lo, ldb_lo, incB_lo, sftB, C_hi, ldc_hi, skip_scalA, skip_scalB); break;
    case 7: scaling_launch<backend, T, 7>(stream, op_A, op_B, m, n, k, A, lda, A_lo, lda_lo, incA_lo, sftA, B, ldb, B_lo, ldb_lo, incB_lo, sftB, C_hi, ldc_hi, skip_scalA, skip_scalB); break;
    case 8: scaling_launch<backend, T, 8>(stream, op_A, op_B, m, n, k, A, lda, A_lo, lda_lo, incA_lo, sftA, B, ldb, B_lo, ldb_lo, incB_lo, sftB, C_hi, ldc_hi, skip_scalA, skip_scalB); break;
    case 9: scaling_launch<backend, T, 9>(stream, op_A, op_B, m, n, k, A, lda, A_lo, lda_lo, incA_lo, sftA, B, ldb, B_lo, ldb_lo, incB_lo, sftB, C_hi, ldc_hi, skip_scalA, skip_scalB); break;
    case 10: scaling_launch<backend, T, 10>(stream, op_A, op_B, m, n, k, A, lda, A_lo, lda_lo, incA_lo, sftA, B, ldb, B_lo, ldb_lo, incB_lo, sftB, C_hi, ldc_hi, skip_scalA, skip_scalB); break;
    case 11: scaling_launch<backend, T, 11>(stream, op_A, op_B, m, n, k, A, lda, A_lo, lda_lo, incA_lo, sftA, B, ldb, B_lo, ldb_lo, incB_lo, sftB, C_hi, ldc_hi, skip_scalA, skip_scalB); break;
    case 12: scaling_launch<backend, T, 12>(stream, op_A, op_B, m, n, k, A, lda, A_lo, lda_lo, incA_lo, sftA, B, ldb, B_lo, ldb_lo, incB_lo, sftB, C_hi, ldc_hi, skip_scalA, skip_scalB); break;
    case 13: scaling_launch<backend, T, 13>(stream, op_A, op_B, m, n, k, A, lda, A_lo, lda_lo, incA_lo, sftA, B, ldb, B_lo, ldb_lo, incB_lo, sftB, C_hi, ldc_hi, skip_scalA, skip_scalB); break;
    case 14: scaling_launch<backend, T, 14>(stream, op_A, op_B, m, n, k, A, lda, A_lo, lda_lo, incA_lo, sftA, B, ldb, B_lo, ldb_lo, incB_lo, sftB, C_hi, ldc_hi, skip_scalA, skip_scalB); break;
    case 15: scaling_launch<backend, T, 15>(stream, op_A, op_B, m, n, k, A, lda, A_lo, lda_lo, incA_lo, sftA, B, ldb, B_lo, ldb_lo, incB_lo, sftB, C_hi, ldc_hi, skip_scalA, skip_scalB); break;
    case 16: scaling_launch<backend, T, 16>(stream, op_A, op_B, m, n, k, A, lda, A_lo, lda_lo, incA_lo, sftA, B, ldb, B_lo, ldb_lo, incB_lo, sftB, C_hi, ldc_hi, skip_scalA, skip_scalB); break;
    case 17: scaling_launch<backend, T, 17>(stream, op_A, op_B, m, n, k, A, lda, A_lo, lda_lo, incA_lo, sftA, B, ldb, B_lo, ldb_lo, incB_lo, sftB, C_hi, ldc_hi, skip_scalA, skip_scalB); break;
    case 18: scaling_launch<backend, T, 18>(stream, op_A, op_B, m, n, k, A, lda, A_lo, lda_lo, incA_lo, sftA, B, ldb, B_lo, ldb_lo, incB_lo, sftB, C_hi, ldc_hi, skip_scalA, skip_scalB); break;
    case 19: scaling_launch<backend, T, 19>(stream, op_A, op_B, m, n, k, A, lda, A_lo, lda_lo, incA_lo, sftA, B, ldb, B_lo, ldb_lo, incB_lo, sftB, C_hi, ldc_hi, skip_scalA, skip_scalB); break;
    case 20: scaling_launch<backend, T, 20>(stream, op_A, op_B, m, n, k, A, lda, A_lo, lda_lo, incA_lo, sftA, B, ldb, B_lo, ldb_lo, incB_lo, sftB, C_hi, ldc_hi, skip_scalA, skip_scalB); break;
    default: break;
    }
}

} // namespace accu
} // namespace real
