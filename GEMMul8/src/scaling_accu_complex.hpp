#pragma once

namespace complex {
namespace accu {

__device__ __forceinline__ int8_t sub_ru_8bit(int8_t a, int8_t b) { return int8_t(a - b); }
__device__ __forceinline__ __nv_fp8_e4m3 sub_ru_8bit(__nv_fp8_e4m3 a, __nv_fp8_e4m3 b) {
    const __half r_minus_i = __hsub(__half(a), __half(b));
    return fp8_e4m3_ru(r_minus_i);
}

//------------------------------
// Extract first several-bit of abs(A^T)
//------------------------------
template <gemmul8::Backend backend, typename T>
__global__ void extract_A_lo_kernel(
    const unsigned m,                          // size(A,1)
    const unsigned k,                          // size(A,2)
    const T *const __restrict__ A,             // input (lda * k)
    const size_t lda,                          // leading dimension
    low_t<backend> *const __restrict__ A_lo_1, // output (lda_lo * (m+pad))
    low_t<backend> *const __restrict__ A_lo_2, // output (lda_lo * (m+pad))
    low_t<backend> *const __restrict__ A_lo_3, // output (lda_lo * (m+pad))
    const size_t lda_lo,                       // leading dimension
    int16_t *const __restrict__ sftA           // exponent of shift values
) {
    __shared__ lowx2_t<backend> tile[TILE_DIM][TILE_DIM + 1];

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

    const lowx2_t<backend> out = tile[threadIdx.x][threadIdx.y];

    A_lo_1[out_col * lda_lo + out_row] = out.x;                     // Re(A) <= 2^6 (INT8) or 2^8 (FP8)
    A_lo_2[out_col * lda_lo + out_row] = out.y;                     // Im(A) <= 2^6 (INT8) or 2^8 (FP8)
    A_lo_3[out_col * lda_lo + out_row] = sub_ru_8bit(out.x, out.y); // Re(A) - Im(A) <= 2^6 (INT8)
}

//------------------------------
// Extract first several-bit of abs(B)
//------------------------------
template <gemmul8::Backend backend, typename T>
__global__ void extract_B_lo_kernel(
    const unsigned k,                            // size(B,1)
    const T *const __restrict__ B,               // input (ldb * n)
    const size_t ldb,                            // leading dimension
    lowx4_t<backend> *const __restrict__ B_lo_1, // output (ldb_lo / 4 * n)
    lowx4_t<backend> *const __restrict__ B_lo_2, // output (ldb_lo / 4 * n)
    lowx4_t<backend> *const __restrict__ B_lo_3, // output (ldb_lo / 4 * n)
    const size_t ldb_lo,                         // leading dimension ldb_lo / 4
    int16_t *const __restrict__ sftB             // exponent of shift values
) {
    using U = underlying_t<T>;
    __shared__ U shm[32];
    const auto col_idx             = blockIdx.x;
    const T *const __restrict__ in = B + col_idx * ldb;
    const U amax                   = find_amax<T>(in, k, shm);
    const int sft                  = maxUFP<backend> - Tilogb<U>(amax); // sftA(j)*|A(ij)| < 2^6 (INT8) or 2^8 (FP8)
    if (threadIdx.x == 0) {
        sftB[col_idx] = int16_t(sft);
    }

    const size_t inc                           = col_idx * ldb_lo;
    lowx4_t<backend> *const __restrict__ out_1 = B_lo_1 + inc;
    lowx4_t<backend> *const __restrict__ out_2 = B_lo_2 + inc;
    lowx4_t<backend> *const __restrict__ out_3 = B_lo_3 + inc;
    const unsigned kmax                        = k >> 2;
    unsigned i                                 = threadIdx.x;

    for (; i < kmax; i += blockDim.x) {
        unsigned idx = i << 2;

        T in0                 = in[idx];
        lowx2_t<backend> out0 = upperBound_lo<backend, T>(in0, sft); // <= 2^6 (INT8), 2^8 (FP8)

        T in1                 = in[idx + 1];
        lowx2_t<backend> out1 = upperBound_lo<backend, T>(in1, sft); // <= 2^6 (INT8), 2^8 (FP8)

        T in2                 = in[idx + 2];
        lowx2_t<backend> out2 = upperBound_lo<backend, T>(in2, sft); // <= 2^6 (INT8), 2^8 (FP8)

        T in3                 = in[idx + 3];
        lowx2_t<backend> out3 = upperBound_lo<backend, T>(in3, sft); // <= 2^6 (INT8), 2^8 (FP8)

        out_1[i] = concat(out0.x, out1.x, out2.x, out3.x); // Re(B) <= 2^6
        out_2[i] = concat(out0.y, out1.y, out2.y, out3.y); // Im(B) <= 2^6
        out_3[i] = concat(sub_ru_8bit(out0.x, out0.y),
                          sub_ru_8bit(out1.x, out1.y),
                          sub_ru_8bit(out2.x, out2.y),
                          sub_ru_8bit(out3.x, out3.y)); // Re(B) - Im(B) <= 2^6
    }
    for (; i < ldb_lo; i += blockDim.x) {
        unsigned idx = i << 2;

        T in0                 = (idx < k) ? in[idx] : Tconst<T>::zero();
        lowx2_t<backend> out0 = upperBound_lo<backend, T>(in0, sft); // <= 2^6 (INT8), 2^8 (FP8)

        T in1                 = (idx + 1 < k) ? in[idx + 1] : Tconst<T>::zero();
        lowx2_t<backend> out1 = upperBound_lo<backend, T>(in1, sft); // <= 2^6 (INT8), 2^8 (FP8)

        T in2                 = (idx + 2 < k) ? in[idx + 2] : Tconst<T>::zero();
        lowx2_t<backend> out2 = upperBound_lo<backend, T>(in2, sft); // <= 2^6 (INT8), 2^8 (FP8)

        T in3                 = (idx + 3 < k) ? in[idx + 3] : Tconst<T>::zero();
        lowx2_t<backend> out3 = upperBound_lo<backend, T>(in3, sft); // <= 2^6 (INT8), 2^8 (FP8)

        out_1[i] = concat(out0.x, out1.x, out2.x, out3.x); // Re(B) <= 2^6
        out_2[i] = concat(out0.y, out1.y, out2.y, out3.y); // Im(B) <= 2^6
        out_3[i] = concat(sub_ru_8bit(out0.x, out0.y),
                          sub_ru_8bit(out1.x, out1.y),
                          sub_ru_8bit(out2.x, out2.y),
                          sub_ru_8bit(out3.x, out3.y)); // Re(B) - Im(B) <= 2^6
    }
}

//------------------------------
// Determine row-wise shift values after A*B
//------------------------------
// INT8
template <int num_moduli>
__global__ void compute_sft_rowwise_kernel(
    const unsigned m,                         // size(C_hi,1)
    const unsigned n,                         // size(C_hi,2)
    const int32_t *const __restrict__ C_hi_1, // (Re(A)-Im(A)) * (Re(B)-Im(B))
    const int32_t *const __restrict__ C_hi_2, // Re(A)*Im(B) + Im(A)*Re(B)
    const size_t ldc_hi,                      // leading dimension
    int16_t *const __restrict__ sft           // exponent of shift values
) {
    __shared__ int32_t samax[TILE_DIM][TILE_DIM + 1];

    unsigned row_idx   = blockIdx.x * TILE_DIM + threadIdx.x;
    const int32_t amax = find_max_tile_complex(m, n, row_idx, C_hi_1, C_hi_2, ldc_hi, samax);

    row_idx = blockIdx.x * TILE_DIM + threadIdx.y;
    if (row_idx < m && threadIdx.x == 0) {
        int32_t sft_tmp = sft[row_idx];
        sft_tmp += real::accu::compute_sft<num_moduli>(amax);
        sft[row_idx] = int16_t(-sft_tmp);
    }
}

// FP8
template <int num_moduli>
__global__ void compute_sft_rowwise_kernel(
    const unsigned m,                       // size(C_hi,1)
    const unsigned n,                       // size(C_hi,2)
    const unsigned k,                       // inner dimension
    const float *const __restrict__ C_hi_1, // (Re(A)-Im(A)) * (Re(B)-Im(B))
    const float *const __restrict__ C_hi_2, // Re(A)*Im(B)
    const float *const __restrict__ C_hi_3, // Im(A)*Re(B)
    const size_t ldc_hi,                    // leading dimension
    int16_t *const __restrict__ sft         // exponent of shift values
) {
    __shared__ float samax[TILE_DIM][TILE_DIM + 1];

    unsigned row_idx = blockIdx.x * TILE_DIM + threadIdx.x;
    const float amax = find_max_tile_complex(m, n, k, row_idx, C_hi_1, C_hi_2, C_hi_3, ldc_hi, samax);

    row_idx = blockIdx.x * TILE_DIM + threadIdx.y;
    if (row_idx < m && threadIdx.x == 0) {
        int32_t sft_tmp = sft[row_idx];
        sft_tmp += real::accu::compute_sft<num_moduli>(amax);
        sft[row_idx] = int16_t(-sft_tmp);
    }
}

//------------------------------
// Determine column-wise shift values after A*B
//------------------------------
// INT8
template <int num_moduli>
__global__ void compute_sft_colwise_kernel(
    const unsigned m,                         // size(C_hi,1)
    const int32_t *const __restrict__ C_hi_1, // (Re(A)-Im(A)) * (Re(B)-Im(B))
    const int32_t *const __restrict__ C_hi_2, // Re(A)*Im(B) + Im(A)*Re(B)
    const size_t ldc_hi,                      // leading dimension
    int16_t *const __restrict__ sft           // exponent of shift values
) {
    __shared__ int32_t shm[32];
    const auto col_idx  = blockIdx.x;
    const size_t in_idx = col_idx * ldc_hi;
    const int32_t amax  = find_max_complex(C_hi_1 + in_idx, C_hi_2 + in_idx, m, shm);

    if (threadIdx.x == 0) {
        int32_t sft_tmp = sft[col_idx];
        sft_tmp += real::accu::compute_sft<num_moduli>(amax);
        sft[col_idx] = int16_t(-sft_tmp);
    }
}

// FP8
template <int num_moduli>
__global__ void compute_sft_colwise_kernel(
    const unsigned m,                       // size(C_hi,1)
    const unsigned k,                       // inner dimension
    const float *const __restrict__ C_hi_1, // (Re(A)-Im(A)) * (Re(B)-Im(B))
    const float *const __restrict__ C_hi_2, // Re(A)*Im(B)
    const float *const __restrict__ C_hi_3, // Im(A)*Re(B)
    const size_t ldc_hi,                    // leading dimension
    int16_t *const __restrict__ sft         // exponent of shift values
) {
    __shared__ float shm[32];
    const auto col_idx  = blockIdx.x;
    const size_t in_idx = col_idx * ldc_hi;
    const float amax    = find_max_complex(k, C_hi_1 + in_idx, C_hi_2 + in_idx, C_hi_3 + in_idx, m, shm);

    if (threadIdx.x == 0) {
        int32_t sft_tmp = sft[col_idx];
        sft_tmp += real::accu::compute_sft<num_moduli>(amax);
        sft[col_idx] = int16_t(-sft_tmp);
    }
}

//------------------------------
// Convert trunc(B*diag(2^sftB)) to B_lo
//------------------------------
template <gemmul8::Backend backend, typename T, int num_moduli, bool CONJ = false>
__global__ void scalingB_kernel(
    const unsigned k,                            // size(B,1)
    const size_t incB_lo,                        // ldb_lo / 4 * n
    const T *const __restrict__ B,               // input (ldb * n)
    const size_t ldb,                            // leading dimension
    lowx4_t<backend> *const __restrict__ B_lo_1, // output (ldb_lo / 4 * n)
    lowx4_t<backend> *const __restrict__ B_lo_2, // output (ldb_lo / 4 * n)
    lowx4_t<backend> *const __restrict__ B_lo_3, // output (ldb_lo / 4 * n)
    const size_t ldb_lo,                         // leading dimension ldb_lo / 4
    int16_t *const __restrict__ sftB             // exponent of shift values
) {
    const auto col_idx = blockIdx.x;
    int sft            = -sftB[col_idx];

    const T *const __restrict__ in       = B + col_idx * ldb;
    const size_t out_idx                 = col_idx * ldb_lo;
    lowx4_t<backend> *__restrict__ out_1 = B_lo_1 + out_idx;
    lowx4_t<backend> *__restrict__ out_2 = B_lo_2 + out_idx;
    lowx4_t<backend> *__restrict__ out_3 = B_lo_3 + out_idx;
    complex::fast::scalingB_device<backend, T, num_moduli, CONJ>(k, incB_lo, in, out_1, out_2, out_3, ldb_lo, sft);
}

//------------------------------
// Launcher!!
//------------------------------
template <gemmul8::Backend backend, typename T>
__forceinline__ void extract_lo_launch(
    const cudaStream_t stream,         //
    const cublasOperation_t op_A,      // CUBLAS_OP_N, CUBLAS_OP_T, or CUBLAS_OP_C
    const cublasOperation_t op_B,      // CUBLAS_OP_N, CUBLAS_OP_T, or CUBLAS_OP_C
    const size_t m,                    // Number of rows of C
    const size_t n,                    // Number of columns of C
    const size_t k,                    // Inner dimension
    const T *const A,                  // input
    const size_t lda,                  // leading dimension
    low_t<backend> *const *const A_lo, // output (lda_lo * (m+pad))
    const size_t lda_lo,               // leading dimension
    int16_t *const sftA,               // exponent of shift values for rows of A
    const T *const B,                  // input
    const size_t ldb,                  // leading dimension
    low_t<backend> *const *const B_lo, // output (ldb_lo * n)
    const size_t ldb_lo,               // leading dimension
    int16_t *const sftB,               // exponent of shift values for cols of B
    const bool skip_scalA,             // false (unskip scaling_A) or true (skip scaling_A)
    const bool skip_scalB              // false (unskip scaling_B) or true (skip scaling_B)
) {
    if (!skip_scalA) {
        if (op_A == CUBLAS_OP_N) {
            // m*k -> k*m
            constexpr dim3 threads1(TILE_DIM, TILE_DIM);
            dim3 grid((m + (TILE_DIM - 1)) / TILE_DIM, (lda_lo + (TILE_DIM - 1)) / TILE_DIM);
            real::accu::compute_sft_extract_kernel<backend, T><<<grid.x, threads1, 0, stream>>>(m, k, A, lda, sftA);
            extract_A_lo_kernel<backend, T><<<grid, threads1, 0, stream>>>(m, k, A, lda, A_lo[0], A_lo[1], A_lo[2], lda_lo, sftA);
        } else {
            // k*m -> k*m
            extract_B_lo_kernel<backend, T><<<m, threads_scaling, 0, stream>>>(k, A, lda,
                                                                               reinterpret_cast<lowx4_t<backend> *>(A_lo[0]),
                                                                               reinterpret_cast<lowx4_t<backend> *>(A_lo[1]),
                                                                               reinterpret_cast<lowx4_t<backend> *>(A_lo[2]),
                                                                               lda_lo >> 2, sftA);
        }
    }

    if (!skip_scalB) {
        if (op_B == CUBLAS_OP_N) {
            // k*n -> k*n
            extract_B_lo_kernel<backend, T><<<n, threads_scaling, 0, stream>>>(k, B, ldb,
                                                                               reinterpret_cast<lowx4_t<backend> *>(B_lo[0]),
                                                                               reinterpret_cast<lowx4_t<backend> *>(B_lo[1]),
                                                                               reinterpret_cast<lowx4_t<backend> *>(B_lo[2]),
                                                                               ldb_lo >> 2, sftB);
        } else {
            // n*k -> k*n
            constexpr dim3 threads1(TILE_DIM, TILE_DIM);
            dim3 grid((n + (TILE_DIM - 1)) / TILE_DIM, (ldb_lo + (TILE_DIM - 1)) / TILE_DIM);
            real::accu::compute_sft_extract_kernel<backend, T><<<grid.x, threads1, 0, stream>>>(n, k, B, ldb, sftB);
            extract_A_lo_kernel<backend, T><<<grid, threads1, 0, stream>>>(n, k, B, ldb, B_lo[0], B_lo[1], B_lo[2], ldb_lo, sftB);
        }
    }
}

//------------------------------
// Launcher!!
//------------------------------
template <gemmul8::Backend backend, typename T, int num_moduli>
__forceinline__ void scaling_launch(
    const cudaStream_t stream,         //
    const cublasOperation_t op_A,      // CUBLAS_OP_N, CUBLAS_OP_T, or CUBLAS_OP_C
    const cublasOperation_t op_B,      // CUBLAS_OP_N, CUBLAS_OP_T, or CUBLAS_OP_C
    const size_t m,                    // Number of rows of C
    const size_t n,                    // Number of columns of C
    const size_t k,                    // Inner dimension
    const T *const A,                  // input
    const size_t lda,                  // leading dimension
    low_t<backend> *const *const A_lo, // output (lda_lo * (m+pad))
    const size_t lda_lo,               // leading dimension
    const size_t incA_lo,              // increment between the A_lo
    int16_t *const sftA,               // exponent of shift values for rows of A
    const T *const B,                  // input
    const size_t ldb,                  // leading dimension
    low_t<backend> *const *const B_lo, // output (ldb_lo * n)
    const size_t ldb_lo,               // leading dimension
    const size_t incB_lo,              // increment between the B_lo
    int16_t *const sftB,               // exponent of shift values for cols of B
    hi_t<backend> *const *const C_hi,  // tmp (ldc_hi * n)
    const size_t ldc_hi,               // ((m + 15) >> 4) << 4
    const bool skip_scalA,             // false (unskip scaling_A) or true (skip scaling_A)
    const bool skip_scalB              // false (unskip scaling_B) or true (skip scaling_B)
) {
    if (!skip_scalA) {
        if (op_A == CUBLAS_OP_N) {
            // m*k -> k*m
            dim3 grid((m + (TILE_DIM - 1)) / TILE_DIM, (lda_lo + (TILE_DIM - 1)) / TILE_DIM);
            constexpr dim3 threads1(TILE_DIM, TILE_DIM);
            if constexpr (backend == gemmul8::Backend::INT8) {
                compute_sft_rowwise_kernel<num_moduli><<<grid.x, threads1, 0, stream>>>(m, n, C_hi[0], C_hi[1], ldc_hi, sftA);
            } else {
                compute_sft_rowwise_kernel<num_moduli><<<grid.x, threads1, 0, stream>>>(m, n, k, C_hi[0], C_hi[1], C_hi[2], ldc_hi, sftA);
            }
            complex::fast::scalingA_kernel<backend, T, num_moduli><<<grid, threads1, 0, stream>>>(m, k, incA_lo, A, lda, A_lo[0], A_lo[1], A_lo[2], lda_lo, sftA);
        } else {
            // k*m -> k*m
            constexpr dim3 threads1(TILE_DIM, TILE_DIM);
            if constexpr (backend == gemmul8::Backend::INT8) {
                compute_sft_rowwise_kernel<num_moduli><<<(m + (TILE_DIM - 1)) / TILE_DIM, threads1, 0, stream>>>(m, n, C_hi[0], C_hi[1], ldc_hi, sftA);
            } else {
                compute_sft_rowwise_kernel<num_moduli><<<(m + (TILE_DIM - 1)) / TILE_DIM, threads1, 0, stream>>>(m, n, k, C_hi[0], C_hi[1], C_hi[2], ldc_hi, sftA);
            }
            if (op_A == CUBLAS_OP_T) {
                scalingB_kernel<backend, T, num_moduli, false><<<m, threads_scaling, 0, stream>>>(k, incA_lo >> 2, A, lda,
                                                                                                  reinterpret_cast<lowx4_t<backend> *>(A_lo[0]),
                                                                                                  reinterpret_cast<lowx4_t<backend> *>(A_lo[1]),
                                                                                                  reinterpret_cast<lowx4_t<backend> *>(A_lo[2]),
                                                                                                  lda_lo >> 2, sftA);
            } else {
                scalingB_kernel<backend, T, num_moduli, true><<<m, threads_scaling, 0, stream>>>(k, incA_lo >> 2, A, lda,
                                                                                                 reinterpret_cast<lowx4_t<backend> *>(A_lo[0]),
                                                                                                 reinterpret_cast<lowx4_t<backend> *>(A_lo[1]),
                                                                                                 reinterpret_cast<lowx4_t<backend> *>(A_lo[2]),
                                                                                                 lda_lo >> 2, sftA);
            }
        }
    }

    if (!skip_scalB) {
        if constexpr (backend == gemmul8::Backend::INT8) {
            compute_sft_colwise_kernel<num_moduli><<<n, threads_scaling, 0, stream>>>(m, C_hi[0], C_hi[1], ldc_hi, sftB);
        } else {
            compute_sft_colwise_kernel<num_moduli><<<n, threads_scaling, 0, stream>>>(m, k, C_hi[0], C_hi[1], C_hi[2], ldc_hi, sftB);
        }
        if (op_B == CUBLAS_OP_N) {
            // k*n -> k*n
            scalingB_kernel<backend, T, num_moduli><<<n, threads_scaling, 0, stream>>>(k, incB_lo >> 2, B, ldb,
                                                                                       reinterpret_cast<lowx4_t<backend> *>(B_lo[0]),
                                                                                       reinterpret_cast<lowx4_t<backend> *>(B_lo[1]),
                                                                                       reinterpret_cast<lowx4_t<backend> *>(B_lo[2]),
                                                                                       ldb_lo >> 2, sftB);
        } else {
            // n*k -> k*n
            dim3 grid((n + (TILE_DIM - 1)) / TILE_DIM, (ldb_lo + (TILE_DIM - 1)) / TILE_DIM);
            constexpr dim3 threads2(TILE_DIM, TILE_DIM);
            if (op_B == CUBLAS_OP_T) {
                complex::fast::scalingA_kernel<backend, T, num_moduli, false><<<grid, threads2, 0, stream>>>(n, k, incB_lo, B, ldb, B_lo[0], B_lo[1], B_lo[2], ldb_lo, sftB);
            } else {
                complex::fast::scalingA_kernel<backend, T, num_moduli, true><<<grid, threads2, 0, stream>>>(n, k, incB_lo, B, ldb, B_lo[0], B_lo[1], B_lo[2], ldb_lo, sftB);
            }
        }
    }
}

//------------------------------
// Interface!!
//------------------------------
template <gemmul8::Backend backend, typename T>
__inline__ void scaling(
    const cudaStream_t stream,              //
    Handle_t &handle,                       // Handle to the cuBLAS/cuBLASLt library context
    const cublasOperation_t op_A,           // CUBLAS_OP_N, CUBLAS_OP_T, or CUBLAS_OP_C
    const cublasOperation_t op_B,           // CUBLAS_OP_N, CUBLAS_OP_T, or CUBLAS_OP_C
    const size_t m,                         // Number of rows of C
    const size_t n,                         // Number of columns of C
    const size_t k,                         // Inner dimension
    const unsigned num_moduli,              // #moduli
    const T *const A,                       // input
    const size_t lda,                       // leading dimension
    low_t<backend> *const *const A_lo,      // output (lda_lo * m)
    low_t<backend> *const *const A_lo_high, // work/output (lda_lo * (m+pad))
    const size_t lda_lo,                    // leading dimension
    const size_t incA_lo,                   // increment between the A_lo
    int16_t *const sftA,                    // exponent of shift values for rows of A
    const T *const B,                       // input
    const size_t ldb,                       // leading dimension
    low_t<backend> *const *const B_lo,      // output (ldb_lo * n)
    low_t<backend> *const *const B_lo_high, // work/output (ldb_lo * n)
    const size_t ldb_lo,                    // leading dimension
    const size_t incB_lo,                   // increment between the B_lo
    int16_t *const sftB,                    // exponent of shift values for cols of B
    hi_t<backend> *const *const C_hi,       // tmp (ldc_hi * n)
    const size_t ldc_hi,                    // ((m + 15) >> 4) << 4
    const bool skip_scalA,                  // false (unskip scaling_A) or true (skip scaling_A)
    const bool skip_scalB                   // false (unskip scaling_B) or true (skip scaling_B)
) {
    // Extract first several-bit from A and B
    extract_lo_launch<backend, T>(stream, op_A, op_B, m, n, k,
                                  A, lda, A_lo_high, lda_lo, sftA,
                                  B, ldb, B_lo_high, ldb_lo, sftB,
                                  skip_scalA, skip_scalB);

    // Im(A)*Re(B)
    constexpr hi_t<backend> one  = 1;
    constexpr hi_t<backend> zero = 0;
    if constexpr (backend == gemmul8::Backend::INT8) {
        // C_hi[0] = (Re(A)-Im(A)) * (Re(B)-Im(B))
        // C_hi[1] = Re(A)*Im(B) + Im(A)*Re(B)
        gemm_low_prec_i8x3(stream, handle, ldc_hi, n, lda_lo,
                           &one, &one, &one,
                           A_lo_high[0], A_lo_high[1], A_lo_high[2], lda_lo,
                           B_lo_high[1], B_lo_high[0], B_lo_high[2], ldb_lo,
                           &zero, &one, &zero,
                           C_hi[1], C_hi[1], C_hi[0], ldc_hi);
    } else {
        // C_hi[0] = (Re(A)-Im(A)) * (Re(B)-Im(B))
        // C_hi[1] = Re(A)*Im(B)
        // C_hi[2] = Im(A)*Re(B)
        gemm_low_prec_f8x3(stream, handle, ldc_hi, n, lda_lo,
                           &one, &one, &one,
                           A_lo_high[0], A_lo_high[1], A_lo_high[2], lda_lo,
                           B_lo_high[1], B_lo_high[0], B_lo_high[2], ldb_lo,
                           &zero, &zero, &zero,
                           C_hi[1], C_hi[2], C_hi[0], ldc_hi);
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
} // namespace complex
