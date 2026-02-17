#pragma once

namespace real {
namespace fast {

template <gemmul8::Backend backend, int num_moduli> __forceinline__ __device__ int compute_sft(double amax, double vecnrm) {
    const int exponent    = Tilogb<double>(vecnrm);
    const float vecnrmf   = __double2float_ru(scalbn(vecnrm, -exponent));
    const float log2vsum  =  __fadd_ru(__log2f(vecnrmf), exponent);
    const float log2vnrm  = __fmul_ru(0x1.0000060000000p-1F, log2vsum);
    constexpr float log2P = table::log2P<backend, num_moduli>; 
    const float exp1      = __fsub_rd(__fsub_rd(log2P, 1.5f), fmaxf(1.0f, log2vnrm));
    return __float2int_rd(exp1) - Tilogb<float>(amax);
}

template <gemmul8::Backend backend, int num_moduli> __forceinline__ __device__ int compute_sft(float amax, float vecnrm) {
    const float log2vsum  = __log2f(vecnrm);
    const float log2vnrm  = __fmul_ru(0x1.0000060000000p-1F, log2vsum);
    constexpr float log2P = table::log2P<backend, num_moduli>; 
    const float exp1      = __fsub_rd(__fsub_rd(log2P, 1.5f), fmaxf(1.0f, log2vnrm));
    return __float2int_rd(exp1) - Tilogb<float>(amax);
}

//------------------------------
// Determine row-wise shift values
//------------------------------
template <gemmul8::Backend backend, typename T, int num_moduli>
__global__ void compute_sftA_kernel(
    const unsigned m,                // size(A,1)
    const unsigned k,                // size(A,2)
    const T *const __restrict__ A,   // input (lda * k)
    const size_t lda,                // leading dimension
    int16_t *const __restrict__ sftA // exponent of shift values
) {
    using U = underlying_t<T>;
    __shared__ U samax[TILE_DIM][TILE_DIM + 1];
    __shared__ U ssum[TILE_DIM][TILE_DIM + 1];

    unsigned row_idx = blockIdx.x * TILE_DIM + threadIdx.x;

    U sum;
    U amax = find_amax_and_nrm_tile<T>(m, k, row_idx, A, lda, samax, ssum, sum);

    row_idx = blockIdx.x * TILE_DIM + threadIdx.y;
    if (row_idx < m && threadIdx.x == 0) {
        int sft       = compute_sft<backend, num_moduli>(amax, sum);
        sftA[row_idx] = int16_t(-sft);
    }
}

//------------------------------
// Convert trunc(diag(2^sftA)*A) to A_lo
//------------------------------
template <gemmul8::Backend backend, typename T, int num_moduli>
__global__ void scalingA_kernel(
    const unsigned m,                        // size(A,1)
    const unsigned k,                        // size(A,2)
    const size_t incA_lo,                    // lda_lo * (m+pad)
    const T *const __restrict__ A,           // input (lda * k)
    const size_t lda,                        // leading dimension
    low_t<backend> *const __restrict__ A_lo, // output (lda_lo * (m+pad))
    const size_t lda_lo,                     // leading dimension
    int16_t *const __restrict__ sftA         // exponent of shift values (m+pad)
) {
    using ValT = decltype(trunc_scalbn<backend, T, num_moduli>::run(T{}, 0));

    __shared__ ValT tile[TILE_DIM][TILE_DIM + 1];

    const auto rowBase = blockIdx.x * TILE_DIM;
    const auto colBase = blockIdx.y * TILE_DIM;

    const auto in_row              = rowBase + threadIdx.x;
    const auto in_col              = colBase + threadIdx.y;
    const int sft                  = (in_row < m) ? -sftA[in_row] : 0;
    const T Atmp                   = (in_row < m && in_col < k) ? A[in_col * lda + in_row] : Tconst<T>::zero();
    tile[threadIdx.y][threadIdx.x] = trunc_scalbn<backend, T, num_moduli>::run(Atmp, sft);
    __syncthreads();

    const ValT in = tile[threadIdx.x][threadIdx.y];

    const auto out_col = rowBase + threadIdx.y;
    const auto out_row = colBase + threadIdx.x;
    if (out_col >= m || out_row >= lda_lo) return;

    low_t<backend> *__restrict__ out = A_lo + out_col * lda_lo + out_row;

    ModUnroll<num_moduli, ValT>::run(out, incA_lo, in);
}

//------------------------------
// Convert trunc(B*diag(2^sftB)) to B_lo
//------------------------------
template <gemmul8::Backend backend, typename T, int num_moduli>
__forceinline__ __device__ void scalingB_device(
    const unsigned k,                   // size(B,1)
    const size_t incB_lo,               // ldb_lo / 4 * n
    const T *const __restrict__ in,     // input (ldb * n)
    lowx4_t<backend> *__restrict__ out, // output (ldb_lo / 4 * n)
    const size_t ldb_lo,                // leading dimension ldb_lo / 4
    const int sft                       // shift value
) {
    using ValT = decltype(trunc_scalbn<backend, T, num_moduli>::run(T{}, 0));

    const unsigned kmax = k >> 2;
    unsigned i          = threadIdx.x;

    for (; i < kmax; i += blockDim.x) {
        unsigned idx = i << 2;

        T in0 = in[idx];
        T in1 = in[idx + 1];
        T in2 = in[idx + 2];
        T in3 = in[idx + 3];

        ValT v0 = trunc_scalbn<backend, T, num_moduli>::run(in0, sft);
        ValT v1 = trunc_scalbn<backend, T, num_moduli>::run(in1, sft);
        ValT v2 = trunc_scalbn<backend, T, num_moduli>::run(in2, sft);
        ValT v3 = trunc_scalbn<backend, T, num_moduli>::run(in3, sft);

        ModUnroll<num_moduli, ValT>::run(out + i, incB_lo, v0, v1, v2, v3);
    }
    for (; i < ldb_lo; i += blockDim.x) {
        unsigned idx = i << 2;

        T in0 = (idx < k) ? in[idx] : Tconst<T>::zero();
        T in1 = (idx + 1 < k) ? in[idx + 1] : Tconst<T>::zero();
        T in2 = (idx + 2 < k) ? in[idx + 2] : Tconst<T>::zero();
        T in3 = (idx + 3 < k) ? in[idx + 3] : Tconst<T>::zero();

        ValT v0 = trunc_scalbn<backend, T, num_moduli>::run(in0, sft);
        ValT v1 = trunc_scalbn<backend, T, num_moduli>::run(in1, sft);
        ValT v2 = trunc_scalbn<backend, T, num_moduli>::run(in2, sft);
        ValT v3 = trunc_scalbn<backend, T, num_moduli>::run(in3, sft);

        ModUnroll<num_moduli, ValT>::run(out + i, incB_lo, v0, v1, v2, v3);
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
    __shared__ T shm[64];
    const auto col_idx             = blockIdx.x;
    const T *const __restrict__ in = B + col_idx * ldb;
    T vecnrm;
    const T amax  = find_amax_and_nrm<T>(in, k, shm, vecnrm);
    const int sft = compute_sft<backend, num_moduli>(amax, vecnrm);
    if (threadIdx.x == 0) {
        sftB[col_idx] = int16_t(-sft);
    }

    lowx4_t<backend> *__restrict__ out = B_lo + col_idx * ldb_lo;
    scalingB_device<backend, T, num_moduli>(k, incB_lo, in, out, ldb_lo, sft);
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
    const bool skip_scalA,        // false (unskip scaling_A) or true (skip scaling_A)
    const bool skip_scalB         // false (unskip scaling_B) or true (skip scaling_B)
) {
    if (!skip_scalA) {
        if (op_A == CUBLAS_OP_N) {
            // m*k -> k*m
            dim3 grid((m + (TILE_DIM - 1)) / TILE_DIM, (lda_lo + (TILE_DIM - 1)) / TILE_DIM);
            constexpr dim3 threads1(TILE_DIM, TILE_DIM);
            compute_sftA_kernel<backend, T, num_moduli><<<grid.x, threads1, 0, stream>>>(m, k, A, lda, sftA);
            scalingA_kernel<backend, T, num_moduli><<<grid, threads1, 0, stream>>>(m, k, incA_lo, A, lda, A_lo, lda_lo, sftA);
        } else {
            // k*m -> k*m
            scalingB_kernel<backend, T, num_moduli><<<m, threads_scaling, 0, stream>>>(k, incA_lo >> 2, A, lda, reinterpret_cast<lowx4_t<backend> *>(A_lo), lda_lo >> 2, sftA);
        }
    }

    if (!skip_scalB) {
        if (op_B == CUBLAS_OP_N) {
            // k*n -> k*n
            scalingB_kernel<backend, T, num_moduli><<<n, threads_scaling, 0, stream>>>(k, incB_lo >> 2, B, ldb, reinterpret_cast<lowx4_t<backend> *>(B_lo), ldb_lo >> 2, sftB);
        } else {
            // n*k -> k*n
            dim3 grid((n + (TILE_DIM - 1)) / TILE_DIM, (ldb_lo + (TILE_DIM - 1)) / TILE_DIM);
            constexpr dim3 threads1(TILE_DIM, TILE_DIM);
            compute_sftA_kernel<backend, T, num_moduli><<<grid.x, threads1, 0, stream>>>(n, k, B, ldb, sftB);
            scalingA_kernel<backend, T, num_moduli><<<grid, threads1, 0, stream>>>(n, k, incB_lo, B, ldb, B_lo, ldb_lo, sftB);
        }
    }
}

//------------------------------
// Interface!!
//------------------------------
template <gemmul8::Backend backend, typename T>
__inline__ void scaling(
    const cudaStream_t stream,    //
    const cublasOperation_t op_A, // CUBLAS_OP_N or CUBLAS_OP_T
    const cublasOperation_t op_B, // CUBLAS_OP_N or CUBLAS_OP_T
    const size_t m,               // Number of rows of C
    const size_t n,               // Number of columns of C
    const size_t k,               // Inner dimension
    const unsigned num_moduli,    // #moduli
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
    const bool skip_scalA,        // false (unskip scaling_A) or true (skip scaling_A)
    const bool skip_scalB         // false (unskip scaling_B) or true (skip scaling_B)
) {
    switch (num_moduli) {
    case 2: scaling_launch<backend, T, 2>(stream, op_A, op_B, m, n, k, A, lda, A_lo, lda_lo, incA_lo, sftA, B, ldb, B_lo, ldb_lo, incB_lo, sftB, skip_scalA, skip_scalB); break;
    case 3: scaling_launch<backend, T, 3>(stream, op_A, op_B, m, n, k, A, lda, A_lo, lda_lo, incA_lo, sftA, B, ldb, B_lo, ldb_lo, incB_lo, sftB, skip_scalA, skip_scalB); break;
    case 4: scaling_launch<backend, T, 4>(stream, op_A, op_B, m, n, k, A, lda, A_lo, lda_lo, incA_lo, sftA, B, ldb, B_lo, ldb_lo, incB_lo, sftB, skip_scalA, skip_scalB); break;
    case 5: scaling_launch<backend, T, 5>(stream, op_A, op_B, m, n, k, A, lda, A_lo, lda_lo, incA_lo, sftA, B, ldb, B_lo, ldb_lo, incB_lo, sftB, skip_scalA, skip_scalB); break;
    case 6: scaling_launch<backend, T, 6>(stream, op_A, op_B, m, n, k, A, lda, A_lo, lda_lo, incA_lo, sftA, B, ldb, B_lo, ldb_lo, incB_lo, sftB, skip_scalA, skip_scalB); break;
    case 7: scaling_launch<backend, T, 7>(stream, op_A, op_B, m, n, k, A, lda, A_lo, lda_lo, incA_lo, sftA, B, ldb, B_lo, ldb_lo, incB_lo, sftB, skip_scalA, skip_scalB); break;
    case 8: scaling_launch<backend, T, 8>(stream, op_A, op_B, m, n, k, A, lda, A_lo, lda_lo, incA_lo, sftA, B, ldb, B_lo, ldb_lo, incB_lo, sftB, skip_scalA, skip_scalB); break;
    case 9: scaling_launch<backend, T, 9>(stream, op_A, op_B, m, n, k, A, lda, A_lo, lda_lo, incA_lo, sftA, B, ldb, B_lo, ldb_lo, incB_lo, sftB, skip_scalA, skip_scalB); break;
    case 10: scaling_launch<backend, T, 10>(stream, op_A, op_B, m, n, k, A, lda, A_lo, lda_lo, incA_lo, sftA, B, ldb, B_lo, ldb_lo, incB_lo, sftB, skip_scalA, skip_scalB); break;
    case 11: scaling_launch<backend, T, 11>(stream, op_A, op_B, m, n, k, A, lda, A_lo, lda_lo, incA_lo, sftA, B, ldb, B_lo, ldb_lo, incB_lo, sftB, skip_scalA, skip_scalB); break;
    case 12: scaling_launch<backend, T, 12>(stream, op_A, op_B, m, n, k, A, lda, A_lo, lda_lo, incA_lo, sftA, B, ldb, B_lo, ldb_lo, incB_lo, sftB, skip_scalA, skip_scalB); break;
    case 13: scaling_launch<backend, T, 13>(stream, op_A, op_B, m, n, k, A, lda, A_lo, lda_lo, incA_lo, sftA, B, ldb, B_lo, ldb_lo, incB_lo, sftB, skip_scalA, skip_scalB); break;
    case 14: scaling_launch<backend, T, 14>(stream, op_A, op_B, m, n, k, A, lda, A_lo, lda_lo, incA_lo, sftA, B, ldb, B_lo, ldb_lo, incB_lo, sftB, skip_scalA, skip_scalB); break;
    case 15: scaling_launch<backend, T, 15>(stream, op_A, op_B, m, n, k, A, lda, A_lo, lda_lo, incA_lo, sftA, B, ldb, B_lo, ldb_lo, incB_lo, sftB, skip_scalA, skip_scalB); break;
    case 16: scaling_launch<backend, T, 16>(stream, op_A, op_B, m, n, k, A, lda, A_lo, lda_lo, incA_lo, sftA, B, ldb, B_lo, ldb_lo, incB_lo, sftB, skip_scalA, skip_scalB); break;
    case 17: scaling_launch<backend, T, 17>(stream, op_A, op_B, m, n, k, A, lda, A_lo, lda_lo, incA_lo, sftA, B, ldb, B_lo, ldb_lo, incB_lo, sftB, skip_scalA, skip_scalB); break;
    case 18: scaling_launch<backend, T, 18>(stream, op_A, op_B, m, n, k, A, lda, A_lo, lda_lo, incA_lo, sftA, B, ldb, B_lo, ldb_lo, incB_lo, sftB, skip_scalA, skip_scalB); break;
    case 19: scaling_launch<backend, T, 19>(stream, op_A, op_B, m, n, k, A, lda, A_lo, lda_lo, incA_lo, sftA, B, ldb, B_lo, ldb_lo, incB_lo, sftB, skip_scalA, skip_scalB); break;
    case 20: scaling_launch<backend, T, 20>(stream, op_A, op_B, m, n, k, A, lda, A_lo, lda_lo, incA_lo, sftA, B, ldb, B_lo, ldb_lo, incB_lo, sftB, skip_scalA, skip_scalB); break;
    default: break;
    }
}

} // namespace fast
} // namespace real
