#pragma once

template <typename T> __forceinline__ __device__ T reduction_max(T max, T *smem) {
    // inner-warp reduction
    max = inner_warp_max<T>(max);

    // inner-threadblock reduction
    if ((threadIdx.x & 31) == 0) smem[threadIdx.x >> 5] = max; // smem[warp-id] = max in warp
    __syncthreads();

    if (threadIdx.x < 32) {
        if (threadIdx.x < (blockDim.x >> 5)) max = smem[threadIdx.x];
        max = inner_warp_max<T>(max);
        if (threadIdx.x == 0) smem[0] = max;
    }
    __syncthreads();

    return smem[0];
}

//------------------------------
// for accurate mode
//------------------------------

// Column-wise absmax of input matrix
template <typename T> __forceinline__ __device__ underlying_t<T> find_amax(
    const T *const __restrict__ ptr, //
    const unsigned length,           //
    underlying_t<T> *samax           // shared memory (workspace)
) {
    // max in thread
    using U = underlying_t<T>;
    U amax  = Tconst<U>::zero();
    for (unsigned i = threadIdx.x; i < length; i += blockDim.x) {
        const T tmp = Tabs<T>(ptr[i]);
        amax        = Tmax<T>(tmp, amax);
    }
    return reduction_max<U>(amax, samax);
}

// Row-wise absmax of input matrix
template <typename T> __forceinline__ __device__ underlying_t<T> find_amax_tile(
    const unsigned m, const unsigned k,   // size(A)
    const unsigned row_idx,               //
    const T *const __restrict__ A,        // input (lda * k)
    const size_t lda,                     // leading dimension
    underlying_t<T> samax[][TILE_DIM + 1] // shared memory (workspace)
) {
    using U = underlying_t<T>;
    U amax  = Tconst<U>::zero();
    if (row_idx < m) {
        const T *row_ptr = A + row_idx;
        for (unsigned col = threadIdx.y; col < k; col += blockDim.y) {
            const T tmp = Tabs<T>(row_ptr[col * lda]);
            amax        = Tmax<T>(tmp, amax);
        }
    }
    samax[threadIdx.y][threadIdx.x] = amax;
    __syncthreads();

    amax = samax[threadIdx.x][threadIdx.y];
    amax = inner_warp_max<U, TILE_DIM>(amax);
    return amax;
}

// Column-wise max of C_hi (gemmul8::Backend::INT8, real)
__forceinline__ __device__ int32_t find_max(
    const int32_t *const __restrict__ C_hi, // A*B
    const unsigned length,                  // number of rows
    int32_t *samax                          // shared memory (workspace)
) {
    // max in thread
    int32_t amax = 0;
    for (unsigned i = threadIdx.x; i < length; i += blockDim.x) {
        const int32_t tmp = C_hi[i];
        amax              = max(tmp, amax);
    }
    return reduction_max<int32_t>(amax, samax);
}

// Column-wise max of C_hi (gemmul8::Backend::FP8, real)
__forceinline__ __device__ float find_max(
    const unsigned k,
    const float *const __restrict__ C_hi, // A*B
    const unsigned length,                // number of rows
    float *samax                          // shared memory (workspace)
) {
    // max in thread
    const float ku = (k + 1) * 0x1.0000000000000p-24F;
    float amax     = 0.0F;
    for (unsigned i = threadIdx.x; i < length; i += blockDim.x) {
        const float tmp = C_hi[i];
        amax            = max(__fmaf_ru(ku, tmp, tmp), amax);
    }
    return reduction_max<float>(amax, samax);
}

// Column-wise max of C_hi (gemmul8::Backend::INT8, complex)
__forceinline__ __device__ int32_t find_max_complex(
    const int32_t *const __restrict__ C_hi_1, // (Re(A)-Im(A)) * (Re(B)-Im(B))
    const int32_t *const __restrict__ C_hi_2, // Re(A)*Im(B) + Im(A)*Re(B)
    const unsigned length,                    // number of rows
    int32_t *samax                            // shared memory (workspace)
) {
    // max in thread
    int32_t amax = 0;
    for (unsigned i = threadIdx.x; i < length; i += blockDim.x) {
        const int32_t tmp1 = C_hi_1[i];   // (Re(A)-Im(A)) * (Re(B)-Im(B))
        const int32_t tmp2 = C_hi_2[i];   // Re(A)*Im(B) + Im(A)*Re(B)
        const int32_t tmp3 = tmp1 + tmp2; // Re(A)*Re(B) + Im(A)*Im(B)
        amax               = max(max(tmp3, tmp2), amax);
    }
    return reduction_max<int32_t>(amax, samax);
}

// Column-wise max of C_hi (gemmul8::Backend::FP8, complex)
__forceinline__ __device__ float find_max_complex(
    const unsigned k,                       // inner dimension
    const float *const __restrict__ C_hi_1, // (Re(A)-Im(A)) * (Re(B)-Im(B))
    const float *const __restrict__ C_hi_2, // Re(A)*Im(B)
    const float *const __restrict__ C_hi_3, // Im(A)*Re(B)
    const unsigned length,                  // number of rows
    float *samax                            // shared memory (workspace)
) {
    // max in thread
    const float ku = (k + 1) * 0x1.0000000000000p-24F;
    float amax     = 0.0F;
    for (unsigned i = threadIdx.x; i < length; i += blockDim.x) {
        const float ArBi           = C_hi_2[i];                            // Re(A)*Im(B)
        const float ArBi_up        = __fmaf_ru(ku, ArBi, ArBi);            // upper bound of Re(A)*Im(B)
        const float AiBr           = C_hi_3[i];                            // Im(A)*Re(B)
        const float AiBr_up        = __fmaf_ru(ku, AiBr, AiBr);            // upper bound of Im(A)*Re(B)
        const float ArBi_plus_AiBr = __fadd_ru(ArBi_up, AiBr_up);          // upper bound of Re(A)*Im(B) + Im(A)*Re(B)
        const float AriBri         = C_hi_1[i];                            // (Re(A)-Im(A)) * (Re(B)-Im(B))
        const float AriBri_up      = __fmaf_ru(ku, AriBri, AriBri);        // upper bound of (Re(A)-Im(A)) * (Re(B)-Im(B))
        const float ArBr_plus_AiBi = __fadd_ru(AriBri_up, ArBi_plus_AiBr); // upper bound of Re(A)*Re(B) + Im(A)*Im(B)
        amax                       = max(max(ArBr_plus_AiBi, ArBi_plus_AiBr), amax);
    }
    return reduction_max<float>(amax, samax);
}

// Row-wise max of C_hi (gemmul8::Backend::INT8, real)
__forceinline__ __device__ int32_t find_max_tile(
    const unsigned m, const unsigned n,     // size(C)
    const unsigned row_idx,                 //
    const int32_t *const __restrict__ C_hi, // A*B
    const size_t ldc,                       // leading dimension
    int32_t samax[][TILE_DIM + 1]           // shared memory (workspace)
) {
    int32_t amax = 0;
    if (row_idx < m) {
        for (unsigned col = threadIdx.y; col < n; col += blockDim.y) {
            const int32_t tmp = C_hi[col * ldc + row_idx];
            amax              = max(tmp, amax);
        }
    }
    samax[threadIdx.y][threadIdx.x] = amax;
    __syncthreads();

    amax = samax[threadIdx.x][threadIdx.y];
    amax = inner_warp_max<int32_t, TILE_DIM>(amax);
    return amax;
}

// Row-wise max of C_hi (gemmul8::Backend::FP8, real)
__forceinline__ __device__ float find_max_tile(
    const unsigned m, const unsigned n,   // size(C)
    const unsigned k,                     // inner dimension
    const unsigned row_idx,               //
    const float *const __restrict__ C_hi, // A*B
    const size_t ldc,                     // leading dimension
    float samax[][TILE_DIM + 1]           // shared memory (workspace)
) {
    const float ku = (k + 1) * 0x1.0000000000000p-24F;
    float amax     = 0.0F;
    if (row_idx < m) {
        for (unsigned col = threadIdx.y; col < n; col += blockDim.y) {
            const float tmp = C_hi[col * ldc + row_idx];
            amax            = max(__fmaf_ru(ku, tmp, tmp), amax);
        }
    }
    samax[threadIdx.y][threadIdx.x] = amax;
    __syncthreads();

    amax = samax[threadIdx.x][threadIdx.y];
    amax = inner_warp_max<float, TILE_DIM>(amax);
    return amax;
}

// Row-wise max of C_hi (gemmul8::Backend::INT8, complex)
__forceinline__ __device__ int32_t find_max_tile_complex(
    const unsigned m, const unsigned n,       // size(C)
    const unsigned row_idx,                   //
    const int32_t *const __restrict__ C_hi_1, // (Re(A)-Im(A)) * (Re(B)-Im(B))
    const int32_t *const __restrict__ C_hi_2, // Re(A)*Im(B) + Im(A)*Re(B)
    const size_t lda,                         // leading dimension
    int32_t samax[][TILE_DIM + 1]             // shared memory (workspace)
) {
    int32_t amax = 0;
    if (row_idx < m) {
        for (unsigned col = threadIdx.y; col < n; col += blockDim.y) {
            const int32_t tmp1 = C_hi_1[col * lda + row_idx];
            const int32_t tmp2 = C_hi_2[col * lda + row_idx]; // Re(A)*Im(B) + Im(A)*Re(B)^
            const int32_t tmp3 = tmp1 + tmp2;                 // Re(A)*Re(B) + Im(A)*Im(B)
            amax               = max(max(tmp3, tmp2), amax);
        }
    }
    samax[threadIdx.y][threadIdx.x] = amax;
    __syncthreads();

    amax = samax[threadIdx.x][threadIdx.y];
    amax = inner_warp_max<int32_t, TILE_DIM>(amax);

    return amax;
}

// Row-wise max of C_hi (gemmul8::Backend::FP8, complex)
__forceinline__ __device__ float find_max_tile_complex(
    const unsigned m, const unsigned n,     // size(C)
    const unsigned k,                       // inner dimension
    const unsigned row_idx,                 //
    const float *const __restrict__ C_hi_1, // (Re(A)-Im(A)) * (Re(B)-Im(B))
    const float *const __restrict__ C_hi_2, // Re(A)*Im(B)
    const float *const __restrict__ C_hi_3, // Im(A)*Re(B)
    const size_t lda,                       // leading dimension
    float samax[][TILE_DIM + 1]             // shared memory (workspace)
) {
    const float ku = (k + 1) * 0x1.0000000000000p-24F;
    float amax     = 0.0F;
    if (row_idx < m) {
        for (unsigned col = threadIdx.y; col < n; col += blockDim.y) {
            const size_t i             = col * lda + row_idx;
            const float ArBi           = C_hi_2[i];                            // Re(A)*Im(B)
            const float ArBi_up        = __fmaf_ru(ku, ArBi, ArBi);            // upper bound of Re(A)*Im(B)
            const float AiBr           = C_hi_3[i];                            // Im(A)*Re(B)
            const float AiBr_up        = __fmaf_ru(ku, AiBr, AiBr);            // upper bound of Im(A)*Re(B)
            const float ArBi_plus_AiBr = __fadd_ru(ArBi_up, AiBr_up);          // upper bound of Re(A)*Im(B) + Im(A)*Re(B)
            const float AriBri         = C_hi_1[i];                            // (Re(A)-Im(A)) * (Re(B)-Im(B))
            const float AriBri_up      = __fmaf_ru(ku, AriBri, AriBri);        // upper bound of (Re(A)-Im(A)) * (Re(B)-Im(B))
            const float ArBr_plus_AiBi = __fadd_ru(AriBri_up, ArBi_plus_AiBr); // upper bound of Re(A)*Re(B) + Im(A)*Im(B)
            amax                       = max(max(ArBr_plus_AiBi, ArBi_plus_AiBr), amax);
        }
    }
    samax[threadIdx.y][threadIdx.x] = amax;
    __syncthreads();

    amax = samax[threadIdx.x][threadIdx.y];
    amax = inner_warp_max<float, TILE_DIM>(amax);

    return amax;
}

//------------------------------
// for fast mode
//------------------------------

// Column-wise amax & sum of squares of input matrix
template <typename T> __forceinline__ __device__ underlying_t<T> find_amax_and_nrm(
    const T *const __restrict__ ptr, //
    const unsigned length,           //
    underlying_t<T> *shm,            // shared memory (workspace)
    underlying_t<T> &vecnrm          // 2-norm^2
) {
    using U = underlying_t<T>;

    U *samax = shm;
    U *ssum  = shm + 32;

    // max in thread
    U amax = Tconst<U>::zero();
    U sum  = Tconst<U>::zero();
    for (unsigned i = threadIdx.x; i < length; i += blockDim.x) {
        T tmp = Tabs<T>(ptr[i]);
        amax  = Tmax<T>(tmp, amax);
        sum   = Tsqr_add_ru<T>(tmp, sum); // round-up mode
    }

    // inner-warp reduction
    amax = inner_warp_max(amax);
    sum  = inner_warp_sum(sum);

    // inner-threadblock reduction
    if ((threadIdx.x & 31) == 0) {
        samax[threadIdx.x >> 5] = amax; // samax[warp-id] = max in warp
        ssum[threadIdx.x >> 5]  = sum;  // ssum[warp-id] = sum in warp
    }
    __syncthreads();

    sum = Tconst<U>::zero();
    if (threadIdx.x < 32) {
        if (threadIdx.x < (blockDim.x >> 5)) {
            amax = samax[threadIdx.x];
            sum  = ssum[threadIdx.x];
        }
        amax = inner_warp_max(amax);
        sum  = inner_warp_sum(sum);
        if (threadIdx.x == 0) {
            samax[0] = amax;
            ssum[0]  = sum;
        }
    }
    __syncthreads();

    vecnrm = ssum[0];
    return samax[0];
}

// Row-wise amax & sum of squares of input matrix
template <typename T> __forceinline__ __device__ underlying_t<T> find_amax_and_nrm_tile(
    const unsigned m, const unsigned k,    // size(A)
    const unsigned row_idx,                //
    const T *const __restrict__ A,         // input (lda * k)
    const size_t lda,                      // leading dimension
    underlying_t<T> samax[][TILE_DIM + 1], // shared memory (workspace)
    underlying_t<T> ssum[][TILE_DIM + 1],  // shared memory (workspace)
    underlying_t<T> &vecnrm                // 2-norm^2
) {
    using U = underlying_t<T>;

    U amax = Tconst<U>::zero();
    U sum  = Tconst<U>::zero();
    if (row_idx < m) {
        const T *row_ptr = A + row_idx;
        for (unsigned col = threadIdx.y; col < k; col += blockDim.y) {
            const T tmp = Tabs<T>(row_ptr[col * lda]);
            amax        = Tmax<T>(tmp, amax);
            sum         = Tsqr_add_ru<T>(tmp, sum); // round-up mode
        }
    }
    samax[threadIdx.y][threadIdx.x] = amax;
    ssum[threadIdx.y][threadIdx.x]  = sum;
    __syncthreads();

    sum    = ssum[threadIdx.x][threadIdx.y];
    vecnrm = inner_warp_sum<U, TILE_DIM>(sum);

    amax = samax[threadIdx.x][threadIdx.y];
    amax = inner_warp_max<U, TILE_DIM>(amax);

    return amax;
}
