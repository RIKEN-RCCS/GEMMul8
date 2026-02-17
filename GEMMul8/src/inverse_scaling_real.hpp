#pragma once

namespace real {

//------------------------------
// accumulation for INT8
//------------------------------
template <int num_moduli, int IDX = 0>
__forceinline__ __device__ void accumulator_double(double &C64f, const int8_t *const __restrict__ Ctmp, const size_t incCtmp) {
    if constexpr (IDX < num_moduli) {
        const double C8i_d = double(Ctmp[IDX * incCtmp]);
        const double qPi   = table::INT8::qPi_double_v<num_moduli, IDX>;
        C64f               = fma(qPi, C8i_d, C64f);
        accumulator_double<num_moduli, IDX + 1>(C64f, Ctmp, incCtmp);
    }
}

template <int num_moduli, int IDX = 0>
__forceinline__ __device__ void accumulator_double2(double2 &C64f, const int8_t *const __restrict__ Ctmp, const size_t incCtmp) {
    if constexpr (IDX < num_moduli) {
        const double C8i_d = double(Ctmp[IDX * incCtmp]);
        const double2 qPi  = table::INT8::qPi_double2_v<num_moduli, IDX>;
        C64f.x             = fma(qPi.x, C8i_d, C64f.x); // error-free
        C64f.y             = fma(qPi.y, C8i_d, C64f.y); // non-error-free
        accumulator_double2<num_moduli, IDX + 1>(C64f, Ctmp, incCtmp);
    }
}

//------------------------------
// accumulation for FP8
//------------------------------
template <int num_moduli, int IDX = 0>
__forceinline__ __device__ void accumulator_double(double &C64f, const int16_t *const __restrict__ Ctmp, const size_t incCtmp) {
    if constexpr (IDX < num_moduli) {
        const double C16i_d = double(Ctmp[IDX * incCtmp]);
        const double qPi    = table::FP8::qPi_double_v<num_moduli, IDX>;
        C64f                = fma(qPi, C16i_d, C64f);
        accumulator_double<num_moduli, IDX + 1>(C64f, Ctmp, incCtmp);
    }
}

template <int num_moduli, int IDX = 0>
__forceinline__ __device__ void accumulator_double2(double2 &C64f, const int16_t *const __restrict__ Ctmp, const size_t incCtmp) {
    if constexpr (IDX < num_moduli) {
        const double C16i_d = double(Ctmp[IDX * incCtmp]);
        const double2 qPi   = table::FP8::qPi_double2_v<num_moduli, IDX>;
        C64f.x              = fma(qPi.x, C16i_d, C64f.x); // error-free
        C64f.y              = fma(qPi.y, C16i_d, C64f.y); // non-error-free
        accumulator_double2<num_moduli, IDX + 1>(C64f, Ctmp, incCtmp);
    }
}

//------------------------------
// final reduction & Undo scaling
//------------------------------
template <int num_moduli, typename TC, typename TP, typename TCtmp>
__forceinline__ __device__ TC invscal_device(
    const size_t incCtmp,                 // increment
    const TCtmp *const __restrict__ Ctmp, // input
    const TP P,                           // prod(moduli)
    const double invP,                    // 1/prod(moduli)
    const int sft                         // exponent of shift values
) {
    if constexpr (std::is_same_v<TP, double>) {

        // sum(qi*Pi*Ctmp[i])
        double C64f = 0;
        accumulator_double<num_moduli>(C64f, Ctmp, incCtmp);

        const double quot    = rint(invP * C64f);                          // round(C64f/P)
        const double CRT_ans = fma(P, quot, C64f);                         // C64f - P*round(C64f/P)
        const TC C           = Tscalbn<TC>(static_cast<TC>(CRT_ans), sft); // undo scaling
        return C;

    } else if constexpr (std::is_same_v<TP, double2>) {

        // sum(qi*Pi*Ctmp[i])
        double2 C64f{};
        accumulator_double2<num_moduli>(C64f, Ctmp, incCtmp);

        const double quot    = rint(invP * C64f.x);                             // round(C64f/P)
        const double CRT_ans = fma(P.y, quot, fma(P.x, quot, C64f.x) + C64f.y); // C64f - P*round(C64f/P)
        const TC C           = Tscalbn<TC>(static_cast<TC>(CRT_ans), sft);      // undo scaling
        return C;

    } else {
        return Tconst<TC>::zero();
    }
}

//------------------------------
// C := alpha*AB + beta*C
//------------------------------
template <int num_moduli, typename TC, typename TP, typename TCtmp>
__global__ void invscal_kernel_general(
    const TC alpha, const TC beta,          //
    const size_t m,                         // size(C64f,1)
    const size_t sizeC,                     // m*n
    const size_t incCtmp,                   // increment
    const TCtmp *const __restrict__ Ctmp,   // input
    const size_t ldctmp,                    // leading dim of Ctmp
    TC *const __restrict__ C,               // output
    const size_t ldc,                       // leading dimension
    const TP P,                             // -prod(moduli)
    const double invP,                      // 1/prod(moduli)
    const int16_t *const __restrict__ sftA, // exponent of shift values for rows of A
    const int16_t *const __restrict__ sftB  // exponent of shift values for cols of B
) {
    const size_t idx = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (idx >= sizeC) return;
    const auto col     = idx / m;
    const auto row     = idx - col * m;
    const auto mem_idx = col * ldctmp + row;
    const TC AB        = invscal_device<num_moduli, TC, TP, TCtmp>(incCtmp, Ctmp + mem_idx, P, invP, sftA[row] + sftB[col]);

    const auto idxC = col * ldc + row;
    C[idxC]         = Taxpby_scal<TC>(alpha, AB, beta, C[idxC]);
}

template <int num_moduli, typename TC, typename TP, typename TCtmp>
__global__ void invscal_kernel_general_deviceScalar(
    const TC *alpha, const TC *beta,        //
    const size_t m,                         // size(C64f,1)
    const size_t sizeC,                     // m*n
    const size_t incCtmp,                   // increment
    const TCtmp *const __restrict__ Ctmp,   // input
    const size_t ldctmp,                    // leading dim of Ctmp
    TC *const __restrict__ C,               // output
    const size_t ldc,                       // leading dimension
    const TP P,                             // -prod(moduli)
    const double invP,                      // 1/prod(moduli)
    const int16_t *const __restrict__ sftA, // exponent of shift values for rows of A
    const int16_t *const __restrict__ sftB  // exponent of shift values for cols of B
) {
    const size_t idx = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (idx >= sizeC) return;
    const auto col     = idx / m;
    const auto row     = idx - col * m;
    const auto mem_idx = col * ldctmp + row;
    const TC AB        = invscal_device<num_moduli, TC, TP, TCtmp>(incCtmp, Ctmp + mem_idx, P, invP, sftA[row] + sftB[col]);

    const auto idxC = col * ldc + row;
    C[idxC]         = Taxpby_scal<TC>(*alpha, AB, *beta, C[idxC]);
}

//------------------------------
// C := alpha*AB + beta*C
//------------------------------
template <int num_moduli, typename TC, typename TP, typename TCtmp, int ALPHA, int BETA>
__global__ void invscal_kernel_special(
    const size_t m,                         // size(C64f,1)
    const size_t sizeC,                     // m*n
    const size_t incCtmp,                   // increment
    const TCtmp *const __restrict__ Ctmp,   // input
    const size_t ldctmp,                    // leading dim of Ctmp
    TC *const __restrict__ C,               // output
    const size_t ldc,                       // leading dimension
    const TP P,                             // -prod(moduli)
    const double invP,                      // 1/prod(moduli)
    const int16_t *const __restrict__ sftA, // exponent of shift values for rows of A
    const int16_t *const __restrict__ sftB  // exponent of shift values for cols of B
) {
    const size_t idx = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (idx >= sizeC) return;
    const auto col     = idx / m;
    const auto row     = idx - col * m;
    const auto mem_idx = col * ldctmp + row;
    const TC AB        = invscal_device<num_moduli, TC, TP, TCtmp>(incCtmp, Ctmp + mem_idx, P, invP, sftA[row] + sftB[col]);

    const auto idxC = col * ldc + row;
    if constexpr (ALPHA == 1 && BETA == 0) {
        // C[idxC] =  1*AB + 0*C[idxC] = AB
        C[idxC] = AB;

    } else if constexpr (ALPHA == 1 && BETA == 1) {
        // C[idxC] =  1*AB + 1*C[idxC] = C[idxC] + AB
        C[idxC] += AB;

    } else if constexpr (ALPHA == -1 && BETA == 0) {
        // C[idxC] = -1*AB + 0*C[idxC] = -AB
        C[idxC] = -AB;

    } else if constexpr (ALPHA == -1 && BETA == 1) {
        // C[idxC] = -1*AB + 1*C[idxC] = C[idxC] - AB
        C[idxC] -= AB;
    }
}

//------------------------------
// Launcher!!!
//------------------------------
template <gemmul8::Backend backend, typename T, int num_moduli>
__forceinline__ void inverse_scaling_launch(
    const cudaStream_t &stream,       //
    const size_t m, const size_t n,   // size(C)
    const mid_t<backend> *const Ctmp, // input
    const size_t ldctmp,              // leading dim of Ctmp
    const size_t incCtmp,             // increment
    T *const C,                       // output
    const size_t ldc,                 // leading dimension
    const int16_t *const sftA,        // exponent of shift values for rows of A
    const int16_t *const sftB,        // exponent of shift values for cols of B
    const T *alpha, const T *beta     //
) {
    const size_t sizeC        = m * n;
    const size_t grid_invscal = (sizeC + threads_invscal - 1) / threads_invscal;
    using TP                  = std::conditional_t<(std::is_same_v<underlying_t<T>, float> || (num_moduli <= threshold<backend>::P_is_double)), double, double2>;
    const double invP         = table::get_invP<backend>(num_moduli);
    const TP P                = table::get_P<backend, TP>(num_moduli);

    cudaPointerAttributes attr{};
    cudaPointerGetAttributes(&attr, alpha);
    const bool is_device = (attr.type != cudaMemoryTypeUnregistered) && (attr.type != cudaMemoryTypeHost);

    if (is_device) {
        invscal_kernel_general_deviceScalar<num_moduli, T, TP, mid_t<backend>><<<grid_invscal, threads_invscal, 0, stream>>>(alpha, beta, m, sizeC, incCtmp, Ctmp, ldctmp, C, ldc, P, invP, sftA, sftB);
    } else {
        if (*alpha == Tconst<T>::one()) {
            if (*beta == Tconst<T>::zero()) {
                invscal_kernel_special<num_moduli, T, TP, mid_t<backend>, 1, 0><<<grid_invscal, threads_invscal, 0, stream>>>(m, sizeC, incCtmp, Ctmp, ldctmp, C, ldc, P, invP, sftA, sftB);
                return;
            } else if (*beta == Tconst<T>::one()) {
                invscal_kernel_special<num_moduli, T, TP, mid_t<backend>, 1, 1><<<grid_invscal, threads_invscal, 0, stream>>>(m, sizeC, incCtmp, Ctmp, ldctmp, C, ldc, P, invP, sftA, sftB);
                return;
            }
        } else if (*alpha == Tconst<T>::mone()) {
            if (*beta == Tconst<T>::zero()) {
                invscal_kernel_special<num_moduli, T, TP, mid_t<backend>, -1, 0><<<grid_invscal, threads_invscal, 0, stream>>>(m, sizeC, incCtmp, Ctmp, ldctmp, C, ldc, P, invP, sftA, sftB);
                return;
            } else if (*beta == Tconst<T>::one()) {
                invscal_kernel_special<num_moduli, T, TP, mid_t<backend>, -1, 1><<<grid_invscal, threads_invscal, 0, stream>>>(m, sizeC, incCtmp, Ctmp, ldctmp, C, ldc, P, invP, sftA, sftB);
                return;
            }
        }
        invscal_kernel_general<num_moduli, T, TP, mid_t<backend>><<<grid_invscal, threads_invscal, 0, stream>>>(*alpha, *beta, m, sizeC, incCtmp, Ctmp, ldctmp, C, ldc, P, invP, sftA, sftB);
    }
}

//------------------------------
// Interface!!
//------------------------------
template <gemmul8::Backend backend, typename T>
__inline__ void inverse_scaling(
    const cudaStream_t &stream,       //
    const unsigned num_moduli,        // number of moduli
    const size_t m, const size_t n,   // size(C)
    const mid_t<backend> *const Ctmp, // input
    const size_t ldctmp,              // leading dim of Ctmp
    const size_t incCtmp,             // increment
    T *const C,                       // output
    const size_t ldc,                 // leading dimension
    const int16_t *const sftA,        // exponent of shift values for rows of A
    const int16_t *const sftB,        // exponent of shift values for cols of B
    const T *alpha, const T *beta     //
) {
    switch (num_moduli) {
    case 2: inverse_scaling_launch<backend, T, 2>(stream, m, n, Ctmp, ldctmp, incCtmp, C, ldc, sftA, sftB, alpha, beta); break;
    case 3: inverse_scaling_launch<backend, T, 3>(stream, m, n, Ctmp, ldctmp, incCtmp, C, ldc, sftA, sftB, alpha, beta); break;
    case 4: inverse_scaling_launch<backend, T, 4>(stream, m, n, Ctmp, ldctmp, incCtmp, C, ldc, sftA, sftB, alpha, beta); break;
    case 5: inverse_scaling_launch<backend, T, 5>(stream, m, n, Ctmp, ldctmp, incCtmp, C, ldc, sftA, sftB, alpha, beta); break;
    case 6: inverse_scaling_launch<backend, T, 6>(stream, m, n, Ctmp, ldctmp, incCtmp, C, ldc, sftA, sftB, alpha, beta); break;
    case 7: inverse_scaling_launch<backend, T, 7>(stream, m, n, Ctmp, ldctmp, incCtmp, C, ldc, sftA, sftB, alpha, beta); break;
    case 8: inverse_scaling_launch<backend, T, 8>(stream, m, n, Ctmp, ldctmp, incCtmp, C, ldc, sftA, sftB, alpha, beta); break;
    case 9: inverse_scaling_launch<backend, T, 9>(stream, m, n, Ctmp, ldctmp, incCtmp, C, ldc, sftA, sftB, alpha, beta); break;
    case 10: inverse_scaling_launch<backend, T, 10>(stream, m, n, Ctmp, ldctmp, incCtmp, C, ldc, sftA, sftB, alpha, beta); break;
    case 11: inverse_scaling_launch<backend, T, 11>(stream, m, n, Ctmp, ldctmp, incCtmp, C, ldc, sftA, sftB, alpha, beta); break;
    case 12: inverse_scaling_launch<backend, T, 12>(stream, m, n, Ctmp, ldctmp, incCtmp, C, ldc, sftA, sftB, alpha, beta); break;
    case 13: inverse_scaling_launch<backend, T, 13>(stream, m, n, Ctmp, ldctmp, incCtmp, C, ldc, sftA, sftB, alpha, beta); break;
    case 14: inverse_scaling_launch<backend, T, 14>(stream, m, n, Ctmp, ldctmp, incCtmp, C, ldc, sftA, sftB, alpha, beta); break;
    case 15: inverse_scaling_launch<backend, T, 15>(stream, m, n, Ctmp, ldctmp, incCtmp, C, ldc, sftA, sftB, alpha, beta); break;
    case 16: inverse_scaling_launch<backend, T, 16>(stream, m, n, Ctmp, ldctmp, incCtmp, C, ldc, sftA, sftB, alpha, beta); break;
    case 17: inverse_scaling_launch<backend, T, 17>(stream, m, n, Ctmp, ldctmp, incCtmp, C, ldc, sftA, sftB, alpha, beta); break;
    case 18: inverse_scaling_launch<backend, T, 18>(stream, m, n, Ctmp, ldctmp, incCtmp, C, ldc, sftA, sftB, alpha, beta); break;
    case 19: inverse_scaling_launch<backend, T, 19>(stream, m, n, Ctmp, ldctmp, incCtmp, C, ldc, sftA, sftB, alpha, beta); break;
    case 20: inverse_scaling_launch<backend, T, 20>(stream, m, n, Ctmp, ldctmp, incCtmp, C, ldc, sftA, sftB, alpha, beta); break;
    default: break;
    }
}

} // namespace real
