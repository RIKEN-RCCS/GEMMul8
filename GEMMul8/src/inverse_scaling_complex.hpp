#pragma once

namespace complex {

//------------------------------
// accumulation for INT8
//------------------------------
template <int num_moduli, int IDX = 0>
__forceinline__ __device__ void accumulator_double(double2 &C64f, const char2 *const __restrict__ Ctmp, const size_t incCtmp) {
    if constexpr (IDX < num_moduli) {
        const char2 C8i_tmp = Ctmp[IDX * incCtmp];
        const double2 C8i_d = {double(C8i_tmp.x), double(C8i_tmp.y)};
        const double qPi    = table::INT8::qPi_double_v<num_moduli, IDX>;
        C64f.x              = fma(qPi, C8i_d.x, C64f.x);
        C64f.y              = fma(qPi, C8i_d.y, C64f.y);
        accumulator_double<num_moduli, IDX + 1>(C64f, Ctmp, incCtmp);
    }
}

template <int num_moduli, int IDX = 0>
__forceinline__ __device__ void accumulator_double2(double2 &C64f_hi, double2 &C64f_lo, const char2 *const __restrict__ Ctmp, const size_t incCtmp) {
    if constexpr (IDX < num_moduli) {
        const char2 C8i_tmp = Ctmp[IDX * incCtmp];
        const double2 C8i_d = {double(C8i_tmp.x), double(C8i_tmp.y)};
        const double2 qPi   = table::INT8::qPi_double2_v<num_moduli, IDX>;
        C64f_hi.x           = fma(qPi.x, C8i_d.x, C64f_hi.x); // error-free
        C64f_hi.y           = fma(qPi.x, C8i_d.y, C64f_hi.y); // error-free
        C64f_lo.x           = fma(qPi.y, C8i_d.x, C64f_lo.x); // non-error-free
        C64f_lo.y           = fma(qPi.y, C8i_d.y, C64f_lo.y); // non-error-free
        accumulator_double2<num_moduli, IDX + 1>(C64f_hi, C64f_lo, Ctmp, incCtmp);
    }
}

//------------------------------
// accumulation for FP8
//------------------------------
template <int num_moduli, int IDX = 0>
__forceinline__ __device__ void accumulator_double(double2 &C64f, const short2 *const __restrict__ Ctmp, const size_t incCtmp) {
    if constexpr (IDX < num_moduli) {
        const short2 C16i_tmp = Ctmp[IDX * incCtmp];
        const double2 C16i_d  = {double(C16i_tmp.x), double(C16i_tmp.y)};
        const double qPi      = table::FP8::qPi_double_v<num_moduli, IDX>;
        C64f.x                = fma(qPi, C16i_d.x, C64f.x);
        C64f.y                = fma(qPi, C16i_d.y, C64f.y);
        accumulator_double<num_moduli, IDX + 1>(C64f, Ctmp, incCtmp);
    }
}

template <int num_moduli, int IDX = 0>
__forceinline__ __device__ void accumulator_double2(double2 &C64f_hi, double2 &C64f_lo, const short2 *const __restrict__ Ctmp, const size_t incCtmp) {
    if constexpr (IDX < num_moduli) {
        const short2 C16i_tmp = Ctmp[IDX * incCtmp];
        const double2 C16i_d  = {double(C16i_tmp.x), double(C16i_tmp.y)};
        const double2 qPi     = table::FP8::qPi_double2_v<num_moduli, IDX>;
        C64f_hi.x             = fma(qPi.x, C16i_d.x, C64f_hi.x); // error-free
        C64f_hi.y             = fma(qPi.x, C16i_d.y, C64f_hi.y); // error-free
        C64f_lo.x             = fma(qPi.y, C16i_d.x, C64f_lo.x); // non-error-free
        C64f_lo.y             = fma(qPi.y, C16i_d.y, C64f_lo.y); // non-error-free
        accumulator_double2<num_moduli, IDX + 1>(C64f_hi, C64f_lo, Ctmp, incCtmp);
    }
}

//------------------------------
// final reduction & Undo scaling
//------------------------------
template <int num_moduli, typename TC, typename TP, typename TCtmp>
__forceinline__ __device__ TC invscal_device(
    const size_t incC_lo,                 // increment
    const TCtmp *const __restrict__ C_lo, // input
    const TP P,                           // prod(moduli)
    const double invP,                    // 1/prod(moduli)
    const int sft                         // exponent of shift values
) {
    if constexpr (std::is_same_v<TP, double>) {

        // sum(qi*Pi*C_lo[i])
        double2 C64f{};
        accumulator_double<num_moduli>(C64f, C_lo, incC_lo);

        // round(C64f/P)
        double2 quot;
        quot.x = rint(invP * C64f.x);
        quot.y = rint(invP * C64f.y);

        // C64f - P*round(C64f/P)
        double2 CRT_ans;
        CRT_ans.x = fma(P, quot.x, C64f.x);
        CRT_ans.y = fma(P, quot.y, C64f.y);

        // undo scaling
        using U = underlying_t<TC>;
        TC C;
        C.x = Tscalbn<U>(static_cast<U>(CRT_ans.x), sft);
        C.y = Tscalbn<U>(static_cast<U>(CRT_ans.y), sft);

        return C;

    } else if constexpr (std::is_same_v<TP, double2>) {

        // sum(qi*Pi*C_lo[i])
        double2 C64f_hi{}, C64f_lo{};
        accumulator_double2<num_moduli>(C64f_hi, C64f_lo, C_lo, incC_lo);

        // round(C64f/P)
        double2 quot;
        quot.x = rint(invP * C64f_hi.x);
        quot.y = rint(invP * C64f_hi.y);

        // C64f - P*round(C64f/P)
        double2 CRT_ans;
        CRT_ans.x = fma(P.y, quot.x, fma(P.x, quot.x, C64f_hi.x) + C64f_lo.x);
        CRT_ans.y = fma(P.y, quot.y, fma(P.x, quot.y, C64f_hi.y) + C64f_lo.y);

        // undo scaling
        using U = underlying_t<TC>;
        TC C;
        C.x = Tscalbn<U>(static_cast<U>(CRT_ans.x), sft);
        C.y = Tscalbn<U>(static_cast<U>(CRT_ans.y), sft);

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
    TC AB              = invscal_device<num_moduli, TC, TP, TCtmp>(incCtmp, Ctmp + mem_idx, P, invP, sftA[row] + sftB[col]);

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
    TC AB              = invscal_device<num_moduli, TC, TP, TCtmp>(incCtmp, Ctmp + mem_idx, P, invP, sftA[row] + sftB[col]);

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
    const TCtmp *const __restrict__ C_lo,   // input
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
    TC AB              = invscal_device<num_moduli, TC, TP, TCtmp>(incCtmp, C_lo + mem_idx, P, invP, sftA[row] + sftB[col]);

    const auto idxC = col * ldc + row;
    if constexpr (ALPHA == 1 && BETA == 0) {
        // C[idxC] =  1*AB + 0*C[idxC] = AB
        C[idxC] = AB;

    } else if constexpr (ALPHA == 1 && BETA == 1) {
        // C[idxC] =  1*AB + 1*C[idxC] = C[idxC] + AB
        TC Ctmp = C[idxC];
        Ctmp.x += AB.x;
        Ctmp.y += AB.y;
        C[idxC] = Ctmp;

    } else if constexpr (ALPHA == -1 && BETA == 0) {
        // C[idxC] = -1*AB + 0*C[idxC] = -AB
        AB.x    = -AB.x;
        AB.y    = -AB.y;
        C[idxC] = AB;

    } else if constexpr (ALPHA == -1 && BETA == 1) {
        // C[idxC] = -1*AB + 1*C[idxC] = C[idxC] - AB
        TC Ctmp = C[idxC];
        Ctmp.x -= AB.x;
        Ctmp.y -= AB.y;
        C[idxC] = Ctmp;
    }
}

//------------------------------
// Launcher!!!
//------------------------------
template <gemmul8::Backend backend, typename T, int num_moduli>
__inline__ void inverse_scaling_launch(
    const cudaStream_t &stream,         //
    const size_t m, const size_t n,     // size(C)
    const midx2_t<backend> *const Ctmp, // input
    const size_t ldctmp,                // leading dim of Ctmp
    const size_t incCtmp,               // increment
    T *const C,                         // output
    const size_t ldc,                   // leading dimension
    const int16_t *const sftA,          // exponent of shift values for rows of A
    const int16_t *const sftB,          // exponent of shift values for cols of B
    const T *alpha, const T *beta       //
) {
    using U                   = underlying_t<T>;
    const size_t sizeC        = m * n;
    const size_t grid_invscal = (sizeC + threads_invscal - 1) / threads_invscal;
    using TP                  = std::conditional_t<(std::is_same_v<underlying_t<T>, float> || (num_moduli <= threshold<backend>::P_is_double)), double, double2>;
    const double invP         = table::get_invP<backend>(num_moduli);
    const TP P                = table::get_P<backend, TP>(num_moduli);

    cudaPointerAttributes attr{};
    cudaPointerGetAttributes(&attr, alpha);
    const bool is_device = (attr.type != cudaMemoryTypeUnregistered) && (attr.type != cudaMemoryTypeHost);

    if (is_device) {
        invscal_kernel_general_deviceScalar<num_moduli, T, TP, midx2_t<backend>><<<grid_invscal, threads_invscal, 0, stream>>>(alpha, beta, m, sizeC, incCtmp, Ctmp, ldctmp, C, ldc, P, invP, sftA, sftB);
    } else {
        const T alpha_h = *alpha;
        const T beta_h  = *beta;
        if (alpha_h.x == Tconst<U>::one() && alpha_h.y == Tconst<U>::zero()) {
            if (beta_h.x == Tconst<U>::zero() && beta_h.y == Tconst<U>::zero()) {
                invscal_kernel_special<num_moduli, T, TP, midx2_t<backend>, 1, 0><<<grid_invscal, threads_invscal, 0, stream>>>(m, sizeC, incCtmp, Ctmp, ldctmp, C, ldc, P, invP, sftA, sftB);
                return;
            } else if (beta_h.x == Tconst<U>::one() && beta_h.y == Tconst<U>::zero()) {
                invscal_kernel_special<num_moduli, T, TP, midx2_t<backend>, 1, 1><<<grid_invscal, threads_invscal, 0, stream>>>(m, sizeC, incCtmp, Ctmp, ldctmp, C, ldc, P, invP, sftA, sftB);
                return;
            }
        } else if (alpha_h.x == Tconst<U>::mone() && alpha_h.y == Tconst<U>::zero()) {
            if (beta_h.x == Tconst<U>::zero() && beta_h.y == Tconst<U>::zero()) {
                invscal_kernel_special<num_moduli, T, TP, midx2_t<backend>, -1, 0><<<grid_invscal, threads_invscal, 0, stream>>>(m, sizeC, incCtmp, Ctmp, ldctmp, C, ldc, P, invP, sftA, sftB);
                return;
            } else if (beta_h.x == Tconst<U>::one() && beta_h.y == Tconst<U>::zero()) {
                invscal_kernel_special<num_moduli, T, TP, midx2_t<backend>, -1, 1><<<grid_invscal, threads_invscal, 0, stream>>>(m, sizeC, incCtmp, Ctmp, ldctmp, C, ldc, P, invP, sftA, sftB);
                return;
            }
        }
        invscal_kernel_general<num_moduli, T, TP, midx2_t<backend>><<<grid_invscal, threads_invscal, 0, stream>>>(alpha_h, beta_h, m, sizeC, incCtmp, Ctmp, ldctmp, C, ldc, P, invP, sftA, sftB);
    }
}

//------------------------------
// Interface!!
//------------------------------
template <gemmul8::Backend backend, typename T>
__inline__ void inverse_scaling(
    const cudaStream_t &stream,         //
    const unsigned num_moduli,          // number of moduli
    const size_t m, const size_t n,     // size(C)
    const midx2_t<backend> *const Ctmp, // input
    const size_t ldctmp,                // leading dim of Ctmp
    const size_t incCtmp,               // increment
    T *const C,                         // output
    const size_t ldc,                   // leading dimension
    const int16_t *const sftA,          // exponent of shift values for rows of A
    const int16_t *const sftB,          // exponent of shift values for cols of B
    const T *alpha, const T *beta       //
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
} // namespace complex
