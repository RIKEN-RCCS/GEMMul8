#pragma once

namespace complex {

//------------------------------
// Calculate required work size
//------------------------------
template <gemmul8::Backend backend>
inline size_t workSize(
    size_t m, size_t n, size_t k, unsigned num_moduli,
    bool enable_skip_scalA, bool enable_skip_scalB,
    size_t *workSizeA, size_t *workSizeB //
) {
    const size_t lda_lo  = padding(k);
    const size_t cola_lo = padding(m);
    const size_t ldb_lo  = lda_lo;
    const size_t colb_lo = n;
    const size_t ldc_hi  = cola_lo;
    const size_t colc_hi = colb_lo;

    const size_t sizeA     = lda_lo * cola_lo;
    const size_t sizeB     = ldb_lo * colb_lo;
    const size_t sizeC     = ldc_hi * colc_hi;
    const size_t size_vecA = cola_lo;
    const size_t size_vecB = padding(n);

    const unsigned num_mat        = table::num_mat<backend>(num_moduli);
    const unsigned num_A_lo       = (num_mat + ((enable_skip_scalA) ? 1 : 0)) * 3; // +1 for skip_scalA in accurate mode
    const unsigned num_B_lo       = (num_mat + ((enable_skip_scalB) ? 1 : 0)) * 3; // +1 for skip_scalB in accurate mode
    const unsigned num_C_mid      = num_moduli;
    const unsigned num_C_hi       = ((backend == gemmul8::Backend::INT8) ? 1 : 3) * 3;
    const size_t lwork_nativeGemm = size_t(1) << 25; // 32 MiB

    size_t total_size_A = 255;
    size_t total_size_B = 255;
    size_t total_size_C = 255;
    total_size_A += sizeof(low_t<backend>) * sizeA * num_A_lo;
    total_size_A += sizeof(int16_t) * size_vecA;
    total_size_B += sizeof(low_t<backend>) * sizeB * num_B_lo;
    total_size_B += sizeof(int16_t) * size_vecB;
    total_size_C += sizeof(midx2_t<backend>) * sizeC * (num_C_mid - 1) + std::max(lwork_nativeGemm, sizeof(midx2_t<backend>) * sizeC);
    total_size_C += sizeof(hi_t<backend>) * sizeC * num_C_hi;

    if (workSizeA != nullptr) *workSizeA = total_size_A;
    if (workSizeB != nullptr) *workSizeB = total_size_B;
    return total_size_A + total_size_B + total_size_C;
}

//------------------------------
// GEMM emulation using INT8 Tensor Cores
//------------------------------
template <typename T, gemmul8::Backend backend>
inline std::vector<double> gemm(
    Handle_t handle, cublasOperation_t op_A, cublasOperation_t op_B,
    size_t m, size_t n, size_t k,
    const T *alpha, const T *const A, size_t lda, const T *const B, size_t ldb,
    const T *beta, T *const C, size_t ldc,
    unsigned num_moduli, bool fastmode,
    void *const work, void *const workA, void *const workB,
    bool enable_skip_scalA, bool enable_skip_scalB,
    bool skip_scalA, bool skip_scalB,
    cudaStream_t stream //
) {
    //------------------------------
    // Timer
    //------------------------------
    std::chrono::system_clock::time_point time_stamp;
    std::vector<double> timer(4, 0.0);

    //------------------------------
    // Set constants
    //------------------------------
    const size_t lda_lo          = padding(k);
    const size_t ldb_lo          = lda_lo;
    const size_t ldc_hi          = padding(m);
    const size_t sizeA           = lda_lo * ldc_hi;
    const size_t sizeB           = ldb_lo * n;
    const size_t sizeC           = ldc_hi * n;
    const size_t sizeC_4         = sizeC >> 2;
    const size_t size_vecA       = ldc_hi;
    const size_t size_vecB       = padding(n);
    const bool skipA             = skip_scalA && enable_skip_scalA;
    const bool skipB             = skip_scalB && enable_skip_scalB;
    constexpr hi_t<backend> one  = 1;
    constexpr hi_t<backend> zero = 0;

    //------------------------------
    // Set constant memory
    //------------------------------
    table::upload_constants<backend>(stream);

    //------------------------------
    // Set workspace
    //------------------------------
    const size_t lwork_nativeGemm = size_t(1) << 25; // 32 MiB
    const size_t offsetA          = sizeA * table::num_mat<backend>(num_moduli);
    const size_t offsetB          = sizeB * table::num_mat<backend>(num_moduli);
    const size_t offsetC          = sizeC * ((backend == gemmul8::Backend::INT8) ? 1 : 3);
    void *work_aligned            = align256(work);
    void *workA_aligned           = align256(workA);
    void *workB_aligned           = align256(workB);
    low_t<backend> *const A_lo_1  = reinterpret_cast<low_t<backend> *>((workA_aligned) ? workA_aligned : work_aligned);
    low_t<backend> *const A_lo_2  = A_lo_1 + offsetA;
    low_t<backend> *const A_lo_3  = A_lo_2 + offsetA;
    int16_t *const sftA           = reinterpret_cast<int16_t *>(A_lo_3 + offsetA + ((enable_skip_scalA) ? (sizeA * 3) : 0));
    low_t<backend> *const B_lo_1  = reinterpret_cast<low_t<backend> *>((workB_aligned) ? workB_aligned : ((workA_aligned) ? work_aligned : (sftA + size_vecA)));
    low_t<backend> *const B_lo_2  = B_lo_1 + offsetB;
    low_t<backend> *const B_lo_3  = B_lo_2 + offsetB;
    int16_t *const sftB           = reinterpret_cast<int16_t *>(B_lo_3 + offsetB + ((enable_skip_scalB) ? (sizeB * 3) : 0));
    midx2_t<backend> *const C_mid = reinterpret_cast<midx2_t<backend> *>((workB_aligned) ? ((workA_aligned) ? work_aligned : (sftA + size_vecA)) : (sftB + size_vecB));
    void *work_nativeGemm         = reinterpret_cast<void *>(C_mid + (num_moduli - 1) * sizeC);
    hi_t<backend> *const C_hi_1   = reinterpret_cast<hi_t<backend> *>(static_cast<int8_t *>(work_nativeGemm) + std::max(lwork_nativeGemm, sizeof(midx2_t<backend>) * sizeC));
    hi_t<backend> *const C_hi_2   = C_hi_1 + offsetC;
    hi_t<backend> *const C_hi_3   = C_hi_2 + offsetC;
    low_t<backend> *const A_lo[3] = {A_lo_1, A_lo_2, A_lo_3};
    low_t<backend> *const B_lo[3] = {B_lo_1, B_lo_2, B_lo_3};
    hi_t<backend> *const C_hi[3]  = {C_hi_1, C_hi_2, C_hi_3};

    //------------------------------
    // Set handle
    //------------------------------
    set_handle<backend>(stream, handle, ldc_hi, n, lda_lo, ldb_lo, ldb_lo, ldc_hi, work_nativeGemm, lwork_nativeGemm);

    //------------------------------
    // Scaling
    // A =: diag(2^sftA) * A', A' is integer
    // B =: B' * diag(2^sftB), B' is integer
    // Then, calculating mod for all moduli
    // A_lo := A' - p[i]*round(A'/p[i])  (A_lo is INT8 or 2*FP8)
    // B_lo := B' - p[i]*round(B'/p[i])  (B_lo is INT8 or 2*FP8)
    //------------------------------
    timing(stream, time_stamp);
    if (!(skipA && skipB)) {
        // When both scalingA & scalingB are skipped, this is skiped.
        if (fastmode) {
            fast::scaling<backend, T>(stream, op_A, op_B, m, n, k, num_moduli,
                                      A, lda, A_lo, lda_lo, sizeA, sftA,
                                      B, ldb, B_lo, ldb_lo, sizeB, sftB,
                                      skipA, skipB);
        } else {
            low_t<backend> *const A_lo_high_1  = ((enable_skip_scalA) ? (A_lo_3 + offsetA) : A_lo_1);
            low_t<backend> *const B_lo_high_1  = ((enable_skip_scalB) ? (B_lo_3 + offsetB) : B_lo_1);
            low_t<backend> *const A_lo_high[3] = {A_lo_high_1, A_lo_high_1 + sizeA, A_lo_high_1 + sizeA * 2};
            low_t<backend> *const B_lo_high[3] = {B_lo_high_1, B_lo_high_1 + sizeB, B_lo_high_1 + sizeB * 2};
            accu::scaling<backend, T>(stream, handle, op_A, op_B, m, n, k, num_moduli,
                                      A, lda, A_lo, A_lo_high, lda_lo, sizeA, sftA,
                                      B, ldb, B_lo, B_lo_high, ldb_lo, sizeB, sftB,
                                      C_hi, ldc_hi,
                                      skipA, skipB);
        }
    }
    timing(stream, time_stamp, timer[0]);

    low_t<backend> *Ar = A_lo[0], *Ai = A_lo[1], *Ari = A_lo[2];
    low_t<backend> *Br = B_lo[0], *Bi = B_lo[1], *Bri = B_lo[2];
    hi_t<backend> *Cr = C_hi[0], *Ci = C_hi[1], *Cri = C_hi[2];
    for (unsigned i = 0; i < num_moduli; ++i) {
        //-----------------------------
        // Error-free matrix multiplication
        // C_hi := A_lo*B_lo
        //------------------------------
        if constexpr (backend == gemmul8::Backend::INT8) {
            gemm_low_prec_i8x3(
                stream, handle, ldc_hi, n, lda_lo,
                &one, &one, &one, Ar, Ai, Ari, lda_lo, Br, Bi, Bri, ldb_lo,
                &zero, &zero, &zero, Cr, Ci, Cri, ldc_hi);
            Ar += sizeA;
            Ai += sizeA;
            Ari += sizeA;
            Br += sizeB;
            Bi += sizeB;
            Bri += sizeB;
        } else {
            if (i < table::not_Karatsuba) {
                low_t<backend> *Ax9[9] = {Ar, Ar + sizeA, Ar + sizeA, Ai, Ai + sizeA, Ai + sizeA, Ari, Ari + sizeA, Ari + sizeA};
                low_t<backend> *Bx9[9] = {Br + sizeB, Br, Br + sizeB, Bi + sizeB, Bi, Bi + sizeB, Bri + sizeB, Bri, Bri + sizeB};
                hi_t<backend> *Cx9[9]  = {Cr, Cr + sizeC, Cr + 2 * sizeC, Ci, Ci + sizeC, Ci + 2 * sizeC, Cri, Cri + sizeC, Cri + 2 * sizeC};
                gemm_low_prec_f8x9(stream, handle, ldc_hi, n, lda_lo, &one, Ax9, lda_lo, Bx9, ldb_lo, &zero, Cx9, ldc_hi);
                Ar += 2 * sizeA;
                Ai += 2 * sizeA;
                Ari += 2 * sizeA;
                Br += 2 * sizeB;
                Bi += 2 * sizeB;
                Bri += 2 * sizeB;
            } else {
                low_t<backend> *Ax9[9] = {Ar, Ar + sizeA, Ar + 2 * sizeA, Ai, Ai + sizeA, Ai + 2 * sizeA, Ari, Ari + sizeA, Ari + 2 * sizeA};
                low_t<backend> *Bx9[9] = {Br, Br + sizeB, Br + 2 * sizeB, Bi, Bi + sizeB, Bi + 2 * sizeB, Bri, Bri + sizeB, Bri + 2 * sizeB};
                hi_t<backend> *Cx9[9]  = {Cr, Cr + sizeC, Cr + 2 * sizeC, Ci, Ci + sizeC, Ci + 2 * sizeC, Cri, Cri + sizeC, Cri + 2 * sizeC};
                gemm_low_prec_f8x9(stream, handle, ldc_hi, n, lda_lo, &one, Ax9, lda_lo, Bx9, ldb_lo, &zero, Cx9, ldc_hi);
                Ar += 3 * sizeA;
                Ai += 3 * sizeA;
                Ari += 3 * sizeA;
                Br += 3 * sizeB;
                Bi += 3 * sizeB;
                Bri += 3 * sizeB;
            }
        }
        timing(stream, time_stamp, timer[1]);

        //------------------------------
        // Calculating mod
        // C_mid[i] := C_hi - round(C_hi/p[i])*p[i]
        //------------------------------
        conv_hi2mid<backend>(stream, i, sizeC_4, C_hi, C_mid + i * sizeC);
        timing(stream, time_stamp, timer[2]);
    }

    //------------------------------
    // Accumulation and Inverse scaling
    // C64f = sum(qi*Pi*C_mid[i]),
    //  where
    //      Pi := P/p[i],
    //      P  := prod(p[all]),
    //      mod(qi*Pi, p[i]) \equiv 1.
    // C := C64f - round(C64f/P)*P
    // C := diag(2^sftA) * C * diag(2^sftB)
    //------------------------------
    inverse_scaling<backend, T>(stream, num_moduli, m, n, C_mid, ldc_hi, sizeC, C, ldc, sftA, sftB, alpha, beta);
    timing(stream, time_stamp, timer[3]);

    //------------------------------
    // Clean up Lt memory
    //------------------------------
    cleanup_handle(handle);
    return timer;
}

} // namespace complex
