#pragma once

//------------------------------
// CUDA grid and thread configuration
//------------------------------
inline constexpr size_t threads_scaling     = 256;
inline constexpr size_t threads_conv_hi2mid = 256;
inline constexpr size_t threads_invscal     = 128;
inline constexpr int TILE_DIM               = 32; // better than 16 for A100, GH200

//------------------------------
// Iteration threshold for modular reduction
// Used to decide mod implementation based on num_moduli
//------------------------------
template <gemmul8::Backend backend = gemmul8::Backend::INT8> struct threshold;
template <> struct threshold<gemmul8::Backend::INT8> {
    static constexpr int P_is_double = 6;
    static constexpr int S           = 7;
    static constexpr int M           = 15;
    static constexpr int L           = 25;
};
template <> struct threshold<gemmul8::Backend::FP8> {
    static constexpr int P_is_double = 5;
    static constexpr int S           = 5;
    static constexpr int M           = 12;
    static constexpr int L           = 20;
};

//------------------------------
// Pad size to multiple of 32 (for alignment)
//------------------------------
static __forceinline__ __host__ __device__ size_t padding(const size_t n) { return 256 * ((n + 255) / 256); }

inline void *align256(void *p) {
    constexpr std::uintptr_t A = 256;
    std::uintptr_t x           = reinterpret_cast<std::uintptr_t>(p);
    x                          = (x + (A - 1)) & ~(A - 1);
    return reinterpret_cast<void *>(x);
}

//------------------------------
// Start timing measurement
//------------------------------
static inline void timing(cudaStream_t &stream, std::chrono::system_clock::time_point &time_stamp) {
    cudaStreamSynchronize(stream);
    time_stamp = std::chrono::system_clock::now();
}

//------------------------------
// Stop timing and accumulate elapsed time (ns)
//------------------------------
static inline void timing(cudaStream_t &stream, std::chrono::system_clock::time_point &time_stamp, double &timer) {
    cudaStreamSynchronize(stream);
    std::chrono::system_clock::time_point time_now = std::chrono::system_clock::now();
    timer += std::chrono::duration_cast<std::chrono::nanoseconds>(time_now - time_stamp).count();
    time_stamp = time_now;
}
