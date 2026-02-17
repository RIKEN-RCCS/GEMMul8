#pragma once

#define CHECK_CUDA(x)                                                 \
    do {                                                              \
        cudaError_t _e = (x);                                         \
        if (_e != cudaSuccess) {                                      \
            std::fprintf(stderr, "CUDA error %s:%d: %s\n",            \
                         __FILE__, __LINE__, cudaGetErrorString(_e)); \
            std::abort();                                             \
        }                                                             \
    } while (0)

#define CHECK_CUBLAS(x)                                             \
    do {                                                            \
        cublasStatus_t _s = (x);                                    \
        if (_s != CUBLAS_STATUS_SUCCESS) {                          \
            std::fprintf(stderr, "cuBLAS error %s:%d: status=%d\n", \
                         __FILE__, __LINE__, (int)_s);              \
            std::abort();                                           \
        }                                                           \
    } while (0)

#define CHECK_CUBLASLT(x)                                             \
    do {                                                              \
        cublasStatus_t _s = (x);                                      \
        if (_s != CUBLAS_STATUS_SUCCESS) {                            \
            std::fprintf(stderr, "cuBLASLt error %s:%d: status=%d\n", \
                         __FILE__, __LINE__, (int)_s);                \
            std::abort();                                             \
        }                                                             \
    } while (0)

inline constexpr unsigned warmup          = 30;
inline constexpr unsigned mainloop        = 30;
inline constexpr unsigned long long seedA = 12345;
inline constexpr unsigned long long seedB = 54321;
inline constexpr gemmul8::Backend backend = gemmul8::Backend::FP8;
std::vector<double> phi_list{0.0, 0.5, 1.0, 2.0, 4.0};
template <typename T> inline constexpr unsigned NUM_MODULI_MIN        = 3;
template <typename T> inline constexpr unsigned NUM_MODULI_MAX        = 10;
template <> inline constexpr unsigned NUM_MODULI_MIN<double>          = 9;
template <> inline constexpr unsigned NUM_MODULI_MAX<double>          = 16;
template <> inline constexpr unsigned NUM_MODULI_MIN<cuDoubleComplex> = 9;
template <> inline constexpr unsigned NUM_MODULI_MAX<cuDoubleComplex> = 16;
std::vector<size_t> size_list{1024, 2048, 4096, 8192, 16384};

std::string getDeviceName() {
    cudaDeviceProp deviceProp;
    CHECK_CUDA(cudaGetDeviceProperties(&deviceProp, 0));
    std::string deviceName = deviceProp.name;
    for (char &c : deviceName) {
        if (c == ' ' || c == '/' || c == '\\') c = '_';
    }
    return deviceName;
}

std::string getCurrentDateTime(std::chrono::system_clock::time_point &now) {
    now                  = std::chrono::system_clock::now();
    std::time_t now_time = std::chrono::system_clock::to_time_t(now);
    std::stringstream ss;
    ss << std::put_time(std::localtime(&now_time), "%Y-%m-%d_%H-%M-%S");
    return ss.str();
}

constexpr char ascii_upper(char c) {
    return ('a' <= c && c <= 'z') ? static_cast<char>(c - 'a' + 'A') : c;
}

template <typename T> struct gemmTraits;

template <> struct gemmTraits<float> {
    static constexpr char prefix     = 's';
    static constexpr auto gemm       = &cublasSgemm;
    static constexpr bool is_complex = false;
    using ACCU_TYPE                  = double;

    static constexpr float one() { return 1.0f; }
    static constexpr float zero() { return 0.0f; }
    static constexpr char prefix_upper() { return ascii_upper(prefix); }
};

template <> struct gemmTraits<double> {
    static constexpr char prefix     = 'd';
    static constexpr auto gemm       = &cublasDgemm;
    static constexpr bool is_complex = false;
    using ACCU_TYPE                  = double2;

    static constexpr double one() { return 1.0; }
    static constexpr double zero() { return 0.0; }
    static constexpr char prefix_upper() { return ascii_upper(prefix); }
};

template <> struct gemmTraits<cuFloatComplex> {
    static constexpr char prefix     = 'c';
    static constexpr auto gemm       = &cublasCgemm;
    static constexpr bool is_complex = true;
    using ACCU_TYPE                  = cuDoubleComplex;

    static constexpr cuFloatComplex one() { return cuFloatComplex{1.0f, 0.0f}; }
    static constexpr cuFloatComplex zero() { return cuFloatComplex{0.0f, 0.0f}; }
    static constexpr char prefix_upper() { return ascii_upper(prefix); }
};

template <> struct gemmTraits<cuDoubleComplex> {
    static constexpr char prefix     = 'z';
    static constexpr auto gemm       = &cublasZgemm;
    static constexpr bool is_complex = true;
    using ACCU_TYPE                  = eval::dd::double2_complex;

    static constexpr cuDoubleComplex one() { return cuDoubleComplex{1.0, 0.0}; }
    static constexpr cuDoubleComplex zero() { return cuDoubleComplex{0.0, 0.0}; }
    static constexpr char prefix_upper() { return ascii_upper(prefix); }
};

template <typename T>
size_t max_size() {
    size_t total_memory           = GPU_MEM_MB * 1e6;
    size_t n                      = 1024;
    using accu_t                  = typename gemmTraits<T>::ACCU_TYPE;
    const unsigned num_moduli_max = NUM_MODULI_MAX<T>;
    while (1) {
        size_t size_mat   = n * n;
        size_t lwork      = gemmul8::workSize<gemmTraits<T>::is_complex, backend>(n, n, n, num_moduli_max);
        size_t total_work = 0;
        total_work += 3 * size_mat * sizeof(T);
        total_work += std::max(lwork, size_mat * sizeof(accu_t));
        total_work += 256 * sizeof(accu_t); // for safety
        if (total_work > total_memory) {
            n -= 1024;
            break;
        }
        n += 1024;
    }
    return std::min(*max_element(begin(size_list), end(size_list)), n);
}
