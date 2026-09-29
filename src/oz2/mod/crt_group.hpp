#pragma once
#include "../common/crt_group_plan.hpp"
#include "../undo_scaling/crt_tree.hpp"
#include "reconstruct_fp8.hpp"
#include "complex_2m.hpp"

namespace gemmul8::mod {

template <Backend B, unsigned N, bool Complex,
          unsigned First, unsigned End,
          bool Herk = false, unsigned Vector = 1>
struct crt_group_constants {
    static constexpr unsigned first = First;
    static constexpr bool center    = Herk || End == First + 1;
    static constexpr bool raw       = B == Backend::FP8 && !center &&
                                      (!Complex || Vector == 1 || End == First + 2);

    using Moduli = common::crt_moduli<B, N, Complex>;
    using Leaf   = undo_scaling::crt::leaf_constants<Moduli, First, End>;

    struct data {
        int32_t re[End - First][raw ? 3 : 1]{}, im[End - First][raw ? 3 : 1]{};
        uint64_t re_bound = 0, im_bound = 0;
    };

    template <unsigned I>
    static constexpr void coefficient(data &d) {
        constexpr uint32_t p = Moduli::get(I), q = Leaf::q;
        constexpr int64_t e     = Leaf::coefficients.coefficient[I - First];
        constexpr auto centered = [](int64_t v) {
            v %= q;
            if (v < 0) v += q;
            return int32_t(v > q / 2 ? v - q : v);
        };
        int32_t a = int32_t(e), b = 0;
        if constexpr (Complex) {
            constexpr unsigned index = common::table::active_index<B, N, I, true>;
            constexpr int64_t hr     = complex_2m_traits<B, index>::half_root;
            a                        = centered(e * ((p + 1) / 2));
            b                        = centered(e * hr);
        }

        const auto append = [&](unsigned j, int32_t ar, int32_t ai, uint64_t re_limit, uint64_t im_limit) {
            d.re[I - First][j] = ar;
            d.im[I - First][j] = ai;
            d.re_bound += uint64_t(ar < 0 ? -int64_t(ar) : ar) * re_limit;
            d.im_bound += uint64_t(ai < 0 ? -int64_t(ai) : ai) * im_limit;
        };

        if constexpr (raw) {
            constexpr auto scheme = common::fp8_plan::scheme<int32_t(p)>;
            for (unsigned j = 0; j < scheme.products; ++j) {
                append(j, centered(int64_t(a) * scheme.coefficient[j]),
                       centered(int64_t(b) * scheme.coefficient[j]),
                       (1ULL << 24) * (Complex ? 2U : 1U), 1ULL << 25);
            }
        } else {
            append(0, a, b, center ? (Complex ? p - 1 : p / 2) : (Complex ? 4U : 2U) * p,
                   center ? p - 1 : 4U * p);
        }

        if constexpr (I + 1 < End) coefficient<I + 1>(d);
    }

    static constexpr data coefficients = [] { data d{}; coefficient<First>(d); return d; }();
    static constexpr uint64_t re_bias  = ((coefficients.re_bound + Leaf::q - 1) / Leaf::q) * Leaf::q;
    static constexpr uint64_t im_bias  = ((coefficients.im_bound + Leaf::q - 1) / Leaf::q) * Leaf::q;

    static constexpr bool bounded(uint64_t maximum) {
        if (maximum > uint64_t(INT64_MAX)) return false;
        if constexpr (Leaf::power_of_two) return true;
        constexpr uint64_t scale = 1ULL << (Leaf::q_bits - 1);
        constexpr uint64_t error = (1ULL << (Leaf::q_bits + 31)) % Leaf::q;
        const uint64_t high      = maximum / scale;
        return high <= UINT32_MAX && ((scale - 1 + ((high * error + UINT32_MAX) >> 32)) < Leaf::q);
    }

    static_assert(bounded(re_bias + coefficients.re_bound) && bounded(im_bias + coefficients.im_bound));
};

template <Backend B, unsigned N, bool Complex,
          unsigned First, unsigned End, bool Herk = false>
struct crt_group_input {
    using Hi       = common::hi_t<B>;
    using Opposite = std::conditional_t<B == Backend::INT8, uint32_t, uint64_t>;
    const Hi *ptr;
    uint32_t index, sizeC;
    Opposite opposite = 0;

    template <unsigned L, class V> __device__ __forceinline__ static auto component(V v) {
        if constexpr (L == 0) return v.x;
        else if constexpr (L == 1) return v.y;
        else if constexpr (L == 2) return v.z;
        else return v.w;
    }

    static constexpr bool center = Herk || (End == (First + 1));

    template <unsigned Mod>
    __device__ __forceinline__ static int32_t reduce(int32_t a) {
        if constexpr (center) return mod_small<B, Mod, Complex>(a);
        else {
            constexpr int32_t p = common::table::moduli<B, Mod, Complex>;
            if constexpr ((p & (p - 1)) == 0) {
                return a & (p - 1);
            } else {
                return int32_t(uint32_t(a) - uint32_t(p) * uint32_t(__mulhi(a, common::table::p_inv_32_v<B, Mod, Complex>)));
            }
        }
    }

    template <unsigned I, unsigned Vector = 1, unsigned Lane = 0, bool Raw = false>
    __device__ __forceinline__ int32_t residue(const Hi *p, size_t i) const {
        constexpr unsigned mod = common::table::active_index<B, N, I, Complex>;
        if constexpr (Raw) {
            if constexpr (Vector == 1) return int32_t(p[i]);
            else {
                using V = std::conditional_t<Vector == 2, float2, float4>;
                return int32_t(component<Lane>(reinterpret_cast<const V *>(p)[i / Vector]));
            }
        } else if constexpr (Vector == 1) {
            if constexpr (B == Backend::INT8) {
                return reduce<mod>(p[i]);
            } else {
                return load_fp8_product<mod, Complex, !center>(p, i, sizeC);
            }
        } else {
            using V = std::conditional_t<B == Backend::INT8,
                                         std::conditional_t<Vector == 2, int2, int4>,
                                         std::conditional_t<Vector == 2, float2, float4>>;

            const auto a = reinterpret_cast<const V *>(p)[i / Vector];
            if constexpr (B == Backend::INT8) {
                return reduce<mod>(component<Lane>(a));
            } else {
                float b = 0, c = 0;
                if constexpr (fp8_product_count<mod, Complex> > 1) {
                    b = component<Lane>(reinterpret_cast<const V *>(p + sizeC)[i / Vector]);
                }
                if constexpr (fp8_product_count<mod, Complex> > 2) {
                    c = component<Lane>(reinterpret_cast<const V *>(p + 2 * sizeC)[i / Vector]);
                }
                return mod_f32x3_2_i32<mod, !center, Complex>(component<Lane>(a), b, c);
            }
        }
    }

    static constexpr unsigned products(unsigned i) {
        if constexpr (B == Backend::INT8) return 1;
        else return common::fp8_plan::get(common::crt_modulus<B, Complex>(N, i)).products;
    }

    static constexpr unsigned offset(unsigned i) {
        unsigned r = 0;
        for (unsigned j = First; j < i; ++j) {
            r += products(j) * (Complex && !Herk ? 2U : 1U);
        }
        return r;
    }

    template <unsigned I> static constexpr unsigned offset_v   = offset(I);
    template <unsigned I> static constexpr unsigned products_v = products(I);

    template <unsigned I, unsigned Vector = 1, unsigned Lane = 0, bool Raw = false, unsigned Plane = 0>
    __device__ __forceinline__ int2 load() const {
        constexpr unsigned off = offset_v<I> + Plane, count = products_v<I>;
        const Hi *p        = ptr + size_t(sizeC) * off;
        const int32_t plus = residue<I, Vector, Lane, Raw>(p, index);
        if constexpr (!Complex) return {plus, 0};
        else if constexpr (Herk) {
            constexpr unsigned width = B == Backend::INT8 ? 8U : 16U;
            using S                  = std::conditional_t<B == Backend::INT8, int8_t, int16_t>;
            return {plus, int32_t(S(opposite >> (width * (I - First))))};
        } else {
            return {plus, residue<I, Vector, Lane, Raw>(ptr + size_t(sizeC) * (off + count), index)};
        }
    }

    template <unsigned I = First>
    __device__ __forceinline__ Opposite pack_opposite(size_t i) const {
        constexpr unsigned width = B == Backend::INT8 ? 8U : 16U;
        constexpr unsigned off   = offset_v<I>;
        const auto r             = residue<I>(ptr + sizeC * off, i);
        Opposite v               = Opposite(uint32_t(r) & ((1U << width) - 1U)) << (width * (I - First));
        if constexpr (I + 1 < End) {
            v |= pack_opposite<I + 1>(i);
        }
        return v;
    }
};

template <class Constants, unsigned I, unsigned End, bool Complex,
          unsigned Vector, unsigned Lane, unsigned Plane = 0, class Input, class Acc>
__device__ __forceinline__ void crt_group_accumulate(const Input &in, Acc &re, Acc &im) {
    const auto v         = in.template load<I, Vector, Lane, Constants::raw, Plane>();
    constexpr unsigned j = I - Constants::first;
    constexpr int32_t a  = Constants::coefficients.re[j][Plane];
    if constexpr (Complex) {
        constexpr int32_t b = Constants::coefficients.im[j][Plane];
        re += Acc(v.x + v.y) * a;
        im += Acc(v.y - v.x) * b;
    } else {
        re += Acc(v.x) * a;
    }
    if constexpr (Constants::raw && Plane + 1 < Input::template products_v<I>) {
        crt_group_accumulate<Constants, I, End, Complex, Vector, Lane, Plane + 1>(in, re, im);
    } else if constexpr (I + 1 < End) {
        crt_group_accumulate<Constants, I + 1, End, Complex, Vector, Lane>(in, re, im);
    }
}

template <Backend B, unsigned N, bool Complex,
          unsigned First, unsigned End,
          bool Herk = false, bool Flip = false,
          unsigned Vector = 1, unsigned Lane = 0>
__device__ __forceinline__ auto crt_group_value(const crt_group_input<B, N, Complex, First, End, Herk> &in) {
    using K = crt_group_constants<B, N, Complex, First, End, Herk, Vector>;
    using L = typename K::Leaf;
    if constexpr (!Complex && End == First + 1) {
        constexpr int32_t p = int32_t(L::q);
        const int32_t a     = in.template load<First, Vector, Lane>().x;
        return uint32_t(a < 0 ? a + p : a);
    }
    using Acc = std::conditional_t<(K::re_bias + K::coefficients.re_bound <= INT32_MAX &&
                                    K::im_bias + K::coefficients.im_bound <= INT32_MAX),
                                   int32_t, int64_t>;
    Acc re = Acc(K::re_bias), im = Acc(K::im_bias);
    crt_group_accumulate<K, First, End, Complex, Vector, Lane>(in, re, im);
    const uint32_t a = L::template reduce<K::re_bias + K::coefficients.re_bound>(uint64_t(re));
    if constexpr (!Complex) {
        return a;
    } else {
        uint32_t b = L::template reduce<K::im_bias + K::coefficients.im_bound>(uint64_t(im));
        if constexpr (Flip) {
            b = b ? L::q - b : 0;
        }
        return uint2{a, b};
    }
}

#if defined(__CUDA_ARCH__)
    #if (__CUDA_ARCH__ == 900) || (__CUDA_ARCH__ == 1000) || (__CUDA_ARCH__ == 1030)
        #define GEMMUL8_CRT_MIN_BLOCKS 8
    #elif (__CUDA_ARCH__ == 860) || (__CUDA_ARCH__ == 870) || (__CUDA_ARCH__ == 890) || \
        (__CUDA_ARCH__ == 1100) || (1200 <= __CUDA_ARCH__ && __CUDA_ARCH__ < 1300) ||   \
        (__CUDA_ARCH__ == 1070)
        #define GEMMUL8_CRT_MIN_BLOCKS 4
    #else
        #define GEMMUL8_CRT_MIN_BLOCKS 1
    #endif
#else
    #define GEMMUL8_CRT_MIN_BLOCKS 1
#endif

template <Backend B, unsigned N, bool Complex, unsigned First, unsigned End,
          cublasFillMode_t Uplo, bool Herk = false, bool Flip = false>
__global__ __launch_bounds__(256, Herk ? 1 : GEMMUL8_CRT_MIN_BLOCKS) void mod_hi2group_kernel(
    unsigned n, unsigned ldc, unsigned sizeC,
    const common::hi_t<B> *__restrict__ hi,
    std::conditional_t<Complex, uint2, uint32_t> *__restrict__ out //
) {
    using Input = crt_group_input<B, N, Complex, First, End, Herk>;
    if constexpr (Herk) {
        if constexpr (Uplo == CUBLAS_FILL_MODE_UPPER) {
            if (blockIdx.x > blockIdx.y) return;
        } else {
            if (blockIdx.x < blockIdx.y) return;
        }

        __shared__ typename Input::Opposite opposite[32][33];
        Input input{hi, 0, sizeC};

#pragma unroll
        for (unsigned j = 0; j < 32; j += 8) {
            const unsigned row                     = blockIdx.y * 32 + threadIdx.x;
            const unsigned col                     = blockIdx.x * 32 + threadIdx.y + j;
            opposite[threadIdx.y + j][threadIdx.x] = row < n && col < n ? input.pack_opposite(col * ldc + row) : 0;
        }
        __syncthreads();

#pragma unroll 1
        for (unsigned j = 0; j < 32; j += 8) {
            const unsigned row = blockIdx.x * 32 + threadIdx.x;
            const unsigned col = blockIdx.y * 32 + threadIdx.y + j;

            if (row >= n || col >= n) continue;
            if constexpr (Uplo == CUBLAS_FILL_MODE_UPPER) {
                if (row > col) continue;
            } else {
                if (row < col) continue;
            }

            input.index      = col * ldc + row;
            input.opposite   = opposite[threadIdx.x][threadIdx.y + j];
            out[input.index] = crt_group_value<B, N, Complex, First, End, true, Flip>(input);
        }

    } else {
        unsigned index;
        if constexpr (Uplo == CUBLAS_FILL_MODE_FULL) {
            constexpr unsigned vector = Complex ? 2U : 4U;

            index = (blockIdx.x * blockDim.x + threadIdx.x) * vector;
            if (index >= sizeC) return;

            const Input input{hi, index, sizeC};
            const auto a = crt_group_value<B, N, Complex, First, End, false, false, vector, 0>(input);
            const auto b = crt_group_value<B, N, Complex, First, End, false, false, vector, 1>(input);
            if constexpr (Complex) {
                reinterpret_cast<uint4 *>(out)[index / 2] = {a.x, a.y, b.x, b.y};
            } else {
                const auto c = crt_group_value<B, N, Complex, First, End, false, false, vector, 2>(input);
                const auto d = crt_group_value<B, N, Complex, First, End, false, false, vector, 3>(input);

                reinterpret_cast<uint4 *>(out)[index / 4] = {a, b, c, d};
            }

        } else {
            const unsigned row = blockIdx.x * 32 + threadIdx.x;
            const unsigned col = blockIdx.y * 4 + threadIdx.y;

            if (row >= ldc || col >= n) return;
            if constexpr (Uplo == CUBLAS_FILL_MODE_UPPER) {
                if (row > col) return;
            } else {
                if (row < col) return;
            }

            index      = col * ldc + row;
            out[index] = crt_group_value<B, N, Complex, First, End>(Input{hi, index, sizeC});
        }
    }
}

template <Backend B, unsigned N, bool Complex,
          unsigned First, unsigned End,
          cublasFillMode_t Uplo, bool Herk = false, bool Flip = false>
inline void mod_hi2group(
    cudaStream_t stream,
    unsigned n, size_t ldc,
    const common::hi_t<B> *hi,
    void *out //
) {
    const size_t sizeC = ldc * n;
    assert(common::use_grouped_crt(sizeC));

    const dim3 threads = Herk                            ? dim3(32, 8)
                         : Uplo == CUBLAS_FILL_MODE_FULL ? dim3(256)
                                                         : dim3(32, 4);

    const dim3 grid = Herk                            ? dim3((n + 31) / 32, (n + 31) / 32)
                      : Uplo == CUBLAS_FILL_MODE_FULL ? dim3((sizeC / (Complex ? 2U : 4U) + 255) / 256)
                                                      : dim3((ldc + 31) / 32, (n + 3) / 4);

    using Out = std::conditional_t<Complex, uint2, uint32_t>;
    mod_hi2group_kernel<B, N, Complex, First, End, Uplo, Herk, Flip>
        <<<grid, threads, 0, stream>>>(n, ldc, sizeC, hi, static_cast<Out *>(out));
}

} // namespace gemmul8::mod
