#pragma once
#include <cstdint>
#include <cstring>
#include <type_traits>

namespace gemmul8::undo_scaling::crt {

template <unsigned N> struct uintx { uint32_t word[N]{}; };

template <unsigned N>
using packed_uint = std::conditional_t<(N <= 2), uint64_t, unsigned __int128>;

template <unsigned N>
__host__ __device__ __forceinline__ constexpr packed_uint<N> pack(const uintx<N> &a) {
    packed_uint<N> out = 0;
#pragma unroll
    for (unsigned i = 0; i < N; ++i) out |= packed_uint<N>(a.word[i]) << (32 * i);
    return out;
}

template <unsigned N>
__host__ __device__ __forceinline__ constexpr void unpack(packed_uint<N> x, uintx<N> &a) {
#pragma unroll
    for (unsigned i = 0; i < N; ++i) a.word[i] = uint32_t(x >> (32 * i));
}

template <unsigned TO, unsigned FROM>
__host__ __device__ __forceinline__ constexpr uintx<TO> resize(const uintx<FROM> &a) {
    uintx<TO> r{};
#pragma unroll
    for (unsigned i = 0; i < TO && i < FROM; ++i) r.word[i] = a.word[i];
    return r;
}

template <unsigned N>
__host__ __device__ __forceinline__ constexpr bool ge(const uintx<N> &a, const uintx<N> &b) {
    if constexpr (N >= 2 && N <= 4) {
        return pack(a) >= pack(b);
    }
    bool result = true;
#pragma unroll
    for (unsigned i = 0; i < N; ++i) {
        result = a.word[i] == b.word[i] ? result : a.word[i] > b.word[i];
    }
    return result;
}

template <unsigned N>
__host__ __device__ __forceinline__ constexpr uint32_t sub(uintx<N> &a, const uintx<N> &b) {
    if constexpr (N >= 2 && N <= 4) {
        const auto x = pack(a), y = pack(b);
        unpack<N>(x - y, a);
        return uint32_t(x < y);
    }
    uint32_t borrow = 0;
#pragma unroll
    for (unsigned i = 0; i < N; ++i) {
        const uint32_t x = a.word[i], y = b.word[i];
        a.word[i] = x - y - borrow;
        borrow    = uint32_t(x < y || (x == y && borrow));
    }
    return borrow;
}

template <unsigned N>
__host__ __device__ __forceinline__ constexpr uint32_t add(uintx<N> &a, const uintx<N> &b) {
    if constexpr (N >= 2 && N <= 4) {
        const auto x = pack(a), y = pack(b), z = x + y;
        unpack<N>(z, a);
        if constexpr (N == 3) {
            return uint32_t(z >> 96);
        } else {
            return uint32_t(z < x);
        }
    }
    uint32_t carry = 0;
#pragma unroll
    for (unsigned i = 0; i < N; ++i) {
        const uint64_t t = uint64_t(a.word[i]) + b.word[i] + carry;
        a.word[i]        = uint32_t(t);
        carry            = uint32_t(t >> 32);
    }
    return carry;
}

template <unsigned N>
__host__ __device__ __forceinline__ constexpr uintx<N> select(bool take_a, const uintx<N> &a, const uintx<N> &b) {
    uintx<N> r{};
#pragma unroll
    for (unsigned i = 0; i < N; ++i) r.word[i] = take_a ? a.word[i] : b.word[i];
    return r;
}

template <unsigned N>
__host__ __device__ __forceinline__ constexpr uintx<N> sub_mod(uintx<N> a, const uintx<N> &b, const uintx<N> &p) {
    const uint32_t mask = 0U - sub(a, b);
    uintx<N> correction{};
#pragma unroll
    for (unsigned i = 0; i < N; ++i) correction.word[i] = p.word[i] & mask;
    add(a, correction);
    return a;
}

template <unsigned N>
__host__ __device__ __forceinline__ constexpr uintx<N> half(uintx<N> a, uint32_t high = 0) {
#pragma unroll
    for (unsigned j = N; j > 0; --j) {
        const unsigned i    = j - 1;
        const uint32_t next = a.word[i] & 1U;
        a.word[i]           = (a.word[i] >> 1) | (high << 31);
        high                = next;
    }
    return a;
}

template <unsigned N>
__host__ __device__ __forceinline__ constexpr unsigned bits(const uintx<N> &a) {
    unsigned out = 0;
    for (unsigned i = 0; i < N; ++i) {
        uint32_t x = a.word[i];
        unsigned k = 0;
        while (x != 0) {
            ++k;
            x >>= 1;
        }
        if (k != 0) out = 32U * i + k;
    }
    return out;
}

__host__ __device__ __forceinline__ constexpr uint32_t negative_inverse32(uint32_t odd) {
    uint32_t x = 1;
    for (unsigned i = 0; i < 5; ++i) x *= 2U - odd * x;
    return 0U - x;
}

template <unsigned OUT, unsigned A, unsigned B>
__host__ __device__ __forceinline__ constexpr uintx<OUT> multiply(const uintx<A> &a, const uintx<B> &b) {
    uintx<OUT> out{};
#pragma unroll
    for (unsigned i = 0; i < A && i < OUT; ++i) {
        uint32_t carry = 0;
#pragma unroll
        for (unsigned j = 0; j < B && i + j < OUT; ++j) {
            const uint64_t t = uint64_t(a.word[i]) * b.word[j] + out.word[i + j] + carry;
            out.word[i + j]  = uint32_t(t);
            carry            = uint32_t(t >> 32);
        }
        if (i + B < OUT) out.word[i + B] = carry;
    }
    return out;
}

template <unsigned N, bool HALF_RADIX = false>
__device__ __forceinline__ uintx<N> montgomery(
    const uintx<N> &a,
    const uintx<N> &k,
    const uintx<N> &p,
    const uint32_t neg_inv //
) {
    if constexpr (N == 1) {
        const uint32_t m = a.word[0] * (k.word[0] * (0U - neg_inv));
        const uint32_t h = __umulhi(a.word[0], k.word[0]);
        const uint32_t l = __umulhi(m, p.word[0]);
        return uintx<1>{{h - l + (h < l ? p.word[0] : 0U)}};
    } else {
        constexpr unsigned size = 2 * N + (HALF_RADIX ? 0 : 1);
        auto t                  = multiply<size>(a, k);
#pragma unroll
        for (unsigned i = 0; i < N; ++i) {
            const uint32_t m = t.word[i] * neg_inv;
            uint32_t carry   = 0;
#pragma unroll
            for (unsigned j = 0; j < N; ++j) {
                const uint64_t z = uint64_t(m) * p.word[j] + t.word[i + j] + carry;
                t.word[i + j]    = uint32_t(z);
                carry            = uint32_t(z >> 32);
            }
#pragma unroll
            for (unsigned j = i + N; j < size; ++j) {
                const uint64_t z = uint64_t(t.word[j]) + carry;
                t.word[j]        = uint32_t(z);
                carry            = uint32_t(z >> 32);
            }
        }
        uintx<N> out{};
#pragma unroll
        for (unsigned i = 0; i < N; ++i) out.word[i] = t.word[N + i];
        uintx<N> reduced      = out;
        const uint32_t borrow = sub(reduced, p);
        if constexpr (HALF_RADIX) {
            return select(borrow == 0, reduced, out);
        } else {
            return select(t.word[2 * N] != 0 || borrow == 0, reduced, out);
        }
    }
}

template <unsigned TOP, unsigned N>
__device__ __forceinline__ double magnitude_to_double(const uintx<N> &a, uint64_t sign) {
    if constexpr (TOP <= 2) {
        uint64_t value = a.word[0];
        if constexpr (TOP == 2) {
            value |= uint64_t(a.word[1]) << 32;
        }
        const uint64_t bits = uint64_t(__double_as_longlong(__ull2double_rn(value)));
        return __longlong_as_double(static_cast<long long>(sign | bits));
    } else {
        const uint32_t hi = a.word[TOP - 1];
        if (hi == 0) {
            return magnitude_to_double<TOP - 1>(a, sign);
        }

        const unsigned shift = unsigned(__clz(hi));
        uint64_t window      = ((uint64_t(hi) << 32) | a.word[TOP - 2]) << shift;
        const uint32_t third = a.word[TOP - 3];
        window |= uint64_t(third) >> (32 - shift);
        bool sticky = (third << shift) != 0;

#pragma unroll
        for (unsigned i = 0; i + 3 < TOP; ++i) sticky |= a.word[i] != 0;

        const double rounded = __ull2double_rn(window | uint64_t(sticky));
        const uint64_t bits  = uint64_t(__double_as_longlong(rounded));
        const uint64_t scale = uint64_t(32U * TOP - 64U - shift) << 52;
        return __longlong_as_double(static_cast<long long>(sign | (bits + scale)));
    }
}

template <unsigned N>
__device__ __forceinline__ double centered_double(uintx<N> a, const uintx<N> &p) {
    const auto h            = half(p);
    const bool positive     = ge(h, a);
    auto negative_magnitude = p;
    sub(negative_magnitude, a);
    a = select(positive, a, negative_magnitude);
    return magnitude_to_double<N>(a, positive ? 0ULL : (1ULL << 63));
}

using constant_uint = uintx<7>;

template <class Moduli>
constexpr constant_uint product(unsigned first, unsigned end) {
    constant_uint p{{1}};
    for (unsigned j = first; j < end; ++j) {
        uint64_t carry = 0;
        for (unsigned i = 0; i < 7; ++i) {
            const uint64_t x = uint64_t(p.word[i]) * Moduli::get(j) + carry;
            p.word[i]        = uint32_t(x);
            carry            = x >> 32;
        }
    }
    return p;
}

template <class Moduli>
constexpr unsigned split_at(unsigned first, unsigned end) {
    unsigned best = first, best_score = 256;
    for (unsigned i = first; i < end;) {
        uint64_t q = 1;
        do {
            q *= Moduli::get(i++);
        } while (i < end && q * Moduli::get(i) <= UINT32_MAX);
        if (i == end) break;
        const unsigned l     = bits(product<Moduli>(first, i));
        const unsigned r     = bits(product<Moduli>(i, end));
        const unsigned score = l > r ? l - r : r - l;
        if (score < best_score) {
            best       = i;
            best_score = score;
        }
    }
    return best;
}

using constant_wide = unsigned __int128;

template <unsigned N>
constexpr constant_wide constant_value(const uintx<N> &a) {
    constant_wide out = 0;
    for (unsigned i = N; i > 0; --i) out = (out << 32) | a.word[i - 1];
    return out;
}

constexpr constant_wide constant_inverse(constant_wide a, constant_wide p) {
    __int128 t = 0, next_t = 1;
    constant_wide r = p, next_r = a % p;
    while (next_r != 0) {
        const constant_wide q     = r / next_r;
        const __int128 tmp_t      = t - __int128(q) * next_t;
        t                         = next_t;
        next_t                    = tmp_t;
        const constant_wide tmp_r = r - q * next_r;
        r                         = next_r;
        next_r                    = tmp_r;
    }
    return constant_wide(t < 0 ? t + __int128(p) : t);
}

template <class Moduli, unsigned FIRST, unsigned END> struct leaf_constants {
    static constexpr unsigned count  = END - FIRST;
    static constexpr uint32_t q      = product<Moduli>(FIRST, END).word[0];
    static constexpr unsigned q_bits = bits(uintx<1>{{q}});

    struct data {
        int32_t coefficient[count]{};
        uint64_t bound = 0;
    };

    static constexpr data coefficients = [] {
        data out{};
        for (unsigned j = 0; j < count; ++j) {
            const uint32_t p   = Moduli::get(FIRST + j);
            const uint32_t m   = q / p;
            const uint64_t e   = uint64_t(m) * uint32_t(constant_inverse(m, p));
            const int64_t c    = e > q / 2 ? int64_t(e) - q : int64_t(e);
            out.coefficient[j] = int32_t(c);
            out.bound += uint64_t(c < 0 ? -c : c) * (p / 2);
        }
        return out;
    }();

    static constexpr uint64_t bias       = ((coefficients.bound + q - 1) / q) * q;
    static constexpr uint64_t maximum    = bias + coefficients.bound;
    static constexpr bool power_of_two   = (q & (q - 1)) == 0;
    static constexpr uint32_t reciprocal = power_of_two ? 0U : uint32_t((1ULL << (q_bits + 31)) / q);

    static constexpr uint64_t remainder_bound =
        uint64_t(q) + (1ULL << (q_bits - 1)) - 2 + (maximum >> (q_bits - 1));

    template <uint64_t Maximum = maximum>
    __device__ __forceinline__ static uint32_t reduce(uint64_t x) {
        if constexpr (power_of_two) {
            return uint32_t(x) & (q - 1);
        } else {
            constexpr uint64_t bound = uint64_t(q) + (1ULL << (q_bits - 1)) - 2 + (Maximum >> (q_bits - 1));
            using R                  = std::conditional_t<(bound <= UINT32_MAX), uint32_t, uint64_t>;
            const uint32_t high      = uint32_t(x >> (q_bits - 1));
            const uint32_t quotient  = __umulhi(high, reciprocal);
            const R r                = R(x) - R(quotient) * R(q);
            return uint32_t(r >= q ? r - q : r);
        }
    }
};

template <unsigned N, bool COMPLEX> struct value { uintx<N> re; };
template <unsigned N> struct value<N, true> { uintx<N> re, im; };

template <class Moduli, unsigned FIRST, unsigned END> struct node;

template <class Left, class Right> struct merge_constants {
    static constexpr unsigned words   = Left::words > Right::words ? Left::words : Right::words;
    static constexpr auto p           = resize<words>(Right::modulus);
    static constexpr uint32_t neg_inv = negative_inverse32(p.word[0]);

    static constexpr auto inverse = constant_inverse(constant_value(Left::modulus), constant_value(p));

    static constexpr uint64_t shoup_mu = [] {
        if constexpr (words <= 2) {
            return uint64_t((inverse << (32 * words)) / constant_value(p));
        } else {
            return uint64_t(0);
        }
    }();

    static constexpr auto k = [] {
        constexpr auto modulus = constant_value(p);
        auto out               = inverse;
        for (unsigned i = 0; i < 32 * words; ++i) {
            out = (out * 2) % modulus;
        }
        uintx<words> packed{};
        for (unsigned i = 0; i < words; ++i) {
            packed.word[i] = uint32_t(out >> (32 * i));
        }
        return packed;
    }();

    static constexpr constant_wide bias_value =
        ((constant_value(Left::modulus) + constant_value(p) - 2) / constant_value(p)) * constant_value(p);
    static constexpr bool bounded_difference =
        ((bias_value + constant_value(p) - 1) >> (32 * words - 1)) <= 1;
    static constexpr auto difference_bias = [] {
        uintx<words> out{};
        for (unsigned i = 0; i < words; ++i) out.word[i] = uint32_t(bias_value >> (32 * i));
        return out;
    }();

    template <unsigned OUT>
    __device__ __forceinline__ static uintx<OUT> combine(
        const uintx<Left::words> &a,
        const uintx<Right::words> &b //
    ) {
        constexpr auto right_modulus = p;
        constexpr auto coefficient   = k;
        constexpr auto left_modulus  = Left::modulus;
        auto difference              = resize<words>(b);
        if constexpr (bounded_difference) {
            constexpr auto bias = difference_bias;
            add(difference, bias);
        }
        const uint32_t borrow = sub(difference, resize<words>(a));

        uintx<words> t{};
        if constexpr (words <= 2 && p.word[words - 1] < 0x80000000U) {
            // Shoup: 0 <= x*c - floor(x*mu/R)*p < 2*p < R.
            using U             = std::conditional_t<words == 1, uint32_t, uint64_t>;
            using V             = std::conditional_t<(Right::bit_count < 32), uint32_t, U>;
            constexpr U mu      = U(shoup_mu);
            constexpr V c       = V(inverse);
            constexpr V modulus = V(pack(right_modulus));
            const U x           = U(pack(difference));
            U q;
            if constexpr (words == 1) {
                q = __umulhi(x, mu);
            } else {
                q = __umul64hi(x, mu);
            }
            const V r = V(x) * c - V(q) * modulus;
            unpack<words>(r >= modulus ? r - modulus : r, t);
        } else {
            t = montgomery<words, (p.word[words - 1] < 0x80000000U)>(
                difference, coefficient, right_modulus, neg_inv);
        }

        if constexpr (!bounded_difference) {
            uintx<words> correction{};
#pragma unroll
            for (unsigned i = 0; i < words; ++i) {
                correction.word[i] = coefficient.word[i] & (0U - borrow);
            }
            t = sub_mod(t, correction, right_modulus);
        }

        auto out = multiply<OUT>(left_modulus, resize<Right::words>(t));
        add(out, resize<OUT>(a));
        return out;
    }
};

template <class Moduli, unsigned FIRST, unsigned END> struct node {
    static constexpr auto full_modulus  = product<Moduli>(FIRST, END);
    static constexpr unsigned bit_count = bits(full_modulus);
    static constexpr unsigned words     = (bit_count + 31) / 32;
    static constexpr auto modulus       = resize<words>(full_modulus);
    static constexpr unsigned split     = words > 1 ? split_at<Moduli>(FIRST, END) : FIRST;

    template <unsigned I, class Input>
    __device__ __forceinline__ static void leaf_accumulate(const Input &in, int64_t &re, int64_t &im) {
        using C                       = leaf_constants<Moduli, FIRST, END>;
        const auto residue            = in.template load<I>();
        constexpr int32_t coefficient = C::coefficients.coefficient[I - FIRST];
        re += int64_t(residue.x) * coefficient;
        if constexpr (Input::complex) {
            im += int64_t(residue.y) * coefficient;
        }
        if constexpr (I + 1 < END) {
            leaf_accumulate<I + 1>(in, re, im);
        }
    }

    template <class Input>
    __device__ __forceinline__ static value<words, Input::complex> reconstruct(const Input &in) {
        value<words, Input::complex> out{};
        if constexpr (END == FIRST + 1) {
            const auto residue  = in.template load<FIRST>();
            constexpr int32_t p = int32_t(modulus.word[0]);
            if constexpr (std::is_unsigned_v<decltype(residue.x)>) {
                out.re.word[0] = residue.x;
            } else {
                out.re.word[0] = uint32_t(residue.x < 0 ? residue.x + p : residue.x);
            }
            if constexpr (Input::complex) {
                if constexpr (std::is_unsigned_v<decltype(residue.y)>) {
                    out.im.word[0] = residue.y;
                } else {
                    out.im.word[0] = uint32_t(residue.y < 0 ? residue.y + p : residue.y);
                }
            }
        } else if constexpr (words == 1) {
            using C    = leaf_constants<Moduli, FIRST, END>;
            int64_t re = int64_t(C::bias), im = int64_t(C::bias);
            leaf_accumulate<FIRST>(in, re, im);
            out.re.word[0] = C::reduce(uint64_t(re));
            if constexpr (Input::complex) {
                out.im.word[0] = C::reduce(uint64_t(im));
            }
        } else {
            using A = node<Moduli, FIRST, split>;
            using B = node<Moduli, split, END>;

            constexpr bool swap = (B::modulus.word[0] & 1U) == 0;
            using L             = std::conditional_t<swap, B, A>;
            using R             = std::conditional_t<swap, A, B>;
            using M             = merge_constants<L, R>;

            const auto left  = L::reconstruct(in);
            const auto right = R::reconstruct(in);
            out.re           = M::template combine<words>(left.re, right.re);
            if constexpr (Input::complex) {
                out.im = M::template combine<words>(left.im, right.im);
            }
        }
        return out;
    }
};

} // namespace gemmul8::undo_scaling::crt
