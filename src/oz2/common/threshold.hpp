#pragma once
#include "include.hpp"

namespace gemmul8::common {

//------------------------------
// Iteration threshold for modular reduction
// Used to decide mod implementation based on num_moduli
//------------------------------
template <Backend BACKEND = Backend::INT8, bool COMPLEX = false> struct threshold;
template <> struct threshold<Backend::INT8, false> {
    static constexpr int P_is_double = 6;
    static constexpr int S           = 7;
    static constexpr int M           = 15;
};
template <> struct threshold<Backend::INT8, true> {
    static constexpr int P_is_double = 6;
    static constexpr int S           = 7;
    static constexpr int M           = 16;
};
template <> struct threshold<Backend::FP8, false> {
    static constexpr int P_is_double = 4;
    static constexpr int S           = 5;
    static constexpr int M           = 11;
};
template <> struct threshold<Backend::FP8, true> {
    static constexpr int P_is_double = 5;
    static constexpr int S           = 6;
    static constexpr int M           = 11;
};

} // namespace gemmul8::common
