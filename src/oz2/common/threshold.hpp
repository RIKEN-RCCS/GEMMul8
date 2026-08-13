#pragma once
#include "include.hpp"

namespace gemmul8::common {

//------------------------------
// Iteration threshold for modular reduction
// Used to decide mod implementation based on num_moduli
//------------------------------
template <Backend BACKEND = Backend::INT8> struct threshold;
template <> struct threshold<Backend::INT8> {
    static constexpr int P_is_double = 6;
    static constexpr int S           = 7;
    static constexpr int M           = 15;
    static constexpr int L           = 22; // not used
};
template <> struct threshold<Backend::FP8> {
    static constexpr int P_is_double = 5;
    static constexpr int S           = 5;
    static constexpr int M           = 11;
    static constexpr int L           = 17; // not used
};

} // namespace gemmul8::common
