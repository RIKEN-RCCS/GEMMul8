#pragma once
#include "include.hpp"

#if !defined(GEMMul8_PROFILE)
    #define GEMMul8_PROFILE 0
#endif

namespace gemmul8::common {

#if GEMMul8_PROFILE

template <unsigned NUM_MODULI>
struct Timer {
    static constexpr unsigned max_groups = NUM_MODULI;
    static constexpr unsigned num_events = 2 * max_groups + 3;

    cudaEvent_t events[num_events]{};
    unsigned num_groups = 0;

    Timer() {
    #pragma unroll
        for (unsigned i = 0; i < num_events; ++i) {
            cudaEventCreate(&events[i]);
        }
    }

    Timer(const Timer &)            = delete;
    Timer &operator=(const Timer &) = delete;

    ~Timer() {
    #pragma unroll
        for (unsigned i = 0; i < num_events; ++i) {
            if (events[i]) cudaEventDestroy(events[i]);
        }
    }

    static constexpr unsigned scaling_begin() { return 0; }
    static constexpr unsigned scaling_end() { return 1; }

    static constexpr unsigned mm_end(const unsigned g) {
        return 2 + 2 * g;
    }

    static constexpr unsigned mod_hi2mid_end(const unsigned g) {
        return 3 + 2 * g;
    }

    static constexpr unsigned undo_scaling_end() {
        return 2 * max_groups + 2;
    }

    unsigned begin_group() {
        const unsigned g = num_groups;
        ++num_groups;
        return g;
    }

    unsigned mm_begin(const unsigned g) const {
        return (g == 0) ? scaling_end() : mod_hi2mid_end(g - 1);
    }

    static constexpr unsigned mod_begin(const unsigned g) {
        return mm_end(g);
    }

    unsigned undo_scaling_begin() const {
        return (num_groups == 0) ? scaling_end() : mod_hi2mid_end(num_groups - 1);
    }

    void record(const unsigned idx, cudaStream_t stream) {
        cudaEventRecord(events[idx], stream);
    }

    static double elapsed_s(cudaEvent_t begin, cudaEvent_t end) {
        float ms = 0.0f;
        cudaEventElapsedTime(&ms, begin, end);
        return static_cast<double>(ms) * 1.0e-3;
    }

    void collect(std::vector<double> &timer) {
        static_assert(NUM_MODULI > 0, "NUM_MODULI must be positive.");

        cudaEventSynchronize(events[undo_scaling_end()]);

        timer[0] = elapsed_s(events[scaling_begin()], events[scaling_end()]);
        timer[1] = 0.0;
        timer[2] = 0.0;

        for (unsigned g = 0; g < num_groups; ++g) {
            timer[1] += elapsed_s(events[mm_begin(g)], events[mm_end(g)]);
            timer[2] += elapsed_s(events[mod_begin(g)], events[mod_hi2mid_end(g)]);
        }

        timer[3] = elapsed_s(events[undo_scaling_begin()], events[undo_scaling_end()]);
    }
};

#else

template <unsigned NUM_MODULI>
struct Timer {
    static constexpr unsigned max_groups = NUM_MODULI;
    static constexpr unsigned num_events = 2 * max_groups + 3;

    Timer()                         = default;
    Timer(const Timer &)            = delete;
    Timer &operator=(const Timer &) = delete;
    ~Timer()                        = default;

    static constexpr unsigned scaling_begin() { return 0; }
    static constexpr unsigned scaling_end() { return 0; }
    static constexpr unsigned mm_end(unsigned) { return 0; }
    static constexpr unsigned mod_hi2mid_end(unsigned) { return 0; }
    static constexpr unsigned undo_scaling_end() { return 0; }

    unsigned begin_group() { return 0; }
    unsigned mm_begin(unsigned) const { return 0; }
    static constexpr unsigned mod_begin(unsigned) { return 0; }
    unsigned undo_scaling_begin() const { return 0; }

    void record(unsigned, cudaStream_t) {}
    static double elapsed_s(cudaEvent_t, cudaEvent_t) { return 0.0; }

    void collect(std::vector<double> &timer) {
        static_assert(NUM_MODULI > 0, "NUM_MODULI must be positive.");
        timer[0] = timer[1] = timer[2] = timer[3] = 0.0;
    }
};

#endif

} // namespace gemmul8::common
