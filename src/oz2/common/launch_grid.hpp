#pragma once
#include "common.hpp"

namespace gemmul8::common {

template <bool EXTENDED>
__device__ __forceinline__ unsigned block_y() {
    if constexpr (EXTENDED) return blockIdx.y + gridDim.y * blockIdx.z;
    else return blockIdx.y;
}

template <bool EXTENDED>
inline dim3 make_grid_y(dim3 grid) {
    if constexpr (EXTENDED) {
        grid.z = (grid.y - 1U) / 65535U + 1U;
        grid.y = (grid.y - 1U) / grid.z + 1U;
    }
    return grid;
}

template <class Launch>
inline void launch_grid_y(dim3 grid, Launch &&launch) {
    if (grid.y <= 65535U) launch.template operator()<false>(grid);
    else launch.template operator()<true>(make_grid_y<true>(grid));
}

} // namespace gemmul8::common
