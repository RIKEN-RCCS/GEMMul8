#pragma once
#include "matmult.hpp"
#include "../mod/crt_group.hpp"
#include "../common/timer.hpp"

namespace gemmul8::oz2::core {

template <Func F, Backend B, unsigned N, bool Complex,
          cublasFillMode_t UA, cublasFillMode_t UB, cublasFillMode_t UC,
          class Layout, class Multiply, size_t... G>
inline const common::hi_t<B> *grouped_products_impl(
    cudaStream_t stream, common::Handle_t &handle,
    cublasOperation_t opA, cublasOperation_t opB,
    size_t ldc, unsigned n, size_t k,
    int8_t *work,
    common::Timer<N> &timer,
    Multiply &multiply,
    std::index_sequence<G...> //
) {
    using Plan               = common::crt_group_plan<B, N, Complex>;
    using Hi                 = common::hi_t<B>;
    constexpr size_t blas    = size_t(32) << 20;
    const size_t sizeC       = ldc * n;
    const size_t group_bytes = sizeof(uint32_t) * (Complex ? 2U : 1U) * sizeC;

    const Hi *last   = nullptr;
    const auto group = [&]<size_t J>() {
        constexpr unsigned first = Plan::groups.first[J];
        constexpr unsigned end   = Plan::groups.first[J + 1];
        int8_t *base             = work + J * group_bytes;
        const size_t pointers    = Layout::pointer_count(N, first, end - first);
        const size_t ptr_gap     = Layout::pointer_bytes(N, first, end - first);
        constexpr bool fused     = J + 1 == Plan::count && F != Func::herk && F != Func::herkx;
        auto *hi                 = reinterpret_cast<Hi *>(base + std::max(group_bytes, ptr_gap + blas));
        handle.Aarray            = pointers ? reinterpret_cast<void **>(base) : nullptr;
        handle.Barray            = pointers ? reinterpret_cast<void **>(base) + pointers : nullptr;
        handle.Carray            = pointers ? reinterpret_cast<void **>(base) + 2 * pointers : nullptr;
        handle.workspace         = base + ptr_gap;

        for (unsigned i = first; i < end;) {
            const unsigned count = (handle.arch == 121 || (handle.arch == 90 && n > 2048))
                                       ? 1U
                                       : limit_batch_for_k<B, Complex>(i, end - i, k, N);
            configure_matprod_k_blocking<B, Complex>(handle, i, k);
            auto C_hi = common::make_matptr<Hi, Complex, Complex ? 2U : 1U>(
                hi + sizeC * Layout::planes(N, first, i - first), sizeC * Layout::products(N, i));
            if constexpr (F == Func::herk) {
                C_hi.ptr1 = nullptr;
            }
            multiply(i, count, C_hi);
            i += count;
        }
        const unsigned gid = timer.begin_group();
        timer.record(timer.mm_end(gid), stream);

        if constexpr (fused) {
            last = hi;
            timer.record(timer.mod_hi2mid_end(gid), stream);
            return;
        }

        const auto reduce = [&]<cublasFillMode_t U>() {
            if constexpr (F == Func::herk) {
                if (opA == CUBLAS_OP_N) {
                    mod::mod_hi2group<B, N, Complex, first, end, U, true, true>(stream, n, ldc, hi, base);
                } else {
                    mod::mod_hi2group<B, N, Complex, first, end, U, true, false>(stream, n, ldc, hi, base);
                }
            } else {
                mod::mod_hi2group<B, N, Complex, first, end, U>(stream, n, ldc, hi, base);
            }
        };

        if constexpr (F == Func::syr2k || F == Func::her2k) {
            reduce.template operator()<CUBLAS_FILL_MODE_FULL>();
        } else if constexpr (F == Func::trtrmm) {
            if (opA == CUBLAS_OP_N) {
                if (opB == CUBLAS_OP_N) reduce.template operator()<(UA == UB ? UA : UC)>();
                else reduce.template operator()<(UA == flip_uplo<UB> ? UA : UC)>();
            } else {
                if (opB == CUBLAS_OP_N) reduce.template operator()<(flip_uplo<UA> == UB ? flip_uplo<UA> : UC)>();
                else reduce.template operator()<(UA == UB ? flip_uplo<UA> : UC)>();
            }
        } else {
            reduce.template operator()<UC>();
        }

        timer.record(timer.mod_hi2mid_end(gid), stream);
    };

    (group.template operator()<G>(), ...);
    return last;
}

template <Func F, Backend B, unsigned N, bool Complex,
          cublasFillMode_t UA, cublasFillMode_t UB, cublasFillMode_t UC,
          class Layout, class Multiply>
inline const common::hi_t<B> *grouped_products(
    cudaStream_t stream, common::Handle_t &handle,
    cublasOperation_t opA, cublasOperation_t opB,
    size_t ldc, unsigned n, size_t k,
    int8_t *work,
    common::Timer<N> &timer,
    Multiply multiply //
) {
    return grouped_products_impl<F, B, N, Complex, UA, UB, UC, Layout>(
        stream, handle, opA, opB, ldc, n, k, work, timer, multiply,
        std::make_index_sequence<common::crt_group_plan<B, N, Complex>::count>{});
}

} // namespace gemmul8::oz2::core
