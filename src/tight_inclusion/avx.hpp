#pragma once

#include <tight_inclusion/interval.hpp>
#include "types.hpp"

namespace ticcd {
    void convert_tuv_to_array(
        const Interval3 &itv,
        Array8 &t_up,
        Array8 &t_dw,
        Array8 &u_up,
        Array8 &u_dw,
        Array8 &v_up,
        Array8 &v_dw);

    Array8 function_vf(
        const Scalar &v_t0,
        const Scalar &f0_t0,
        const Scalar &f1_t0,
        const Scalar &f2_t0,
        const Scalar &v_t1,
        const Scalar &f0_t1,
        const Scalar &f1_t1,
        const Scalar &f2_t1,
        const Array8 &t_up,
        const Array8 &t_dw,
        const Array8 &u_up,
        const Array8 &u_dw,
        const Array8 &v_up,
        const Array8 &v_dw);

    Array8 function_ee(
        const Scalar &ea0_t0,
        const Scalar &ea1_t0,
        const Scalar &eb0_t0,
        const Scalar &eb1_t0,
        const Scalar &ea0_t1,
        const Scalar &ea1_t1,
        const Scalar &eb0_t1,
        const Scalar &eb1_t1,
        const Array8 &t_up,
        const Array8 &t_dw,
        const Array8 &u_up,
        const Array8 &u_dw,
        const Array8 &v_up,
        const Array8 &v_dw);

} // namespace ticcd
