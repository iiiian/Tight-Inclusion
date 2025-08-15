#include "avx.hpp"

namespace ticcd {

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
        const Array8 &v_dw)
    {
        const Array8 ea0 = (ea0_t1 - ea0_t0) * t_up / t_dw + ea0_t0;
        const Array8 ea1 = (ea1_t1 - ea1_t0) * t_up / t_dw + ea1_t0;
        const Array8 eb0 = (eb0_t1 - eb0_t0) * t_up / t_dw + eb0_t0;
        const Array8 eb1 = (eb1_t1 - eb1_t0) * t_up / t_dw + eb1_t0;

        const Array8 va = (ea1 - ea0) * u_up / u_dw + ea0;
        const Array8 vb = (eb1 - eb0) * v_up / v_dw + eb0;

        return va - vb;
    }

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
        const Array8 &v_dw)
    {
        const Array8 v = (v_t1 - v_t0) * t_up / t_dw + v_t0;
        const Array8 f0 = (f0_t1 - f0_t0) * t_up / t_dw + f0_t0;
        const Array8 f1 = (f1_t1 - f1_t0) * t_up / t_dw + f1_t0;
        const Array8 f2 = (f2_t1 - f2_t0) * t_up / t_dw + f2_t0;
        const Array8 pt =
            (f1 - f0) * u_up / u_dw + (f2 - f0) * v_up / v_dw + f0;

        return v - pt;
    }

    void convert_tuv_to_array(
        const Interval3 &itv,
        Array8 &t_up,
        Array8 &t_dw,
        Array8 &u_up,
        Array8 &u_dw,
        Array8 &v_up,
        Array8 &v_dw)
    {
        // t order: 0,0,0,0,1,1,1,1
        // u order: 0,0,1,1,0,0,1,1
        // v order: 0,1,0,1,0,1,0,1
        const Scalar t0_up = itv[0].lower.numerator;
        const Scalar t0_dw = itv[0].lower.denominator();
        const Scalar t1_up = itv[0].upper.numerator;
        const Scalar t1_dw = itv[0].upper.denominator();
        const Scalar u0_up = itv[1].lower.numerator;
        const Scalar u0_dw = itv[1].lower.denominator();
        const Scalar u1_up = itv[1].upper.numerator;
        const Scalar u1_dw = itv[1].upper.denominator();
        const Scalar v0_up = itv[2].lower.numerator;
        const Scalar v0_dw = itv[2].lower.denominator();
        const Scalar v1_up = itv[2].upper.numerator;
        const Scalar v1_dw = itv[2].upper.denominator();
        t_up = {t0_up, t0_up, t0_up, t0_up, t1_up, t1_up, t1_up, t1_up};
        t_dw = {t0_dw, t0_dw, t0_dw, t0_dw, t1_dw, t1_dw, t1_dw, t1_dw};
        u_up = {u0_up, u0_up, u1_up, u1_up, u0_up, u0_up, u1_up, u1_up};
        u_dw = {u0_dw, u0_dw, u1_dw, u1_dw, u0_dw, u0_dw, u1_dw, u1_dw};
        v_up = {v0_up, v1_up, v0_up, v1_up, v0_up, v1_up, v0_up, v1_up};
        v_dw = {v0_dw, v1_dw, v0_dw, v1_dw, v0_dw, v1_dw, v0_dw, v1_dw};
    }

} // namespace ticcd
