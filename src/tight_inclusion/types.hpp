#pragma once

#include <tight_inclusion/config.hpp>

#include <Eigen/Core>

namespace ticcd {
#ifdef TIGHT_INCLUSION_WITH_DOUBLE_PRECISION
    typedef double Scalar;
#else
    typedef float Scalar;
#endif
    typedef Eigen::Matrix<Scalar, 3, 1> Vector3;
    typedef Eigen::Array<Scalar, 2, 1> Array2;
    typedef Eigen::Array<Scalar, 3, 1> Array3;
    typedef Eigen::Array<Scalar, 6, 1> Array6;
    typedef Eigen::Array<Scalar, 8, 1> Array8;

    struct CCDResult {
        Array2 t;
        Array2 u;
        Array2 v;
        Scalar tolerance;
        bool use_small_ms;
        Array2 small_ms_t;
    };
} // namespace ticcd
