#include "ccd.hpp"

#include <tight_inclusion/interval_root_finder.hpp>
#include <tight_inclusion/timer.hpp>
#include <tight_inclusion/logger.hpp>

#include <vector>
#include <optional>
#include <cassert>

namespace ticcd {

#ifdef TIGHT_INCLUSION_USE_MAX_ABS_TOL
    static constexpr Scalar CCD_MAX_TIME_TOL = 1e-3;
    static constexpr Scalar CCD_MAX_COORD_TOL = 1e-2;
#else
    static constexpr Scalar CCD_MAX_TIME_TOL =
        std::numeric_limits<double>::infinity();
    static constexpr Scalar CCD_MAX_COORD_TOL =
        std::numeric_limits<double>::infinity();
#endif

    inline std::array<Vector3, 2> bbd_4_pts(
        const Vector3 &p0,
        const Vector3 &p1,
        const Vector3 &p2,
        const Vector3 &p3)
    {
        return {
            {p0.cwiseMin(p1).cwiseMin(p2).cwiseMin(p3),
             p0.cwiseMax(p1).cwiseMax(p2).cwiseMax(p3)}};
    }

    // calculate maximum x, y and z diff
    Scalar get_max_axis_diff(
        const std::array<Vector3, 2> &b1, const std::array<Vector3, 2> &b2)
    {
        return std::max({
            (b1[1] - b1[0]).maxCoeff(),
            (b2[1] - b2[0]).maxCoeff(),
            (b2[0] - b1[1]).cwiseAbs().maxCoeff(),
            (b1[0] - b2[1]).cwiseAbs().maxCoeff(),
        });
    }

    inline Scalar max_linf_4(
        const Vector3 &p1,
        const Vector3 &p2,
        const Vector3 &p3,
        const Vector3 &p4,
        const Vector3 &p1e,
        const Vector3 &p2e,
        const Vector3 &p3e,
        const Vector3 &p4e)
    {
        return std::max(
            {(p1e - p1).lpNorm<Eigen::Infinity>(),
             (p2e - p2).lpNorm<Eigen::Infinity>(),
             (p3e - p3).lpNorm<Eigen::Infinity>(),
             (p4e - p4).lpNorm<Eigen::Infinity>()});
    }

    /// @brief Clamp a/b to [-∞, max_val]
    /// @param a numerator
    /// @param b denominator
    /// @param max_val
    /// @return a/b if b != 0, max_val if b == 0
    inline Scalar clamp_div(Scalar a, Scalar b, Scalar max_val)
    {
        if (b == 0) {
            return max_val;
        } else {
            return std::min(a / b, max_val);
        }
    }

    Array3 compute_vertex_face_tolerances(
        const Vector3 &v_t0,
        const Vector3 &f0_t0,
        const Vector3 &f1_t0,
        const Vector3 &f2_t0,
        const Vector3 &v_t1,
        const Vector3 &f0_t1,
        const Vector3 &f1_t1,
        const Vector3 &f2_t1,
        Scalar distance_tolerance)
    {
        const Vector3 p000 = v_t0 - f0_t0;
        const Vector3 p001 = v_t0 - f2_t0;
        const Vector3 p011 = v_t0 - (f1_t0 + f2_t0 - f0_t0);
        const Vector3 p010 = v_t0 - f1_t0;
        const Vector3 p100 = v_t1 - f0_t1;
        const Vector3 p101 = v_t1 - f2_t1;
        const Vector3 p111 = v_t1 - (f1_t1 + f2_t1 - f0_t1);
        const Vector3 p110 = v_t1 - f1_t1;

        const Scalar dl =
            3 * max_linf_4(p000, p001, p011, p010, p100, p101, p111, p110);
        const Scalar edge0_length =
            3 * max_linf_4(p000, p100, p101, p001, p010, p110, p111, p011);
        const Scalar edge1_length =
            3 * max_linf_4(p000, p100, p110, p010, p001, p101, p111, p011);

        return Array3(
            clamp_div(distance_tolerance, dl, CCD_MAX_TIME_TOL),
            clamp_div(distance_tolerance, edge0_length, CCD_MAX_COORD_TOL),
            clamp_div(distance_tolerance, edge1_length, CCD_MAX_COORD_TOL));
    }

    Array3 compute_edge_edge_tolerances(
        const Vector3 &ea0_t0,
        const Vector3 &ea1_t0,
        const Vector3 &eb0_t0,
        const Vector3 &eb1_t0,
        const Vector3 &ea0_t1,
        const Vector3 &ea1_t1,
        const Vector3 &eb0_t1,
        const Vector3 &eb1_t1,
        Scalar distance_tolerance)
    {

        const Vector3 p000 = ea0_t0 - eb0_t0;
        const Vector3 p001 = ea0_t0 - eb1_t0;
        const Vector3 p010 = ea1_t0 - eb0_t0;
        const Vector3 p011 = ea1_t0 - eb1_t0;
        const Vector3 p100 = ea0_t1 - eb0_t1;
        const Vector3 p101 = ea0_t1 - eb1_t1;
        const Vector3 p110 = ea1_t1 - eb0_t1;
        const Vector3 p111 = ea1_t1 - eb1_t1;

        const Scalar dl =
            3 * max_linf_4(p000, p001, p011, p010, p100, p101, p111, p110);
        const Scalar edge0_length =
            3 * max_linf_4(p000, p100, p101, p001, p010, p110, p111, p011);
        const Scalar edge1_length =
            3 * max_linf_4(p000, p100, p110, p010, p001, p101, p111, p011);

        return Array3(
            clamp_div(distance_tolerance, dl, CCD_MAX_TIME_TOL),
            clamp_div(distance_tolerance, edge0_length, CCD_MAX_COORD_TOL),
            clamp_div(distance_tolerance, edge1_length, CCD_MAX_COORD_TOL));
    }

    template <bool is_vertex_face>
    std::optional<Collision>
    CCD(const Vector3 &a_t0,
        const Vector3 &b_t0,
        const Vector3 &c_t0,
        const Vector3 &d_t0,
        const Vector3 &a_t1,
        const Vector3 &b_t1,
        const Vector3 &c_t1,
        const Vector3 &d_t1,
        const Array3 &err,
        Scalar ms,
        Scalar tolerance,
        Scalar t_max,
        long max_itr,
        bool no_zero_toi)
    {

        Array3 tol;
        if constexpr (is_vertex_face) {
            tol = compute_vertex_face_tolerances(
                a_t0, b_t0, c_t0, d_t0, a_t1, b_t1, c_t1, d_t1, tolerance);
        } else {
            tol = compute_edge_edge_tolerances(
                a_t0, b_t0, c_t0, d_t0, a_t1, b_t1, c_t1, d_t1, tolerance);
        }

        //////////////////////////////////////////////////////////
        // this should be the error of the whole mesh.
        // but if err = {-1, -1, -1}, this is a special value indicating
        // we need to compute err ourself.
        Array3 true_err;
        if (err(0) == -1.0f && err(1) == -1.0f && err(2) == -1.0f) {
            bool has_minimal_separation = (ms > 0.0f);
            true_err = get_numerical_error(
                std::vector<Vector3>{
                    {a_t0, b_t0, c_t0, d_t0, a_t1, b_t1, c_t1, d_t1}},
                is_vertex_face, has_minimal_separation);
        } else {
            true_err = err;
        }

        if (!no_zero_toi) {
            return interval_root_finder_BFS<is_vertex_face>(
                a_t0, b_t0, c_t0, d_t0, a_t1, b_t1, c_t1, d_t1, tol, tolerance,
                true_err, ms, t_max, max_itr);
        }

        // strategies for dealing with zero toi:
        // 1. if reach max_iter, shrink t_max.
        // 2. if ms is too large, shrink ms.
        // 3. if tolerance is too large, shrink tolerance.
        constexpr int MAX_NO_ZERO_TOI_ITER = 4;
        for (int i = 0; i < MAX_NO_ZERO_TOI_ITER; ++i) {
            assert(t_max >= 0.0f && t_max <= 1.0f);

            auto collision = interval_root_finder_BFS<is_vertex_face>(
                a_t0, b_t0, c_t0, d_t0, a_t1, b_t1, c_t1, d_t1, tol, tolerance,
                true_err, ms, t_max, max_itr);

            bool is_zero_toi = (collision && collision->t(0) == 0.0f);
            if (!is_zero_toi) {
                return collision;
            }

            // case 1
            if (collision->tolerance != tolerance) {
                logger().debug("toi refine strategy 1: shrink t max");
                t_max = collision->t(1);
            }
            //case 2
            else if (10.0f * tolerance < ms) {
                logger().debug("toi refine strategy 2: shrink ms");
                ms *= 0.5f;
            }
            //case 3
            else {
                logger().debug("toi refine strategy 3: shrink tolerance");
                tolerance *= 0.5f;

                if constexpr (is_vertex_face) {
                    tol = compute_vertex_face_tolerances(
                        a_t0, b_t0, c_t0, d_t0, a_t1, b_t1, c_t1, d_t1,
                        tolerance);
                } else {
                    tol = compute_edge_edge_tolerances(
                        a_t0, b_t0, c_t0, d_t0, a_t1, b_t1, c_t1, d_t1,
                        tolerance);
                }
            }
        }

        // if after max no zero toi iterations the toi is still zero,
        // ignore this collision.
        logger().debug("toi refine fail, ignore potential collision");
        return std::nullopt;
    }

    std::optional<Collision> edgeEdgeCCD(
        const Vector3 &ea0_t0,
        const Vector3 &ea1_t0,
        const Vector3 &eb0_t0,
        const Vector3 &eb1_t0,
        const Vector3 &ea0_t1,
        const Vector3 &ea1_t1,
        const Vector3 &eb0_t1,
        const Vector3 &eb1_t1,
        const Array3 &err,
        Scalar ms,
        Scalar tolerance,
        Scalar t_max,
        long max_itr,
        bool no_zero_toi)
    {
        return CCD</*is_vertex_face=*/false>(
            ea0_t0, ea1_t0, eb0_t0, eb1_t0, ea0_t1, ea1_t1, eb0_t1, eb1_t1, err,
            ms, tolerance, t_max, max_itr, no_zero_toi);
    }

    std::optional<Collision> vertexFaceCCD(
        const Vector3 &v_t0,
        const Vector3 &f0_t0,
        const Vector3 &f1_t0,
        const Vector3 &f2_t0,
        const Vector3 &v_t1,
        const Vector3 &f0_t1,
        const Vector3 &f1_t1,
        const Vector3 &f2_t1,
        const Array3 &err,
        Scalar ms,
        Scalar tolerance,
        Scalar t_max,
        long max_itr,
        bool no_zero_toi)
    {
        return CCD</*is_vertex_face=*/true>(
            v_t0, f0_t0, f1_t0, f2_t0, v_t1, f0_t1, f1_t1, f2_t1, err, ms,
            tolerance, t_max, max_itr, no_zero_toi);
    }

} // namespace ticcd
