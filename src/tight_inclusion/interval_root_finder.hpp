// A root finder using interval arithmetic.
#pragma once

#include <tight_inclusion/types.hpp>
#include <tight_inclusion/interval.hpp>
#include <tight_inclusion/config.hpp>
#include <tight_inclusion/types.hpp>

#include <vector>
#include <optional>

namespace ticcd {

    /// @brief Perform interval root finding for CCD using BFS.
    /// @param[in] a_t0 Vertex a at t=0
    /// @param[in] b_t0 Vertex b at t=0
    /// @param[in] c_t0 Vertex c at t=0
    /// @param[in] d_t0 Vertex d at t=0
    /// @param[in] a_t1 Vertex a at t=1
    /// @param[in] b_t1 Vertex b at t=1
    /// @param[in] c_t1 Vertex c at t=1
    /// @param[in] d_t1 Vertex d at t=1
    /// @param[in] tol The tolerance of the interval.
    /// @param[in] co_domain_tolerance The tolerance of the co-domain.
    /// @param[in] err The maximum error on each axis when calculating the vertices, err, aka, filter.
    /// @param[in] ms The minimum separation.
    /// @param[in] max_time The maximum time to check.
    /// @param[in] max_itr The maximum number of iterations.
    /// @tparam is_vertex_face Whether to check vertex-face or edge-edge collision.
    /// @return True if there is a root (collision), false otherwise.
    template <bool is_vertex_face>
    std::optional<CCDResult> interval_root_finder_BFS(
        const Vector3 &a_t0,
        const Vector3 &b_t0,
        const Vector3 &c_t0,
        const Vector3 &d_t0,
        const Vector3 &a_t1,
        const Vector3 &b_t1,
        const Vector3 &c_t1,
        const Vector3 &d_t1,
        const Array3 &tol,
        Scalar co_domain_tolerance,
        const Array3 &err,
        Scalar ms,
        Scalar max_time,
        long max_itr);

    /// @brief Perform interval root finding for edge-edge CCD using BFS.
    /// @param[in] ea0_t0 The start position of the first vertex of the first edge.
    /// @param[in] ea1_t0 The start position of the second vertex of the first edge.
    /// @param[in] eb0_t0 The start position of the first vertex of the second edge.
    /// @param[in] eb1_t0 The start position of the second vertex of the second edge.
    /// @param[in] ea0_t1 The end position of the first vertex of the first edge.
    /// @param[in] ea1_t1 The end position of the second vertex of the first edge.
    /// @param[in] eb0_t1 The end position of the first vertex of the second edge.
    /// @param[in] eb1_t1 The end position of the second vertex of the second edge.
    /// @param[in] tol The tolerance of the interval.
    /// @param[in] co_domain_tolerance The tolerance of the co-domain.
    /// @param[in] err The maximum error on each axis when calculating the vertices, err, aka, filter.
    /// @param[in] ms The minimum separation.
    /// @param[in] max_time The maximum time to check.
    /// @param[in] max_itr The maximum number of iterations.
    /// @return True if there is a root (collision), false otherwise.
    std::optional<CCDResult> edge_edge_interval_root_finder_BFS(
        const Vector3 &ea0_t0,
        const Vector3 &ea1_t0,
        const Vector3 &eb0_t0,
        const Vector3 &eb1_t0,
        const Vector3 &ea0_t1,
        const Vector3 &ea1_t1,
        const Vector3 &eb0_t1,
        const Vector3 &eb1_t1,
        const Array3 &tol,
        Scalar co_domain_tolerance,
        const Array3 &err,
        Scalar ms,
        Scalar max_time,
        long max_itr);

    /// @brief Perform interval root finding for vertex-face CCD using BFS.
    /// @param[in] v_t0  The start position of the vertex.
    /// @param[in] f0_t0 The start position of the first vertex of the face.
    /// @param[in] f1_t0 The start position of the second vertex of the face.
    /// @param[in] f2_t0 The start position of the third vertex of the face.
    /// @param[in] v_t1  The end position of the vertex.
    /// @param[in] f0_t1 The end position of the first vertex of the face.
    /// @param[in] f1_t1 The end position of the second vertex of the face.
    /// @param[in] f2_t1 The end position of the third vertex of the face.
    /// @param[in] tol The tolerance of the interval.
    /// @param[in] co_domain_tolerance The tolerance of the co-domain.
    /// @param[in] err The maximum error on each axis when calculating the vertices, err, aka, filter.
    /// @param[in] ms The minimum separation.
    /// @param[in] max_time The maximum time to check.
    /// @param[in] max_itr The maximum number of iterations.
    /// @return True if there is a root (collision), false otherwise.
    std::optional<CCDResult> vertex_face_interval_root_finder_BFS(
        const Vector3 &v_t0,
        const Vector3 &f0_t0,
        const Vector3 &f1_t0,
        const Vector3 &f2_t0,
        const Vector3 &v_t1,
        const Vector3 &f0_t1,
        const Vector3 &f1_t1,
        const Vector3 &f2_t1,
        const Array3 &tol,
        Scalar co_domain_tolerance,
        // this is the maximum error on each axis when calculating the vertices, err, aka, filter
        Array3 &err,
        Scalar ms,
        Scalar max_time,
        long max_itr);

    // calculate the sign of f. dim is the dimension we are evaluating.
    template <typename T>
    T function_f_ee(
        const NumCCD &tpara,
        const NumCCD &upara,
        const NumCCD &vpara,
        const int dim,
        const Vector3 &ea0_t0,
        const Vector3 &ea1_t0,
        const Vector3 &eb0_t0,
        const Vector3 &eb1_t0,
        const Vector3 &ea0_t1,
        const Vector3 &ea1_t1,
        const Vector3 &eb0_t1,
        const Vector3 &eb1_t1);

    template <typename T>
    T function_f_vf(
        const NumCCD &tpara,
        const NumCCD &upara,
        const NumCCD &vpara,
        int dim,
        const Vector3 &v_t0,
        const Vector3 &f0_t0,
        const Vector3 &f1_t0,
        const Vector3 &f2_t0,
        const Vector3 &v_t1,
        const Vector3 &f0_t1,
        const Vector3 &f1_t1,
        const Vector3 &f2_t1);

    void print_times();

    // get the filter of ccd. the inputs are the vertices of the bounding box of the simulation scene
    Array3 get_numerical_error(
        const std::vector<Vector3> &vertices,
        bool is_vertex_face,
        bool using_minimum_separation);

    // get the filter of ccd. the inputs are the vertices of the bounding box of the simulation scene
    Array3 get_numerical_error(
        const Vector3 &abs_max,
        bool is_vertex_face,
        bool using_minimum_separation);

} // namespace ticcd
