// Time-of-impact computation for rigid bodies with angular trajectories.
#pragma once

#include <tight_inclusion/interval.hpp>
#include <tight_inclusion/types.hpp>
#include <tight_inclusion/config.hpp>

#include <optional>

namespace ticcd {
    static constexpr bool DEFAULT_NO_ZERO_TOI = false;
    static constexpr Scalar DEFAULT_CCD_DISTANCE_TOL = 1e-6;

    /// @brief This function can give you the answer of continuous collision detection with minimum
    /// separation, and the earliest collision time if collision happens.
    ///
    /// @param[in] ea0_t0 The start position of the first vertex of the first edge.
    /// @param[in] ea1_t0 The start position of the second vertex of the first edge.
    /// @param[in] eb0_t0 The start position of the first vertex of the second edge.
    /// @param[in] eb1_t0 The start position of the second vertex of the second edge.
    /// @param[in] ea0_t1 The end position of the first vertex of the first edge.
    /// @param[in] ea1_t1 The end position of the second vertex of the first edge.
    /// @param[in] eb0_t1 The end position of the first vertex of the second edge.
    /// @param[in] eb1_t1 The end position of the second vertex of the second edge.
    /// @param[in] err The filters calculated using the bounding box of the simulation scene.
    ///                If you are checking a single query without a scene, please set it as {-1,-1,-1}.
    /// @param[in] ms The minimum separation. should set: ms < max(abs(x),1), ms < max(abs(y),1), ms < max(abs(z),1) of the QUERY (NOT THE SCENE!).
    /// @param[in] tolerance A user - input solving precision. We suggest to use 1e-6.
    /// @param[in] t_max The upper bound of the time interval [0,t_max] to be checked. 0<=t_max<=1
    /// @param[in] max_itr A user-defined value to terminate the algorithm earlier, and return a result under current
    ///                    precision. please set max_itr either a big number like 1e7, or -1 which means it will not be terminated
    ///                    earlier and the precision will be user-defined precision -- tolerance.
    /// @param[in] no_zero_toi Refine further if a zero toi is produced (assumes not initially in contact).
    /// @return True if there is a collision, false otherwise.
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
        bool no_zero_toi = DEFAULT_NO_ZERO_TOI);

    /// This function can give you the answer of continuous collision detection with minimum
    /// separation, and the earliest collision time if collision happens.
    ///
    /// @param[in] v_t0  The start position of the vertex.
    /// @param[in] f0_t0 The start position of the first vertex of the face.
    /// @param[in] f1_t0 The start position of the second vertex of the face.
    /// @param[in] f2_t0 The start position of the third vertex of the face.
    /// @param[in] v_t1  The end position of the vertex.
    /// @param[in] f0_t1 The end position of the first vertex of the face.
    /// @param[in] f1_t1 The end position of the second vertex of the face.
    /// @param[in] f2_t1 The end position of the third vertex of the face.
    /// @param[in] err The filters calculated using the bounding box of the simulation scene.
    ///                If you are checking a single query without a scene, please set it as {-1,-1,-1}.
    /// @param[in] ms The minimum separation. should set: ms < max(abs(x),1), ms < max(abs(y),1), ms < max(abs(z),1) of the QUERY (NOT THE SCENE!).
    /// @param[in] tolerance A user - input solving precision. We suggest to use 1e-6.
    /// @param[in] t_max The upper bound of the time interval [0,t_max] to be checked. 0<=t_max<=1
    /// @param[in] max_itr A user-defined value to terminate the algorithm earlier, and return a result under current
    ///                    precision. please set max_itr either a big number like 1e7, or -1 which means it will not be terminated
    ///                    earlier and the precision will be user-defined precision -- tolerance.
    /// @param[in] no_zero_toi Refine further if a zero toi is produced (assumes not initially in contact).
    /// @return True if there is a collision, false otherwise.
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
        bool no_zero_toi = DEFAULT_NO_ZERO_TOI);

    Array3 compute_vertex_face_tolerances(
        const Vector3 &v_t0,
        const Vector3 &f0_t0,
        const Vector3 &f1_t0,
        const Vector3 &f2_t0,
        const Vector3 &v_t1,
        const Vector3 &f0_t1,
        const Vector3 &f1_t1,
        const Vector3 &f2_t1,
        Scalar distance_tolerance = DEFAULT_CCD_DISTANCE_TOL);

    Array3 compute_edge_edge_tolerances(
        const Vector3 &ea0_t0,
        const Vector3 &ea1_t0,
        const Vector3 &eb0_t0,
        const Vector3 &eb1_t0,
        const Vector3 &ea0_t1,
        const Vector3 &ea1_t1,
        const Vector3 &eb0_t1,
        const Vector3 &eb1_t1,
        Scalar distance_tolerance = DEFAULT_CCD_DISTANCE_TOL);

    long return_queue_size();

} // namespace ticcd
