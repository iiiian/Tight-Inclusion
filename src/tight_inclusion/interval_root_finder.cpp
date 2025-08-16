// A root finder using interval arithmetic.

#include <tight_inclusion/interval_root_finder.hpp>

#include <tight_inclusion/types.hpp>
#include <tight_inclusion/avx.hpp>
#include <tight_inclusion/logger.hpp>

#include <optional>
#include <queue>
#include <vector>
#include <algorithm>

namespace ticcd {
    template <bool is_vertex_face>
    bool eval_unit_bbox_1d(
        const Vector3 &a_t0,
        const Vector3 &b_t0,
        const Vector3 &c_t0,
        const Vector3 &d_t0,
        const Vector3 &a_t1,
        const Vector3 &b_t1,
        const Vector3 &c_t1,
        const Vector3 &d_t1,
        Scalar eps,
        Scalar ms,
        int dim,
        bool &bbox_in_eps,
        Scalar &tol)
    {
        Scalar minv;
        Scalar maxv;
        if constexpr (is_vertex_face) {
            Array6 A;
            A(0) = a_t0(dim);
            A(1) = a_t0(dim);
            A(2) = a_t0(dim);
            A(3) = a_t1(dim);
            A(4) = a_t1(dim);
            A(5) = a_t1(dim);

            Array6 B;
            B(0) = b_t0(dim);
            B(1) = c_t0(dim);
            B(2) = d_t0(dim);
            B(3) = b_t1(dim);
            B(4) = c_t1(dim);
            B(5) = d_t1(dim);

            const Array6 D = A - B;
            minv = D.minCoeff();
            maxv = D.maxCoeff();
        } else {
            Array8 A;
            A(0) = a_t0(dim);
            A(1) = a_t0(dim);
            A(2) = b_t0(dim);
            A(3) = b_t0(dim);
            A(4) = a_t1(dim);
            A(5) = a_t1(dim);
            A(6) = b_t1(dim);
            A(7) = b_t1(dim);

            Array8 B;
            B(0) = c_t0(dim);
            B(1) = d_t0(dim);
            B(2) = c_t0(dim);
            B(3) = d_t0(dim);
            B(4) = c_t1(dim);
            B(5) = d_t1(dim);
            B(6) = c_t1(dim);
            B(7) = d_t1(dim);

            const Array8 D = A - B;
            minv = D.minCoeff();
            maxv = D.maxCoeff();
        }

        tol = maxv - minv; // this is the real tolerance
        bbox_in_eps = false;
        const Scalar eps_and_ms = eps + ms;

        if (minv > eps_and_ms || maxv < -eps_and_ms) {
            return false;
        }

        if (minv >= -eps_and_ms && maxv <= eps_and_ms) {
            bbox_in_eps = true;
        }

        return true;
    }

    // ** this version can return the true x or y or z tolerance of the co-domain **
    // eps is the interval [-eps,eps] we need to check
    // if [-eps,eps] overlap, return true
    // bbox_in_eps tell us if the box is totally in eps box
    // ms is the minimum seperation
    template <bool is_vertex_face>
    bool eval_bbox_1d(
        Array8 &t_up,
        Array8 &t_dw,
        Array8 &u_up,
        Array8 &u_dw,
        Array8 &v_up,
        Array8 &v_dw,
        const Vector3 &a_t0,
        const Vector3 &b_t0,
        const Vector3 &c_t0,
        const Vector3 &d_t0,
        const Vector3 &a_t1,
        const Vector3 &b_t1,
        const Vector3 &c_t1,
        const Vector3 &d_t1,
        Scalar eps,
        Scalar ms,
        int dim,
        bool &bbox_in_eps,
        Scalar &tol)
    {
        Array8 vs;
        if constexpr (is_vertex_face) {
            vs = function_vf(
                a_t0(dim), b_t0(dim), c_t0(dim), d_t0(dim), a_t1(dim),
                b_t1(dim), c_t1(dim), d_t1(dim), t_up, t_dw, u_up, u_dw, v_up,
                v_dw);
        } else {
            vs = function_ee(
                a_t0(dim), b_t0(dim), c_t0(dim), d_t0(dim), a_t1(dim),
                b_t1(dim), c_t1(dim), d_t1(dim), t_up, t_dw, u_up, u_dw, v_up,
                v_dw);
        }

        Scalar minv = vs.minCoeff();
        Scalar maxv = vs.maxCoeff();

        tol = maxv - minv; // this is the real tolerance
        bbox_in_eps = false;
        const Scalar eps_and_ms = eps + ms;

        if (minv > eps_and_ms || maxv < -eps_and_ms) {
            return false;
        }

        if (minv >= -eps_and_ms && maxv <= eps_and_ms) {
            bbox_in_eps = true;
        }

        return true;
    }

    // ** this version can return the true tolerance of the co-domain **
    // give the result of if the hex overlaps the input eps box around origin
    // use vectorized hex-vertex-solving function for acceleration
    // box_in_eps shows if this hex is totally inside box. if so, no need to do further bisection
    template <bool is_vertex_face, bool is_unit_tuv>
    bool origin_in_bbox_eval(
        const Interval3 &tuv,
        const Vector3 &a_t0,
        const Vector3 &b_t0,
        const Vector3 &c_t0,
        const Vector3 &d_t0,
        const Vector3 &a_t1,
        const Vector3 &b_t1,
        const Vector3 &c_t1,
        const Vector3 &d_t1,
        const Array3 &eps,
        Scalar ms,
        bool &bbox_in_eps,
        Array3 &tolerance)
    {
        bool xyz_bbox_in_eps[3];
        if constexpr (is_unit_tuv) {
            for (int dim = 0; dim < 3; dim++) {
                if (!eval_unit_bbox_1d<is_vertex_face>(
                        a_t0, b_t0, c_t0, d_t0, a_t1, b_t1, c_t1, d_t1,
                        eps(dim), ms, dim, xyz_bbox_in_eps[dim],
                        tolerance(dim))) {
                    return false;
                }
            }
        } else {
            Array8 t_up, t_dw, u_up, u_dw, v_up, v_dw;
            convert_tuv_to_array(tuv, t_up, t_dw, u_up, u_dw, v_up, v_dw);

            for (int dim = 0; dim < 3; dim++) {
                if (!eval_bbox_1d<is_vertex_face>(
                        t_up, t_dw, u_up, u_dw, v_up, v_dw, a_t0, b_t0, c_t0,
                        d_t0, a_t1, b_t1, c_t1, d_t1, eps(dim), ms, dim,
                        xyz_bbox_in_eps[dim], tolerance(dim))) {
                    return false;
                }
            }
        }

        bbox_in_eps =
            xyz_bbox_in_eps[0] && xyz_bbox_in_eps[1] && xyz_bbox_in_eps[2];
        return true;
    }

    // find the largest width/tol dimension that is greater than its tolerance
    int find_next_split(const Array3 &widths, const Array3 &tols)
    {
        Array3 tmp =
            (widths > tols)
                .select(
                    widths / tols, -std::numeric_limits<Scalar>::infinity());
        int max_index;
        tmp.maxCoeff(&max_index);
        return max_index;
    }

    bool split_and_push(
        const Interval3 &tuv,
        int split_i,
        std::function<void(const Interval3 &)> push,
        bool is_vertex_face,
        Scalar t_upper_bound = 1)
    {
        std::pair<Interval, Interval> halves = tuv[split_i].bisect();
        if (halves.first.lower >= halves.first.upper
            || halves.second.lower >= halves.second.upper) {
            logger().error("overflow occured when splitting intervals!");
            return true;
        }

        Interval3 tmp = tuv;

        if (split_i == 0) {
            if (t_upper_bound == 1
                || halves.second.overlaps(0, t_upper_bound)) {
                tmp[split_i] = halves.second;
                push(tmp);
            }
            if (t_upper_bound == 1 || halves.first.overlaps(0, t_upper_bound)) {
                tmp[split_i] = halves.first;
                push(tmp);
            }
        } else if (!is_vertex_face) { // edge uv
            tmp[split_i] = halves.second;
            push(tmp);
            tmp[split_i] = halves.first;
            push(tmp);
        } else {
            assert(is_vertex_face && split_i != 0);
            // u + v ≤ 1
            if (split_i == 1) {
                const Interval &v = tuv[2];
                if (NumCCD::is_sum_leq_1(halves.second.lower, v.lower)) {
                    tmp[split_i] = halves.second;
                    push(tmp);
                }
                if (NumCCD::is_sum_leq_1(halves.first.lower, v.lower)) {
                    tmp[split_i] = halves.first;
                    push(tmp);
                }
            } else if (split_i == 2) {
                const Interval &u = tuv[1];
                if (NumCCD::is_sum_leq_1(u.lower, halves.second.lower)) {
                    tmp[split_i] = halves.second;
                    push(tmp);
                }
                if (NumCCD::is_sum_leq_1(u.lower, halves.first.lower)) {
                    tmp[split_i] = halves.first;
                    push(tmp);
                }
            }
        }
        return false; // no overflow
    }

    // this version cannot give the impact time at t=1, although this collision can
    // be detected at t=0 of the next time step, but still may cause problems in
    // line-search based physical simulation
    template <bool is_vertex_face>
    bool interval_root_finder_DFS(
        const Vector3 &a_t0,
        const Vector3 &b_t0,
        const Vector3 &c_t0,
        const Vector3 &d_t0,
        const Vector3 &a_t1,
        const Vector3 &b_t1,
        const Vector3 &c_t1,
        const Vector3 &d_t1,
        const Array3 &tol,
        const Array3 &err,
        const Scalar ms,
        Scalar &toi)
    {
        auto cmp_time = [](const Interval3 &i1, const Interval3 &i2) {
            return i1[0].lower >= i2[0].lower;
        };

        // build interval set [0,1]x[0,1]x[0,1]
        const Interval zero_to_one = Interval(NumCCD(0, 0), NumCCD(1, 0));
        Interval3 iset = {{zero_to_one, zero_to_one, zero_to_one}};

        // Stack of intervals and the last split dimension
        std::priority_queue<
            Interval3, std::vector<Interval3>, decltype(cmp_time)>
            istack(cmp_time);
        istack.push(std::move(iset));

        Array3 err_and_ms = err + ms;

        int refine = 0; // interval split count. aka tree depth

        toi = std::numeric_limits<Scalar>::infinity();
        NumCCD TOI(1, 0);

        bool collision = false;
        int rnbr = 0;
        while (!istack.empty()) {
            Interval3 current = istack.top();
            istack.pop();

            // if(rnbr>0&&less_than( current[0].first,TOI)){
            //     continue;
            // }

            //TOI should always be no larger than current
            if (current[0].lower >= TOI) {
                continue;
            }

            refine++;

            bool zero_in, box_in;
            Array3 tol_placeholder;
            {
                TIGHT_INCLUSION_SCOPED_TIMER(time_predicates);
                zero_in = origin_in_bbox_eval<is_vertex_face, false>(
                    current, a_t0, b_t0, c_t0, d_t0, a_t1, b_t1, c_t1, d_t1,
                    err, ms, box_in, tol_placeholder);
            }

            // #ifdef TIGHT_INCLUSION_WITH_RATIONAL // this is defined in the begining of this file
            // zero_in = origin_in_function_bounding_box_rational<is_vertex_face>(
            //     current, a_t0, b_t0, c_t0, d_t0, a_t1, b_t1, c_t1, d_t1);
            // #endif

            if (!zero_in) {
                continue;
            }

            Array3 widths;
            {
                TIGHT_INCLUSION_SCOPED_TIMER(time_width);
                widths = width(current);
            }

            if (box_in || (widths <= tol).all()) {
                TOI = current[0].lower;
                collision = true;
                rnbr++;
                // continue;
                toi = TOI.value();
                return true;
            }

            // find the next dimension to split
            int split_i = find_next_split(widths, tol);

            bool overflowed = split_and_push(
                current, split_i,
                [&](const Interval3 &i) { istack.emplace(i); }, is_vertex_face);
            if (overflowed) {
                logger().error("overflow occured when splitting intervals!");
                return true;
            }
        }
        if (collision)
            toi = TOI.value();
        return collision;
    }

    NumCCD get_toi(const Interval3 &tuv) { return tuv[0].lower; }

    template <bool is_vertex_face>
    bool interval_root_finder_BFS(
        const Vector3 &a_t0,
        const Vector3 &b_t0,
        const Vector3 &c_t0,
        const Vector3 &d_t0,
        const Vector3 &a_t1,
        const Vector3 &b_t1,
        const Vector3 &c_t1,
        const Vector3 &d_t1,
        const Interval3 &iset,
        const Array3 &tol,
        Scalar co_domain_tolerance,
        const Array3 &err,
        Scalar ms,
        Scalar max_time,
        long max_iter,
        bool is_unit_interval,
        Scalar &toi,
        Scalar &output_tolerance)
    {
        // stack of (interval, tree level) pair.
        // will be sorted by level first then toi lower bound second
        std::vector<std::pair<Interval3, int>> istack;
        istack.emplace_back(iset, 0);

        int iter_count = 0;
        int tree_level = 0;

        // a root interval and its corresponding toi might be skipped if an earlier root interval might exists
        std::optional<Interval3> skipped_candidate;

        // the first root interval at each tree level is the earliest candidate
        bool is_earliest_candidate = true;
        Interval3 earliest_candidate;
        Array3 earliest_bbox_eval_tolerance;

        for (int stack_idx = 0; stack_idx < istack.size(); ++stack_idx) {
            // reach new tree level
            if (istack[stack_idx].second != tree_level) {
                auto cmp = [](const std::pair<Interval3, int> &i1,
                              const std::pair<Interval3, int> &i2) {
                    return i1.first[0].lower < i2.first[0].lower;
                };
                std::sort(
                    istack.data() + stack_idx, istack.data() + istack.size(),
                    cmp);

                ++tree_level;
                is_earliest_candidate = true;
            }

            // current intervals
            Interval3 current = istack[stack_idx].first;

            // if this box is later than skipped toi, skip.
            if (skipped_candidate
                && get_toi(current) >= get_toi(*skipped_candidate)) {
                continue;
            }

            iter_count++;

            // true if bbox eval intersect function root region
            bool bbox_in_zero;
            // true if bbox eval is inside function root region
            bool bbox_in_eps;
            // the xyz width of bbox eval
            Array3 bbox_eval_tolerance;
            if (is_unit_interval) {
                bbox_in_zero = origin_in_bbox_eval<is_vertex_face, true>(
                    current, a_t0, b_t0, c_t0, d_t0, a_t1, b_t1, c_t1, d_t1,
                    err, ms, bbox_in_eps, bbox_eval_tolerance);
                is_unit_interval = false;
            } else {
                bbox_in_zero = origin_in_bbox_eval<is_vertex_face, false>(
                    current, a_t0, b_t0, c_t0, d_t0, a_t1, b_t1, c_t1, d_t1,
                    err, ms, bbox_in_eps, bbox_eval_tolerance);
            }

            // current interval does not contain function root. skip.
            if (!bbox_in_zero) {
                continue;
            }

            // current interval might contain function root.
            // check if interval is small enough to avoid false positive.

            bool is_co_domain_tol_small_enough =
                (bbox_eval_tolerance <= co_domain_tolerance).all();
            bool is_interval_small_enough =
                is_co_domain_tol_small_enough || bbox_in_eps;

            // earlist root interval that is small enough. terminate.
            if (is_earliest_candidate && is_interval_small_enough) {
                // return time interval lower bound as toi
                toi = get_toi(current).value();
                output_tolerance = co_domain_tolerance;
                return true;
            }

            // root interval that is small enough but not the earliest,
            // record then skip.
            if (is_interval_small_enough) {
                if (skipped_candidate) {
                    if (get_toi(*skipped_candidate) > get_toi(current)) {
                        skipped_candidate = current;
                    }
                }
                skipped_candidate = current;
                continue;
            }

            // root interval that is not small enough but is the earliest,
            // record if we have max_iter enabled.
            if (is_earliest_candidate) {
                is_earliest_candidate = false;

                if (max_iter > 0) {
                    earliest_candidate = current;
                    earliest_bbox_eval_tolerance = bbox_eval_tolerance;
                }
            }

            // max iter enabled and reach max iter. return earliest candidate and terminate.
            if (max_iter > 0 && iter_count > max_iter) {
                toi = get_toi(earliest_candidate).value();
                output_tolerance = std::max(
                    co_domain_tolerance,
                    earliest_bbox_eval_tolerance.maxCoeff());
                return true;
            }

            // find the next dimension to split
            Array3 widths = width(current);
            int split_dimension = find_next_split(widths, tol);
            auto push = [&istack, &tree_level](const Interval3 &i) {
                istack.emplace_back(i, tree_level + 1);
            };
            bool overflow = split_and_push(
                current, split_dimension, push, is_vertex_face, max_time);
            if (overflow) {
                logger().error("overflow occured when splitting intervals!");
                return true;
            }
        }

        if (skipped_candidate) {
            toi = get_toi(*skipped_candidate).value();
            output_tolerance = co_domain_tolerance;
            return true;
        }

        // no collision, set toi and output_tolerance to infinity
        toi = std::numeric_limits<Scalar>::infinity();
        output_tolerance = std::numeric_limits<Scalar>::infinity();
        return false;
    }

    template <bool is_vertex_face>
    bool interval_root_finder_BFS(
        const Vector3 &a_t0,
        const Vector3 &b_t0,
        const Vector3 &c_t0,
        const Vector3 &d_t0,
        const Vector3 &a_t1,
        const Vector3 &b_t1,
        const Vector3 &c_t1,
        const Vector3 &d_t1,
        const Array3 &tol,
        const Scalar co_domain_tolerance,
        const Array3 &err,
        const Scalar ms,
        const Scalar max_time,
        const long max_itr,
        Scalar &toi,
        Scalar &output_tolerance)
    {
        // build interval set [0,t_max]x[0,1]x[0,1]
        const Interval zero_to_one = Interval(NumCCD(0, 0), NumCCD(1, 0));
        Interval3 iset = {{
            // Interval(NumCCD(0, 0), NumCCD(max_time)),
            zero_to_one,
            zero_to_one,
            zero_to_one,
        }};

        return interval_root_finder_BFS<is_vertex_face>(
            a_t0, b_t0, c_t0, d_t0, a_t1, b_t1, c_t1, d_t1, iset, tol,
            co_domain_tolerance, err, ms, max_time, max_itr, true, toi,
            output_tolerance);
    }

    Array3 get_numerical_error(
        const std::vector<Vector3> &vertices,
        const bool is_vertex_face,
        const bool using_minimum_separation)
    {
        Scalar eefilter;
        Scalar vffilter;
        if (!using_minimum_separation) {
#ifdef TIGHT_INCLUSION_WITH_DOUBLE_PRECISION
            eefilter = 6.217248937900877e-15;
            vffilter = 6.661338147750939e-15;
#else
            eefilter = 3.337861e-06;
            vffilter = 3.576279e-06;
#endif
        } else { // using minimum separation
#ifdef TIGHT_INCLUSION_WITH_DOUBLE_PRECISION
            eefilter = 7.105427357601002e-15;
            vffilter = 7.549516567451064e-15;
#else
            eefilter = 3.814698e-06;
            vffilter = 4.053116e-06;
#endif
        }

        Vector3 max = vertices[0].cwiseAbs();
        for (int i = 1; i < vertices.size(); i++) {
            max = max.cwiseMax(vertices[i].cwiseAbs());
        }
        Vector3 delta = max.cwiseMin(1);
        Scalar filter = is_vertex_face ? vffilter : eefilter;
        return filter * delta.array().pow(3);
    }

    bool edge_edge_interval_root_finder_DFS(
        const Vector3 &ea0_t0,
        const Vector3 &ea1_t0,
        const Vector3 &eb0_t0,
        const Vector3 &eb1_t0,
        const Vector3 &ea0_t1,
        const Vector3 &ea1_t1,
        const Vector3 &eb0_t1,
        const Vector3 &eb1_t1,
        const Array3 &tol,
        const Array3 &err,
        const Scalar ms,
        Scalar &toi)
    {
        return interval_root_finder_DFS<false>(
            ea0_t0, ea1_t0, eb0_t0, eb1_t0, ea0_t1, ea1_t1, eb0_t1, eb1_t1, tol,
            err, ms, toi);
    }

    bool vertex_face_interval_root_finder_DFS(
        const Vector3 &v_t0,
        const Vector3 &f0_t0,
        const Vector3 &f1_t0,
        const Vector3 &f2_t0,
        const Vector3 &v_t1,
        const Vector3 &f0_t1,
        const Vector3 &f1_t1,
        const Vector3 &f2_t1,
        const Array3 &tol,
        const Array3 &err,
        const Scalar ms,
        Scalar &toi)
    {
        return interval_root_finder_DFS<true>(
            v_t0, f0_t0, f1_t0, f2_t0, v_t1, f0_t1, f1_t1, f2_t1, tol, err, ms,
            toi);
    }

    bool edge_edge_interval_root_finder_BFS(
        const Vector3 &ea0_t0,
        const Vector3 &ea1_t0,
        const Vector3 &eb0_t0,
        const Vector3 &eb1_t0,
        const Vector3 &ea0_t1,
        const Vector3 &ea1_t1,
        const Vector3 &eb0_t1,
        const Vector3 &eb1_t1,
        const Array3 &tol,
        const Scalar co_domain_tolerance,
        // this is the maximum error on each axis when calculating the vertices, err, aka, filter
        const Array3 &err,
        const Scalar ms,
        const Scalar max_time,
        const long max_itr,
        Scalar &toi,
        Scalar &output_tolerance)
    {
        return interval_root_finder_BFS<false>(
            ea0_t0, ea1_t0, eb0_t0, eb1_t0, ea0_t1, ea1_t1, eb0_t1, eb1_t1, tol,
            co_domain_tolerance, err, ms, max_time, max_itr, toi,
            output_tolerance);
    }

    bool vertex_face_interval_root_finder_BFS(
        const Vector3 &v_t0,
        const Vector3 &f0_t0,
        const Vector3 &f1_t0,
        const Vector3 &f2_t0,
        const Vector3 &v_t1,
        const Vector3 &f0_t1,
        const Vector3 &f1_t1,
        const Vector3 &f2_t1,
        const Array3 &tol,
        const Scalar co_domain_tolerance,
        // this is the maximum error on each axis when calculating the vertices, err, aka, filter
        const Array3 &err,
        const Scalar ms,
        const Scalar max_time,
        const long max_itr,
        Scalar &toi,
        Scalar &output_tolerance)
    {
        return interval_root_finder_BFS<true>(
            v_t0, f0_t0, f1_t0, f2_t0, v_t1, f0_t1, f1_t1, f2_t1, tol,
            co_domain_tolerance, err, ms, max_time, max_itr, toi,
            output_tolerance);
    }

    // ------------------------------------------------------------------------
    // Template instantiation
    // ------------------------------------------------------------------------

    // clang-format off
    template bool interval_root_finder_DFS<false>(const Vector3 &,const Vector3 &,const Vector3 &,const Vector3 &,const Vector3 &,const Vector3 &,const Vector3 &,const Vector3 &,const Array3 &,const Array3 &,const Scalar,Scalar &);
    template bool interval_root_finder_DFS<true>(const Vector3 &,const Vector3 &,const Vector3 &,const Vector3 &,const Vector3 &,const Vector3 &,const Vector3 &,const Vector3 &,const Array3 &,const Array3 &,const Scalar,Scalar &);
    template bool interval_root_finder_BFS<false>(const Vector3 &,const Vector3 &,const Vector3 &,const Vector3 &,const Vector3 &,const Vector3 &,const Vector3 &,const Vector3 &,const Array3 &,const Scalar,const Array3 &,const Scalar,const Scalar,const long,Scalar &,Scalar &);
    template bool interval_root_finder_BFS<true>(const Vector3 &,const Vector3 &,const Vector3 &,const Vector3 &,const Vector3 &,const Vector3 &,const Vector3 &,const Vector3 &,const Array3 &,const Scalar,const Array3 &,const Scalar,const Scalar,const long,Scalar &,Scalar &);
    // clang-format on

} // namespace ticcd
