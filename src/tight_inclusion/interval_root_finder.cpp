// A root finder using interval arithmetic.
#include "interval_root_finder.hpp"
#include "types.hpp"

#include <tight_inclusion/timer.hpp>
#include <tight_inclusion/avx.hpp>
#include <tight_inclusion/logger.hpp>

#include <queue>

namespace ticcd {
    double time_predicates = 0, time_width = 0, time_bisect = 0,
           time_eval_origin_1D = 0, time_eval_origin_tuv = 0,
           time_vertex_solving = 0;

    // ** this version can return the true x or y or z tolerance of the co-domain **
    // eps is the interval [-eps,eps] we need to check
    // if [-eps,eps] overlap, return true
    // bbox_in_eps tell us if the box is totally in eps box
    // ms is the minimum seperation
    template <bool is_vertex_face>
    bool evaluate_bbox_one_dimension_vector(
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
        const int dim,
        const Scalar eps,
        bool &bbox_in_eps,
        const Scalar ms = 0,
        Scalar *tol = nullptr)
    {
        TIGHT_INCLUSION_SCOPED_TIMER(time_vertex_solving);

        Array8 vs;
        if constexpr (is_vertex_face) {
            vs = function_vf(
                a_t0[dim], b_t0[dim], c_t0[dim], d_t0[dim], a_t1[dim],
                b_t1[dim], c_t1[dim], d_t1[dim], t_up, t_dw, u_up, u_dw, v_up,
                v_dw);
        } else {
            vs = function_ee(
                a_t0[dim], b_t0[dim], c_t0[dim], d_t0[dim], a_t1[dim],
                b_t1[dim], c_t1[dim], d_t1[dim], t_up, t_dw, u_up, u_dw, v_up,
                v_dw);
        }

        Scalar minv = vs.minCoeff();
        Scalar maxv = vs.maxCoeff();

        if (tol != nullptr) {
            *tol = maxv - minv; // this is the real tolerance
        }

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
    template <bool is_vertex_face>
    bool origin_in_function_bounding_box_vector(
        const Interval3 &paras,
        const Vector3 &a_t0,
        const Vector3 &b_t0,
        const Vector3 &c_t0,
        const Vector3 &d_t0,
        const Vector3 &a_t1,
        const Vector3 &b_t1,
        const Vector3 &c_t1,
        const Vector3 &d_t1,
        const Array3 &eps,
        bool &box_in_eps,
        const Scalar ms = 0,
        Array3 *tolerance = nullptr)
    {
        box_in_eps = false;

        Array8 t_up, t_dw, u_up, u_dw, v_up, v_dw;
        {
            TIGHT_INCLUSION_SCOPED_TIMER(time_eval_origin_tuv);
            convert_tuv_to_array(paras, t_up, t_dw, u_up, u_dw, v_up, v_dw);
        }

        bool box_in[3];
        for (int i = 0; i < 3; i++) {
            TIGHT_INCLUSION_SCOPED_TIMER(time_eval_origin_1D);
            Scalar *tol = tolerance == nullptr ? nullptr : &((*tolerance)[i]);
            if (!evaluate_bbox_one_dimension_vector<is_vertex_face>(
                    t_up, t_dw, u_up, u_dw, v_up, v_dw, a_t0, b_t0, c_t0, d_t0,
                    a_t1, b_t1, c_t1, d_t1, i, eps[i], box_in[i], ms, tol)) {
                return false;
            }
        }

        if (box_in[0] && box_in[1] && box_in[2]) {
            box_in_eps = true;
        }

        return true;
    }

    // find the largest width/tol dimension that is greater than its tolerance
    int find_next_split(const Array3 &widths, const Array3 &tols)
    {
        // assert((widths > tols).any());
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

        //build interval set [0,1]x[0,1]x[0,1]
        const Interval zero_to_one = Interval(NumCCD(0, 0), NumCCD(1, 0));
        Interval3 iset = {{zero_to_one, zero_to_one, zero_to_one}};

        // Stack of intervals and the last split dimension
        // std::stack<std::pair<Interval3,int>> istack;
        std::priority_queue<
            Interval3, std::vector<Interval3>, decltype(cmp_time)>
            istack(cmp_time);
        istack.emplace(iset);

        Array3 err_and_ms = err + ms;

        int refine = 0;

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
            {
                TIGHT_INCLUSION_SCOPED_TIMER(time_predicates);
                zero_in =
                    origin_in_function_bounding_box_vector<is_vertex_face>(
                        current, a_t0, b_t0, c_t0, d_t0, a_t1, b_t1, c_t1, d_t1,
                        err_and_ms, box_in);
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

    // when check_t_overlap = false, check [0,1]x[0,1]x[0,1]; otherwise, check [0, t_max]x[0,1]x[0,1]
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
        const Scalar co_domain_tolerance,
        const Array3 &err,
        const Scalar ms,
        const Scalar max_time,
        const long max_itr,
        Scalar &toi,
        Scalar &output_tolerance)
    {
        long queue_size = 0;
        // if max_itr <0, output_tolerance= co_domain_tolerance;
        // else, output_tolearancewill be the precision after iteration time > max_itr
        output_tolerance = co_domain_tolerance;

        // this is used to catch the tolerance for each level
        Scalar temp_output_tolerance = co_domain_tolerance;

        // check the tree level by level instead of going deep
        // (if level 1 != level 2, return level 1 >= level 2; else, return time1 >= time2)
        auto cmp = [](const std::pair<Interval3, int> &i1,
                      const std::pair<Interval3, int> &i2) {
            if (i1.second != i2.second) {
                return i1.second >= i2.second;
            } else {
                return i1.first[0].lower > i2.first[0].lower;
            }
        };

        // Stack of intervals and the last split dimension.
        // Sorted by interval level first then by toi lower bound.

        std::priority_queue<
            std::pair<Interval3, int>, std::vector<std::pair<Interval3, int>>,
            decltype(cmp)>
            istack(cmp);
        istack.emplace(iset, -1);

        // current intervals
        Interval3 current;
        int refine = 0;

        // set TOI to 4. this is to record the impact time of this level
        NumCCD TOI(4, 0);
        // this is to record the element that already small enough or contained in eps-box
        NumCCD TOI_SKIP = TOI;
        bool use_skip = false;  // this is to record if TOI_SKIP is used.
        int current_level = -2; // in the begining, current_level != level
        // level < tolerance. only true, we can return when we find one overlaps eps box and smaller than tolerance or eps-box
        bool this_level_less_tol = true;
        bool find_level_root = false;
        while (!istack.empty()) {
#ifdef TIGHT_INCLUSION_CHECK_QUEUE_SIZE
            if (istack.size() > queue_size) {
                queue_size = istack.size();
            }
#endif
#ifdef TIGHT_INCLUSION_LIMIT_QUEUE_SIZE
            if (istack.size() > MAX_QSIZE) {
                return true;
            }
#endif

            current = istack.top().first;
            int level = istack.top().second;
            istack.pop();

            // if this box is later than TOI_SKIP in time, we can skip this one.
            // TOI_SKIP is only updated when the box is small enough or totally contained in eps-box
            if (current[0].lower >= TOI_SKIP) {
                continue;
            }
            // before check a new level, set this_level_less_tol=true
            if (current_level != level) {
                current_level = level;
                this_level_less_tol = true;
                find_level_root = false;
            }

            refine++;
            bool zero_in, box_in;
            Array3 true_tol;
            {
                TIGHT_INCLUSION_SCOPED_TIMER(time_predicates);
                // #ifdef TIGHT_INCLUSION_WITH_RATIONAL // this is defined in the begining of this file
                // Array3 ms_3d = Array3::Constant(ms);
                // zero_in = origin_in_function_bounding_box_rational_return_tolerance<is_vertex_face>(
                //     current, a_t0, b_t0, c_t0, d_t0, a_t1, b_t1, c_t1, d_t1,
                //     ms_3d, box_in, true_tol);
                // #else
                zero_in =
                    origin_in_function_bounding_box_vector<is_vertex_face>(
                        current, a_t0, b_t0, c_t0, d_t0, a_t1, b_t1, c_t1, d_t1,
                        err, box_in, ms, &true_tol);
                // #endif
            }

            if (!zero_in)
                continue;

            Array3 widths;
            {
                TIGHT_INCLUSION_SCOPED_TIMER(time_width);
                widths = width(current);
            }

            bool tol_condition = (true_tol <= co_domain_tolerance).all();

            // terminate early if following condition is true:
            // a. no previous interval with smaller toi contains the root box
            // b. current interval is small enough to pass tolerance test
            bool is_interval_small_enough = tol_condition || box_in;
            if (this_level_less_tol && is_interval_small_enough) {
                TOI = current[0].lower;
                toi = TOI.value();
                return true;
            }

            if (!tol_condition) {
                this_level_less_tol = false;
            }

            if (max_itr > 0) { // if max_itr <= 0 ⟹ unlimited iterations
                // current_tolerance=std::max(
                // std::max(std::max(current_tolerance,true_tol[0]),true_tol[1]),true_tol[2]
                // );
                if (!find_level_root) {
                    TOI = current[0].lower;
                    // collision=true;
                    // continue;

                    // if the real tolerance is larger than input, use the real one;
                    // if the real tolerance is smaller than input, use input
                    temp_output_tolerance = std::max(
                        {true_tol[0], true_tol[1], true_tol[2],
                         co_domain_tolerance});
                    // this ensures always find the earlist root
                    find_level_root = true;
                }
                if (refine > max_itr) {
                    toi = TOI.value();
                    output_tolerance = temp_output_tolerance;

                    return true;
                }
                // get the time of impact down here
            }

            // if this box is small enough, or inside of eps-box, then just continue,
            // but we need to record the collision time
            if (is_interval_small_enough) {
                if (current[0].lower < TOI_SKIP) {
                    TOI_SKIP = current[0].lower;
                }
                use_skip = true;
                continue;
            }

            // find the next dimension to split
            int split_i = find_next_split(widths, tol);

            bool overflow = split_and_push(
                current, split_i,
                [&](const Interval3 &i) { istack.emplace(i, level + 1); },
                is_vertex_face, max_time);
            if (overflow) {
                logger().error("overflow occured when splitting intervals!");
                return true;
            }
        }

        if (use_skip) {
            toi = TOI_SKIP.value();
            return true;
        }

        toi = std::numeric_limits<Scalar>::infinity(); //set toi as infinity
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
            co_domain_tolerance, err, ms, max_time, max_itr, toi,
            output_tolerance);
    }

    void print_times()
    {
        logger().trace("[time] origin predicates, {}", time_predicates);
        logger().trace("[time] width, {}", time_width);
        logger().trace("[time] bisect, {}", time_bisect);
        logger().trace(
            "[time] origin part1(evaluate 1 dimension), {}",
            time_eval_origin_1D);
        logger().trace(
            "[time] origin part2(convert tuv), {}", time_eval_origin_tuv);
        logger().trace(
            "[time] time of call the vertex solving function, {}",
            time_vertex_solving);
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
