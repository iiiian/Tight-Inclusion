#include "config.hpp"

#include <tight_inclusion/ccd.hpp>
#include <tight_inclusion/timer.hpp>
#include <tight_inclusion/logger.hpp>
#include <tight_inclusion/config.hpp>

#ifdef TIGHT_INCLUSION_WITH_SAMPLE_QUERIES
#include <ccd_io/read_ccd_queries.hpp>
#endif

#include <pbar.hpp>
#include <CLI/CLI.hpp>

#include <vector>
#include <filesystem>
namespace fs = std::filesystem;

using namespace ticcd;

#ifdef TIGHT_INCLUSION_WITH_SAMPLE_QUERIES

static const std::string root_path(CCD_IO_SAMPLE_QUERIES_DIR);

static const std::vector<std::string> simulation_folders = {{
    "chain",
    "cow-heads",
    "golf-ball",
    "mat-twist",
}};

static const std::vector<std::string> handcrafted_folders = {{
    "erleben-sliding-spike",
    "erleben-spike-wedge",
    "erleben-sliding-wedge",
    "erleben-wedge-crack",
    "erleben-spike-crack",
    "erleben-wedges",
    "erleben-cube-cliff-edges",
    "erleben-spike-hole",
    "erleben-cube-internal-edges",
    "erleben-spikes",
    "unit-tests",
}};

void check_sample_queries(
    const bool is_edge_edge,
    const bool is_simulation_data,
    const double minimum_seperation,
    const double tolerance,
    const long max_itr,
    const bool print_progress = true)
{
    using Matrix8x3 = Eigen::Matrix<double, 8, 3, Eigen::RowMajor>;

    const Eigen::Array3d err(-1, -1, -1);
    constexpr double t_max = 1;
    constexpr bool no_zero_toi = false;

    int total_positives = 0;
    int total_false_positives = 0;
    int total_false_negatives = 0;

    Timer timer;
    double total_time_in_micro_sec = 0.0;

    const auto folders =
        is_simulation_data ? simulation_folders : handcrafted_folders;
    const std::string sub_folder = is_edge_edge ? "edge-edge" : "vertex-face";

    size_t total_number_of_queries = 0;
    pbar::pbar bar(-1);
    if (print_progress) {
        logger().trace("Determining total number of queries...");
        for (const std::string &folder : folders) {
            fs::path dir = fs::path(root_path) / folder / sub_folder;
            for (const auto &csv : fs::directory_iterator(dir)) {
                std::ifstream in_stream(csv.path().string());
                unsigned int line_count = std::count_if(
                    std::istreambuf_iterator<char>{in_stream}, {},
                    [](char c) { return c == '\n'; });
                assert(line_count % 8 == 0);
                total_number_of_queries += line_count / 8;
            }
        }
        logger().trace("Total number of queries: {}", total_number_of_queries);

        bar = pbar::pbar(total_number_of_queries);
        bar.enable_recalc_console_width(1); // check console width every tick
        bar.init();
    }

    for (const std::string &folder : folders) {
        fs::path dir = fs::path(root_path) / folder / sub_folder;
        for (const auto &csv : fs::directory_iterator(dir)) {

            const std::vector<ccd_io::CCDQuery> queries =
                ccd_io::read_ccd_queries(csv.path().string());

            for (int i = 0; i < queries.size(); i++) {
                Eigen::Map<const Matrix8x3> V(&queries[i].vertices[0][0]);
                const bool expected_result = queries[i].ground_truth;
                total_positives += expected_result;

                // Output of CCD
                std::optional<Collision> collision;
                double toi;
                double u;
                double v;
                double output_tolerance = tolerance;

                timer.start();
                if (is_edge_edge) {
                    collision = edgeEdgeCCD(
                        V.row(0), V.row(1), V.row(2), V.row(3), V.row(4),
                        V.row(5), V.row(6), V.row(7), err, minimum_seperation,
                        tolerance, t_max, max_itr, no_zero_toi);
                } else {
                    collision = vertexFaceCCD(
                        V.row(0), V.row(1), V.row(2), V.row(3), V.row(4),
                        V.row(5), V.row(6), V.row(7), err, minimum_seperation,
                        tolerance, t_max, max_itr, no_zero_toi);
                    // result = rational::vertexFaceCCD(
                    //     V.row(0), V.row(1), V.row(2), V.row(3), V.row(4),
                    //     V.row(5), V.row(6), V.row(7), err, minimum_seperation,
                    //     toi);
                }
                timer.stop();
                total_time_in_micro_sec += timer.getElapsedTimeInMicroSec();

                if (collision.has_value() != expected_result) {
                    if (collision) {
                        total_false_positives++;
                    } else {
                        total_false_negatives++;

                        logger().error(
                            "False negative encountered in file \"{}\" (query #{})!",
                            csv.path().string(), i);
                        for (int j = 0; j < 8; j++) {
                            logger().debug(
                                "V{}: {:.17f} {:.17f} {:.17f}", //
                                j, V(j, 0), V(j, 1), V(j, 2));
                        }
                        logger().debug("Is edge: {}", is_edge_edge);
                        assert(false);
                        // exit(1);
                    }
                }

                if (print_progress) {
                    ++bar; // total_number_of_queries already computed
                } else {
                    ++total_number_of_queries;
                }
            }
        }
    }

    logger().info("total # of queries:   {:7d}", total_number_of_queries);
    logger().info("total positives:      {:7d}", total_positives);
    logger().info("# of false positives: {:7d}", total_false_positives);
    logger().info("# of false negatives: {:7d}", total_false_negatives);
    logger().info(
        "total time:           {:g} s", total_time_in_micro_sec / 1e6);
    logger().info(
        "average time:         {:g} μs",
        total_time_in_micro_sec / double(total_number_of_queries));
}

void check_all_sample_queries(const bool print_progress = true)
{
    constexpr double ms = 0.0;
    constexpr double tolerance = 1e-6;
    constexpr long max_itr = 1e6;

    for (const bool is_simulation : std::array<bool, 2>{{true, false}}) {
        for (const bool is_edge_edge : std::array<bool, 2>{{false, true}}) {
            fmt::print(
                "\nRunning {} {} data:\n",
                is_simulation ? "simulation" : "handcrafted",
                is_edge_edge ? "edge-edge" : "vertex-face");
            check_sample_queries(
                is_edge_edge, is_simulation, ms, tolerance, max_itr,
                print_progress);
        }
    }
}
#endif

void check_single_case()
{
    const Vector3 ea0_t0(0.1, 0.1, 0.1);
    const Vector3 ea1_t0(0, 0, 1);
    const Vector3 ea0_t1(1, 0, 1);
    const Vector3 ea1_t1(0, 1, 1);
    const Vector3 eb0_t0(0.1, 0.1, 0.1);
    const Vector3 eb1_t0(0, 0, 0);
    const Vector3 eb0_t1(0, 1, 0);
    const Vector3 eb1_t1(1, 0, 0);

    const Array3 err(-1, -1, -1);
    constexpr double ms = 1e-8;
    constexpr double tolerance = 1e-6;
    constexpr double t_max = 1;
    constexpr long max_itr = 1e6;

    double toi, u, v, output_tolerance;
    auto collision = edgeEdgeCCD(
        ea0_t0, ea1_t0, ea0_t1, ea1_t1, eb0_t0, eb1_t0, eb0_t1, eb1_t1, err, ms,
        tolerance, t_max, max_itr);

    logger().info("Double CCD result: {}", collision.has_value());
#ifdef TIGHT_INCLUSION_CHECK_QUEUE_SIZE
    logger().info("queue size max {}", return_queue_size());
#endif
}

int main(int argc, char *argv[])
{
    logger().set_level(spdlog::level::trace);

    CLI::App app{"Tight Inclusion CCD"};

    bool print_progress = true;
    app.add_flag("--pbar,!--no-pbar", print_progress, "Hide progress bar");

    CLI11_PARSE(app, argc, argv);
    logger().debug("Progress bar: {}", print_progress);

    logger().debug(
        "Using {} precision floating point numbers",
#ifdef TIGHT_INCLUSION_WITH_DOUBLE_PRECISION
        "double"
#else
        "single"
#endif
    );

#ifdef TIGHT_INCLUSION_WITH_SAMPLE_QUERIES
    check_all_sample_queries(print_progress);
#else
    check_single_case();
#endif

    return 0;
}
