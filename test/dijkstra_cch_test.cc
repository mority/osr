#ifdef _WIN32
#include "windows.h"
#endif

#include "gtest/gtest.h"

#include <filesystem>
#include <random>

#include "cista/mmap.h"

#include "utl/parallel_for.h"

#include "fmt/core.h"

#include "osr/extract/extract.h"
#include "osr/geojson.h"
#include "osr/location.h"
#include "osr/lookup.h"
#include "osr/routing/cch/cch.h"
#include "osr/routing/dijkstra.h"
#include "osr/routing/profile.h"
#include "osr/routing/profiles/car.h"
#include "osr/routing/route.h"
#include "osr/types.h"
#include "osr/ways.h"

namespace fs = std::filesystem;
using namespace osr;

namespace {

constexpr auto const kUseMultithreading = true;
constexpr auto const kPrintDebugGeojson = false;
constexpr auto const kMaxMatchDistance = 100;
constexpr auto const kMaxAllowedPathDifferenceRatio = 0.5;

void load(std::string_view raw_data, std::string_view data_dir) {
  if (fs::exists(raw_data)) {
    auto const p = fs::path{data_dir};
    auto ec = std::error_code{};
    fs::remove_all(p, ec);
    fs::create_directories(p, ec);
    osr::extract(false, raw_data, data_dir, fs::path{});
  }
}

void run(ways const& w,
         lookup const& l,
         unsigned const n_samples,
         unsigned const max_cost,
         direction const dir) {

  auto const from_tos = [&]() {
    auto prng = std::mt19937{};
    auto distr =
        std::uniform_int_distribution<std::uint32_t>{0, w.n_nodes() - 1};
    auto from_tos = std::vector<std::pair<node_idx_t, node_idx_t>>{};
    for (auto i = 0U; i != n_samples; ++i) {
      from_tos.emplace_back(distr(prng), distr(prng));
    }
    return from_tos;
  }();

  auto n_congruent = std::atomic<unsigned>{0U};
  auto n_empty_matches = std::atomic<unsigned>{0U};
  auto reference_times = std::vector<std::chrono::steady_clock::duration>{};
  auto experiment_times = std::vector<std::chrono::steady_clock::duration>{};

  auto m = std::mutex{};

  auto const single_run = [&](std::pair<node_idx_t, node_idx_t> const from_to) {
    auto const from_node = from_to.first;
    auto const from_loc = location{w.get_node_pos(from_node)};
    auto const to_node = from_to.second;
    auto const to_loc = location{w.get_node_pos(to_node)};

    // Keeps only the candidates that touch `n`. Candidates cannot be erased
    // in place, so the survivors are appended to a fresh match.
    auto const node_pinned_matches =
        [&](location const& loc, node_idx_t const n, bool const reverse) {
          auto all = match_result{};
          l.match<car>(car::parameters{}, loc, reverse, dir, kMaxMatchDistance,
                       nullptr, all);
          auto const m = all[match_idx_t{0U}];
          auto pinned = match_result{};
          pinned.start(m.lvl_);
          for (auto j = std::size_t{0U}; j != m.size(); ++j) {
            if (m.nodes_[j].left_.node_ == n || m.nodes_[j].right_.node_ == n) {
              pinned.add(m.dist_to_way_[j], m.way_[j], m.nodes_[j]);
            }
          }
          pinned.finish();
          return pinned;
        };
    auto const from_matches = node_pinned_matches(from_loc, from_node, false);
    auto const to_matches = node_pinned_matches(to_loc, to_node, true);
    auto const from_matches_span = from_matches[match_idx_t{0U}];
    auto const to_matches_span = to_matches[match_idx_t{0U}];
    if (from_matches_span.empty() || to_matches_span.empty()) {
      ++n_empty_matches;
    }

    auto const reference_start = std::chrono::steady_clock::now();
    auto const reference = [&]() {
      try {
        return route(car::parameters{}, w, l, search_profile::kCar, from_loc,
                     to_loc, from_matches_span, to_matches_span, max_cost, dir,
                     nullptr, nullptr, nullptr, routing_algorithm::kDijkstra);
      } catch (std::exception const& ex) {
        fmt::println("dijkstra exception: {}", ex.what());
        throw ex;
      }
    }();
    auto const reference_time =
        std::chrono::steady_clock::now() - reference_start;

    auto const experiment_start = std::chrono::steady_clock::now();
    auto const experiment = [&]() {
      try {
        return route(car::parameters{}, w, l, search_profile::kCar, from_loc,
                     to_loc, from_matches_span, to_matches_span, max_cost, dir,
                     nullptr, nullptr, nullptr, routing_algorithm::kCCH);
      } catch (std::exception const& ex) {
        fmt::println("cch exception: {}", ex.what());
        throw ex;
      }
    }();
    auto const experiment_time =
        std::chrono::steady_clock::now() - experiment_start;

    if (reference.has_value() != experiment.has_value() ||
        (reference && experiment &&
         (reference->cost_ != experiment->cost_ /*||
          std::abs(reference->dist_ - experiment->dist_) / reference->dist_ >
              kMaxAllowedPathDifferenceRatio*/))) {
      auto const print_result = [&](std::string_view name, auto const& p,
                                    auto const& t) {
        fmt::println(
            "{:10}: {:11} --> {:11} | {} | time: "
            "{}:{:0>3}:{:0>3} s",
            name, w.node_to_osm_[from_node], w.node_to_osm_[to_node],
            p ? fmt::format("cost: {:5} | dist: {:>10.2f}", p->cost_, p->dist_)
              : "no result",
            std::chrono::duration_cast<std::chrono::seconds>(t).count(),
            std::chrono::duration_cast<std::chrono::milliseconds>(t).count() %
                1000,
            std::chrono::duration_cast<std::chrono::microseconds>(t).count() %
                1000);
        if (p.has_value() && kPrintDebugGeojson) {
          fmt::println("{}\n", to_featurecollection(w, p));
        }
      };

      print_result("dijkstra", reference, reference_time);
      print_result("cch", experiment, experiment_time);

    } else {
      ++n_congruent;
    }

    if (!from_matches.empty() && !to_matches.empty()) {
      auto const guard = std::lock_guard{m};
      reference_times.emplace_back(reference_time);
      experiment_times.emplace_back(experiment_time);
    }
  };

  if (kUseMultithreading) {
    utl::parallel_for(from_tos, single_run);
  } else {
    std::for_each(begin(from_tos), end(from_tos), single_run);
  }

  auto const non_empty_congruent = n_congruent - n_empty_matches;
  auto const non_empty_samples = n_samples - n_empty_matches;

  EXPECT_EQ(non_empty_samples, non_empty_congruent);

  fmt::println("congruent on non-empty: {}/{} ({:3.1f}%)", non_empty_congruent,
               non_empty_samples,
               (static_cast<double>(non_empty_congruent) /
                static_cast<double>(non_empty_samples)) *
                   100);
  if (non_empty_congruent == non_empty_samples) {
    fmt::println(
        "speedup on non-empty: {:.2f}",
        static_cast<double>(
            std::reduce(begin(reference_times), end(reference_times)).count()) /
            static_cast<double>(
                std::reduce(begin(experiment_times), end(experiment_times))
                    .count()));
  }
}

// Every case runs the same comparison against the same oracle; only the graph,
// the sample count, the cost bound and the direction differ.
void run_case(char const* raw_data,
              char const* data_dir,
              unsigned const num_samples,
              unsigned const max_cost,
              direction const dir) {
  if (!fs::exists(raw_data) && !fs::exists(data_dir)) {
    GTEST_SKIP() << raw_data << " not found";
  }

  load(raw_data, data_dir);
  if (!cch::exists(data_dir)) {
    GTEST_SKIP() << "no cch in " << data_dir;
  }
  auto const w = osr::ways{data_dir, cista::mmap::protection::READ};
  auto const l = osr::lookup{w, data_dir, cista::mmap::protection::READ};

  run(w, l, num_samples, max_cost, dir);
}

}  // namespace

// Forward only, deliberately. `route_cch` applies `dir` to the first and last
// mile but the hierarchy search itself takes no direction and always walks the
// forward graph, so a backward point to point query is not implemented --
// `cch_supported` routes it to the Dijkstra instead. A backward case here would
// therefore compare the Dijkstra against itself and pass without asserting
// anything about the CCH, which is what the disabled case below did.
TEST(dijkstra_cch, monaco_fwd) {
  run_case("test/monaco.osm.pbf", "test/monaco", 10000U, 2 * 3600U,
           direction::kForward);
}

TEST(dijkstra_cch, hamburg) {
  run_case("test/hamburg.osm.pbf", "test/hamburg", 5000U, 3 * 3600U,
           direction::kForward);
}

TEST(dijkstra_cch, switzerland) {
  run_case("test/switzerland.osm.pbf", "test/switzerland", 1000U, 5 * 3600U,
           direction::kForward);
}

TEST(dijkstra_cch, DISABLED_germany) {
  run_case("test/germany.osm.pbf", "test/germany", 50U, 12 * 3600U,
           direction::kForward);
}

// Kept disabled, and not worth enabling. `cch_supported` sends a backward
// point to point query to the Dijkstra, so both arms here are the same
// algorithm and the case passes without asserting anything. Wiring `dir` into
// `route_cch` would not change that verdict either: with no time dependence a
// backward one to one query is the forward query with its endpoints swapped,
// so it makes nothing computable that `monaco_fwd` does not already cover.
// (One to many is a different matter -- there, backward is many to one, which
// is a genuinely distinct query and is enabled.)
TEST(dijkstra_cch, DISABLED_monaco_bwd) {
  run_case("test/monaco.osm.pbf", "test/monaco", 10000U, 2 * 3600U,
           direction::kBackward);
}
