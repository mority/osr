#ifdef _WIN32
#include "windows.h"
#endif

#include <atomic>
#include <chrono>
#include <filesystem>
#include <iostream>
#include <mutex>
#include <random>

#include "cista/mmap.h"

#include "fmt/core.h"
#include "fmt/std.h"

#include "conf/options_parser.h"

#include "utl/parallel_for.h"

#include "osr/location.h"
#include "osr/lookup.h"
#include "osr/routing/cch/cch.h"
#include "osr/routing/parameters.h"
#include "osr/routing/profile.h"
#include "osr/routing/profiles/car.h"
#include "osr/routing/route.h"
#include "osr/types.h"
#include "osr/ways.h"

namespace fs = std::filesystem;
using namespace osr;

struct settings : public conf::configuration {
  settings() : configuration{"Options"} {
    param(data_, "data,d", "osr data directory");
    param(n_queries_, "queries,n", "number of queries");
    param(threads_, "threads,t", "number of routing threads");
    param(max_cost_, "max_cost,c", "maximum cost of a route [s]");
    param(max_match_distance_, "match,m", "maximum matching distance [m]");
  }

  fs::path data_{"osr"};
  unsigned n_queries_{1000U};
  unsigned threads_{8U};
  unsigned max_cost_{12U * 3600U};
  unsigned max_match_distance_{100U};
};

// Both algorithms have to see exactly the same candidates, otherwise a
// difference in the result does not say anything about the search itself.
match_result pinned_match(ways const& w,
                          lookup const& l,
                          location const& loc,
                          node_idx_t const n,
                          bool const reverse,
                          direction const dir,
                          double const max_match_distance) {
  auto all = match_result{};
  l.match<car>(car::parameters{}, loc, reverse, dir, max_match_distance,
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
}

int main(int ac, char const** av) {
  auto opt = settings{};
  auto parser = conf::options_parser{{&opt}};
  parser.read_command_line_args(ac, av);
  if (parser.help()) {
    parser.print_help(std::cout);
    return 0;
  }
  parser.print_used(std::cout);

  if (!fs::is_directory(opt.data_)) {
    fmt::println("data directory {} not found", opt.data_);
    return 1;
  }
  if (!cch::exists(opt.data_)) {
    fmt::println("no cch.bin in {} - run osr-cch first", opt.data_);
    return 1;
  }

  auto const w = ways{opt.data_, cista::mmap::protection::READ};
  auto const l = lookup{w, opt.data_, cista::mmap::protection::READ};
  auto const params = profile_parameters{car::parameters{}};

  auto prng = std::mt19937{42U};
  auto distr =
      std::uniform_int_distribution<std::uint32_t>{0U, w.n_nodes() - 1U};
  auto from_tos = std::vector<std::pair<node_idx_t, node_idx_t>>{};
  for (auto i = 0U; i != opt.n_queries_; ++i) {
    from_tos.emplace_back(node_idx_t{distr(prng)}, node_idx_t{distr(prng)});
  }

  auto n_congruent = std::atomic<unsigned>{0U};
  auto n_compared = std::atomic<unsigned>{0U};
  auto reference_time = std::atomic<std::int64_t>{0};
  auto experiment_time = std::atomic<std::int64_t>{0};
  auto mutex = std::mutex{};

  auto const single_run = [&](std::pair<node_idx_t, node_idx_t> const& ft) {
    auto const from_loc = location{w.get_node_pos(ft.first)};
    auto const to_loc = location{w.get_node_pos(ft.second)};
    auto const from_m = pinned_match(w, l, from_loc, ft.first, false,
                                     direction::kForward,
                                     opt.max_match_distance_);
    auto const to_m = pinned_match(w, l, to_loc, ft.second, true,
                                   direction::kForward,
                                   opt.max_match_distance_);
    auto const from_match = from_m[match_idx_t{0U}];
    auto const to_match = to_m[match_idx_t{0U}];
    if (from_match.empty() || to_match.empty()) {
      return;
    }

    auto const t0 = std::chrono::steady_clock::now();
    auto const reference =
        route(params, w, l, search_profile::kCar, from_loc, to_loc, from_match,
              to_match, opt.max_cost_, direction::kForward, nullptr, nullptr,
              nullptr, routing_algorithm::kDijkstra);
    auto const t1 = std::chrono::steady_clock::now();
    auto const experiment =
        route(params, w, l, search_profile::kCar, from_loc, to_loc, from_match,
              to_match, opt.max_cost_, direction::kForward, nullptr, nullptr,
              nullptr, routing_algorithm::kCCH);
    auto const t2 = std::chrono::steady_clock::now();

    reference_time += (t1 - t0).count();
    experiment_time += (t2 - t1).count();
    ++n_compared;

    if (reference.has_value() == experiment.has_value() &&
        (!reference.has_value() ||
         reference->cost_ == experiment->cost_)) {
      ++n_congruent;
      return;
    }

    auto const guard = std::lock_guard{mutex};
    fmt::println("MISMATCH {} -> {}: dijkstra {} | cch {}",
                 w.node_to_osm_[ft.first], w.node_to_osm_[ft.second],
                 reference ? fmt::format("cost {} dist {:.1f}",
                                         reference->cost_, reference->dist_)
                           : "no result",
                 experiment ? fmt::format("cost {} dist {:.1f}",
                                          experiment->cost_, experiment->dist_)
                            : "no result");
  };

  utl::parallel_for(from_tos, single_run, utl::noop_progress_update{},
                    utl::parallel_error_strategy::QUIT_EXEC, opt.threads_);

  fmt::println("{}/{} congruent", n_congruent.load(), n_compared.load());
  fmt::println("dijkstra: {:.2f} ms/query, cch: {:.2f} ms/query, speedup {:.2f}",
               static_cast<double>(reference_time.load()) /
                   std::max(1U, n_compared.load()) / 1e6,
               static_cast<double>(experiment_time.load()) /
                   std::max(1U, n_compared.load()) / 1e6,
               static_cast<double>(reference_time.load()) /
                   std::max(std::int64_t{1}, experiment_time.load()));
  return n_congruent.load() == n_compared.load() ? 0 : 1;
}
