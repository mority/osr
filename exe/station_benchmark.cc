// One-to-many to a real target set: every German (and nearby) long distance
// railway station, from a car-accessible start drawn uniformly at random inside
// Germany.
//
// This is the shape of query a journey planner actually issues for a first mile
// by car, and it differs from the synthetic one-to-many benchmarks in two ways
// that matter: the targets are a fixed, spatially spread set rather than a box
// around the source, and the caller reads all of them.
//
// One arm per invocation so that the algorithms cannot perturb each other's
// caches, with the same seed in each so all three answer identical queries.
// Costs are written per (source, target) so the arms can be diffed: a speedup
// only counts if the answers agree.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "fmt/core.h"
#include "fmt/std.h"

#include "conf/options_parser.h"

#include "geo/box.h"

#include "osr/location.h"
#include "osr/lookup.h"
#include "osr/routing/algorithms.h"
#include "osr/routing/profiles/car.h"
#include "osr/routing/route.h"
#include "osr/ways.h"

namespace fs = std::filesystem;
using namespace osr;

namespace {

struct settings : public conf::configuration {
  settings() : configuration{"Options"} {
    param(data_, "data,d", "osr data directory");
    param(stations_, "stations,s", "csv: city,station,lat,lon");
    param(bbox_, "bbox", "source sampling box min_lat,min_lon,max_lat,max_lon");
    param(n_sources_, "sources,n", "number of random sources");
    param(max_cost_, "max_cost,m", "cost limit in seconds");
    param(algo_, "algorithm,a", "dijkstra or cch");
    param(match_dist_, "match_distance", "max matching distance, metres");
    param(out_, "out,o", "write per-query costs here");
    param(random_targets_, "random_targets",
          "instead of the station csv, draw this many targets from a box "
          "around each source -- the shape an ODM offset query has");
    param(target_radius_km_, "target_radius", "half-width of that box, km");
    param(seed_, "seed", "seed for the source sample");
  }
  fs::path data_{"osr"}, stations_{"stations.csv"}, out_{};
  std::string bbox_{"47.27,5.86,55.15,15.06"}, algo_{"cch"};
  unsigned n_sources_{100U}, max_cost_{3600U}, seed_{20260902U};
  double match_dist_{250.0};
  unsigned random_targets_{0U};
  double target_radius_km_{60.0};
};

std::vector<std::pair<std::string, location>> read_stations(fs::path const& p) {
  auto out = std::vector<std::pair<std::string, location>>{};
  auto f = std::ifstream{p};
  auto line = std::string{};
  std::getline(f, line);  // header
  while (std::getline(f, line)) {
    if (line.empty()) {
      continue;
    }
    auto ss = std::istringstream{line};
    auto city = std::string{}, name = std::string{}, lat = std::string{},
         lon = std::string{};
    if (!std::getline(ss, city, ',') || !std::getline(ss, name, ',') ||
        !std::getline(ss, lat, ',') || !std::getline(ss, lon, ',')) {
      continue;
    }
    out.emplace_back(
        name, location{geo::latlng{std::stod(lat), std::stod(lon)}, level_t{}});
  }
  return out;
}

using clk = std::chrono::steady_clock;
double ms(clk::time_point a, clk::time_point b) {
  return std::chrono::duration<double, std::milli>(b - a).count();
}
double pct(std::vector<double> v, double p) {
  if (v.empty()) {
    return 0.0;
  }
  std::sort(begin(v), end(v));
  return v[static_cast<std::size_t>(
      std::clamp(p / 100.0 * static_cast<double>(v.size() - 1U), 0.0,
                 static_cast<double>(v.size() - 1U)))];
}

}  // namespace

int main(int ac, char const** av) {
  auto opt = settings{};
  auto parser = conf::options_parser{{&opt}};
  parser.read_command_line_args(ac, av);
  if (parser.help()) {
    parser.print_help(std::cout);
    return 0;
  }
  parser.print_used(std::cout);

  auto min_lat = 0.0, min_lon = 0.0, max_lat = 0.0, max_lon = 0.0;
  if (std::sscanf(opt.bbox_.c_str(), "%lf,%lf,%lf,%lf", &min_lat, &min_lon,
                  &max_lat, &max_lon) != 4) {
    fmt::println("bad bbox {}", opt.bbox_);
    return 1;
  }
  auto const algo = opt.algo_ == "dijkstra" ? routing_algorithm::kDijkstra
                                            : routing_algorithm::kCCH;

  auto const stations = read_stations(opt.stations_);
  if (stations.empty()) {
    fmt::println("no stations read from {}", opt.stations_);
    return 1;
  }
  auto const w = ways{opt.data_, cista::mmap::protection::READ};
  auto const l = lookup{w, opt.data_, cista::mmap::protection::READ};
  auto const params = profile_parameters{car::parameters{}};
  auto const max_cost = static_cast<cost_t>(opt.max_cost_);

  auto targets = std::vector<location>{};
  targets.reserve(stations.size());
  for (auto const& [n, loc] : stations) {
    targets.push_back(loc);
  }

  // The target set is fixed, so it is matched once and reused by every source.
  // Matching is identical work in all three arms and is reported separately so
  // it cannot be mistaken for a difference between the algorithms.
  auto const t_match0 = clk::now();
  auto to_match = match_result{};
  for (auto const& t : targets) {
    l.match(params, t, true, direction::kForward, opt.match_dist_, nullptr,
            search_profile::kCar, {}, to_match);
  }
  auto const match_ms = ms(t_match0, clk::now());

  auto prng = std::mt19937{opt.seed_};
  auto distr = std::uniform_int_distribution<std::uint32_t>{
      0U, static_cast<std::uint32_t>(w.n_nodes() - 1U)};

  auto out = std::ofstream{};
  if (!opt.out_.empty()) {
    out.open(opt.out_);
    out << "source,target,cost\n";
  }

  auto route_ms = std::vector<double>{};
  auto reached = std::vector<double>{};
  auto n_done = 0U;
  for (auto tries = 0U; n_done != opt.n_sources_ && tries != 100000000U;
       ++tries) {
    auto const n = node_idx_t{distr(prng)};
    auto const p = w.get_node_pos(n);
    if (p.lat() < min_lat || p.lat() > max_lat || p.lng() < min_lon ||
        p.lng() > max_lon) {
      continue;
    }
    auto const from = location{p, level_t{}};
    auto from_match = match_result{};
    l.match(params, from, false, direction::kForward, opt.match_dist_, nullptr,
            search_profile::kCar, {}, from_match);
    if (from_match[match_idx_t{0U}].empty()) {
      continue;  // not on the car network
    }

    // With random targets the set is per source: sampled from a box around it,
    // which is the shape an ODM offset query has. Matching stays outside the
    // timer, exactly as with the fixed station set, so the arms differ only in
    // the routing.
    auto per_src_targets = targets;
    auto per_src_match = match_result{};
    if (opt.random_targets_ != 0U) {
      auto const dlat = opt.target_radius_km_ / 111.2;
      auto const dlon =
          dlat / std::max(0.2, std::cos(p.lat() * 3.14159265 / 180.0));
      auto box = geo::box{};
      box.extend(geo::latlng{p.lat() - dlat, p.lng() - dlon});
      box.extend(geo::latlng{p.lat() + dlat, p.lng() + dlon});
      auto cand = std::vector<geo::latlng>{};
      l.find(box, [&](way_idx_t const way) {
        for (auto const nn : w.r_->way_nodes_[way]) {
          cand.push_back(w.get_node_pos(nn));
        }
      });
      if (cand.size() < opt.random_targets_) {
        continue;
      }
      std::shuffle(begin(cand), end(cand), prng);
      per_src_targets.clear();
      per_src_match = match_result{};
      for (auto k = 0U; k != opt.random_targets_; ++k) {
        per_src_targets.push_back(location{cand[k], level_t{}});
      }
      for (auto const& t : per_src_targets) {
        l.match(params, t, true, direction::kForward, opt.match_dist_, nullptr,
                search_profile::kCar, {}, per_src_match);
      }
    }
    auto const& use_targets =
        opt.random_targets_ != 0U ? per_src_targets : targets;
    auto const& use_match = opt.random_targets_ != 0U ? per_src_match : to_match;

    auto const t0 = clk::now();
    auto const paths = route(params, w, l, search_profile::kCar, from, use_targets,
                             from_match[match_idx_t{0U}], use_match, max_cost,
                             direction::kForward, nullptr, nullptr, nullptr,
                             [](path const&) { return false; }, std::nullopt,
                             algo);
    route_ms.push_back(ms(t0, clk::now()));

    auto n_reached = 0U;
    for (auto i = std::size_t{0U}; i != paths.size(); ++i) {
      auto const c = paths[i].has_value() ? paths[i]->cost_ : kInfeasible;
      if (paths[i].has_value()) {
        ++n_reached;
      }
      if (out.is_open()) {
        out << n_done << ',' << i << ',' << c << '\n';
      }
    }
    reached.push_back(static_cast<double>(n_reached));
    ++n_done;
  }

  // The first query faults in the hierarchy and the metric from disk -- on the
  // planet that is tens of gigabytes and costs two orders of magnitude more
  // than a warm query. It is discarded rather than allowed to dominate the mean.
  auto cold_ms = 0.0;
  if (route_ms.size() > 1U) {
    cold_ms = route_ms.front();
    route_ms.erase(begin(route_ms));
  }

  fmt::println(
      "STATIONS {} | targets {} | sources {} | match_ms {:.1f} | "
      "cold_first_ms {:.1f} | reached_median {:.0f}",
      opt.algo_, targets.size(), n_done, match_ms, cold_ms, pct(reached, 50));
  fmt::println(
      "ROUTE {} | n {} | p25 {:.2f} | p50 {:.2f} | p75 {:.2f} | p90 {:.2f} | "
      "p99 {:.2f} | mean {:.2f} | total_s {:.1f}",
      opt.algo_, route_ms.size(), pct(route_ms, 25), pct(route_ms, 50),
      pct(route_ms, 75), pct(route_ms, 90), pct(route_ms, 99),
      route_ms.empty() ? 0.0
                       : std::accumulate(begin(route_ms), end(route_ms), 0.0) /
                             static_cast<double>(route_ms.size()),
      std::accumulate(begin(route_ms), end(route_ms), 0.0) / 1000.0);
  return 0;
}
