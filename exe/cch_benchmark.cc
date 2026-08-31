// Random-pair query benchmark for the CCH against a plain Dijkstra on the same
// graph. One algorithm per process so that the reported RSS is attributable.
//
// Percentiles are reported separately for queries that found a path and those
// that did not: at continental scale most random pairs are unreachable within
// the cost limit, and those queries have a very different cost profile, so
// pooling them hides the behaviour of both.
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <mutex>
#include <random>
#include <string>
#include <vector>

#include "fmt/core.h"

#include "conf/options_parser.h"

#include "utl/parallel_for.h"

#include "osr/location.h"
#include "osr/lookup.h"
#include "osr/routing/algorithms.h"
#include "osr/routing/parameters.h"
#include "osr/routing/profiles/car.h"
#include "osr/routing/route.h"
#include "osr/ways.h"

namespace fs = std::filesystem;
using namespace osr;

namespace {

struct settings : public conf::configuration {
  settings() : configuration{"Options"} {
    param(data_, "data,d", "osr data directory");
    param(algorithm_, "algorithm,a", "routing algorithm: cch or dijkstra");
    param(n_queries_, "queries,n", "number of random queries");
    param(threads_, "threads,t", "query threads");
    param(max_cost_, "max_cost,m", "cost limit in seconds");
    param(seed_, "seed", "seed for the random query set");
    param(bbox_, "bbox", "restrict both endpoints to a lat/lon box: "
                         "min_lat,min_lon,max_lat,max_lon");
    param(region_, "region", "label for the sampled region, used in output");
    param(sample_only_, "sample_only",
          "report bbox acceptance for the sampled query set and exit, without "
          "running any query");
  }

  fs::path data_{"osr"};
  std::string algorithm_{"cch"};
  unsigned n_queries_{1000U};
  unsigned threads_{8U};
  unsigned max_cost_{43200U};
  unsigned seed_{42U};
  std::string bbox_{};
  std::string region_{};
  bool sample_only_{false};
};

std::size_t read_status(char const* key) {
  auto in = std::ifstream{"/proc/self/status"};
  auto line = std::string{};
  while (std::getline(in, line)) {
    if (line.rfind(key, 0) == 0) {
      auto kb = std::size_t{0U};
      std::sscanf(line.c_str() + std::string_view{key}.size(), " %zu", &kb);
      return kb * 1024U;
    }
  }
  return 0U;
}

std::size_t peak_rss() { return read_status("VmHWM:"); }
std::size_t rss() { return read_status("VmRSS:"); }

void reset_peak() {
  auto out = std::ofstream{"/proc/self/clear_refs"};
  out << "5\n";
}

double mb(std::size_t const b) { return b / 1024.0 / 1024.0; }

// Matches `loc` but keeps only the candidates that resolve to `n`, so that the
// query really starts at the sampled node instead of wherever the matching
// would have put it.
match_result pinned(ways const& w,
                    lookup const& l,
                    location const& loc,
                    node_idx_t const n,
                    bool const reverse) {
  auto all = match_result{};
  l.match<car>(car::parameters{}, loc, reverse, direction::kForward, 100.0,
               nullptr, all);
  auto const m = all[match_idx_t{0U}];
  auto out = match_result{};
  out.start(m.lvl_);
  for (auto j = std::size_t{0U}; j != m.size(); ++j) {
    if (m.nodes_[j].left_.node_ == n || m.nodes_[j].right_.node_ == n) {
      out.add(m.dist_to_way_[j], m.way_[j], m.nodes_[j]);
    }
  }
  out.finish();
  return out;
}

void print_percentiles(std::string_view const tag,
                       std::string const& name,
                       std::string const& algorithm,
                       std::vector<double> const& sorted_ms) {
  if (sorted_ms.empty()) {
    fmt::println("PCTL{} {} | {} | n 0", tag, name, algorithm);
    return;
  }
  auto const q = [&](double const p) {
    return sorted_ms[static_cast<std::size_t>(
        std::clamp(p / 100.0 * static_cast<double>(sorted_ms.size() - 1U), 0.0,
                   static_cast<double>(sorted_ms.size() - 1U)))];
  };
  fmt::println("PCTL{} {} | {} | n {} | p25 {:.3f} | p50 {:.3f} | p75 {:.3f} | "
               "p90 {:.3f} | p99 {:.3f} | max {:.3f}",
               tag, name, algorithm, sorted_ms.size(), q(25), q(50), q(75),
               q(90), q(99), sorted_ms.back());
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

  if (!fs::is_directory(opt.data_)) {
    fmt::println("data directory {} not found", opt.data_.string());
    return 1;
  }
  if (opt.algorithm_ != "cch" && opt.algorithm_ != "dijkstra") {
    fmt::println("unknown algorithm {}, expected cch or dijkstra",
                 opt.algorithm_);
    return 1;
  }

  auto const name = opt.region_.empty()
                        ? opt.data_.filename().string()
                        : opt.data_.filename().string() + "@" + opt.region_;
  auto const algo = opt.algorithm_ == "cch" ? routing_algorithm::kCCH
                                            : routing_algorithm::kDijkstra;
  auto const max_cost = static_cast<cost_t>(opt.max_cost_);

  auto const rss_start = rss();
  auto const w = ways{opt.data_, cista::mmap::protection::READ};
  auto const l = lookup{w, opt.data_, cista::mmap::protection::READ};
  auto const params = profile_parameters{car::parameters{}};
  auto const rss_graph = rss();

  // Optional geographic restriction. Sampling the same box on two different
  // extracts is what makes them comparable: the query population is then held
  // fixed and only the size of the hierarchy underneath it differs.
  auto min_lat = -90.0, min_lon = -180.0, max_lat = 90.0, max_lon = 180.0;
  auto const has_bbox = !opt.bbox_.empty();
  if (has_bbox && std::sscanf(opt.bbox_.c_str(), "%lf,%lf,%lf,%lf", &min_lat,
                              &min_lon, &max_lat, &max_lon) != 4) {
    fmt::println("bad bbox {}, expected min_lat,min_lon,max_lat,max_lon",
                 opt.bbox_);
    return 1;
  }

  auto prng = std::mt19937{opt.seed_};
  auto distr =
      std::uniform_int_distribution<std::uint32_t>{0U, w.n_nodes() - 1U};

  auto const in_bbox = [&](node_idx_t const n) {
    auto const p = w.get_node_pos(n);
    return p.lat() >= min_lat && p.lat() <= max_lat && p.lng() >= min_lon &&
           p.lng() <= max_lon;
  };

  // Rejection sampling. The cap only trips for a box that holds (almost) no
  // nodes, which is a mistyped box rather than a result worth reporting.
  auto n_rejected = std::uint64_t{0U};
  auto const sample_node = [&](node_idx_t& out) {
    for (auto tries = 0U; tries != 1000000U; ++tries) {
      auto const n = node_idx_t{distr(prng)};
      if (!has_bbox || in_bbox(n)) {
        out = n;
        return true;
      }
      ++n_rejected;
    }
    return false;
  };

  auto from_tos = std::vector<std::pair<node_idx_t, node_idx_t>>{};
  from_tos.reserve(opt.n_queries_);
  for (auto i = 0U; i != opt.n_queries_; ++i) {
    auto a = node_idx_t{}, b = node_idx_t{};
    if (!sample_node(a) || !sample_node(b)) {
      fmt::println("bbox {} holds too few nodes to sample", opt.bbox_);
      return 1;
    }
    from_tos.emplace_back(a, b);
  }
  if (has_bbox) {
    // acceptance rate doubles as the share of the network inside the box
    fmt::println("BBOX {} | {} | box {} | accepted {} | rejected {} | "
                 "acceptance {:.4f}",
                 name, opt.algorithm_, opt.bbox_, 2U * opt.n_queries_,
                 n_rejected,
                 static_cast<double>(2U * opt.n_queries_) /
                     static_cast<double>(2U * opt.n_queries_ + n_rejected));
  }
  if (opt.sample_only_) {
    return 0;
  }

  struct sample {
    double total_ms_{0.0};
    bool found_{false};
  };
  auto samples = std::vector<sample>{};
  auto samples_mutex = std::mutex{};

  auto const one = [&](std::pair<node_idx_t, node_idx_t> const& ft,
                       std::atomic<std::int64_t>* time,
                       std::atomic<unsigned>* n_done,
                       std::atomic<unsigned>* n_found) {
    auto const from_loc = location{w.get_node_pos(ft.first)};
    auto const to_loc = location{w.get_node_pos(ft.second)};
    auto const fm = pinned(w, l, from_loc, ft.first, false);
    auto const tm = pinned(w, l, to_loc, ft.second, true);
    if (fm[match_idx_t{0U}].empty() || tm[match_idx_t{0U}].empty()) {
      return;  // not reachable by car at all, not part of the measurement
    }
    auto const t0 = std::chrono::steady_clock::now();
    auto const p = route(params, w, l, search_profile::kCar, from_loc, to_loc,
                         fm[match_idx_t{0U}], tm[match_idx_t{0U}], max_cost,
                         direction::kForward, nullptr, nullptr, nullptr, algo);
    auto const t1 = std::chrono::steady_clock::now();
    if (time == nullptr) {
      return;  // warm-up
    }
    *time += (t1 - t0).count();
    ++*n_done;
    *n_found += (p.has_value() ? 1U : 0U);
    auto const ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    auto const guard = std::lock_guard{samples_mutex};
    samples.emplace_back(sample{ms, p.has_value()});
  };

  // maps the cch and its metric and faults in the graph pages
  for (auto i = 0U; i != std::min<std::size_t>(20U, from_tos.size()); ++i) {
    one(from_tos[i], nullptr, nullptr, nullptr);
  }
  auto const rss_warm = rss();
  reset_peak();

  auto time = std::atomic<std::int64_t>{0};
  auto n_done = std::atomic<unsigned>{0U};
  auto n_found = std::atomic<unsigned>{0U};
  auto const t0 = std::chrono::steady_clock::now();
  utl::parallel_for(
      from_tos, [&](auto const& ft) { one(ft, &time, &n_done, &n_found); },
      utl::noop_progress_update{}, utl::parallel_error_strategy::QUIT_EXEC,
      opt.threads_);
  auto const wall =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - t0)
          .count();

  fmt::println("RESULT {} | {} | queries {} found {} | per_query_ms {:.3f} | "
               "wall_s {:.2f}",
               name, opt.algorithm_, n_done.load(), n_found.load(),
               static_cast<double>(time.load()) /
                   std::max(1U, n_done.load()) / 1e6,
               wall);

  std::sort(begin(samples), end(samples),
            [](sample const& a, sample const& b) {
              return a.total_ms_ < b.total_ms_;
            });
  auto all_ms = std::vector<double>{};
  auto found_ms = std::vector<double>{};
  auto not_found_ms = std::vector<double>{};
  for (auto const& s : samples) {
    all_ms.push_back(s.total_ms_);
    (s.found_ ? found_ms : not_found_ms).push_back(s.total_ms_);
  }
  print_percentiles("", name, opt.algorithm_, all_ms);
  print_percentiles("_FOUND", name, opt.algorithm_, found_ms);
  print_percentiles("_NOTFOUND", name, opt.algorithm_, not_found_ms);

  fmt::println("RESULT {} | {} | rss_start {:.1f} | rss_graph {:.1f} | "
               "rss_warm {:.1f} | peak_query {:.1f} | peak_total {:.1f}",
               name, opt.algorithm_, mb(rss_start), mb(rss_graph),
               mb(rss_warm), mb(peak_rss()), mb(rss()));
}
