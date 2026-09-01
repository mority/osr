// Decides whether RPHAST can beat the bounded Dijkstra for car first/last mile
// to public transit stations, before any of it is implemented.
//
// The two quantities that settle it:
//
//   |T'|      the RPHAST target selection closure -- the upward closure of the
//             station set over the reverse of the downward graph. The scanning
//             phase of an RPHAST query visits exactly these, so |T'| is what a
//             query costs once the targets are known.
//   |labels|  the nodes the cost-bounded Dijkstra labels today. That is the
//             whole ball of radius `max` around the source, independent of how
//             many stations were asked for.
//
// The ratio is the headroom. It is not the speedup: RPHAST scans sequentially
// in rank order while Dijkstra runs a priority queue with random access, so the
// constant factors differ substantially in RPHAST's favour -- but only once the
// restricted graph is laid out contiguously, which is the second step. If |T'|
// is not clearly below the Dijkstra ball, the first step alone cannot win.
//
// Target selection is metric independent: it walks the CCH topology only, so
// one closure stays valid across customizations.
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "fmt/core.h"

#include "conf/options_parser.h"

#include "geo/latlng.h"

#include "osr/location.h"
#include "osr/lookup.h"
#include "osr/routing/algorithms.h"
#include "osr/routing/cch/cch.h"
#include "osr/routing/dijkstra.h"
#include "osr/routing/profiles/car.h"
#include "osr/routing/route.h"
#include "osr/ways.h"

namespace fs = std::filesystem;
using namespace osr;

namespace {

struct settings : public conf::configuration {
  settings() : configuration{"Options"} {
    param(data_, "data,d", "osr data directory");
    param(stations_, "stations,s", "csv of lat,lon station coordinates");
    param(n_sources_, "sources,n", "number of random source coordinates");
    param(budgets_, "budgets,b", "comma separated cost budgets in seconds");
    param(kmh_, "kmh", "assumed average speed, sets the station search radius");
    param(seed_, "seed", "seed for the random source set");
  }

  fs::path data_{"osr"};
  fs::path stations_{"stations.csv"};
  unsigned n_sources_{25U};
  std::string budgets_{"900,1800,3600"};
  double kmh_{72.0};
  unsigned seed_{42U};
};

double median(std::vector<double> v) {
  if (v.empty()) {
    return 0.0;
  }
  std::sort(begin(v), end(v));
  return v[v.size() / 2U];
}

double median(std::vector<std::size_t> const& v) {
  auto d = std::vector<double>{};
  d.reserve(v.size());
  for (auto const x : v) {
    d.push_back(static_cast<double>(x));
  }
  return median(d);
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
  if (!cch::exists(opt.data_)) {
    fmt::println("no cch in {}", opt.data_.string());
    return 1;
  }

  auto const w = ways{opt.data_, cista::mmap::protection::READ};
  auto const l = lookup{w, opt.data_, cista::mmap::protection::READ};
  auto const c = cch::read(opt.data_);
  auto const params = car::parameters{};

  auto budgets = std::vector<cost_t>{};
  {
    auto s = std::string{};
    auto in = std::stringstream{opt.budgets_};
    while (std::getline(in, s, ',')) {
      budgets.push_back(static_cast<cost_t>(std::stoul(s)));
    }
  }

  // ---- stations -----------------------------------------------------------
  auto st_pos = std::vector<geo::latlng>{};
  {
    auto in = std::ifstream{opt.stations_};
    if (!in) {
      fmt::println("cannot open {}", opt.stations_.string());
      return 1;
    }
    auto line = std::string{};
    while (std::getline(in, line)) {
      auto const comma = line.find(',');
      if (comma == std::string::npos) {
        continue;
      }
      st_pos.push_back(geo::latlng{std::stod(line.substr(0U, comma)),
                                   std::stod(line.substr(comma + 1U))});
    }
  }
  fmt::println("stations read: {}", st_pos.size());

  // Match every station once. Doing it per source would dominate the run and
  // says nothing about either algorithm.
  auto st_ranks = std::vector<std::vector<cch_rank_t>>{};
  st_ranks.resize(st_pos.size());
  auto n_matched = 0U;
  {
    auto const t0 = std::chrono::steady_clock::now();
    for (auto i = std::size_t{0U}; i != st_pos.size(); ++i) {
      auto all = match_result{};
      l.match<car>(params, location{st_pos[i], level_t{}}, true,
                   direction::kForward, 100.0, nullptr, all);
      if (all.empty()) {
        continue;
      }
      auto const m = all[match_idx_t{0U}];
      auto& out = st_ranks[i];
      for (auto j = std::size_t{0U}; j != m.size(); ++j) {
        for (auto const n :
             {m.nodes_[j].left_.node_, m.nodes_[j].right_.node_}) {
          if (n == node_idx_t::invalid() || n >= c->rank_.size()) {
            continue;
          }
          auto const r = c->rank_[n];
          if (r != cch_rank_t::invalid() &&
              std::find(begin(out), end(out), r) == end(out)) {
            out.push_back(r);
          }
        }
      }
      if (!out.empty()) {
        ++n_matched;
      }
    }
    auto const secs = std::chrono::duration<double>(
                          std::chrono::steady_clock::now() - t0)
                          .count();
    fmt::println("stations matched to the car network: {} ({:.1f}s)", n_matched,
                 secs);
  }

  // ---- sources ------------------------------------------------------------
  auto prng = std::mt19937{opt.seed_};
  auto distr =
      std::uniform_int_distribution<std::uint32_t>{0U, w.n_nodes() - 1U};
  auto sources = std::vector<node_idx_t>{};
  while (sources.size() != opt.n_sources_) {
    auto const n = node_idx_t{distr(prng)};
    if (c->rank_[n] != cch_rank_t::invalid()) {
      sources.push_back(n);
    }
  }

  // Generation stamps rather than a cleared bitset: the closure touches a small
  // part of a 15M entry array and clearing it per query would dominate.
  auto stamp = std::vector<std::uint32_t>(c->n_ranks(), 0U);
  auto gen = std::uint32_t{0U};
  auto stack = std::vector<cch_rank_t>{};

  auto d = std::make_unique<dijkstra<car>>();

  fmt::println("");
  fmt::println("{:>7} {:>7} {:>8} {:>10} {:>10} {:>8} {:>9} {:>9} {:>9} "
               "{:>8} {:>7}",
               "budget", "radius", "stations", "|T'|", "labels", "ratio",
               "sel_ms", "dij_ms", "rph_ms", "speedup", "n");
  fmt::println("{}", std::string(104U, '-'));

  for (auto const budget : budgets) {
    auto const radius_m = opt.kmh_ * 1000.0 / 3600.0 * budget;

    auto v_targets = std::vector<std::size_t>{};
    auto v_tprime = std::vector<std::size_t>{};
    auto v_labels = std::vector<std::size_t>{};
    auto v_sel_ms = std::vector<double>{};
    auto v_dij_ms = std::vector<double>{};
    auto v_rt_dij_ms = std::vector<double>{};
    auto v_rt_rph_ms = std::vector<double>{};
    auto n_mismatch = 0U, n_rphast_only = 0U, n_dij_only = 0U,
         n_rphast_lower = 0U, n_rphast_higher = 0U;

    for (auto const src : sources) {
      auto const src_pos = w.get_node_pos(src).as_latlng();
      auto to_locs = std::vector<location>{};
      for (auto ii = std::size_t{0U}; ii != st_pos.size(); ++ii) {
        if (!st_ranks[ii].empty() &&
            geo::distance(src_pos, st_pos[ii]) <= radius_m) {
          to_locs.push_back(location{st_pos[ii], level_t{}});
        }
      }

      // --- target selection ------------------------------------------------
      auto const t0 = std::chrono::steady_clock::now();
      ++gen;
      stack.clear();
      auto n_targets = std::size_t{0U};
      for (auto i = std::size_t{0U}; i != st_pos.size(); ++i) {
        if (st_ranks[i].empty() ||
            geo::distance(src_pos, st_pos[i]) > radius_m) {
          continue;
        }
        ++n_targets;
        for (auto const r : st_ranks[i]) {
          if (stamp[to_idx(r)] != gen) {
            stamp[to_idx(r)] = gen;
            stack.push_back(r);
          }
        }
      }
      auto n_tprime = stack.size();
      while (!stack.empty()) {
        auto const u = stack.back();
        stack.pop_back();
        // downward arcs into u come from the slots where u is the lower node
        for (auto s = c->upper_begin(u); s != c->upper_end(u); ++s) {
          if (c->dn_entries(s).empty()) {
            continue;
          }
          auto const h = c->adj_head_[s];
          if (stamp[to_idx(h)] != gen) {
            stamp[to_idx(h)] = gen;
            stack.push_back(h);
            ++n_tprime;
          }
        }
      }
      auto const sel_ms =
          std::chrono::duration<double, std::milli>(
              std::chrono::steady_clock::now() - t0)
              .count();

      // --- bounded dijkstra, what we do today ------------------------------
      auto from_all = match_result{};
      l.match<car>(params, location{src_pos, level_t{}}, false,
                   direction::kForward, 100.0, nullptr, from_all);
      if (from_all.empty()) {
        continue;
      }
      auto const fm = from_all[match_idx_t{0U}];
      if (fm.empty()) {
        continue;
      }

      auto const t1 = std::chrono::steady_clock::now();
      d->reset(budget);
      for (auto i = std::size_t{0U}; i != fm.size(); ++i) {
        auto const start_way = fm.way_[i];
        auto const start_left = fm.left(i);
        auto const start_right = fm.right(i);
        for (auto const* nc : {&start_left, &start_right}) {
          if (!nc->valid() || nc->cost_ >= budget) {
            continue;
          }
          auto const sc = car::way_cost(
              params, *w.r_, w.timezones_, start_way,
              w.r_->way_properties_[start_way],
              flip(direction::kForward, nc->way_dir_),
              static_cast<distance_t>(nc->dist_to_node_), std::nullopt,
              duration_t{0}, direction::kForward);
          if (sc.cost_ == kInfeasible || sc.cost_ >= budget) {
            continue;
          }
          car::resolve_start_node(
              *w.r_, start_way, nc->node_, level_t{}, direction::kForward,
              [&](auto const node) {
                auto lab = car::label{node, sc.cost_};
                lab.track(lab, *w.r_, start_way, node.get_node(), false);
                d->add_start(w, lab, sc.duration_);
              });
        }
      }
      d->run(params, w, *w.r_, budget, std::nullopt, nullptr, nullptr, nullptr,
             direction::kForward);
      auto const dij_ms =
          std::chrono::duration<double, std::milli>(
              std::chrono::steady_clock::now() - t1)
              .count();

      // --- the two one-to-many implementations, same matched candidates ----
      if (!to_locs.empty()) {
        auto to_match = match_result{};
        for (auto const& x : to_locs) {
          l.match<car>(params, x, true, direction::kForward, 100.0, nullptr,
                       to_match);
        }
        auto const from_loc = location{src_pos, level_t{}};

        auto const t2 = std::chrono::steady_clock::now();
        auto const ref = route(params, w, l, search_profile::kCar, from_loc,
                               to_locs, fm, to_match, budget,
                               direction::kForward, nullptr, nullptr, nullptr,
                               [](path const&) { return false; }, std::nullopt,
                               routing_algorithm::kDijkstra);
        auto const t3 = std::chrono::steady_clock::now();
        auto const exp = route(params, w, l, search_profile::kCar, from_loc,
                               to_locs, fm, to_match, budget,
                               direction::kForward, nullptr, nullptr, nullptr,
                               [](path const&) { return false; }, std::nullopt,
                               routing_algorithm::kCCH);
        auto const t4 = std::chrono::steady_clock::now();

        for (auto z = std::size_t{0U}; z != ref.size(); ++z) {
          if (ref[z].has_value() != exp[z].has_value()) {
            ++n_mismatch;
            if (exp[z].has_value()) {
              ++n_rphast_only;
            } else {
              ++n_dij_only;
            }
          } else if (ref[z].has_value() && ref[z]->cost_ != exp[z]->cost_) {
            ++n_mismatch;
            if (exp[z]->cost_ < ref[z]->cost_) {
              ++n_rphast_lower;
            } else {
              ++n_rphast_higher;
            }
          }
        }
        v_rt_dij_ms.push_back(
            std::chrono::duration<double, std::milli>(t3 - t2).count());
        v_rt_rph_ms.push_back(
            std::chrono::duration<double, std::milli>(t4 - t3).count());
      }

      v_targets.push_back(n_targets);
      v_tprime.push_back(n_tprime);
      v_labels.push_back(d->cost_.size());
      v_sel_ms.push_back(sel_ms);
      v_dij_ms.push_back(dij_ms);
    }

    auto const mt = median(v_tprime);
    auto const ml = median(v_labels);
    auto const rd = median(v_rt_dij_ms);
    auto const rr = median(v_rt_rph_ms);
    fmt::println(
        "{:>6}s {:>6.0f}km {:>8.0f} {:>10.0f} {:>10.0f} {:>7.2f}x {:>9.2f} "
        "{:>9.1f} {:>9.1f} {:>7.2f}x {:>7}",
        budget, radius_m / 1000.0, median(v_targets), mt, ml,
        ml > 0.0 ? ml / std::max(1.0, mt) : 0.0, median(v_sel_ms),
        rd, rr, rr > 0.0 ? rd / rr : 0.0, v_tprime.size());
    if (n_mismatch != 0U) {
      fmt::println("  !! {} mismatches: rphast_only {} dij_only {} "
                   "rphast_lower {} rphast_higher {}",
                   n_mismatch, n_rphast_only, n_dij_only, n_rphast_lower,
                   n_rphast_higher);
    }
  }

  fmt::println("");
  fmt::println("ratio = dijkstra labels / |T'|, the structural headroom.");
  fmt::println("dij_ms / rph_ms are the two one-to-many route() calls on the");
  fmt::println("same matched candidates; speedup is their quotient.");
  return 0;
}
