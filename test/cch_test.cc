#ifdef _WIN32
#include "windows.h"
#endif

#include "gtest/gtest.h"

#include <filesystem>
#include <queue>
#include <random>

#include "fmt/core.h"

#include "utl/timer.h"

#include "osr/extract/extract.h"
#include "osr/routing/cch/build.h"
#include "osr/routing/cch/cch.h"
#include "osr/routing/cch/customize.h"
#include "osr/routing/cch/query.h"
#include "osr/routing/cch/turns.h"
#include "osr/routing/dijkstra.h"
#include "osr/routing/profiles/car.h"
#include "osr/ways.h"

namespace fs = std::filesystem;
using namespace osr;

namespace {

void load(std::string_view raw_data, std::string_view data_dir) {
  if (fs::exists(raw_data)) {
    auto const p = fs::path{data_dir};
    auto ec = std::error_code{};
    fs::remove_all(p, ec);
    fs::create_directories(p, ec);
    osr::extract(false, raw_data, data_dir, fs::path{});
  }
}

// Plain Dijkstra on the CCH graph, using *all* arcs in both directions. If
// this matches the reference but the up/down search does not, the hierarchy
// (contraction / order) is broken; if it does not match either, the arc set or
// the customization is wrong.
cost_t cch_flat_dijkstra(ways const& w,
                         cch const& c,
                         cch_metric const& m,
                         node_idx_t const from,
                         port_t const from_port,
                         node_idx_t const to,
                         port_t const to_port,
                         cost_t const max,
                         std::vector<std::pair<std::uint32_t, port_t>>* path =
                             nullptr) {
  auto const& r = *w.r_;
  auto const params = car::parameters{};

  struct state {
    cost_t cost_;
    cch_rank_t rank_;
    port_t port_;
    bool operator>(state const& o) const { return cost_ > o.cost_; }
  };

  auto dist = ankerl::unordered_dense::map<std::uint64_t, cost_t>{};
  auto const key = [](cch_rank_t const x, port_t const p) {
    return (static_cast<std::uint64_t>(to_idx(x)) << 8U) | p;
  };
  auto pq = std::priority_queue<state, std::vector<state>, std::greater<>>{};
  auto pred = ankerl::unordered_dense::map<std::uint64_t, std::uint64_t>{};

  dist[key(c.rank_[from], from_port)] = 0U;
  pq.push(state{0U, c.rank_[from], from_port});

  while (!pq.empty()) {
    auto const s = pq.top();
    pq.pop();
    auto const it = dist.find(key(s.rank_, s.port_));
    if (it == end(dist) || it->second < s.cost_) {
      continue;
    }
    if (c.order_[s.rank_] == to && s.port_ == to_port) {
      if (path != nullptr) {
        auto k = key(s.rank_, s.port_);
        while (true) {
          path->emplace_back(static_cast<std::uint32_t>(k >> 8U),
                             static_cast<port_t>(k & 0xFFU));
          auto const p = pred.find(k);
          if (p == end(pred)) {
            break;
          }
          k = p->second;
        }
        std::reverse(begin(*path), end(*path));
      }
      return s.cost_;
    }
    auto const n = c.order_[s.rank_];

    auto const relax = [&](cch_rank_t const target, port_t const entry_port,
                           port_t const exit_port, cost_t const weight) {
      if (weight == kInfeasible) {
        return;
      }
      auto const turn =
          cch_turn_cost<car>(params, r, w.timezones_, n, s.port_, entry_port);
      if (turn == kInfeasible) {
        return;
      }
      auto const cost = s.cost_ + turn + weight;
      if (cost >= max) {
        return;
      }
      auto const k = key(target, exit_port);
      auto const found = dist.find(k);
      if (found == end(dist) || found->second > cost) {
        dist[k] = cost;
        pred[k] = key(s.rank_, s.port_);
        pq.push(state{cost, target, exit_port});
      }
    };

    for (auto slot = c.upper_begin(s.rank_); slot != c.upper_end(s.rank_);
         ++slot) {
      auto const entries = c.up_entries(slot);
      for (auto k = std::size_t{0U}; k != entries.size(); ++k) {
        relax(c.adj_head_[slot], entries[k].entry_, entries[k].exit_,
              m.up(c.up_ofs_[slot] + k));
      }
    }
    for (auto const slot : c.lower_slots(s.rank_)) {
      auto const entries = c.dn_entries(slot);
      for (auto k = std::size_t{0U}; k != entries.size(); ++k) {
        relax(c.tail(slot), entries[k].entry_, entries[k].exit_,
              m.dn(c.dn_ofs_[slot] + k));
      }
    }
  }
  return kInfeasible;
}

// Walks the reference Dijkstra path and checks that every step is present in
// the CCH graph as an original edge with the same weight.
void check_dijkstra_path(ways const& w,
                         cch const& c,
                         cch_metric const& m,
                         dijkstra<car> const& d,
                         car::node const dest) {
  auto const& r = *w.r_;
  auto const params = car::parameters{};

  auto n = dest;
  while (true) {
    auto const& e = d.cost_.at(n.get_key());
    auto const pred = e.pred(n);
    if (!pred.has_value()) {
      break;
    }
    auto const step = e.cost(n) - d.get_cost(*pred);
    auto const in_port = make_port(pred->way_, pred->dir_);
    auto const out_port = make_port(n.way_, n.dir_);

    if (!c.contains(pred->n_)) {
      fmt::println("  node {} not in cch", w.node_to_osm_[pred->n_]);
      return;
    }

    auto found = false;
    for_each_edge(r, pred->n_,
                  [&](node_idx_t const v, way_idx_t const way,
                      port_t const tail_port, port_t const head_port,
                      distance_t const dist, std::uint16_t, std::uint16_t) {
                    if (found || v != n.n_ || head_port != out_port) {
                      return;
                    }
                    found = true;
                    auto const edge_cost =
                        cch_edge_cost<car>(params, w, way, v, tail_port, dist);
                    auto const turn = cch_turn_cost<car>(
                        params, r, w.timezones_, pred->n_, in_port, tail_port);
                    auto const up = c.rank_[pred->n_] < c.rank_[v];
                    auto const slot =
                        up ? c.find_slot(c.rank_[pred->n_], c.rank_[v])
                           : c.find_slot(c.rank_[v], c.rank_[pred->n_]);
                    auto arc_cost = kInfeasible;
                    if (slot != cch::kNoSlot) {
                      auto const entries =
                          up ? c.up_entries(slot) : c.dn_entries(slot);
                      auto const idx = cch::find_entry(
                          entries, cch_entry{tail_port, head_port});
                      if (idx != std::numeric_limits<cch_entry_idx_t>::max()) {
                        arc_cost = up ? m.up(c.up_ofs_[slot] + idx)
                                      : m.dn(c.dn_ofs_[slot] + idx);
                      }
                    }
                    if (arc_cost > edge_cost || turn + edge_cost != step) {
                      fmt::println(
                          "  bad step {} -> {} way {}: step={} turn={} "
                          "edge={} arc={} slot={}",
                          w.node_to_osm_[pred->n_], w.node_to_osm_[v],
                          w.way_osm_idx_[way], step, turn, edge_cost, arc_cost,
                          slot);
                    }
                  });
    if (!found) {
      fmt::println("  no original edge {} -> {} (port {})",
                   w.node_to_osm_[pred->n_], w.node_to_osm_[n.n_], out_port);
    }
    n = *pred;
  }
}

// Compares the plain node-to-node distances of the CCH against a Dijkstra on
// the same state space, i.e. without the map matching / first mile logic.
void compare_core(ways const& w,
                  cch const& c,
                  cch_metric const& m,
                  unsigned const n_samples,
                  cost_t const max) {
  auto const& r = *w.r_;
  auto const params = car::parameters{};

  auto prng = std::mt19937{42U};
  auto distr =
      std::uniform_int_distribution<std::uint32_t>{0U, c.n_ranks() - 1U};

  auto s = cch_search<car>{};
  auto d = dijkstra<car>{};

  auto n_ok = 0U;
  auto n_checked = 0U;
  auto settled_cch = std::size_t{0U};
  auto settled_dijkstra = std::size_t{0U};
  auto time_cch = std::chrono::steady_clock::duration{};
  auto time_dijkstra = std::chrono::steady_clock::duration{};
  for (auto i = 0U; i != n_samples; ++i) {
    auto const from = c.order_[cch_rank_t{distr(prng)}];
    auto const to = c.order_[cch_rank_t{distr(prng)}];
    if (from == to) {
      continue;
    }
    auto const from_port = port_t{0U};
    auto const to_port = port_t{0U};

    auto const t0 = std::chrono::steady_clock::now();
    d.reset(max);
    d.add_start(w, car::label{car::node{from, port_way_pos(from_port),
                                        port_dir(from_port)},
                              0U});
    d.run(params, w, r, max, std::nullopt, nullptr, nullptr, nullptr,
          direction::kForward);
    auto const expected = d.get_cost(
        car::node{to, port_way_pos(to_port), port_dir(to_port)});
    auto const t1 = std::chrono::steady_clock::now();

    s.clear();
    s.add_start(c, from, from_port, 0U);
    s.add_target(c, to, to_port, 0U);
    s.run(params, w, c, m, max);
    auto const actual = s.found() ? s.best() : kInfeasible;
    auto const t2 = std::chrono::steady_clock::now();

    time_dijkstra += t1 - t0;
    time_cch += t2 - t1;
    settled_cch += s.n_settled_;
    settled_dijkstra += d.cost_.size();

    ++n_checked;
    if (expected == actual) {
      ++n_ok;
    } else if (n_ok + 8U > n_checked) {
      auto p = std::vector<std::pair<std::uint32_t, port_t>>{};
      auto const flat = cch_flat_dijkstra(w, c, m, from, from_port, to,
                                          to_port, max, &p);
      fmt::println("MISMATCH {} ({}) -> {} ({}): dijkstra={} cch={} flat={}",
                   w.node_to_osm_[from], to_idx(c.rank_[from]),
                   w.node_to_osm_[to], to_idx(c.rank_[to]), expected, actual,
                   flat);
      auto ranks = std::string{};
      for (auto const& [rk, pt] : p) {
        ranks += fmt::format("{}/{} ", rk, pt);
      }
      fmt::println("  flat path ranks: {}", ranks);
      check_dijkstra_path(w, c, m,  d,
                          car::node{to, port_way_pos(to_port),
                                    port_dir(to_port)});
    }
  }
  fmt::println(
      "core: {}/{} match | settled: dijkstra {} cch {} | avg query: "
      "dijkstra {}us cch {}us | speedup {:.1f}",
      n_ok, n_checked, settled_dijkstra / std::max(n_checked, 1U),
      settled_cch / std::max(n_checked, 1U),
      std::chrono::duration_cast<std::chrono::microseconds>(time_dijkstra)
              .count() /
          std::max(n_checked, 1U),
      std::chrono::duration_cast<std::chrono::microseconds>(time_cch).count() /
          std::max(n_checked, 1U),
      static_cast<double>(time_dijkstra.count()) /
          static_cast<double>(std::max(time_cch.count(), std::int64_t{1})));
  EXPECT_EQ(n_checked, n_ok);
}

void report(ways const& w, cch const& c, cch_metric const& m) {
  auto max_ways = std::size_t{0U};
  auto n_over_16 = std::size_t{0U};
  for (auto const n : c.order_) {
    auto const d = w.r_->node_ways_[n].size();
    max_ways = std::max(max_ways, static_cast<std::size_t>(d));
    n_over_16 += (d > 16U ? 1U : 0U);
  }
  fmt::println("max node ways: {} ({} nodes with > 16)", max_ways, n_over_16);

  auto n_inf = std::size_t{0U};
  for (auto const x : m.up_) {
    n_inf += (x == kInfeasible ? 1U : 0U);
  }
  for (auto const x : m.dn_) {
    n_inf += (x == kInfeasible ? 1U : 0U);
  }
  fmt::println(
      "cch: {} nodes, {} slots, {} up entries, {} dn entries, {} infeasible, "
      "{} triangles applied, {} missing slots, {} missing entries",
      c.n_ranks(), c.n_slots(), c.up_.size(), c.dn_.size(), n_inf, m.applied_,
      m.missing_slots_, m.missing_entries_);
}

void run_core(std::string_view raw_data,
              std::string_view data_dir,
              unsigned const n_samples,
              cost_t const max) {
  load(raw_data, data_dir);
  auto const w = osr::ways{data_dir, cista::mmap::protection::READ};
  auto const c = cch::read(data_dir);
  auto m = cch_metric{};
  {
    auto const t = utl::scoped_timer{"customize"};
    customize<car>(car::parameters{}, w, *c, m);
  }
  report(w, *c, m);
  compare_core(w, *c, m, n_samples, max);
}

}  // namespace

TEST(cch, monaco_core) {
  if (!fs::exists("test/monaco.osm.pbf") && !fs::exists("test/monaco")) {
    GTEST_SKIP() << "monaco not found";
  }
  run_core("test/monaco.osm.pbf", "test/monaco", 2000U, 2 * 3600U);
}

TEST(cch, hamburg_core) {
  if (!fs::exists("test/hamburg.osm.pbf") && !fs::exists("test/hamburg")) {
    GTEST_SKIP() << "hamburg not found";
  }
  run_core("test/hamburg.osm.pbf", "test/hamburg", 300U, 3 * 3600U);
}

TEST(cch, switzerland_core) {
  if (!fs::exists("test/switzerland.osm.pbf") &&
      !fs::exists("test/switzerland")) {
    GTEST_SKIP() << "switzerland not found";
  }
  run_core("test/switzerland.osm.pbf", "test/switzerland", 100U, 5 * 3600U);
}

TEST(cch, DISABLED_germany_core) {
  if (!fs::exists("test/germany.osm.pbf") && !fs::exists("test/germany")) {
    GTEST_SKIP() << "germany not found";
  }
  run_core("test/germany.osm.pbf", "test/germany", 20U, 12 * 3600U);
}
