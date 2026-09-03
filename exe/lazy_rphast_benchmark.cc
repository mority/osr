// Eager RPHAST vs Lazy RPHAST on the same (source, target set) pairs.
//
// Both answer the same question -- cost from one source to k targets over the
// customized hierarchy -- but they spend the work differently:
//
//   eager  selects T', the upward closure of the whole target set, then sweeps
//          all of it. Cost is O(|T'|) whether one target is read or all of
//          them, and T' always reaches the dense top of the hierarchy.
//   lazy   has no selection and no sweep. Each target walks up the elimination
//          tree until it meets a memoized node and computes the distances back
//          down. Targets share the upper part of the tree, so the marginal
//          cost of a target falls as the memo fills.
//
// The interesting variable is therefore the target count, and the interesting
// output is not just the total but where the crossover sits.
//
// Distances are compared entry by entry: the two must agree exactly, otherwise
// a speedup means nothing.
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <cmath>
#include <numeric>
#include <optional>
#include <random>
#include <string>
#include <vector>

#include "fmt/core.h"
#include "fmt/std.h"

#include "conf/options_parser.h"

#include "geo/box.h"

#include "osr/lookup.h"
#include "osr/routing/cch/cch.h"
#include "osr/routing/cch/customize.h"
#include "osr/routing/cch/lazy_rphast.h"
#include "osr/routing/cch/rphast.h"
#include "osr/routing/profiles/car.h"
#include "osr/ways.h"

namespace fs = std::filesystem;
using namespace osr;

namespace {

struct settings : public conf::configuration {
  settings() : configuration{"Options"} {
    param(data_, "data,d", "osr data directory");
    param(n_sources_, "sources,n", "number of random sources per target count");
    param(targets_, "targets,t", "comma separated target counts");
    param(max_cost_, "max", "cost bound in seconds");
    param(radius_km_, "radius", "target sampling radius around the source, km");
    param(read_frac_, "read",
          "fraction of the targets whose distance is actually read (the "
          "incremental setting lazy RPHAST is built for)");
    param(backward_, "backward",
          "run the reversed (last mile) direction instead of the forward one");
    param(seed_, "seed", "seed for the random query set");
  }

  fs::path data_{"osr"};
  unsigned n_sources_{50U};
  std::string targets_{"10,50,100,200,500,1000,2000"};
  unsigned max_cost_{3600U};
  double radius_km_{40.0};
  double read_frac_{1.0};
  bool backward_{false};
  unsigned seed_{20260902U};
};

std::vector<unsigned> parse_list(std::string const& s) {
  auto out = std::vector<unsigned>{};
  auto pos = std::size_t{0U};
  while (pos <= s.size()) {
    auto const comma = s.find(',', pos);
    auto const tok = s.substr(pos, comma == std::string::npos ? std::string::npos
                                                             : comma - pos);
    if (!tok.empty()) {
      out.push_back(static_cast<unsigned>(std::stoul(tok)));
    }
    if (comma == std::string::npos) {
      break;
    }
    pos = comma + 1U;
  }
  return out;
}

double median(std::vector<double> v) {
  if (v.empty()) {
    return 0.0;
  }
  std::sort(begin(v), end(v));
  return v[v.size() / 2U];
}

double pct(std::vector<double> v, double const p) {
  if (v.empty()) {
    return 0.0;
  }
  std::sort(begin(v), end(v));
  return v[static_cast<std::size_t>(p * static_cast<double>(v.size() - 1U))];
}

using clk = std::chrono::steady_clock;

double ms_since(clk::time_point const t) {
  return std::chrono::duration<double, std::milli>(clk::now() - t).count();
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
    fmt::println("data directory {} not found", opt.data_);
    return 1;
  }
  if (!cch::exists(opt.data_)) {
    fmt::println("no cch in {}", opt.data_);
    return 1;
  }

  auto const w = ways{opt.data_, cista::mmap::protection::READ};
  // The rtree, not rejection sampling: at continental scale the chance of a
  // uniformly drawn node landing in a 40 km box is far too small to sample
  // target sets that way.
  auto const l = lookup{w, opt.data_, cista::mmap::protection::READ};
  auto const c = cch::read(opt.data_);
  auto const params = car::parameters{};

  // Prefer the stored metric; customizing Germany takes seconds and is not
  // what this benchmark is measuring.
  auto owned = cch_metric{};
  auto mapped = std::optional<cista::wrapped<cch_metric>>{};
  auto const* m = static_cast<cch_metric const*>(nullptr);
  if (cch_metric::exists(opt.data_)) {
    mapped = cch_metric::read(opt.data_);
    if ((*mapped)->matches(params_blob(params))) {
      m = mapped->get();
    }
  }
  if (m == nullptr) {
    fmt::println("customizing (no matching stored metric)...");
    customize<car>(params, w, *c, owned);
    m = &owned;
  }

  auto const max_cost = static_cast<cost_t>(opt.max_cost_);
  auto const counts = parse_list(opt.targets_);

  auto eager = rphast<car>{};
  auto lazy = lazy_rphast<car>{};

  auto prng = std::mt19937{opt.seed_};
  auto distr = std::uniform_int_distribution<std::uint32_t>{
      0U, static_cast<std::uint32_t>(w.n_nodes() - 1U)};

  // Only nodes that are in the hierarchy can be sources or targets.
  auto const sample_hierarchy_node = [&](node_idx_t& out) {
    for (auto tries = 0U; tries != 1000000U; ++tries) {
      auto const n = node_idx_t{distr(prng)};
      if (c->contains(n)) {
        out = n;
        return true;
      }
    }
    return false;
  };

  auto const deg_lat = opt.radius_km_ / 111.2;
  auto candidates = std::vector<node_idx_t>{};
  auto seen = hash_set<node_idx_t>{};

  fmt::println(
      "\n{:>8} {:>10} {:>10} {:>9} {:>10} {:>10} {:>12} {:>12}", "targets",
      "eager ms", "lazy ms", "speedup", "eager p90", "lazy p90", "|T'| nodes",
      "memo nodes");
  fmt::println("{}", std::string(96, '-'));

  auto total_mismatch = std::uint64_t{0U};
  auto total_compared = std::uint64_t{0U};

  for (auto const k : counts) {
    auto eager_ms = std::vector<double>{};
    auto lazy_ms = std::vector<double>{};
    auto sel = std::vector<double>{};
    auto memo = std::vector<double>{};

    // Reset the stream per target count so every count sees the same sources.
    prng.seed(opt.seed_);

    for (auto q = 0U; q != opt.n_sources_; ++q) {
      auto src = node_idx_t{};
      if (!sample_hierarchy_node(src)) {
        continue;
      }
      auto const sp = w.get_node_pos(src);
      auto const deg_lon =
          deg_lat / std::max(0.2, std::cos(sp.lat() * 3.14159265 / 180.0));

      // Targets drawn from a box around the source, which is what a first/last
      // mile target set looks like: near the source, not spread over the map.
      auto box = geo::box{};
      box.extend(geo::latlng{sp.lat() - deg_lat, sp.lng() - deg_lon});
      box.extend(geo::latlng{sp.lat() + deg_lat, sp.lng() + deg_lon});
      candidates.clear();
      seen.clear();
      l.find(box, [&](way_idx_t const way) {
        for (auto const n : w.r_->way_nodes_[way]) {
          if (c->contains(n) && seen.emplace(n).second) {
            candidates.push_back(n);
          }
        }
      });
      if (candidates.size() < k) {
        continue;
      }
      std::shuffle(begin(candidates), end(candidates), prng);
      auto targets =
          std::vector<node_idx_t>(begin(candidates), begin(candidates) + k);

      // Both arms read the same prefix. Eager still has to select and sweep
      // the closure of the *whole* target set, because it cannot know in
      // advance which ones will be asked for -- that is precisely the
      // asymmetry lazy RPHAST exploits.
      auto const n_read = std::max(
          std::size_t{1U},
          static_cast<std::size_t>(opt.read_frac_ *
                                   static_cast<double>(targets.size())));

      // ---- eager ---------------------------------------------------------
      auto eager_out = std::vector<cost_t>(targets.size(), kInfeasible);
      auto const t_eager = clk::now();
      eager.clear();
      eager.add_start(*c, src, port_t{0U}, cost_t{0U});
      for (auto const t : targets) {
        eager.add_target(*c, t);
      }
      eager.select(*w.r_, *c, /* pack */ false, opt.backward_);
      eager.run(params, w, *c, *m, max_cost);
      for (auto i = std::size_t{0U}; i != n_read; ++i) {
        eager_out[i] = eager.get(*c, targets[i], port_t{0U});
      }
      eager_ms.push_back(ms_since(t_eager));
      sel.push_back(static_cast<double>(eager.n_selected()));

      // ---- lazy ----------------------------------------------------------
      auto lazy_out = std::vector<cost_t>(targets.size(), kInfeasible);
      auto const t_lazy = clk::now();
      lazy.clear();
      lazy.add_start(*c, src, port_t{0U}, cost_t{0U});
      lazy.run(params, w, *c, *m, max_cost, opt.backward_);
      for (auto i = std::size_t{0U}; i != n_read; ++i) {
        lazy_out[i] = lazy.get(params, w, *c, *m, targets[i], port_t{0U});
      }
      lazy_ms.push_back(ms_since(t_lazy));
      memo.push_back(static_cast<double>(lazy.n_memoized()));

      for (auto i = std::size_t{0U}; i != n_read; ++i) {
        ++total_compared;
        if (eager_out[i] != lazy_out[i]) {
          if (total_mismatch < 5U) {
            fmt::println(
                "MISMATCH src={} target={} eager={} lazy={}", src.v_,
                targets[i].v_, eager_out[i], lazy_out[i]);
          }
          ++total_mismatch;
        }
      }
    }

    auto const e = median(eager_ms);
    auto const l = median(lazy_ms);
    fmt::println("{:>8} {:>10.2f} {:>10.2f} {:>8.2f}x {:>10.2f} {:>10.2f} "
                 "{:>12.0f} {:>12.0f}",
                 k, e, l, l > 0.0 ? e / l : 0.0, pct(eager_ms, 0.9),
                 pct(lazy_ms, 0.9), median(sel), median(memo));
  }

  fmt::println("\ncompared {} distances, {} mismatches", total_compared,
               total_mismatch);
  return total_mismatch == 0U ? 0 : 1;
}
