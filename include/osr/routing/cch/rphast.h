#pragma once

#include <algorithm>
#include <cstdint>
#include <limits>
#include <vector>

#include "osr/routing/cch/cch.h"
#include "osr/routing/cch/customize.h"
#include "osr/routing/cch/query.h"
#include "osr/routing/cch/turns.h"
#include "osr/types.h"
#include "osr/ways.h"

namespace osr {

// RPHAST (restricted PHAST) for one-to-many queries.
//
//   Delling, Goldberg, Werneck: "Faster Batched Shortest Paths in Road
//   Networks", ATMOS 2011.
//
// Three phases. Preprocessing is the CCH we already have. Target selection
// runs once the target set is known and extracts T', the upward closure of the
// targets over the reverse of the downward graph -- every node that can lie on
// a downward path into a target. The query then runs an ordinary upward search
// from the source and sweeps T' once in decreasing rank order, which is where
// the win over a cost bounded Dijkstra comes from: the sweep is proportional to
// |T'| rather than to the ball around the source.
//
// Turn awareness changes the algorithm as stated in the paper. A label here is
// not one distance per node but one per (node, port), because a turn cost is a
// function of the port arrived on and the port left by. The sweep therefore
// relaxes per port, and after pulling from the higher ranked neighbours it has
// to close the self loops at the node -- a loop leaves the node and comes back
// to it on a different port, which is how a forbidden turn gets driven around.
//
// Selection walks the hierarchy topology only and never reads the metric, so
// one T' stays valid across customizations and across metrics.
template <Profile P>
struct rphast {
  static constexpr auto const kNone = std::numeric_limits<std::uint32_t>::max();

  void clear() {
    search_.clear();
    targets_.clear();
  }

  void add_start(cch const& c,
                 node_idx_t const n,
                 port_t const p,
                 cost_t const cost) {
    search_.add_start(c, n, p, cost);
  }

  // Registers a node as a target for the selection phase. Which query target it
  // belongs to and what the last mile costs stay with the caller: selection
  // only needs the set of nodes.
  void add_target(cch const& c, node_idx_t const n) {
    if (c.contains(n)) {
      targets_.push_back(c.rank_[n]);
    }
  }

  std::size_t n_selected() const { return order_.size(); }

  bool no_starts() const { return search_.no_starts(); }

  // ---- phase 2: target selection -----------------------------------------
  void select(ways::routing const& r, cch const& c) {
    // Only the entries touched by the previous selection are stale, so reset
    // those instead of the whole rank indexed array.
    if (pos_.size() != c.n_ranks()) {
      pos_.assign(c.n_ranks(), kNone);
    } else {
      for (auto const rk : order_) {
        pos_[to_idx(rk)] = kNone;
      }
    }

    order_.clear();
    for (auto const rk : targets_) {
      if (pos_[to_idx(rk)] == kNone) {
        pos_[to_idx(rk)] = 0U;  // provisional, real index assigned below
        order_.push_back(rk);
      }
    }

    // Upward closure. A downward arc into `u` comes from a higher ranked node,
    // and those are exactly the slots in which `u` is the lower node. `order_`
    // doubles as the queue; the closure only ever adds higher ranks, so it
    // terminates.
    for (auto i = std::size_t{0U}; i != order_.size(); ++i) {
      auto const u = order_[i];
      for (auto s = c.upper_begin(u); s != c.upper_end(u); ++s) {
        if (c.dn_entries(s).empty()) {
          continue;
        }
        auto const h = c.adj_head_[s];
        if (pos_[to_idx(h)] == kNone) {
          pos_[to_idx(h)] = 0U;
          order_.push_back(h);
        }
      }
    }

    // The sweep needs decreasing rank so that every node is final by the time
    // a lower one pulls from it.
    std::sort(begin(order_), end(order_), std::greater<>{});

    ofs_.clear();
    ofs_.reserve(order_.size() + 1U);
    ofs_.push_back(0U);
    for (auto i = std::size_t{0U}; i != order_.size(); ++i) {
      pos_[to_idx(order_[i])] = static_cast<std::uint32_t>(i);
      auto const ports =
          std::min(n_ports(r, c.order_[order_[i]]), kMaxPorts);
      ofs_.push_back(ofs_.back() + ports);
    }
    cost_.assign(ofs_.empty() ? 0U : ofs_.back(), kInfeasible);
  }

  // ---- phase 3: query ----------------------------------------------------
  bool run(typename P::parameters const& params,
           ways const& w,
           cch const& c,
           cch_metric const& m,
           cost_t const max) {
    auto const& r = *w.r_;
    auto const ok = search_.run_forward(params, w, c, m, max);

    std::fill(begin(cost_), end(cost_), kInfeasible);

    // Seed T' from the upward search space. The meeting node of any source to
    // target path is the highest ranked node on it, so it lies in both the
    // upward search space and -- being on the downward half -- in T'.
    for (auto const& [rank, e] : search_.forward()) {
      if (rank >= pos_.size()) {
        continue;
      }
      auto const idx = pos_[rank];
      if (idx == kNone) {
        continue;
      }
      auto const base = ofs_[idx];
      auto const ports = static_cast<port_t>(ofs_[idx + 1U] - base);
      for (auto p = port_t{0U}; p != ports; ++p) {
        if (e.cost_[p] < cost_[base + p]) {
          cost_[base + p] = e.cost_[p];
        }
      }
    }

    auto const turn = [&](node_idx_t const n, port_t const in,
                          port_t const out) {
      return cch_turn_cost<P>(params, r, w.timezones_, n, in, out);
    };

    for (auto i = std::size_t{0U}; i != order_.size(); ++i) {
      auto const rk = order_[i];
      auto const n = c.order_[rk];
      auto const base = ofs_[i];
      auto const ports = static_cast<port_t>(ofs_[i + 1U] - base);

      // Pull along the downward arcs that end here. Their tails rank higher
      // and the sweep is descending, so every one of them is already final.
      for (auto s = c.upper_begin(rk); s != c.upper_end(rk); ++s) {
        auto const entries = c.dn_entries(s);
        if (entries.empty()) {
          continue;
        }
        auto const h = c.adj_head_[s];
        auto const hidx = pos_[to_idx(h)];
        if (hidx == kNone) {
          continue;  // cannot happen: the closure put every such tail in T'
        }
        auto const hbase = ofs_[hidx];
        auto const hports = static_cast<port_t>(ofs_[hidx + 1U] - hbase);
        auto const hn = c.order_[h];
        auto const ofs = c.dn_ofs_[s];
        for (auto k = std::size_t{0U}; k != entries.size(); ++k) {
          auto const e = entries[k];
          if (e.exit_ >= ports) {
            continue;
          }
          auto const weight = m.dn(ofs + k);
          if (weight == kInfeasible) {
            continue;
          }
          for (auto p = port_t{0U}; p != hports; ++p) {
            auto const ch = cost_[hbase + p];
            if (ch == kInfeasible) {
              continue;
            }
            auto const tc = turn(hn, p, e.entry_);
            if (tc == kInfeasible) {
              continue;
            }
            auto const nc =
                clamp_cost(static_cast<std::uint64_t>(ch) + tc + weight);
            if (nc < max && nc < cost_[base + e.exit_]) {
              cost_[base + e.exit_] = nc;
            }
          }
        }
      }

      close_loops(c, m, rk, n, base, ports, max, turn);
    }

    return ok;
  }

  // Cost from the source to `n` arriving on port `p`, or kInfeasible.
  cost_t get(cch const& c, node_idx_t const n, port_t const p) const {
    if (!c.contains(n)) {
      return kInfeasible;
    }
    auto const rk = to_idx(c.rank_[n]);
    if (rk >= pos_.size()) {
      return kInfeasible;
    }
    auto const idx = pos_[rk];
    if (idx == kNone) {
      return kInfeasible;
    }
    auto const base = ofs_[idx];
    return base + p < ofs_[idx + 1U] ? cost_[base + p] : kInfeasible;
  }

private:
  // A self loop leaves the node on one port and returns on another, so it can
  // improve a port after the downward arcs have been pulled. Settling the ports
  // in increasing cost makes chains of loops converge in one pass, the same way
  // `close_target_loops` does it on the target side of a point to point query.
  template <typename TurnFn>
  void close_loops(cch const& c,
                   cch_metric const& m,
                   cch_rank_t const rk,
                   node_idx_t const n,
                   std::uint32_t const base,
                   port_t const ports,
                   cost_t const max,
                   TurnFn const& turn) {
    auto const loops = c.loop_entries(rk);
    if (loops.empty()) {
      return;
    }

    auto done = std::uint32_t{0U};
    for (auto step = port_t{0U}; step != ports; ++step) {
      auto best = kInfeasible;
      auto at = kMaxPorts;
      for (auto p = port_t{0U}; p != ports; ++p) {
        if ((done & (std::uint32_t{1U} << p)) == 0U && cost_[base + p] < best) {
          best = cost_[base + p];
          at = p;
        }
      }
      if (at == kMaxPorts) {
        break;
      }
      done |= std::uint32_t{1U} << at;

      for (auto k = std::size_t{0U}; k != loops.size(); ++k) {
        auto const lk = loops[k];
        if (lk.exit_ >= ports) {
          continue;
        }
        auto const weight =
            m.loop(static_cast<cch_entry_idx_t>(c.loop_begin(rk) + k));
        if (weight == kInfeasible) {
          continue;
        }
        auto const tc = turn(n, at, lk.entry_);
        if (tc == kInfeasible) {
          continue;
        }
        auto const nc =
            clamp_cost(static_cast<std::uint64_t>(best) + tc + weight);
        if (nc < max && nc < cost_[base + lk.exit_]) {
          cost_[base + lk.exit_] = nc;
        }
      }
    }
  }

  cch_search<P> search_;
  std::vector<cch_rank_t> targets_;

  // T' in decreasing rank order, and rank -> index into it.
  std::vector<cch_rank_t> order_;
  std::vector<std::uint32_t> pos_;

  // Ragged per port costs over `order_`: node i owns [ofs_[i], ofs_[i+1]).
  std::vector<std::uint32_t> ofs_;
  std::vector<cost_t> cost_;
};

}  // namespace osr
