#pragma once

#include <algorithm>
#include <cstdint>
#include <limits>
#include <span>
#include <vector>

#include "osr/routing/cch/cch.h"
#include "osr/routing/cch/customize.h"
#include "osr/routing/cch/query.h"
#include "osr/routing/cch/turns.h"
#include "utl/verify.h"

#include "osr/types.h"
#include "osr/ways.h"

namespace osr {

// RPHAST (restricted PHAST) for one-to-many queries.
//
//   Delling, Goldberg, Werneck: "Faster Batched Shortest Paths in Road
//   Networks", ATMOS 2011.
//
// Three phases. Preprocessing is the CCH we already have. Target selection runs
// once the target set is known and extracts T', the upward closure of the
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
// The work is split by what each piece depends on, so the expensive parts can
// be hoisted out of the per source query:
//
//   select()       depends only on the target set and the hierarchy topology.
//                  Builds T' and packs the restricted downward graph: arcs are
//                  renumbered into T' local indices and laid out contiguously
//                  in sweep order, so the sweep walks memory forwards instead
//                  of chasing the global CSR. Survives customization.
//   materialize()  depends additionally on the metric and the profile
//                  parameters. Copies the arc weights in and precomputes a turn
//                  table per node. Profiling put ~40% of the sweep in
//                  `cch_turn_cost`, which is a pure function of (node, in port,
//                  out port) and so is worth tabulating: a node with P ports and
//                  A incoming arcs costs P*P table entries instead of A*P calls,
//                  and A exceeds P comfortably.
//   run()          depends on the source. Upward search plus the sweep.
template <Profile P>
struct rphast {
  static constexpr auto const kNone = std::numeric_limits<std::uint32_t>::max();

  // One incoming downward arc, tail renumbered into the selected set.
  struct dn_arc {
    std::uint32_t src_;  // index into `order_` of the tail, always < head index
    port_t entry_;  // port left at the tail
    port_t exit_;  // port arrived at the head
    cost_t weight_;
  };

  struct loop_arc {
    port_t entry_;
    port_t exit_;
    cost_t weight_;
  };

  void clear() {
    search_.clear();
    targets_.clear();
  }

  // Drops the sources but keeps the selection, so the next source can reuse it.
  void clear_starts() { search_.clear_starts(); }

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

  bool no_starts() const { return search_.no_starts(); }
  std::size_t n_selected() const { return order_.size(); }
  std::size_t n_arcs() const { return arcs_.size(); }
  bool materialized() const { return materialized_; }
  bool packed() const { return packed_; }

  // ---- phase 2a: target selection, metric independent ---------------------
  //
  // `pack` builds the restricted downward graph as a contiguous array in sweep
  // order. It is off by default and deliberately so: packing walks exactly the
  // slots the sweep walks, so a one-shot query traverses them twice and comes
  // out ~1.5x slower. It pays only when one selection serves many sources --
  // measured at ~2x on the query once selection is amortized -- which needs the
  // target set to be a property of the region rather than of the request.
  // `backward` mirrors the whole construction into the reversed graph, which is
  // what a last-mile query needs: it asks for the cost from each target *to*
  // the source, not from it. Reversing means the upward search runs over `dn_`
  // arcs and the sweep over `up_`, with the turn at a node evaluated in the
  // opposite order. Flipping only the terminal ports is not sufficient and
  // produces answers wrong in both directions.
  void select(ways::routing const& r,
              cch const& c,
              bool const pack = false,
              bool const backward = false) {
    backward_ = backward;
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
        if ((backward_ ? c.up_entries(s) : c.dn_entries(s)).empty()) {
          continue;
        }
        auto const h = c.adj_head_[s];
        if (pos_[to_idx(h)] == kNone) {
          pos_[to_idx(h)] = 0U;
          order_.push_back(h);
        }
      }
    }

    // Decreasing rank, so that every node is final by the time a lower one
    // pulls from it -- and so that a tail always sits earlier in the arrays
    // than its head, which is what makes the sweep a forward walk.
    std::sort(begin(order_), end(order_), std::greater<>{});

    auto const n = order_.size();
    node_.resize(n);
    ports_.resize(n);
    ofs_.clear();
    ofs_.reserve(n + 1U);
    ofs_.push_back(0U);
    turn_ofs_.clear();
    turn_ofs_.reserve(n + 1U);
    turn_ofs_.push_back(0U);

    for (auto i = std::size_t{0U}; i != n; ++i) {
      pos_[to_idx(order_[i])] = static_cast<std::uint32_t>(i);
      auto const nd = c.order_[order_[i]];
      node_[i] = nd;
      auto const p = std::min(n_ports(r, nd), kMaxPorts);
      ports_[i] = p;
      ofs_.push_back(ofs_.back() + p);
      turn_ofs_.push_back(turn_ofs_.back() +
                          static_cast<std::uint32_t>(p) * p);
    }
    cost_.assign(ofs_.back(), kInfeasible);

    packed_ = pack;
    materialized_ = false;
    if (!pack) {
      arcs_.clear();
      arc_widx_.clear();
      lps_.clear();
      lp_widx_.clear();
      turn_.clear();
      return;
    }
    turn_.assign(turn_ofs_.back(), kInfeasible);

    // Pack the restricted downward graph, grouped by head in sweep order.
    arc_ofs_.clear();
    arc_ofs_.reserve(n + 1U);
    arc_ofs_.push_back(0U);
    arcs_.clear();
    arc_widx_.clear();
    lp_ofs_.clear();
    lp_ofs_.reserve(n + 1U);
    lp_ofs_.push_back(0U);
    lps_.clear();
    lp_widx_.clear();

    for (auto i = std::size_t{0U}; i != n; ++i) {
      auto const rk = order_[i];
      auto const head_ports = ports_[i];

      for (auto s = c.upper_begin(rk); s != c.upper_end(rk); ++s) {
        auto const entries = backward_ ? c.up_entries(s) : c.dn_entries(s);
        if (entries.empty()) {
          continue;
        }
        auto const src = pos_[to_idx(c.adj_head_[s])];
        if (src == kNone) {
          continue;  // cannot happen: the closure put every tail in T'
        }
        auto const tail_ports = ports_[src];
        auto const ofs = backward_ ? c.up_ofs_[s] : c.dn_ofs_[s];
        for (auto k = std::size_t{0U}; k != entries.size(); ++k) {
          auto const e = entries[k];
          // Ports are stored by their role in the sweep -- at the tail (the
          // higher ranked node) and at the head -- which swaps between the two
          // directions because a `dn_` arc runs high to low and an `up_` arc
          // low to high.
          auto const tail_port = backward_ ? e.exit_ : e.entry_;
          auto const head_port = backward_ ? e.entry_ : e.exit_;
          // A node with more than `kMaxPorts` ways has ports the label arrays
          // cannot address; those arcs are dropped rather than aliased.
          if (head_port >= head_ports || tail_port >= tail_ports) {
            continue;
          }
          arcs_.push_back(dn_arc{src, tail_port, head_port, kInfeasible});
          arc_widx_.push_back(static_cast<cch_entry_idx_t>(ofs + k));
        }
      }
      arc_ofs_.push_back(static_cast<std::uint32_t>(arcs_.size()));

      auto const loops = c.loop_entries(rk);
      for (auto k = std::size_t{0U}; k != loops.size(); ++k) {
        auto const lk = loops[k];
        if (lk.exit_ >= head_ports || lk.entry_ >= head_ports) {
          continue;
        }
        lps_.push_back(loop_arc{backward_ ? lk.exit_ : lk.entry_,
                                backward_ ? lk.entry_ : lk.exit_, kInfeasible});
        lp_widx_.push_back(static_cast<cch_entry_idx_t>(c.loop_begin(rk) + k));
      }
      lp_ofs_.push_back(static_cast<std::uint32_t>(lps_.size()));
    }
  }

  // ---- phase 2b: bind the metric and the profile parameters ---------------
  // Binds the metric and the profile parameters into the packed arrays. Only
  // meaningful after `select(..., true)`; a caller reusing one selection across
  // customizations calls this again instead of re-selecting.
  void materialize(typename P::parameters const& params,
                   ways const& w,
                   cch_metric const& m) {
    utl::verify(packed_, "rphast: materialize requires a packed selection");
    auto const& r = *w.r_;

    for (auto a = std::size_t{0U}; a != arcs_.size(); ++a) {
      arcs_[a].weight_ = backward_ ? m.up(arc_widx_[a])
                                   : m.dn(arc_widx_[a]);
    }
    for (auto a = std::size_t{0U}; a != lps_.size(); ++a) {
      lps_[a].weight_ = m.loop(lp_widx_[a]);
    }

    // `cch_turn_cost` was ~40% of the sweep, and is a pure function of
    // (node, in port, out port), so it is tabulated once here: P*P entries per
    // node instead of a call per arc and port.
    for (auto i = std::size_t{0U}; i != order_.size(); ++i) {
      auto const nd = node_[i];
      auto const p = ports_[i];
      auto const base = turn_ofs_[i];
      for (auto in = port_t{0U}; in != p; ++in) {
        for (auto out = port_t{0U}; out != p; ++out) {
          turn_[base + static_cast<std::uint32_t>(in) * p + out] =
              cch_turn_cost<P>(params, r, w.timezones_, nd, in, out);
        }
      }
    }
    materialized_ = true;
  }

  // ---- phase 3: query ----------------------------------------------------
  bool run(typename P::parameters const& params,
           ways const& w,
           cch const& c,
           cch_metric const& m,
           cost_t const max) {
    if (packed_ && !materialized_) {
      materialize(params, w, m);
    }

    auto const ok = backward_ ? search_.run_backward(params, w, c, m, max)
                              : search_.run_forward(params, w, c, m, max);

    std::fill(begin(cost_), end(cost_), kInfeasible);

    // Seed T' from the upward search space. The meeting node of any source to
    // target path is the highest ranked node on it, so it lies in both the
    // upward search space and -- being on the downward half -- in T'.
    auto const& seed_from = backward_ ? search_.backward() : search_.forward();
    for (auto const& [rank, e] : seed_from) {
      if (rank >= pos_.size()) {
        continue;
      }
      auto const idx = pos_[rank];
      if (idx == kNone) {
        continue;
      }
      auto const base = ofs_[idx];
      auto const p = ports_[idx];
      for (auto q = port_t{0U}; q != p; ++q) {
        if (e.cost_[q] < cost_[base + q]) {
          cost_[base + q] = e.cost_[q];
        }
      }
    }

    if (packed_) {
      sweep_packed(max);
    } else {
      sweep_hierarchy(params, w, c, m, max);
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
    return idx != kNone && p < ports_[idx] ? cost_[ofs_[idx] + p] : kInfeasible;
  }

private:
  // Default sweep: reads the hierarchy's own CSR directly and evaluates turn
  // costs as it goes. One pass over the selected set, which is what makes it
  // the right choice when the selection is not reused.
  void sweep_hierarchy(typename P::parameters const& params,
                       ways const& w,
                       cch const& c,
                       cch_metric const& m,
                       cost_t const max) {
    auto const& r = *w.r_;
    auto const turn = [&](node_idx_t const n, port_t const in,
                          port_t const out) {
      return cch_turn_cost<P>(params, r, w.timezones_, n, in, out);
    };

    for (auto i = std::size_t{0U}; i != order_.size(); ++i) {
      auto const rk = order_[i];
      auto const base = ofs_[i];
      auto const ports = ports_[i];

      // Pull along the downward arcs that end here. Their tails rank higher
      // and the sweep is descending, so every one of them is already final.
      for (auto s = c.upper_begin(rk); s != c.upper_end(rk); ++s) {
        auto const entries = backward_ ? c.up_entries(s) : c.dn_entries(s);
        if (entries.empty()) {
          continue;
        }
        auto const hidx = pos_[to_idx(c.adj_head_[s])];
        if (hidx == kNone) {
          continue;  // cannot happen: the closure put every tail in T'
        }
        auto const hbase = ofs_[hidx];
        auto const hports = ports_[hidx];
        auto const hn = node_[hidx];
        auto const ofs = backward_ ? c.up_ofs_[s] : c.dn_ofs_[s];
        for (auto k = std::size_t{0U}; k != entries.size(); ++k) {
          auto const e = entries[k];
          auto const tail_port = backward_ ? e.exit_ : e.entry_;
          auto const head_port = backward_ ? e.entry_ : e.exit_;
          if (head_port >= ports || tail_port >= hports) {
            continue;
          }
          auto const weight = backward_ ? m.up(ofs + k) : m.dn(ofs + k);
          if (weight == kInfeasible) {
            continue;
          }
          auto& dst = cost_[base + head_port];
          for (auto q = port_t{0U}; q != hports; ++q) {
            auto const ch = cost_[hbase + q];
            if (ch == kInfeasible) {
              continue;
            }
            // Reversed, the turn at the tail is taken arriving on the arc and
            // leaving by the port the state already holds.
            auto const tc = backward_ ? turn(hn, tail_port, q)
                                      : turn(hn, q, tail_port);
            if (tc == kInfeasible) {
              continue;
            }
            auto const nc =
                clamp_cost(static_cast<std::uint64_t>(ch) + tc + weight);
            if (nc < max && nc < dst) {
              dst = nc;
            }
          }
        }
      }

      close_loops(base, ports, max, mirrored_loops(c, rk),
                  [&](std::size_t const k) {
                    return m.loop(
                        static_cast<cch_entry_idx_t>(c.loop_begin(rk) + k));
                  },
                  [&](port_t const in, port_t const out) {
                    return turn(node_[i], in, out);
                  });
    }
  }

  // Opt-in sweep over the packed restricted graph: arcs contiguous in sweep
  // order, weights inlined, turn costs read from the precomputed table.
  void sweep_packed(cost_t const max) {
    for (auto i = std::size_t{0U}; i != order_.size(); ++i) {
      auto const base = ofs_[i];
      auto const arc_end = arc_ofs_[i + 1U];
      for (auto a = arc_ofs_[i]; a != arc_end; ++a) {
        auto const& arc = arcs_[a];
        if (arc.weight_ == kInfeasible) {
          continue;
        }
        auto const sbase = ofs_[arc.src_];
        auto const sp = ports_[arc.src_];
        auto const tbase_src = turn_ofs_[arc.src_];
        auto& dst = cost_[base + arc.exit_];
        for (auto q = port_t{0U}; q != sp; ++q) {
          auto const ch = cost_[sbase + q];
          if (ch == kInfeasible) {
            continue;
          }
          auto const tc =
              backward_
                  ? turn_[tbase_src +
                          static_cast<std::uint32_t>(arc.entry_) * sp + q]
                  : turn_[tbase_src + static_cast<std::uint32_t>(q) * sp +
                          arc.entry_];
          if (tc == kInfeasible) {
            continue;
          }
          auto const nc =
              clamp_cost(static_cast<std::uint64_t>(ch) + tc + arc.weight_);
          if (nc < max && nc < dst) {
            dst = nc;
          }
        }
      }

      auto const tbase = turn_ofs_[i];
      auto const ports = ports_[i];
      close_loops(base, ports, max,
                  std::span{lps_.data() + lp_ofs_[i],
                            lp_ofs_[i + 1U] - lp_ofs_[i]},
                  [&](std::size_t const k) {
                    return lps_[lp_ofs_[i] + k].weight_;
                  },
                  [&](port_t const in, port_t const out) {
                    return turn_[tbase + static_cast<std::uint32_t>(in) *
                                             ports + out];
                  });
    }
  }

  // A self loop leaves the node on one port and returns on another, so it can
  // improve a port after the downward arcs have been pulled. Settling the ports
  // in increasing cost makes chains of loops converge in one pass, the same way
  // `close_target_loops` does it on the target side of a point to point query.
  // The entries, their weights and the turn costs come from the caller so that
  // both sweeps share the logic.
  template <typename Entries, typename WeightFn, typename TurnFn>
  void close_loops(std::uint32_t const base,
                   port_t const ports,
                   cost_t const max,
                   Entries const& loops,
                   WeightFn const& weight_of,
                   TurnFn const& turn) {
    if (loops.empty()) {
      return;
    }

    auto done = std::uint32_t{0U};
    for (auto step = port_t{0U}; step != ports; ++step) {
      auto best = kInfeasible;
      auto at = kMaxPorts;
      for (auto q = port_t{0U}; q != ports; ++q) {
        if ((done & (std::uint32_t{1U} << q)) == 0U && cost_[base + q] < best) {
          best = cost_[base + q];
          at = q;
        }
      }
      if (at == kMaxPorts) {
        break;
      }
      done |= std::uint32_t{1U} << at;

      for (auto k = std::size_t{0U}; k != loops.size(); ++k) {
        auto const entry = loops[k].entry_;
        auto const exit = loops[k].exit_;
        if (entry >= ports || exit >= ports) {
          continue;
        }
        auto const wgt = weight_of(k);
        if (wgt == kInfeasible) {
          continue;
        }
        auto const tc = backward_ ? turn(entry, at) : turn(at, entry);
        if (tc == kInfeasible) {
          continue;
        }
        auto const nc =
            clamp_cost(static_cast<std::uint64_t>(best) + tc + wgt);
        if (nc < max && nc < cost_[base + exit]) {
          cost_[base + exit] = nc;
        }
      }
    }
  }

  cch_search<P> search_;
  std::vector<cch_rank_t> targets_;

  // T' in decreasing rank order, plus rank -> index and the per node data the
  // sweep needs, all indexed by position so the walk stays local.
  std::vector<cch_rank_t> order_;
  std::vector<std::uint32_t> pos_;
  std::vector<node_idx_t> node_;
  std::vector<port_t> ports_;

  // Ragged per port costs over `order_`: node i owns [ofs_[i], ofs_[i+1]).
  std::vector<std::uint32_t> ofs_;
  std::vector<cost_t> cost_;

  // Restricted downward graph and self loops, both CSR over `order_`.
  std::vector<std::uint32_t> arc_ofs_;
  std::vector<dn_arc> arcs_;
  std::vector<cch_entry_idx_t> arc_widx_;
  std::vector<std::uint32_t> lp_ofs_;
  std::vector<loop_arc> lps_;
  std::vector<cch_entry_idx_t> lp_widx_;

  // Turn cost table: node i owns a ports_[i] x ports_[i] row major block at
  // turn_ofs_[i], indexed [in * ports + out].
  std::vector<std::uint32_t> turn_ofs_;
  std::vector<cost_t> turn_;

  // Loop entries with the ports put in sweep role order, matching how the
  // packed path stores them.
  std::span<loop_arc const> mirrored_loops(cch const& c, cch_rank_t const rk) {
    auto const loops = c.loop_entries(rk);
    scratch_loops_.clear();
    scratch_loops_.reserve(loops.size());
    for (auto const& lk : loops) {
      scratch_loops_.push_back(loop_arc{backward_ ? lk.exit_ : lk.entry_,
                                        backward_ ? lk.entry_ : lk.exit_,
                                        kInfeasible});
    }
    return {scratch_loops_.data(), scratch_loops_.size()};
  }

  std::vector<loop_arc> scratch_loops_;
  bool backward_{false};
  bool packed_{false};
  bool materialized_{false};
};

}  // namespace osr
