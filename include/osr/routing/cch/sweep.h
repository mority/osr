#pragma once

#include <cstdint>

#include <optional>
#include <span>
#include <vector>

#include "osr/routing/cch/cch.h"
#include "osr/routing/cch/customize.h"
#include "osr/types.h"

namespace osr {

// The downward half of a one-to-many query, shared by the two RPHAST variants
// in this directory. Both evaluate the same recurrence
//
//   D[v] = min( D_up[v], min over downward arcs w->v of ( l[w->v] + D[w] ) )
//
// over per (node, port) labels; they only differ in the order they visit the
// nodes in and in where a node's label block lives. That difference is the
// `tail_of` callback, everything below is common.

// The per port cost block of one node: its offset in the cost array, how many
// ports of it are addressable, and the node itself, which the turn cost needs.
struct cch_sweep_node {
  std::uint32_t base_;
  port_t ports_;
  node_idx_t node_;
};

// Relaxes every downward arc that ends at `rk` into the cost block of its head.
// The tails rank higher, so both callers have made them final before this runs
// -- `tail_of` returning `std::nullopt` therefore cannot happen and is only
// guarded against.
//
// Reversed, the downward arcs of the search are the `up_` arcs of the
// hierarchy, and the turn at the tail is taken arriving on the arc and leaving
// by the port the state already holds, so both the entry list and the two port
// roles swap.
template <typename TailFn, typename TurnFn>
void cch_relax_dn_arcs(cch const& c,
                       cch_metric const& m,
                       cch_rank_t const rk,
                       bool const backward,
                       cost_t const max,
                       cch_sweep_node const& head,
                       std::vector<cost_t>& cost,
                       TailFn const& tail_of,
                       TurnFn const& turn) {
  for (auto s = c.upper_begin(rk); s != c.upper_end(rk); ++s) {
    auto const entries = backward ? c.up_entries(s) : c.dn_entries(s);
    if (entries.empty()) {
      continue;
    }
    auto const tail = tail_of(c.adj_head_[s]);
    if (!tail.has_value()) {
      continue;
    }
    auto const ofs = backward ? c.up_ofs_[s] : c.dn_ofs_[s];
    for (auto k = std::size_t{0U}; k != entries.size(); ++k) {
      auto const e = entries[k];
      auto const tail_port = backward ? e.exit_ : e.entry_;
      auto const head_port = backward ? e.entry_ : e.exit_;
      if (head_port >= head.ports_ || tail_port >= tail->ports_) {
        continue;
      }
      auto const weight = backward ? m.up(ofs + k) : m.dn(ofs + k);
      if (weight == kInfeasible) {
        continue;
      }
      auto& dst = cost[head.base_ + head_port];
      for (auto q = port_t{0U}; q != tail->ports_; ++q) {
        auto const ch = cost[tail->base_ + q];
        if (ch == kInfeasible) {
          continue;
        }
        auto const tc = backward ? turn(tail->node_, tail_port, q)
                                 : turn(tail->node_, q, tail_port);
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
}

// Closes the self loops at one node. A loop leaves the node on one port and
// comes back on another, which is how a turn that is forbidden here gets
// driven around, so it can still improve a port after the downward arcs have
// been pulled. Settling the ports cheapest first makes chains of loops
// converge in a single pass, the same way `close_target_loops` does it on the
// target side of a point to point query.
//
// `cost` is the node's own block. `loops` holds the entries in hierarchy
// orientation -- `entry_` is the port the loop leaves the node on -- which a
// reversed search swaps; `weight_of(k)` is the metric weight of `loops[k]`.
template <typename Loops, typename WeightFn, typename TurnFn>
void cch_close_loops(std::span<cost_t> const cost,
                     port_t const ports,
                     cost_t const max,
                     bool const backward,
                     Loops const& loops,
                     WeightFn const& weight_of,
                     TurnFn const& turn) {
  if (loops.empty()) {
    return;
  }
  for_ports_by_cost(
      ports, [&](port_t const p) { return cost[p]; },
      [&](port_t const at, cost_t const best) {
        for (auto k = std::size_t{0U}; k != loops.size(); ++k) {
          auto const entry = backward ? loops[k].exit_ : loops[k].entry_;
          auto const exit = backward ? loops[k].entry_ : loops[k].exit_;
          if (entry >= ports || exit >= ports) {
            continue;
          }
          auto const weight = weight_of(k);
          if (weight == kInfeasible) {
            continue;
          }
          auto const tc = backward ? turn(entry, at) : turn(at, entry);
          if (tc == kInfeasible) {
            continue;
          }
          auto const nc =
              clamp_cost(static_cast<std::uint64_t>(best) + tc + weight);
          if (nc < max && nc < cost[exit]) {
            cost[exit] = nc;
          }
        }
      });
}

}  // namespace osr
