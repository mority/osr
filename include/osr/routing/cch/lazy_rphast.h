#pragma once

#include <cinttypes>

#include <limits>
#include <span>
#include <vector>

#include "osr/routing/cch/cch.h"
#include "osr/routing/cch/customize.h"
#include "osr/routing/cch/query.h"
#include "osr/routing/cch/turns.h"
#include "osr/routing/profile.h"
#include "osr/types.h"
#include "osr/ways.h"

namespace osr {

// Lazy RPHAST for one-to-many queries.
//
//   Strasser, Zeitz: "A Simple Algorithm for Bounded-Hop Shortest Path
//   Queries" (Lazy RPHAST), and the CCH adaptation in
//   Bläsius, Buchhold, Wagner, Zeitz, Zündorf: "Customizable Contraction
//   Hierarchies - A Survey", Algorithm 4.1.
//
// Eager RPHAST (Delling/Goldberg/Werneck, ATMOS 2011) selects T', the upward
// closure of the whole target set, and sweeps every node of it once per source.
// Its cost is O(|T'|) no matter how many of the targets are actually read, and
// on a country sized hierarchy the closure always climbs into the dense top
// levels, so even a handful of targets pays for the whole ascent.
//
// Lazy RPHAST has no selection phase and no sweep. For a target u it walks up
// the elimination tree until it meets a node whose distance is already known,
// then computes the distances back down that path and memoizes them. The
// second and later targets stop as soon as they hit the part of the tree the
// earlier ones already paid for, so the shared upper region is computed once
// and the marginal cost of a target is the length of its own private tail.
//
// The eager variant used to live next to this one, behind an env var. It was
// removed after measuring both on Germany: reading *every* target the two do
// provably identical work -- the memo ends up holding exactly |T'| nodes -- and
// lazy still won on constant factors (0.99x to 1.05x). Reading a quarter of
// them it won by 1.3x to 1.6x, and the advantage grows as the read fraction
// falls. Lazy is bounded above by eager by construction: it can never touch
// more of the hierarchy than the closure eager commits to up front.
//
// The recurrence is the one the eager sweep applies, only evaluated on demand:
//
//   D[v] = min( D_up[v], min over upward arcs v->w of ( l_dn[w->v] + D[w] ) )
//
// where `D_up` is the tentative distance from the source's upward search.
//
// Turn awareness: a label is one cost per (node, port) rather than per node,
// the relaxation charges the turn taken at the tail node, and after pulling
// from the higher ranked neighbours the self loops at the node have to be
// closed so a forbidden turn can be driven around.
template <Profile P>
struct lazy_rphast {
  static constexpr auto const kNone = std::numeric_limits<std::uint32_t>::max();

  void clear() {
    search_.clear();
    reset_memo();
  }

  void add_start(cch const& c,
                 node_idx_t const n,
                 port_t const p,
                 cost_t const cost) {
    search_.add_start(c, n, p, cost);
  }

  bool no_starts() const { return search_.no_starts(); }

  // ---- phase 1: the source ------------------------------------------------
  //
  // Only the upward search. There is deliberately nothing else here: the
  // downward half is what `get` does lazily.
  bool run(typename P::parameters const& params,
           ways const& w,
           cch const& c,
           cch_metric const& m,
           cost_t const max,
           bool const backward = false) {
    backward_ = backward;
    max_ = max;
    auto const ok = search_.run_upward(params, w, c, m, max, backward);
    reset_memo();
    if (pos_.size() != c.n_ranks()) {
      pos_.assign(c.n_ranks(), kNone);
    }
    return ok;
  }

  // ---- phase 2: one target ------------------------------------------------
  //
  // Cost from the source to `n` arriving on port `p`, or kInfeasible.
  cost_t get(typename P::parameters const& params,
             ways const& w,
             cch const& c,
             cch_metric const& m,
             node_idx_t const n,
             port_t const p) {
    if (!c.contains(n)) {
      return kInfeasible;
    }
    auto const rk = c.rank_[n];
    if (to_idx(rk) >= pos_.size()) {
      return kInfeasible;
    }
    resolve(params, w, c, m, rk);
    auto const idx = pos_[to_idx(rk)];
    if (idx == kNone) {
      return kInfeasible;
    }
    return p < n_capped_ports(*w.r_, c.order_[rk]) ? cost_[idx + p]
                                                   : kInfeasible;
  }

private:
  void reset_memo() {
    for (auto const rk : touched_) {
      pos_[to_idx(rk)] = kNone;
    }
    touched_.clear();
    cost_.clear();
  }

  // Computes `D` for `rk` and for every ancestor of it that is not memoized
  // yet.
  //
  // The survey walks the elimination tree parent chain, which is correct
  // because every upward neighbour of a node is one of its ancestors. This
  // does the equivalent as an explicit post-order DFS over the upward arcs
  // instead: it visits the same set, but it does not depend on that structural
  // property holding, so a hierarchy whose upward neighbourhoods are not
  // cliques still gets the right answer rather than a silently too large one.
  void resolve(typename P::parameters const& params,
               ways const& w,
               cch const& c,
               cch_metric const& m,
               cch_rank_t const rk) {
    if (pos_[to_idx(rk)] != kNone) {
      return;
    }

    stack_.clear();
    stack_.push_back(rk);
    while (!stack_.empty()) {
      auto const v = stack_.back();
      if (pos_[to_idx(v)] != kNone) {
        stack_.pop_back();
        continue;
      }

      // Every higher ranked neighbour has to be final before `v` can be, which
      // is the same precondition the descending sweep gets for free.
      auto ready = true;
      for (auto s = c.upper_begin(v); s != c.upper_end(v); ++s) {
        auto const entries = backward_ ? c.up_entries(s) : c.dn_entries(s);
        if (entries.empty()) {
          continue;
        }
        auto const h = c.adj_head_[s];
        if (pos_[to_idx(h)] == kNone) {
          stack_.push_back(h);
          ready = false;
        }
      }
      if (!ready) {
        continue;
      }

      compute(params, w, c, m, v);
      stack_.pop_back();
    }
  }

  void compute(typename P::parameters const& params,
               ways const& w,
               cch const& c,
               cch_metric const& m,
               cch_rank_t const rk) {
    auto const& r = *w.r_;
    auto const node = c.order_[rk];
    auto const ports = n_capped_ports(r, node);
    auto const base = static_cast<std::uint32_t>(cost_.size());
    cost_.resize(cost_.size() + ports, kInfeasible);
    pos_[to_idx(rk)] = base;
    touched_.push_back(rk);

    // Seed with the tentative distance from the source's upward search. For a
    // node the source never reached this is simply absent.
    auto const& up = backward_ ? search_.backward() : search_.forward();
    if (auto const it = up.find(to_idx(rk)); it != end(up)) {
      for (auto q = port_t{0U}; q != ports; ++q) {
        auto const cu = it->second.cost_[q];
        if (cu < cost_[base + q]) {
          cost_[base + q] = cu;
        }
      }
    }

    auto const turn = cch_turn_fn<P>(params, w);

    // Pull from the downward arcs that end here:
    //
    //   D[v] = min( D_up[v], min over downward arcs w->v of ( l[w->v] + D[w] ) )
    //
    // The tails rank higher, so `resolve` made every one of them final before
    // this runs. Reversed, the downward arcs of the search are the `up_` arcs
    // of the hierarchy, and the turn at the tail is taken arriving on the arc
    // and leaving by the port the state already holds, so both the entry list
    // and the two port roles swap.
    for (auto s = c.upper_begin(rk); s != c.upper_end(rk); ++s) {
      auto const entries = backward_ ? c.up_entries(s) : c.dn_entries(s);
      if (entries.empty()) {
        continue;
      }
      auto const h = c.adj_head_[s];
      auto const tail_base = pos_[to_idx(h)];
      if (tail_base == kNone) {
        continue;  // cannot happen: `resolve` made every tail final first
      }
      auto const tail_node = c.order_[h];
      auto const tail_ports = n_capped_ports(r, tail_node);
      auto const ofs = backward_ ? c.up_ofs_[s] : c.dn_ofs_[s];
      for (auto k = std::size_t{0U}; k != entries.size(); ++k) {
        auto const e = entries[k];
        auto const tail_port = backward_ ? e.exit_ : e.entry_;
        auto const head_port = backward_ ? e.entry_ : e.exit_;
        if (head_port >= ports || tail_port >= tail_ports) {
          continue;
        }
        auto const weight = backward_ ? m.up(ofs + k) : m.dn(ofs + k);
        if (weight == kInfeasible) {
          continue;
        }
        auto& dst = cost_[base + head_port];
        for (auto q = port_t{0U}; q != tail_ports; ++q) {
          auto const ch = cost_[tail_base + q];
          if (ch == kInfeasible) {
            continue;
          }
          auto const tc = backward_ ? turn(tail_node, tail_port, q)
                                    : turn(tail_node, q, tail_port);
          if (tc == kInfeasible) {
            continue;
          }
          auto const nc =
              clamp_cost(static_cast<std::uint64_t>(ch) + tc + weight);
          if (nc < max_ && nc < dst) {
            dst = nc;
          }
        }
      }
    }

    // Close the self loops at the node. A loop leaves the node and comes back
    // to it on another port, which is how a turn that is forbidden here gets
    // driven around, so it can still improve a port after the downward arcs
    // have been pulled. Settling the ports cheapest first makes chains of loops
    // converge in a single pass.
    auto const loops = c.loop_entries(rk);
    if (loops.empty()) {
      return;
    }
    for_ports_by_cost(
        ports, [&](port_t const p) { return cost_[base + p]; },
        [&](port_t const at, cost_t const best) {
          for (auto k = std::size_t{0U}; k != loops.size(); ++k) {
            auto const entry = backward_ ? loops[k].exit_ : loops[k].entry_;
            auto const exit = backward_ ? loops[k].entry_ : loops[k].exit_;
            if (entry >= ports || exit >= ports) {
              continue;
            }
            auto const weight = m.loop(c.loop_begin(rk) + k);
            if (weight == kInfeasible) {
              continue;
            }
            auto const tc =
                backward_ ? turn(node, entry, at) : turn(node, at, entry);
            if (tc == kInfeasible) {
              continue;
            }
            auto const nc =
                clamp_cost(static_cast<std::uint64_t>(best) + tc + weight);
            if (nc < max_ && nc < cost_[base + exit]) {
              cost_[base + exit] = nc;
            }
          }
        });
  }

  cch_search<P, /* WithPred */ false> search_;

  // Memo. `pos_` is rank indexed and holds the offset of a node's per port
  // block in `cost_`, or kNone. Only the entries in `touched_` are stale
  // between sources, so the reset is proportional to what was computed rather
  // than to the size of the hierarchy.
  std::vector<std::uint32_t> pos_;
  std::vector<cost_t> cost_;
  std::vector<cch_rank_t> touched_;

  std::vector<cch_rank_t> stack_;

  cost_t max_{kInfeasible};
  bool backward_{false};
};

}  // namespace osr
