#pragma once

#include <cinttypes>

#include <limits>
#include <optional>
#include <span>
#include <vector>

#include "osr/routing/cch/cch.h"
#include "osr/routing/cch/customize.h"
#include "osr/routing/cch/query.h"
#include "osr/routing/cch/sweep.h"
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
// The difference to `rphast` in this directory is where the work goes.
// Eager RPHAST selects T', the upward closure of the whole target set, and
// then sweeps every node of it once per source. Its cost is O(|T'|) no matter
// how many of the targets are actually read, and on a country sized hierarchy
// the closure always climbs into the dense top levels, so even a handful of
// targets pays for the whole ascent.
//
// Lazy RPHAST has no selection phase and no sweep. For a target u it walks up
// the elimination tree until it meets a node whose distance is already known,
// then computes the distances back down that path and memoizes them. The
// second and later targets stop as soon as they hit the part of the tree the
// earlier ones already paid for, so the shared upper region is computed once
// and the marginal cost of a target is the length of its own private tail.
//
// The recurrence is the same one the eager sweep applies, only evaluated on
// demand:
//
//   D[v] = min( D_up[v], min over upward arcs v->w of ( l_dn[w->v] + D[w] ) )
//
// where `D_up` is the tentative distance from the source's upward search.
//
// Turn awareness works exactly as in `rphast`: a label is one cost per
// (node, port) rather than per node, the relaxation charges the turn taken at
// the tail node, and after pulling from the higher ranked neighbours the self
// loops at the node have to be closed so a forbidden turn can be driven
// around.
template <Profile P>
struct lazy_rphast {
  static constexpr auto const kNone = std::numeric_limits<std::uint32_t>::max();

  void clear() {
    search_.clear();
    reset_memo();
  }

  void clear_starts() { search_.clear_starts(); }

  void add_start(cch const& c,
                 node_idx_t const n,
                 port_t const p,
                 cost_t const cost) {
    search_.add_start(c, n, p, cost);
  }

  bool no_starts() const { return search_.no_starts(); }

  // How many nodes the memo holds, i.e. how much of the hierarchy this source
  // actually had to touch. The eager equivalent is `n_selected()`, which is
  // fixed by the target set before the first distance is read.
  std::size_t n_memoized() const { return touched_.size(); }

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
    auto const ok = backward ? search_.run_backward(params, w, c, m, max)
                             : search_.run_forward(params, w, c, m, max);
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
    auto const here = cch_sweep_node{static_cast<std::uint32_t>(cost_.size()),
                                     n_capped_ports(r, c.order_[rk]),
                                     c.order_[rk]};
    cost_.resize(cost_.size() + here.ports_, kInfeasible);
    pos_[to_idx(rk)] = here.base_;
    touched_.push_back(rk);

    // Seed with the tentative distance from the source's upward search. For a
    // node the source never reached this is simply absent.
    auto const& up = backward_ ? search_.backward() : search_.forward();
    if (auto const it = up.find(to_idx(rk)); it != end(up)) {
      for (auto q = port_t{0U}; q != here.ports_; ++q) {
        auto const cu = it->second.cost_[q];
        if (cu < cost_[here.base_ + q]) {
          cost_[here.base_ + q] = cu;
        }
      }
    }

    auto const turn = cch_turn_fn<P>(params, w);

    cch_relax_dn_arcs(
        c, m, rk, backward_, max_, here, cost_,
        [&](cch_rank_t const h) -> std::optional<cch_sweep_node> {
          auto const hidx = pos_[to_idx(h)];
          // kNone cannot happen: `resolve` made every tail final first
          return hidx == kNone
                     ? std::nullopt
                     : std::optional{cch_sweep_node{
                           hidx, n_capped_ports(r, c.order_[h]), c.order_[h]}};
        },
        turn);

    cch_close_loops(
        std::span{cost_}.subspan(here.base_, here.ports_), here.ports_, max_,
        backward_, c.loop_entries(rk),
        [&](std::size_t const k) { return m.loop(c.loop_begin(rk) + k); },
        [&](port_t const in, port_t const out) {
          return turn(here.node_, in, out);
        });
  }

  cch_search<P> search_;

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
