#pragma once

#include <array>
#include <vector>

#include "ankerl/unordered_dense.h"

#include "osr/routing/cch/cch.h"
#include "osr/routing/cch/customize.h"
#include "osr/routing/cch/turns.h"
#include "osr/routing/dial.h"
#include "osr/routing/profile.h"
#include "osr/types.h"
#include "osr/ways.h"

namespace osr {

// One original edge of an unpacked path.
struct cch_path_edge {
  node_idx_t from_, to_;
  port_t from_port_, to_port_;
  cost_t cost_;
};

// Bidirectional up/down search on the customized contraction hierarchy.
//
// A search state is a (node, port) pair, i.e. exactly the state of the way
// aware profiles: the port says which original edge we arrived at the node
// with (forward search) resp. which original edge we will leave the node with
// (backward search). Turn costs and turn restrictions are evaluated whenever
// two arcs meet at a node, which keeps them out of the stored graph.
// `WithPred` says whether a state also records where it was reached from.
// Only a point to point query reconstructs a path; the one to many searches
// read costs and nothing else. The predecessor block is the larger part of a
// state and is paid for every node the search touches, so it is left out
// entirely rather than filled and ignored.
template <WayAwareProfile P, bool WithPred = true>
struct cch_search {
  using node = typename P::node;

  struct label {
    cost_t cost_;
    cch_rank_t rank_;
    port_t port_;
  };

  struct get_bucket {
    cost_t operator()(label const& l) const noexcept { return l.cost_; }
  };

  // How a state was reached. `arc_` is the index of the entry *within* its
  // slot's list, resp. within the node's loop list -- both bounded by
  // `kMaxPorts * kMaxPorts` -- and not an absolute entry index, which would
  // need the full 64 bits for a planet sized hierarchy. The reconstruction
  // knows the slot and the rank at that point and can add the base back.
  struct pred_block {
    std::array<cch_slot_idx_t, kMaxPorts> slot_;
    std::array<std::uint16_t, kMaxPorts> arc_;
    std::array<port_t, kMaxPorts> pred_port_{};
    std::uint32_t up_{};  // one bit per port
    // predecessor is a target state, i.e. the chain ends here (bit per port)
    std::uint32_t seed_pred_{};
  };
  static_assert(kMaxPorts * kMaxPorts <=
                std::numeric_limits<std::uint16_t>::max());

  struct no_pred {};

  struct entry : std::conditional_t<WithPred, pred_block, no_pred> {
    entry() {
      cost_.fill(kInfeasible);
      if constexpr (WithPred) {
        this->slot_.fill(cch::kNoSlot);
      }
    }

    std::array<cost_t, kMaxPorts> cost_;
  };

  struct terminal {
    cch_rank_t rank_;
    port_t port_;
    cost_t cost_;
  };

  static constexpr auto const kNoLoop = kNoEntry;

  // Cost of the last mile per arriving port of a target node. Ports that can
  // only reach the last mile by driving a self loop first (turn around) get
  // the loop cost added and remember the loop to drive.
  struct target_state {
    target_state() {
      cost_.fill(kInfeasible);
      loop_.fill(kNoLoop);
    }

    std::array<cost_t, kMaxPorts> cost_;
    std::array<cch_entry_idx_t, kMaxPorts> loop_;
    std::array<port_t, kMaxPorts> next_{};
  };

  using map_t = ankerl::unordered_dense::map<cch_rank_t::value_t, entry>;

  void clear() {
    starts_.clear();
    targets_.clear();
  }

  // `port` is the edge the start node is entered with, `cost` the cost of the
  // first mile.
  void add_start(cch const& c, node_idx_t const n, port_t const p,
                 cost_t const cost) {
    if (c.contains(n)) {
      starts_.emplace_back(terminal{c.rank_[n], p, cost});
    }
  }

  // `port` is the edge the target node is entered with, `cost` the cost of the
  // last mile. The turn from `port` onto the last mile is free (only checked
  // for restrictions) - this matches what `best_candidate` does for the
  // reference Dijkstra.
  void add_target(cch const& c, node_idx_t const n, port_t const p,
                  cost_t const cost) {
    if (c.contains(n)) {
      targets_.emplace_back(terminal{c.rank_[n], p, cost});
    }
  }

  bool no_starts() const { return starts_.empty(); }

  // Common setup of the three entry points. The bucket counts stay with the
  // caller: a one directional run only fills one of the two queues, and sizing
  // the other one is not free.
  void reset(cost_t const max) {
    max_ = max;
    best_ = kInfeasible;
    meet_rank_ = cch_rank_t::invalid();
    max_reached_ = false;
    n_settled_ = 0U;
    f_.clear();
    b_.clear();
    target_at_.clear();
    pq_f_.clear();
    pq_b_.clear();
  }

  void seed_starts(map_t& map, dial<label, get_bucket>& pq) {
    for (auto const& s : starts_) {
      if (s.cost_ < max_) {
        relax(map, pq, s.rank_, s.port_, s.cost_, cch::kNoSlot, 0U, false,
              0U);
      }
    }
  }

  bool run(typename P::parameters const& params,
           ways const& w,
           cch const& c,
           cch_metric const& m,
           cost_t const max) {
    auto const& r = *w.r_;

    reset(max);
    pq_f_.n_buckets(max);
    pq_b_.n_buckets(max);

    auto const turn = cch_turn_fn<P>(params, w);
    seed_starts(f_, pq_f_);

    // The target states are "arrived at the target node", the last mile does
    // not add a turn. Therefore they are expanded directly instead of being
    // pushed as ordinary backward states.
    for (auto const& t : targets_) {
      if (t.cost_ >= max) {
        continue;
      }
      auto const [it, inserted] = target_at_.try_emplace(to_idx(t.rank_));
      auto& cost = it->second.cost_[t.port_];
      cost = std::min(cost, t.cost_);
    }
    for (auto& [rank, ts] : target_at_) {
      close_target_loops(r, w.timezones_, params, c, m, cch_rank_t{rank}, ts);
      for (auto p = port_t{0U}; p != kMaxPorts; ++p) {
        if (ts.cost_[p] == kInfeasible || ts.cost_[p] >= max) {
          continue;
        }
        expand_target(c, m, cch_rank_t{rank}, p, ts.cost_[p]);
      }
    }

    while (true) {
      auto progressed = false;

      if (!pq_f_.empty() && pq_f_.get_next_bucket() < best_) {
        step<true>(r, c, m, turn);
        progressed = true;
      }
      if (!pq_b_.empty() && pq_b_.get_next_bucket() < best_) {
        step<false>(r, c, m, turn);
        progressed = true;
      }

      if (!progressed) {
        break;
      }
    }

    return !max_reached_;
  }

  // One directional upward search from the starts, run to exhaustion within
  // the cost bound. RPHAST needs the complete upward search space of the
  // source: there is no opposite search to meet, so the bidirectional
  // termination criterion does not apply and `best_` stays infeasible
  // throughout. Targets are not set up at all -- the scanning phase handles the
  // downward half.
  //
  // `Forward == false` mirrors it into the reversed graph, which is what a last
  // mile query needs: the *backward* states are seeded from the starts and
  // `dn_` arcs are relaxed upward, which is the upward search of the reversed
  // graph. A state (node, port) then carries the cost from that node, leaving
  // by that port, to the seeded location.
  template <bool Forward>
  bool run_upward(typename P::parameters const& params,
                  ways const& w,
                  cch const& c,
                  cch_metric const& m,
                  cost_t const max) {
    reset(max);

    auto& pq = Forward ? pq_f_ : pq_b_;
    pq.n_buckets(max);

    auto const turn = cch_turn_fn<P>(params, w);
    seed_starts(Forward ? f_ : b_, pq);

    while (!pq.empty()) {
      step<Forward>(*w.r_, c, m, turn);
    }

    return !max_reached_;
  }

  bool run_upward(typename P::parameters const& params,
                  ways const& w,
                  cch const& c,
                  cch_metric const& m,
                  cost_t const max,
                  bool const backward) {
    return backward ? run_upward<false>(params, w, c, m, max)
                    : run_upward<true>(params, w, c, m, max);
  }

  // Per port costs of the forward search space, keyed by rank.
  map_t const& forward() const { return f_; }

  // Per port costs of the backward search space, keyed by rank.
  map_t const& backward() const { return b_; }

  cost_t best() const { return best_; }
  bool found() const { return meet_rank_ != cch_rank_t::invalid(); }

  // ---------------------------------------------------------------------

  void relax(map_t& map,
             dial<label, get_bucket>& pq,
             cch_rank_t const rank,
             port_t const p,
             cost_t const cost,
             cch_slot_idx_t const slot,
             std::uint16_t const arc,
             bool const up,
             port_t const pred_port,
             bool const seed_pred = false) {
    if (cost >= max_) {
      max_reached_ = true;
      return;
    }
    auto& e = map[to_idx(rank)];
    if (cost < e.cost_[p]) {
      e.cost_[p] = cost;
      if constexpr (WithPred) {
        auto const bit = std::uint32_t{1U} << p;
        e.slot_[p] = slot;
        e.arc_[p] = arc;
        e.pred_port_[p] = pred_port;
        e.up_ = up ? (e.up_ | bit) : (e.up_ & ~bit);
        e.seed_pred_ = seed_pred ? (e.seed_pred_ | bit) : (e.seed_pred_ & ~bit);
      }
      pq.push(label{cost, rank, p});
    }
  }

  // Dijkstra over the ports of a target node: a port that cannot reach the
  // last mile directly may still be able to after driving a self loop.
  void close_target_loops(ways::routing const& r,
                          timezone_cache_t const& timezones,
                          typename P::parameters const& params,
                          cch const& c,
                          cch_metric const& m,
                          cch_rank_t const rank,
                          target_state& ts) {
    auto const loops = c.loop_entries(rank);
    if (loops.empty()) {
      return;
    }
    auto const n = c.order_[rank];
    auto const ports = n_capped_ports(r, n);

    for_ports_by_cost(
        ports, [&](port_t const p) { return ts.cost_[p]; },
        [&](port_t const at, cost_t const best) {
          // which port `p` can reach `at` by driving one loop?
          for (auto k = std::size_t{0U}; k != loops.size(); ++k) {
            if (loops[k].exit_ != at) {
              continue;
            }
            auto const idx = c.loop_begin(rank) + k;
            auto const weight = m.loop(idx);
            if (weight == kInfeasible) {
              continue;
            }
            for (auto p = port_t{0U}; p != ports; ++p) {
              auto const tc = cch_turn_cost<P>(params, r, timezones, n, p,
                                               loops[k].entry_);
              if (tc == kInfeasible) {
                continue;
              }
              auto const cost =
                  clamp_cost(static_cast<std::uint64_t>(best) + tc + weight);
              if (cost < ts.cost_[p]) {
                ts.cost_[p] = cost;
                ts.loop_[p] = idx;
                ts.next_[p] = at;
              }
            }
          }
        });
  }

  // Relaxes all arcs that arrive at `rank` with port `p` without charging a
  // turn: used for the target states where the last mile follows immediately.
  void expand_target(cch const& c,
                     cch_metric const& m,
                     cch_rank_t const rank,
                     port_t const p,
                     cost_t const cost) {
    for (auto s = c.upper_begin(rank); s != c.upper_end(rank); ++s) {
      auto const h = c.adj_head_[s];
      auto const entries = c.dn_entries(s);
      for (auto k = std::size_t{0U}; k != entries.size(); ++k) {
        if (entries[k].exit_ != p) {
          continue;
        }
        auto const weight = m.dn(c.dn_ofs_[s] + k);
        if (weight == kInfeasible) {
          continue;
        }
        relax(b_, pq_b_, h, entries[k].entry_,
              clamp_cost(static_cast<std::uint64_t>(cost) + weight), s,
              static_cast<std::uint16_t>(k), false, p, true);
      }
    }
  }

  template <bool Forward, typename TurnFn>
  void step(ways::routing const& r,
            cch const& c,
            cch_metric const& m,
            TurnFn const& turn) {
    auto& pq = Forward ? pq_f_ : pq_b_;
    auto& map = Forward ? f_ : b_;
    auto& other = Forward ? b_ : f_;

    auto const l = pq.pop();
    auto const& e = map.at(to_idx(l.rank_));
    if (e.cost_[l.port_] < l.cost_) {
      return;
    }
    ++n_settled_;

    auto const n = c.order_[l.rank_];

    // meeting with the opposite search
    if (auto const it = other.find(to_idx(l.rank_)); it != end(other)) {
      auto const ports = n_capped_ports(r, n);
      for (auto p = port_t{0U}; p != ports; ++p) {
        auto const c2 = it->second.cost_[p];
        if (c2 == kInfeasible) {
          continue;
        }
        auto const tc = Forward ? turn(n, l.port_, p) : turn(n, p, l.port_);
        if (tc == kInfeasible) {
          continue;
        }
        auto const total = clamp_cost(static_cast<std::uint64_t>(l.cost_) +
                                      tc + c2);
        if (total < best_) {
          best_ = total;
          meet_rank_ = l.rank_;
          meet_f_port_ = Forward ? l.port_ : p;
          meet_b_port_ = Forward ? p : l.port_;
        }
      }
    }

    // meeting directly at a target node (free turn onto the last mile)
    if constexpr (Forward) {
      if (auto const it = target_at_.find(to_idx(l.rank_));
          it != end(target_at_)) {
        auto const c2 = it->second.cost_[l.port_];
        if (c2 != kInfeasible) {
          auto const total =
              clamp_cost(static_cast<std::uint64_t>(l.cost_) + c2);
          if (total < best_) {
            best_ = total;
            meet_rank_ = l.rank_;
            meet_f_port_ = l.port_;
            meet_b_port_ = kMaxPorts;  // last mile follows directly
          }
        }
      }
    }

    // self loops keep the node but change the port: driving around the block
    // to get a turn that is forbidden at this intersection
    {
      auto const loops = c.loop_entries(l.rank_);
      for (auto k = std::size_t{0U}; k != loops.size(); ++k) {
        auto const idx =
            static_cast<cch_entry_idx_t>(c.loop_begin(l.rank_) + k);
        auto const weight = m.loop(idx);
        if (weight == kInfeasible) {
          continue;
        }
        auto const tc = Forward ? turn(n, l.port_, loops[k].entry_)
                                : turn(n, loops[k].exit_, l.port_);
        if (tc == kInfeasible) {
          continue;
        }
        relax(map, pq, l.rank_,
              Forward ? loops[k].exit_ : loops[k].entry_,
              clamp_cost(static_cast<std::uint64_t>(l.cost_) + tc + weight),
              cch::kLoopSlot, static_cast<std::uint16_t>(k), Forward, l.port_);
      }
    }

    for (auto s = c.upper_begin(l.rank_); s != c.upper_end(l.rank_); ++s) {
      auto const h = c.adj_head_[s];
      auto const entries = Forward ? c.up_entries(s) : c.dn_entries(s);
      auto const ofs = Forward ? c.up_ofs_[s] : c.dn_ofs_[s];
      for (auto k = std::size_t{0U}; k != entries.size(); ++k) {
        auto const weight = (Forward ? m.up(ofs + k) : m.dn(ofs + k));
        if (weight == kInfeasible) {
          continue;
        }
        auto const tc = Forward ? turn(n, l.port_, entries[k].entry_)
                                : turn(n, entries[k].exit_, l.port_);
        if (tc == kInfeasible) {
          continue;
        }
        relax(map, pq, h, Forward ? entries[k].exit_ : entries[k].entry_,
              clamp_cost(static_cast<std::uint64_t>(l.cost_) + tc + weight), s,
              static_cast<std::uint16_t>(k), Forward, l.port_);
      }
    }
  }

  std::vector<terminal> starts_, targets_;

  dial<label, get_bucket> pq_f_{get_bucket{}}, pq_b_{get_bucket{}};
  map_t f_, b_;
  ankerl::unordered_dense::map<cch_rank_t::value_t, target_state> target_at_;

  cost_t max_{kInfeasible};
  cost_t best_{kInfeasible};
  cch_rank_t meet_rank_{cch_rank_t::invalid()};
  port_t meet_f_port_{0U};
  port_t meet_b_port_{0U};
  bool max_reached_{false};
  std::size_t n_settled_{0U};
};

}  // namespace osr
