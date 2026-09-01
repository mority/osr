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
template <WayAwareProfile P>
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

  struct entry {
    entry() {
      cost_.fill(kInfeasible);
      slot_.fill(cch::kNoSlot);
    }

    std::array<cost_t, kMaxPorts> cost_;
    std::array<cch_slot_idx_t, kMaxPorts> slot_;
    std::array<cch_entry_idx_t, kMaxPorts> arc_;
    std::array<port_t, kMaxPorts> pred_port_{};
    std::array<std::uint8_t, kMaxPorts> up_{};
    // predecessor is a target state, i.e. the chain ends here
    std::array<std::uint8_t, kMaxPorts> seed_pred_{};
  };

  struct terminal {
    cch_rank_t rank_;
    port_t port_;
    cost_t cost_;
  };

  static constexpr auto const kNoLoop =
      std::numeric_limits<cch_entry_idx_t>::max();

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

  bool empty() const { return starts_.empty() || targets_.empty(); }

  bool no_starts() const { return starts_.empty(); }

  // Stall-on-demand: a state that can be reached cheaper via a higher ranked
  // node cannot be on an up-down path, so it does not have to be expanded.
  // Measured on Hamburg it only prunes ~8% of the settled states while the
  // check itself costs more than that, so it is off by default.
  static constexpr auto const kStallOnDemand = false;

  bool run(typename P::parameters const& params,
           ways const& w,
           cch const& c,
           cch_metric const& m,
           cost_t const max) {
    auto const& r = *w.r_;

    max_ = max;
    best_ = kInfeasible;
    meet_rank_ = cch_rank_t::invalid();
    max_reached_ = false;
    n_settled_ = 0U;
    f_.clear();
    b_.clear();
    target_at_.clear();
    pq_f_.clear();
    pq_f_.n_buckets(max);
    pq_b_.clear();
    pq_b_.n_buckets(max);

    auto const turn = [&](node_idx_t const n, port_t const in,
                          port_t const out) {
      return cch_turn_cost<P>(params, r, w.timezones_, n, in, out);
    };

    for (auto const& s : starts_) {
      if (s.cost_ < max) {
        relax(f_, pq_f_, s.rank_, s.port_, s.cost_, cch::kNoSlot, 0U, 0U, 0U);
      }
    }

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

  // Forward-only upward search from the starts, run to exhaustion within the
  // cost bound. RPHAST needs the complete upward search space of the source:
  // there is no opposite search to meet, so the bidirectional termination
  // criterion does not apply and `best_` stays infeasible throughout. Targets
  // are not set up at all -- the scanning phase handles the downward half.
  bool run_forward(typename P::parameters const& params,
                   ways const& w,
                   cch const& c,
                   cch_metric const& m,
                   cost_t const max) {
    auto const& r = *w.r_;

    max_ = max;
    best_ = kInfeasible;
    meet_rank_ = cch_rank_t::invalid();
    max_reached_ = false;
    n_settled_ = 0U;
    f_.clear();
    b_.clear();
    target_at_.clear();
    pq_f_.clear();
    pq_f_.n_buckets(max);
    pq_b_.clear();

    auto const turn = [&](node_idx_t const n, port_t const in,
                          port_t const out) {
      return cch_turn_cost<P>(params, r, w.timezones_, n, in, out);
    };

    for (auto const& s : starts_) {
      if (s.cost_ < max) {
        relax(f_, pq_f_, s.rank_, s.port_, s.cost_, cch::kNoSlot, 0U, 0U, 0U);
      }
    }

    while (!pq_f_.empty()) {
      step<true>(r, c, m, turn);
    }

    return !max_reached_;
  }

  // Per port costs of the forward search space, keyed by rank.
  map_t const& forward() const { return f_; }

  cost_t best() const { return best_; }
  bool found() const { return meet_rank_ != cch_rank_t::invalid(); }

  // ---------------------------------------------------------------------

  void relax(map_t& map,
             dial<label, get_bucket>& pq,
             cch_rank_t const rank,
             port_t const p,
             cost_t const cost,
             cch_slot_idx_t const slot,
             cch_entry_idx_t const arc,
             std::uint8_t const up,
             port_t const pred_port,
             std::uint8_t const seed_pred = 0U) {
    if (cost >= max_) {
      max_reached_ = true;
      return;
    }
    auto& e = map[to_idx(rank)];
    if (cost < e.cost_[p]) {
      e.cost_[p] = cost;
      e.slot_[p] = slot;
      e.arc_[p] = arc;
      e.up_[p] = up;
      e.pred_port_[p] = pred_port;
      e.seed_pred_[p] = seed_pred;
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
    auto const ports = std::min(n_ports(r, n), kMaxPorts);

    auto done = std::uint32_t{0U};
    for (auto step = port_t{0U}; step != ports; ++step) {
      auto best = kInfeasible;
      auto at = kMaxPorts;
      for (auto p = port_t{0U}; p != ports; ++p) {
        if ((done & (std::uint32_t{1U} << p)) == 0U && ts.cost_[p] < best) {
          best = ts.cost_[p];
          at = p;
        }
      }
      if (at == kMaxPorts) {
        break;
      }
      done |= std::uint32_t{1U} << at;

      // which port `p` can reach `at` by driving one loop?
      for (auto k = std::size_t{0U}; k != loops.size(); ++k) {
        if (loops[k].exit_ != at) {
          continue;
        }
        auto const idx = static_cast<cch_entry_idx_t>(c.loop_begin(rank) + k);
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
    }
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
              static_cast<cch_entry_idx_t>(k), 0U, p, 1U);
      }
    }
  }

  // Is there a cheaper way to this state coming down from a higher ranked
  // node? Then it is not part of any up-down path and can be skipped.
  template <bool Forward, typename TurnFn>
  bool stalled(ways::routing const& r,
               cch const& c,
               cch_metric const& m,
               map_t const& map,
               label const& l,
               TurnFn const& turn) const {
    for (auto s = c.upper_begin(l.rank_); s != c.upper_end(l.rank_); ++s) {
      auto const h = c.adj_head_[s];
      auto const it = map.find(to_idx(h));
      if (it == end(map)) {
        continue;
      }
      // forward: arcs h -> this node, backward: arcs this node -> h
      auto const entries = Forward ? c.dn_entries(s) : c.up_entries(s);
      auto const ofs = Forward ? c.dn_ofs_[s] : c.up_ofs_[s];
      auto const h_ports =
          std::min(n_ports(r, c.order_[h]), kMaxPorts);
      for (auto k = std::size_t{0U}; k != entries.size(); ++k) {
        if ((Forward ? entries[k].exit_ : entries[k].entry_) != l.port_) {
          continue;
        }
        auto const weight = (Forward ? m.dn(ofs + k) : m.up(ofs + k));
        if (weight == kInfeasible || weight >= l.cost_) {
          continue;
        }
        for (auto p = port_t{0U}; p != h_ports; ++p) {
          auto const c2 = it->second.cost_[p];
          if (c2 == kInfeasible || c2 + weight >= l.cost_) {
            continue;
          }
          auto const tc = Forward ? turn(c.order_[h], p, entries[k].entry_)
                                  : turn(c.order_[h], entries[k].exit_, p);
          if (tc == kInfeasible) {
            continue;
          }
          if (clamp_cost(static_cast<std::uint64_t>(c2) + tc + weight) <
              l.cost_) {
            return true;
          }
        }
      }
    }
    return false;
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

    if constexpr (kStallOnDemand) {
      if (stalled<Forward>(r, c, m, map, l, turn)) {
        return;
      }
    }

    // meeting with the opposite search
    if (auto const it = other.find(to_idx(l.rank_)); it != end(other)) {
      auto const ports = std::min(n_ports(r, n), kMaxPorts);
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
              cch::kLoopSlot, idx, Forward ? std::uint8_t{1U} : std::uint8_t{0U},
              l.port_);
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
              static_cast<cch_entry_idx_t>(k),
              Forward ? std::uint8_t{1U} : std::uint8_t{0U}, l.port_);
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
