#pragma once

#include <cstring>

#include <algorithm>
#include <array>
#include <atomic>
#include <limits>
#include <filesystem>
#include <thread>
#include <vector>

#include "cista/memory_holder.h"
#include "cista/reflection/for_each_field.h"

#include "utl/helpers/algorithm.h"
#include "utl/progress_tracker.h"

#include "osr/routing/cch/cch.h"
#include "osr/routing/cch/parallel.h"
#include "osr/routing/cch/turns.h"
#include "osr/routing/profile.h"
#include "osr/types.h"
#include "osr/ways.h"

namespace osr {

// Identity of a parameter set, used to decide whether a metric can be reused.
//
// The raw object bytes cannot be compared directly: `parameters` structs have
// padding between their members, and that padding is not initialized by
// `parameters{}`, so two parameter sets that are equal field by field can
// still differ byte by byte -- which would silently throw away a matching
// metric and customize again. Concatenating the members skips the padding.
template <typename Params>
std::string params_blob(Params const& p) {
  auto blob = std::string{};
  cista::for_each_field(p, [&](auto const& f) {
    blob.append(reinterpret_cast<char const*>(&f), sizeof(f));
  });
  return blob;
}

// Metric dependent part of the CCH: one weight per arc entry.
//
// The weights are stored as 16 bit values, which halves the largest file the
// hierarchy needs. Two values are reserved: `kWeightInfeasible` for an arc the
// profile cannot use, and `kWeightEscape` for the handful of arcs whose cost
// does not fit -- Germany has two of them in 190 million, at 23.7 hours. Those
// are kept exactly in a sorted side table instead of being saturated, which
// would make an arc look cheaper than it is and let the search return a cost
// that no path can realize.
//
// It can be computed at startup for every profile that is actually used, or
// precomputed once and stored next to the graph. The stored metric is mapped
// read only, so it costs page cache instead of heap and is shared between
// processes. Unlike `cch` it is not checksummed: the file reaches a gigabyte
// on a country sized graph, which would make every open pay for a full scan.
struct cch_metric {
  static constexpr auto const kMode = cista::mode::WITH_STATIC_VERSION;

  using weight_t = std::uint16_t;
  static constexpr auto const kWeightInfeasible =
      std::numeric_limits<weight_t>::max();
  static constexpr auto const kWeightEscape =
      static_cast<weight_t>(kWeightInfeasible - 1U);

  static std::filesystem::path file(std::filesystem::path const&);
  static bool exists(std::filesystem::path const&);
  static cista::wrapped<cch_metric> read(std::filesystem::path const&);
  void write(std::filesystem::path const&) const;

  // A stored metric is only valid for the parameters it was customized with.
  bool matches(std::string_view const params) const {
    return params.size() == params_.size() &&
           std::memcmp(params_.data(), params.data(), params.size()) == 0;
  }

  // The side table is two parallel arrays -- the arc indices, sorted, and
  // their costs -- rather than one array of pairs. A `pair<uint64, uint32>`
  // carries four bytes of tail padding that nothing initializes and cista
  // serializes verbatim, which made two customizations of the same graph
  // produce metric files that differed byte for byte.
  static cost_t decode(vec64<weight_t> const& a,
                       vec<cch_entry_idx_t> const& x_idx,
                       vec<cost_t> const& x_cost,
                       cch_entry_idx_t const i) {
    auto const v = a[i];
    if (v < kWeightEscape) {
      return v;
    }
    if (v == kWeightInfeasible) {
      return kInfeasible;
    }
    auto const it = std::lower_bound(begin(x_idx), end(x_idx), i);
    return x_cost[static_cast<std::uint32_t>(
        std::distance(begin(x_idx), it))];
  }

  cost_t up(cch_entry_idx_t const i) const {
    return decode(up_, up_xi_, up_xc_, i);
  }
  cost_t dn(cch_entry_idx_t const i) const {
    return decode(dn_, dn_xi_, dn_xc_, i);
  }
  cost_t loop(cch_entry_idx_t const i) const {
    return decode(loop_, loop_xi_, loop_xc_, i);
  }

  // Fills one of the arrays from the costs the customization computed.
  static void compress(std::vector<cost_t> const& src,
                       vec64<weight_t>& dst,
                       vec<cch_entry_idx_t>& x_idx,
                       vec<cost_t>& x_cost) {
    dst.clear();
    dst.resize(static_cast<std::uint32_t>(src.size()));
    x_idx.clear();
    x_cost.clear();
    for (auto i = std::size_t{0U}; i != src.size(); ++i) {
      auto const v = src[i];
      if (v == kInfeasible) {
        dst[i] = kWeightInfeasible;
      } else if (v < kWeightEscape) {
        dst[i] = static_cast<weight_t>(v);
      } else {
        dst[i] = kWeightEscape;
        x_idx.emplace_back(static_cast<cch_entry_idx_t>(i));
        x_cost.emplace_back(v);
      }
    }
  }

  vec64<weight_t> up_, dn_, loop_;
  vec<cch_entry_idx_t> up_xi_, dn_xi_, loop_xi_;
  vec<cost_t> up_xc_, dn_xc_, loop_xc_;

  // raw bytes of the profile parameters this metric was customized for
  vec<std::uint8_t> params_;
};

// Weight of the original edge `u --way--> v` for profile P: way cost plus the
// node cost of the target node. The turn cost is *not* included, it depends on
// the edge we arrive at `u` with and is evaluated during the search.
template <WayAwareProfile P>
cost_t cch_edge_cost(typename P::parameters const& params,
                     ways const& w,
                     way_idx_t const way,
                     node_idx_t const target,
                     port_t const tail_port,
                     distance_t const dist) {
  auto const& r = *w.r_;
  auto const nc = P::node_cost(params, r.node_properties_[target]);
  if (nc.cost_ == kInfeasible) {
    return kInfeasible;
  }
  auto const wc =
      P::way_cost(params, r, w.timezones_, way, r.way_properties_[way],
                  port_dir(tail_port), dist, std::nullopt, duration_t{0},
                  direction::kForward);
  if (wc.cost_ == kInfeasible) {
    return kInfeasible;
  }
  return clamp_add(wc, nc).cost_;
}

// Cost of going from arriving port `in` to departing port `out` at a node,
// possibly by driving one or more of its self loops (turn around at a lower
// ranked node). Only nodes that actually have self loops need the closure, for
// all others this is just the plain turn cost.
template <WayAwareProfile P>
struct cch_node_turns {
  void reset(typename P::parameters const& params,
             ways const& w,
             cch const& c,
             auto&& loop_weight,
             cch_rank_t const rank) {
    params_ = &params;
    r_ = w.r_.get();
    timezones_ = &w.timezones_;
    n_ = c.order_[rank];
    ports_ = n_capped_ports(*r_, n_);
    has_loops_ = false;

    auto const loops = c.loop_entries(rank);
    if (loops.empty()) {
      return;
    }

    // Dijkstra between the arriving ports of the node
    utl::fill(dist_, kInfeasible);
    for (auto i = port_t{0U}; i != ports_; ++i) {
      auto* const d = &dist_[i * kMaxPorts];
      d[i] = 0U;
      for_ports_by_cost(
          ports_, [&](port_t const k) { return d[k]; },
          [&](port_t const at, cost_t const best) {
            for (auto e = std::size_t{0U}; e != loops.size(); ++e) {
              auto const weight = loop_weight(c.loop_begin(rank) + e);
              if (weight == kInfeasible || loops[e].exit_ >= ports_) {
                continue;
              }
              auto const turn = plain(at, loops[e].entry_);
              if (turn == kInfeasible) {
                continue;
              }
              auto const cost = clamp_cost(static_cast<std::uint64_t>(best) +
                                           turn + weight);
              if (cost < d[loops[e].exit_]) {
                d[loops[e].exit_] = cost;
                pred_port_[i * kMaxPorts + loops[e].exit_] = at;
                pred_loop_[i * kMaxPorts + loops[e].exit_] =
                    c.loop_begin(rank) + e;
              }
            }
          });
    }
    has_loops_ = true;
  }

  // Cheapest way from `in` to `out`: either the plain turn, or the loops that
  // lead to another arriving port plus the turn from there. Returns the cost
  // and that intermediate port, which is `in` itself if no loop helps.
  std::pair<cost_t, port_t> best_via(port_t const in, port_t const out) const {
    auto best = plain(in, out);
    auto via = in;
    if (has_loops_ && in < ports_) {
      for (auto k = port_t{0U}; k != ports_; ++k) {
        auto const d = dist_[in * kMaxPorts + k];
        if (d == kInfeasible || d == 0U) {
          continue;
        }
        auto const turn = plain(k, out);
        if (turn == kInfeasible) {
          continue;
        }
        auto const total = clamp_cost(static_cast<std::uint64_t>(d) + turn);
        if (total < best) {
          best = total;
          via = k;
        }
      }
    }
    return {best, via};
  }

  // Same, appending the loop entries that have to be driven (in travel order)
  // to `loops`.
  cost_t chain(port_t const in,
               port_t const out,
               std::vector<cch_entry_idx_t>& loops) const {
    auto const [best, via] = best_via(in, out);
    auto const first = loops.size();
    for (auto k = via; k != in;) {
      loops.emplace_back(pred_loop_[in * kMaxPorts + k]);
      k = pred_port_[in * kMaxPorts + k];
    }
    std::reverse(begin(loops) + static_cast<std::ptrdiff_t>(first),
                 end(loops));
    return best;
  }

  cost_t plain(port_t const in, port_t const out) const {
    return cch_turn_cost<P>(*params_, *r_, *timezones_, n_, in, out);
  }

  cost_t operator()(port_t const in, port_t const out) const {
    return best_via(in, out).first;
  }

  typename P::parameters const* params_{nullptr};
  ways::routing const* r_{nullptr};
  timezone_cache_t const* timezones_{nullptr};
  node_idx_t n_{node_idx_t::invalid()};
  port_t ports_{0U};
  bool has_loops_{false};
  std::array<cost_t, kMaxPorts * kMaxPorts> dist_{};
  std::array<port_t, kMaxPorts * kMaxPorts> pred_port_{};
  std::array<cch_entry_idx_t, kMaxPorts * kMaxPorts> pred_loop_{};
};

namespace cch_detail {

// Level of a node in the elimination tree: one more than the deepest of its
// lower neighbours. Nodes on the same level are never adjacent -- if they
// were, one would be a lower neighbour of the other and therefore sit on a
// smaller level -- which is what makes a level a safe unit of parallelism:
// no node of a level reads an arc that another node of the same level writes.
inline vec_map<cch_rank_t, std::uint32_t> compute_levels(cch const& c) {
  auto levels = vec_map<cch_rank_t, std::uint32_t>{};
  levels.resize(c.n_ranks(), 0U);
  // Pushed forward rather than pulled: walking the upper slots of each rank in
  // increasing order reaches every node only after all of its lower neighbours
  // are final, and needs the head of a slot, which is stored.
  for (auto i = cch_rank_t::value_t{0U}; i != c.n_ranks(); ++i) {
    auto const rank = cch_rank_t{i};
    auto const lvl = levels[rank];
    for (auto s = c.upper_begin(rank); s != c.upper_end(rank); ++s) {
      auto& h = levels[c.adj_head_[s]];
      h = std::max(h, lvl + 1U);
    }
  }
  return levels;
}

// Two nodes of the same level can compose the same arc, so the arcs they
// write have to be lowered atomically. The load is a relaxed read and the
// compare exchange only runs when the arc actually improves, so the
// uncontended case costs the same as the plain minimum it replaces.
inline void relax_min(cost_t& weight, cost_t const cost) {
  auto ref = std::atomic_ref<cost_t>{weight};
  auto cur = ref.load(std::memory_order_relaxed);
  while (cost < cur &&
         !ref.compare_exchange_weak(cur, cost, std::memory_order_relaxed,
                                    std::memory_order_relaxed)) {
  }
}

}  // namespace cch_detail

template <WayAwareProfile P>
void customize(typename P::parameters const& params,
               ways const& w,
               cch const& c,
               cch_metric& m,
               unsigned const n_threads = 0U) {
  auto const& r = *w.r_;

  auto arena = cch_arena{n_threads};

  // The customization needs full precision while it composes triangles; the
  // result is compressed into the metric's 16 bit form at the end.
  auto w_up = std::vector<cost_t>(c.up_.size(), kInfeasible);
  auto w_dn = std::vector<cost_t>(c.dn_.size(), kInfeasible);
  auto w_loop = std::vector<cost_t>(c.loop_.size(), kInfeasible);
  auto const blob = params_blob(params);
  m.params_.clear();
  m.params_.resize(static_cast<std::uint32_t>(blob.size()));
  std::memcpy(m.params_.data(), blob.data(), blob.size());

  auto const n = c.n_ranks();

  // ---------------------------------------------------------------------
  // 1. respect the weights of the original edges
  //
  // Every entry is written by exactly one node -- the tail of the original
  // edge it belongs to -- so this phase needs no synchronization.
  // ---------------------------------------------------------------------
  arena.run(n, [&](std::size_t const i) {
        auto const rank = cch_rank_t{static_cast<cch_rank_t::value_t>(i)};
        auto const u = c.order_[rank];
        for_each_edge(
            r, u,
            [&](node_idx_t const v, way_idx_t const way, port_t const tail_port,
                port_t const head_port, distance_t const dist, std::uint16_t,
                std::uint16_t) {
              if (!c.contains(v)) {
                return;
              }
              auto const cost =
                  cch_edge_cost<P>(params, w, way, v, tail_port, dist);
              if (cost == kInfeasible) {
                return;
              }

              auto const ref =
                  c.find_edge(rank, v, cch_entry{tail_port, head_port});
              if (ref.kind_ == cch_arc::kNone) {
                return;
              }
              auto& weight = (ref.kind_ == cch_arc::kUp     ? w_up
                              : ref.kind_ == cch_arc::kDn   ? w_dn
                                                            : w_loop)[ref.idx_];
              weight = std::min(weight, cost);
            });
  });

  // ---------------------------------------------------------------------
  // 2. lower triangles
  //
  // Processed level by level instead of node by node: the contraction order
  // is a valid order, but so is any order that finishes a level before it
  // starts the next one, because a node only ever reads arcs of nodes on
  // smaller levels. Within a level the nodes run in parallel.
  // ---------------------------------------------------------------------
  // The turn closure of a node is read only once it is built, so all slots of
  // a node can share it.
  auto const process_slots = [&](cch_rank_t const rank,
                                 cch_node_turns<P> const& turns,
                                 cch_slot_idx_t const from,
                                 cch_slot_idx_t const to) {
    auto const s_begin = c.upper_begin(rank);
    auto const s_end = c.upper_end(rank);

    for (auto s1 = from; s1 != to; ++s1) {
      auto const x = c.adj_head_[s1];
      auto const in = c.dn_entries(s1);  // x -> via
      if (in.empty()) {
        continue;
      }

      for (auto s2 = s_begin; s2 != s_end; ++s2) {
        auto const y = c.adj_head_[s2];
        auto const out = c.up_entries(s2);  // via -> y
        if (out.empty()) {
          continue;
        }

        // s1 == s2 composes x -> via -> x, i.e. a self loop at x
        auto const is_loop = s1 == s2;
        auto const up = x < y;

        // resolved lazily: only compositions that are actually feasible say
        // something about the completeness of the metric independent arc set
        auto target = cch::kNoSlot;
        auto t_entries = std::span<cch_entry const>{};
        auto* t_weights = static_cast<cost_t*>(nullptr);
        auto resolved = false;

        for (auto e1 = std::size_t{0U}; e1 != in.size(); ++e1) {
          auto const c1 = w_dn[c.dn_ofs_[s1] + e1];
          if (c1 == kInfeasible) {
            continue;
          }
          for (auto e2 = std::size_t{0U}; e2 != out.size(); ++e2) {
            auto const c2 = w_up[c.up_ofs_[s2] + e2];
            if (c2 == kInfeasible) {
              continue;
            }
            auto const turn = turns(in[e1].exit_, out[e2].entry_);
            if (turn == kInfeasible) {
              continue;
            }

            if (!resolved) {
              resolved = true;
              target = is_loop ? cch::kNoSlot
                               : (up ? c.find_slot(x, y) : c.find_slot(y, x));
              // A feasible composition whose arc the hierarchy does not
              // contain cannot happen: the contraction prunes a pair only
              // when the turn is forbidden for every motorized profile, and
              // then it is forbidden for this one too, so the composition
              // would already have been rejected above.
              if (is_loop || target != cch::kNoSlot) {
                t_entries = is_loop ? c.loop_entries(x)
                                    : (up ? c.up_entries(target)
                                          : c.dn_entries(target));
                t_weights = is_loop
                                ? &w_loop[c.loop_begin(x)]
                                : (up ? &w_up[c.up_ofs_[target]]
                                      : &w_dn[c.dn_ofs_[target]]);
              }
            }
            if (t_weights == nullptr) {
              break;
            }

            auto const cost = clamp_cost(static_cast<std::uint64_t>(c1) +
                                         turn + c2);
            auto const idx = cch::find_entry(
                t_entries, cch_entry{in[e1].entry_, out[e2].exit_});
            if (idx == kNoEntry) {
              continue;  // cannot happen, see above
            }
            cch_detail::relax_min(t_weights[idx], cost);
          }
          if (resolved && t_weights == nullptr) {
            break;
          }
        }
      }
    }
  };

  auto const levels = cch_detail::compute_levels(c);
  auto by_level = std::vector<std::vector<cch_rank_t>>{};
  for (auto i = cch_rank_t::value_t{0U}; i != n; ++i) {
    auto const rank = cch_rank_t{i};
    auto const lvl = levels[rank];
    if (by_level.size() <= lvl) {
      by_level.resize(lvl + 1U);
    }
    by_level[lvl].emplace_back(rank);
  }

  auto pt = utl::get_active_progress_tracker_or_activate("osr");
  pt->status("CCH Customization").in_high(n);
  auto done = std::size_t{0U};

  auto const loop_weight = [&](cch_entry_idx_t const i) { return w_loop[i]; };

  auto turns = std::vector<cch_node_turns<P>>{};

  for_each_level(
      arena, by_level, turns,
      [&](cch_node_turns<P>& t, cch_rank_t const rank) {
        t.reset(params, w, c, loop_weight, rank);
      },
      [&](cch_node_turns<P> const&, cch_rank_t const rank) {
        return std::pair{c.upper_begin(rank), c.upper_end(rank)};
      },
      [&](cch_node_turns<P> const& t, cch_rank_t const rank,
          cch_slot_idx_t const s) {
        process_slots(rank, t, s, s + 1U);
      },
      [&](std::size_t const n_nodes) {
        done += n_nodes;
        pt->update_monotonic(done);
      });

  // Freed one at a time: holding all three full precision arrays and the
  // compressed output at once would raise the peak by the size of the result.
  auto const compress_and_free = [](std::vector<cost_t>& src,
                                    vec64<cch_metric::weight_t>& dst,
                                    vec<cch_entry_idx_t>& x_idx,
                                    vec<cost_t>& x_cost) {
    cch_metric::compress(src, dst, x_idx, x_cost);
    src.clear();
    src.shrink_to_fit();
  };
  compress_and_free(w_up, m.up_, m.up_xi_, m.up_xc_);
  compress_and_free(w_dn, m.dn_, m.dn_xi_, m.dn_xc_);
  compress_and_free(w_loop, m.loop_, m.loop_xi_, m.loop_xc_);
}

}  // namespace osr
