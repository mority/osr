#pragma once

#include <vector>

#include "utl/verify.h"

#include "osr/routing/cch/cch.h"
#include "osr/routing/cch/customize.h"
#include "osr/routing/cch/query.h"
#include "osr/routing/cch/turns.h"
#include "osr/ways.h"

namespace osr {

// Recursively resolves one CCH arc into the original edges it represents.
//
// No unpacking information is stored with the metric: the middle node of a
// shortcut is recovered by scanning the lower triangles of the arc for the one
// that reproduces its weight. Lower triangle lists are short, so this is much
// cheaper than keeping a middle node per arc entry around.
template <WayAwareProfile P>
struct cch_unpacker {
  cch_unpacker(typename P::parameters const& params,
               ways const& w,
               cch const& c,
               cch_metric const& m)
      : params_{params}, w_{w}, r_{*w.r_}, c_{c}, m_{m} {}

  void unpack_loop(cch_entry_idx_t const idx,
                   std::vector<cch_path_edge>& out) {
    auto const rank = c_.loop_owner(idx);
    utl::verify(rank != cch_rank_t::invalid(), "cch: unknown self loop {}",
                idx);
    expand(rank, rank, c_.loop_[idx], m_.loop(idx),
           c_.loop_is_edge_.test(idx), out);
  }

  void unpack(cch_slot_idx_t const slot,
              cch_entry_idx_t const idx,
              bool const up,
              std::vector<cch_path_edge>& out) {
    auto const l = c_.tail(slot);
    auto const h = c_.adj_head_[slot];
    expand(up ? l : h, up ? h : l,
           up ? c_.up_entries(slot)[idx] : c_.dn_entries(slot)[idx],
           up ? m_.up(c_.up_ofs_[slot] + idx) : m_.dn(c_.dn_ofs_[slot] + idx),
           up ? c_.up_is_edge_.test(c_.up_ofs_[slot] + idx)
              : c_.dn_is_edge_.test(c_.dn_ofs_[slot] + idx),
           out);
  }

private:
  void expand(cch_rank_t const x,
              cch_rank_t const y,
              cch_entry const e,
              cost_t const weight,
              bool const is_edge,
              std::vector<cch_path_edge>& out) {
    auto const xn = c_.order_[x];
    auto const yn = c_.order_[y];

    // an original edge always wins: it is the shortest representation
    auto found = false;
    if (is_edge) {
      for_each_edge(r_, xn,
                  [&](node_idx_t const v, way_idx_t const way,
                      port_t const tail_port, port_t const head_port,
                      distance_t const dist, std::uint16_t, std::uint16_t) {
                    if (found || v != yn || tail_port != e.entry_ ||
                        head_port != e.exit_) {
                      return;
                    }
                    if (cch_edge_cost<P>(params_, w_, way, v, tail_port,
                                         dist) == weight) {
                      found = true;
                    }
                  });
    }
    if (found) {
      out.emplace_back(cch_path_edge{xn, yn, e.entry_, e.exit_, weight});
      return;
    }

    // otherwise: find the lower triangle that produced the weight
    auto const lx = c_.lower_slots(x);
    auto const ly = c_.lower_slots(y);
    auto i = std::size_t{0U};
    auto j = std::size_t{0U};
    auto loops = std::vector<cch_entry_idx_t>{};
    // The tails are derived, so they are recomputed only when the side they
    // belong to advances instead of on every step of the merge.
    auto v1 = i != lx.size() ? c_.tail(lx[i]) : cch_rank_t::invalid();
    auto v2 = j != ly.size() ? c_.tail(ly[j]) : cch_rank_t::invalid();
    while (i != lx.size() && j != ly.size()) {
      auto const s1 = lx[i];
      auto const s2 = ly[j];
      if (v1 < v2) {
        ++i;
        if (i != lx.size()) {
          v1 = c_.tail(lx[i]);
        }
        continue;
      }
      if (v2 < v1) {
        ++j;
        if (j != ly.size()) {
          v2 = c_.tail(ly[j]);
        }
        continue;
      }
      auto const via = v1;  // == v2: the common lower neighbour
      ++i;
      ++j;
      if (i != lx.size()) {
        v1 = c_.tail(lx[i]);
      }
      if (j != ly.size()) {
        v2 = c_.tail(ly[j]);
      }

      auto const in = c_.dn_entries(s1);  // x -> via
      auto const out_e = c_.up_entries(s2);  // via -> y
      auto turns_ready = false;
      for (auto k1 = std::size_t{0U}; k1 != in.size(); ++k1) {
        if (in[k1].entry_ != e.entry_) {
          continue;
        }
        auto const c1 = m_.dn(c_.dn_ofs_[s1] + k1);
        if (c1 == kInfeasible) {
          continue;
        }
        for (auto k2 = std::size_t{0U}; k2 != out_e.size(); ++k2) {
          if (out_e[k2].exit_ != e.exit_) {
            continue;
          }
          auto const c2 = m_.up(c_.up_ofs_[s2] + k2);
          if (c2 == kInfeasible) {
            continue;
          }
          if (!turns_ready) {
            turns_.reset(
                params_, w_, c_,
                [&](cch_entry_idx_t const i) { return m_.loop(i); }, via);
            turns_ready = true;
          }
          loops.clear();
          auto const tc = turns_.chain(in[k1].exit_, out_e[k2].entry_, loops);
          if (tc == kInfeasible ||
              clamp_cost(static_cast<std::uint64_t>(c1) + tc + c2) != weight) {
            continue;
          }
          unpack(s1, static_cast<cch_entry_idx_t>(k1), false, out);
          for (auto const loop : loops) {
            unpack_loop(loop, out);
          }
          unpack(s2, static_cast<cch_entry_idx_t>(k2), true, out);
          return;
        }
      }
    }

    throw utl::fail("cch: could not unpack arc {} -> {} ({} -> {}, cost {})",
                    to_idx(xn), to_idx(yn), e.entry_, e.exit_, weight);
  }

  typename P::parameters const& params_;
  ways const& w_;
  ways::routing const& r_;
  cch const& c_;
  cch_metric const& m_;
  cch_node_turns<P> turns_;
};

}  // namespace osr
