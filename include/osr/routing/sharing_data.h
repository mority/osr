#pragma once

#include <cstdint>
#include <span>
#include <vector>

#include "geo/latlng.h"

#include "utl/verify.h"

#include "osr/routing/additional_edge.h"
#include "osr/types.h"
#include "osr/ways.h"

namespace osr {

// Colors split the walking state after the vehicle is left - in search order:
// trailing foot for forward searches, initial foot for backward searches - by
// where it was left. Each color keeps its own labels, so a cheaper drop-off
// for one color cannot prune the path of another color.
using color_t = std::uint8_t;
constexpr auto const kNoColor = color_t{0xFFU};
constexpr auto const kMaxColors = std::size_t{kNoColor};

inline void verify_additional_edge_count(
    hash_map<node_idx_t, std::vector<additional_edge>> const& edges,
    node_idx_t::value_t const additional_node_offset) {
  for (auto const& [node, node_edges] : edges) {
    if (to_idx(node) < additional_node_offset) {
      continue;
    }
    utl::verify(node_edges.size() <= kMaxWaysPerNode,
                "additional node {} has {} additional edges, maximum is {}",
                node, node_edges.size(), kMaxWaysPerNode);
  }
}

struct sharing_data {
  geo::latlng get_additional_node_coordinates(node_idx_t const n) const {
    return additional_node_coordinates_.at(to_idx(n) - additional_node_offset_);
  }

  bool is_additional_node(node_idx_t const n) const {
    return n != node_idx_t::invalid() && to_idx(n) >= additional_node_offset_;
  }

  way_pos_t get_way_pos(ways::routing const& w,
                        node_idx_t const node,
                        way_idx_t const way) const {
    if (is_additional_node(node)) {
      if (auto it = additional_edges_.find(node);
          it != end(additional_edges_)) {
        auto const& edges = it->second;
        for (auto const [i, ae] : utl::enumerate(edges)) {
          if (ae.underlying_way_ == way) {
            return static_cast<way_pos_t>(i);
          }
        }
      }
      return 0;
    }
    return w.get_way_pos(node, way);
  }

  bitvec<node_idx_t> const* start_allowed_{};
  bitvec<node_idx_t> const* end_allowed_{};
  bitvec<node_idx_t> const* through_allowed_{};

  bool is_colored() const { return !colors_.empty(); }

  std::size_t n_colors() const { return is_colored() ? colors_.size() : 1U; }

  // Calls fn once per color of node n (fn(kNoColor) if uncolored).
  template <typename Fn>
  void for_each_color(node_idx_t const n, Fn&& fn) const {
    if (!is_colored()) {
      fn(kNoColor);
      return;
    }
    for (auto c = std::size_t{0U}; c != colors_.size(); ++c) {
      if (colors_[c]->test(n)) {
        fn(static_cast<color_t>(c));
      }
    }
  }

  node_idx_t::value_t additional_node_offset_{};
  std::vector<geo::latlng> const& additional_node_coordinates_;
  hash_map<node_idx_t, std::vector<additional_edge>> const& additional_edges_;

  // Optional. colors_[c] holds the nodes where leaving the vehicle yields
  // color c: a subset of end_allowed_ for forward searches, of start_allowed_
  // for backward searches, sized like them (additional nodes included).
  // A node may have several colors (overlapping areas). Empty: one uncolored
  // search, as before.
  std::span<bitvec<node_idx_t> const* const> colors_{};
};

inline bool is_allowed(bitvec<node_idx_t> const* b, node_idx_t const n) {
  return b == nullptr || b->test(n);
}

}  // namespace osr
