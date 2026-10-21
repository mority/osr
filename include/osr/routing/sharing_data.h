#pragma once

#include <cstdint>
#include <vector>

#include "geo/latlng.h"

#include "utl/verify.h"

#include "osr/routing/additional_edge.h"
#include "osr/types.h"
#include "osr/ways.h"

namespace osr {

// Bit i set: class i (see sharing_data::classes_).
using sharing_class_mask_t = std::uint64_t;
constexpr auto const kMaxSharingClasses = 64U;

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

  node_idx_t::value_t additional_node_offset_{};
  std::vector<geo::latlng> const& additional_node_coordinates_;
  hash_map<node_idx_t, std::vector<additional_edge>> const& additional_edges_;

  // Optional, car_sharing only: splits the walk that follows the vehicle in
  // search direction (after the drop-off forward, before the pickup backward)
  // into `n_classes_` classes. A node maps to the classes in which the vehicle
  // may be left there (forward) or taken there (backward); with `classes_`
  // set, nodes without an entry allow no such switch. route_one_to_many()
  // then reports one result per class and destination, at index
  // `class * to.size() + destination`: one search instead of one per class
  // whose searches differ only in that end.
  hash_map<node_idx_t, sharing_class_mask_t> const* classes_{};
  std::uint8_t n_classes_{1U};
};

inline bool is_allowed(bitvec<node_idx_t> const* b, node_idx_t const n) {
  return b == nullptr || b->test(n);
}

}  // namespace osr
