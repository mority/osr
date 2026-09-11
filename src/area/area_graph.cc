#include "osr/area/area_graph.h"

#include <cmath>
#include <bit>

#include "utl/verify.h"

namespace osr {

level_t area_floor(area_levels const& l) {
  if (l.any_) {
    return level_t{0.F};
  }
  if (std::popcount(l.bits_) == 1) {
    return level_t{static_cast<std::uint8_t>(std::countr_zero(l.bits_))};
  }
  return kNoLevel;
}

area_graph::area_graph(node_idx_t::value_t const first_hub,
                       std::vector<area_crossing_input> const& areas,
                       node_lookup const& lookup)
    : first_hub_{first_hub} {
  for (auto ai = std::size_t{0U}; ai != areas.size(); ++ai) {
    auto const& a = areas[ai];
    auto& connector_of = connector_of_.emplace_back();
    if (a.cells_.n_cells() == 0U) {
      continue;
    }
    utl::verify(a.connector_nodes_.size() == a.cells_.n_connectors() &&
                    a.connector_pos_.size() == a.cells_.n_connectors(),
                "area graph: {} has {} connectors, {} nodes, {} positions",
                a.osm_, a.cells_.n_connectors(), a.connector_nodes_.size(),
                a.connector_pos_.size());

    auto const first = hub_pos_.size();
    for (auto const& h : hub_positions(a.cells_, a.connector_pos_)) {
      hub_pos_.push_back(h);
      hub_level_.push_back(a.level_);
      hub_area_.push_back(static_cast<std::uint32_t>(ai));
    }
    auto const hub = [&](std::size_t const cell) {
      return node_idx_t{
          static_cast<node_idx_t::value_t>(first_hub_ + first + cell)};
    };
    auto const half = [](double const x) {
      return static_cast<distance_t>(std::lround(x / 2.0));
    };
    auto const link = [&](node_idx_t const x, node_idx_t const y,
                          distance_t const d) {
      edges_[x].push_back(additional_edge{.to_ = y, .distance_ = d});
      edges_[y].push_back(additional_edge{.to_ = x, .distance_ = d});
    };

    auto const& cost = a.cells_.cost_;
    for (auto i = std::size_t{0U}; i != a.connector_nodes_.size(); ++i) {
      auto const cell = a.cells_.connector_cell_[i];
      if (cell == area_cells::kNoCell) {
        continue;  // reaches nothing inside the area
      }
      ++n_connectors_;
      auto const n = lookup(a.connector_nodes_[i]);
      if (!n.has_value()) {
        ++n_without_routing_node_;
        continue;
      }
      if (connector_of.emplace(*n, i).second) {
        link(*n, hub(cell), half(cost[cell]));
      }
    }
    for (auto x = std::size_t{0U}; x != a.cells_.n_cells(); ++x) {
      for (auto y = x + 1U; y != a.cells_.n_cells(); ++y) {
        if (a.cells_.is_neighbour(static_cast<area_cells::cell_idx_t>(x),
                                  static_cast<area_cells::cell_idx_t>(y))) {
          link(hub(x), hub(y), half(cost[x] + cost[y]));
        }
      }
    }
  }
}

sharing_data area_graph::sharing() const {
  return sharing_data{.start_allowed_ = nullptr,
                      .end_allowed_ = nullptr,
                      .through_allowed_ = nullptr,
                      .additional_node_offset_ = first_hub_,
                      .additional_node_coordinates_ = hub_pos_,
                      .additional_edges_ = edges_,
                      .additional_node_levels_ = &hub_level_};
}

bool area_graph::is_hub(node_idx_t const n) const noexcept {
  return n != node_idx_t::invalid() && to_idx(n) >= first_hub_ &&
         to_idx(n) < first_hub_ + hub_pos_.size();
}

std::size_t area_graph::area_of(node_idx_t const hub) const {
  return hub_area_[to_idx(hub) - first_hub_];
}

level_t area_graph::level_of(node_idx_t const hub) const {
  return hub_level_[to_idx(hub) - first_hub_];
}

std::optional<std::size_t> area_graph::connector_of(std::size_t const area,
                                                    node_idx_t const n) const {
  auto const& m = connector_of_[area];
  auto const it = m.find(n);
  return it == end(m) ? std::nullopt : std::optional{it->second};
}

}  // namespace osr
