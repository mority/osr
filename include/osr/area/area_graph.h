#pragma once

#include <cstdint>
#include <functional>
#include <optional>
#include <string>
#include <vector>

#include "geo/latlng.h"

#include "osr/area/cells.h"
#include "osr/area/walkable.h"
#include "osr/routing/additional_edge.h"
#include "osr/routing/sharing_data.h"
#include "osr/types.h"

namespace osr {

// Layer 4 of area routing (the layers are listed in osr/area/walkable.h): the
// cells as a graph the router walks.
//
// Every cell becomes a hub node, numbered after the dataset's own nodes. A
// connector is joined to the hub of its cell at half the cell's cost, two
// neighbouring hubs at half the sum of both costs, so a route through cells
// a, b, c pays c_a + c_b + c_c. The router gets the graph as sharing_data's
// additional nodes and edges, which route() and every search already carry;
// the foot profile (profiles/foot.h) follows them, and steps onto a hub only
// from the floor its area is on.

// The floor an area's hubs are on. An untagged area is at street level: open
// to routes on the ground floor or on untagged ways, not to one coming up
// from a station below. An area on several levels (stairs mapped as an area)
// is open from any floor.
level_t area_floor(area_levels const&);

// One area as layer 3 left it.
struct area_crossing_input {
  std::string osm_{};  // "way/1", for reporting
  level_t level_{0.F};  // see area_floor
  area_cells cells_{};
  std::vector<std::int64_t> connector_nodes_{};  // OSM node of each connector
  std::vector<geo::latlng> connector_pos_{};
};

// The routing node an OSM node is, if it is one.
using node_lookup =
    std::function<std::optional<node_idx_t>(std::int64_t osm_node)>;

struct area_graph {
  // Hubs are numbered from `first_hub` on - the dataset's node count.
  area_graph(node_idx_t::value_t first_hub,
             std::vector<area_crossing_input> const&,
             node_lookup const&);

  bool empty() const noexcept { return hub_pos_.empty(); }

  // The graph, for route(). Refers into this object.
  sharing_data sharing() const;

  bool is_hub(node_idx_t) const noexcept;

  // Which input area a hub belongs to, and the floor it is on.
  std::size_t area_of(node_idx_t hub) const;
  level_t level_of(node_idx_t hub) const;

  // Which connector of input area `area` a routing node is, if any.
  std::optional<std::size_t> connector_of(std::size_t area,
                                          node_idx_t) const;

  std::size_t n_hubs() const noexcept { return hub_pos_.size(); }

  // Connectors with a cell, and those left out because their OSM node is no
  // routing node - nothing to attach a spoke to.
  std::size_t n_connectors_{0U};
  std::size_t n_without_routing_node_{0U};

private:
  node_idx_t::value_t first_hub_;
  std::vector<geo::latlng> hub_pos_;
  std::vector<level_t> hub_level_;
  std::vector<std::uint32_t> hub_area_;
  hash_map<node_idx_t, std::vector<additional_edge>> edges_;
  std::vector<hash_map<node_idx_t, std::size_t>> connector_of_;
};

}  // namespace osr
