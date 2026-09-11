#pragma once

#include <memory>
#include <string>
#include <vector>

#include "boost/json.hpp"

#include "osr/area/area_graph.h"
#include "osr/area/crossing.h"
#include "osr/routing/path.h"
#include "osr/routing/sharing_data.h"
#include "osr/ways.h"

namespace osr::backend {

// Area routing in the backend: reads the cells written by `osr-area-stats
// --cells-out` - standing in for the stored format the extractor will write -
// into layer 4 (osr/area/area_graph.h) and layer 5 (osr/area/crossing.h).
struct area_routing {
  area_routing(ways const&, boost::json::array const& features);

  bool empty() const noexcept { return graph_ == nullptr || graph_->empty(); }

  // Layer 4, for route(). Refers into this object.
  sharing_data sharing() const { return graph_->sharing(); }
  level_t level_of(node_idx_t const hub) const {
    return graph_->level_of(hub);
  }

  // Layer 5.
  std::vector<area_crossing> crossings(path const& p) const {
    return drawer_->crossings(p);
  }
  std::string const& osm_of(std::size_t const area) const {
    return drawer_->geometry(area).osm_;
  }

  std::size_t n_areas() const noexcept {
    return drawer_ == nullptr ? 0U : drawer_->n_areas();
  }
  std::size_t n_hubs() const noexcept {
    return graph_ == nullptr ? 0U : graph_->n_hubs();
  }
  std::size_t n_connectors_{0U};
  std::size_t n_without_routing_node_{0U};

private:
  std::unique_ptr<area_graph> graph_;
  std::unique_ptr<crossing_drawer> drawer_;  // refers to graph_
};

}  // namespace osr::backend
