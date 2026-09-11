#pragma once

#include <memory>
#include <mutex>
#include <string>
#include <vector>

#include "geo/latlng.h"

#include "osr/area/area_graph.h"
#include "osr/area/geodesic.h"
#include "osr/routing/path.h"
#include "osr/types.h"

namespace osr {

// Layer 5 of area routing (the layers are listed in osr/area/walkable.h):
// drawing a crossing.
//
// The router crosses an area through the hubs of its cells (layer 4), and a
// path built from that runs in straight lines through points nobody walks
// to. Once the route is fixed, each crossing - entering the area at one
// connector, leaving it at another - is drawn as the shortest walkable line
// between the two: exactly what layer 2 computes, so it is reused here.

// What drawing a crossing needs of an area, indexed like the area_graph's
// input.
struct area_geometry {
  std::string osm_{};
  std::vector<std::vector<geo::latlng>> rings_{};  // outline, then holes
  std::vector<std::vector<geo::latlng>> barriers_{};
  std::vector<geo::latlng> connectors_{};  // as in area_crossing_input
};

// One crossing of a routed path: segments first_segment_ .. last_segment_ run
// connector -> hub(s) -> connector. The router paid model_distance_ for it;
// geodesic_ is the shortest walkable line between the two connectors,
// geodesic_distance_ long (both meters), and empty if the connectors are
// not known to the area.
struct area_crossing {
  std::size_t first_segment_{0U}, last_segment_{0U};
  std::size_t area_{0U};
  level_t level_{kNoLevel};
  std::vector<geo::latlng> geodesic_{};
  double model_distance_{0.0};
  double geodesic_distance_{-1.0};
};

struct crossing_drawer {
  // `graph` must outlive the drawer.
  crossing_drawer(area_graph const& graph, std::vector<area_geometry>);

  std::vector<area_crossing> crossings(path const&) const;

  // The path as a walker should see it: each crossing's hub segments
  // replaced by one segment along its shortest walkable line. The path's
  // totals stay what the router computed; the new segment's distance is the
  // line's length.
  path drawn(path const&) const;

  std::size_t n_areas() const noexcept { return areas_.size(); }
  area_geometry const& geometry(std::size_t const area) const {
    return areas_[area]->geometry_;
  }

private:
  struct entry {
    area_geometry geometry_;

    // Layer 2 for this area, built the first time a route crosses it.
    std::once_flag once_;
    std::unique_ptr<area_geodesics> geodesics_;
  };

  area_geodesics const& geodesics(std::size_t area) const;

  area_graph const& graph_;
  std::vector<std::unique_ptr<entry>> areas_;
};

}  // namespace osr
