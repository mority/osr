#pragma once

#include <memory>
#include <utility>
#include <vector>

#include "geo/latlng.h"

namespace osr {

// The medial axis of an area's free space: the lines running down the middle,
// each point of them as far from one wall as from another. It is the other
// way of routing over an area - Valhalla builds its areas this way
// (mjolnir/areabuilder.cc), and van Toll et al. and Geraerts' Explicit
// Corridor Map are navigation meshes on the same structure.
//
// Here it exists to be compared against the cells (osr/area/cells.h) and to be
// drawn. It is built from a segment Voronoi diagram of the free space's own
// edges (Boost.Polygon, exact integer predicates), keeping the edges that lie
// inside; curved ones - between a wall and a corner - are taken as their
// chords. A crossing over it goes onto the skeleton, along it, and off again,
// which is a real walk, so it never comes out shorter than the geodesic.
struct area_skeleton {
  // `rings` and `obstacles` as for area_free_space: the outline, its holes,
  // and the barriers drawn across it.
  area_skeleton(std::vector<std::vector<geo::latlng>> const& rings,
                std::vector<std::vector<geo::latlng>> const& obstacles = {});
  ~area_skeleton();
  area_skeleton(area_skeleton&&) noexcept;
  area_skeleton& operator=(area_skeleton&&) noexcept;

  bool empty() const;

  std::size_t n_vertices() const;
  std::size_t n_edges() const;

  // Every edge, as its two ends - for drawing.
  std::vector<std::pair<geo::latlng, geo::latlng>> edges() const;

  // The way from `a` to `b` over the skeleton, both ends included. Empty if
  // either end cannot reach it, or the two are on separate pieces.
  std::vector<geo::latlng> path(geo::latlng const& a,
                                geo::latlng const& b) const;

  // That way's length in meters, infinite where there is none.
  double distance(geo::latlng const& a, geo::latlng const& b) const;

  struct impl;

private:
  std::unique_ptr<impl> impl_;
};

}  // namespace osr
