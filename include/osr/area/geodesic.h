#pragma once

#include <cstdint>
#include <limits>
#include <vector>

#include "geo/latlng.h"

namespace osr {

// Exact euclidean shortest paths ("geodesics") inside an area polygon.
//
// This is the ground truth the area cost model is measured against: the
// shortest path between two points inside a polygon with holes is a polyline
// that only ever bends at reflex polygon vertices, so running Dijkstra on the
// visibility graph over {connectors} u {polygon vertices} yields the exact
// distance for every connector pair - no discretization, no approximation.
//
// Cost is O(n^3) in the number of polygon vertices (visibility tests dominate),
// which is fine for the areas that occur in OSM (the largest one in the
// Chatelet-Les-Halles dump has 133 vertices) but makes this unsuitable for the
// routing graph itself. It exists to
//   a) produce the per-pair distances the flat per-area cost averages over,
//   b) measure the error that flat cost introduces, per pair and end-to-end,
//   c) later, draw the pretty path across an area once a route is fixed.
struct area_geodesics {
  using vertex_idx_t = std::uint16_t;

  static constexpr auto kNoVertex = std::numeric_limits<vertex_idx_t>::max();
  static constexpr auto kUnreachable = std::numeric_limits<float>::infinity();

  // `rings[0]` is the outer boundary, `rings[1..]` are holes. Rings are closed
  // implicitly: the first point must not be repeated at the end.
  //
  // `connectors` are the points geodesics are computed between (entry/exit
  // nodes). They may sit on the boundary or anywhere inside the polygon;
  // connectors coinciding with a polygon vertex are merged into it.
  //
  // `obstacles` are open polylines that block movement without bounding the
  // area - fences, walls and hedges drawn across it. They are not part of the
  // outline, so they do not affect what counts as inside; they only stop a
  // path crossing them. Their vertices do join the graph, because a path gets
  // round a fence by bending at its end.
  area_geodesics(std::vector<std::vector<geo::latlng>> const& rings,
                 std::vector<geo::latlng> const& connectors,
                 std::vector<std::vector<geo::latlng>> const& obstacles = {});

  std::size_t n_connectors() const noexcept { return n_connectors_; }

  // In-polygon distance in meters between connector `i` and connector `j`,
  // `kUnreachable` if no path inside the polygon exists (connector outside the
  // polygon, or separated from it by holes touching the boundary).
  float distance(std::size_t i, std::size_t j) const;

  // The geodesic itself, including both endpoints. Empty if unreachable.
  // Interior points are always polygon vertices.
  std::vector<geo::latlng> path(std::size_t i, std::size_t j) const;

  // Mean over all unordered connector pairs, ignoring unreachable ones. This
  // is the flat per-area crossing cost. NaN if no pair is reachable.
  float mean_pair_distance() const;

  // Largest |mean_pair_distance() - distance(i, j)| over all reachable pairs,
  // i.e. how far off the flat cost is at its worst in this area.
  float max_pair_error() const;

  bool is_connector_reachable(std::size_t i) const;

  // The visibility predicate this is all built on: is the straight segment
  // a->b contained in the polygon (walls counting as inside)? Exposed so it
  // can be tested directly against an independent implementation, rather than
  // only through the distances it produces.
  static bool is_segment_inside(
      std::vector<std::vector<geo::latlng>> const& rings,
      geo::latlng const& a,
      geo::latlng const& b,
      std::vector<std::vector<geo::latlng>> const& obstacles = {});

private:
  std::size_t n_connectors_{0U};

  // Deduplicated vertex set: connectors first, then polygon vertices.
  std::vector<geo::latlng> vertices_;

  // vertices_ index for each connector (connectors sharing a position with a
  // polygon vertex map onto that vertex).
  std::vector<vertex_idx_t> connector_vertex_;

  // n_connectors_ x n_connectors_, row-major, meters.
  std::vector<float> dist_;

  // n_connectors_ x vertices_.size(), row-major: predecessor on the shortest
  // path from connector `i` to each vertex.
  std::vector<vertex_idx_t> pred_;
};

}  // namespace osr
