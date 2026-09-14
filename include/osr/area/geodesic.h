#pragma once

#include <cstdint>
#include <limits>
#include <memory>
#include <utility>
#include <vector>

#include "geo/latlng.h"

namespace osr {

// Layer 2 of area routing (the layers are listed in osr/area/walkable.h): the
// ground truth everything else is measured against.
//
// Where one can walk in an area is its free space: the outline, minus every
// hole (inner rings, buildings) and every barrier (fences, walls, hedges,
// given a width so they can separate two sides at all). It is built as one
// valid polygon - holes and barriers unioned, then cut out of the outline -
// and snapped to a 0.1 mm grid, so that every geometric decision below is
// exact integer arithmetic: whether a line crosses a wall, touches a corner or
// runs along an edge never depends on a tolerance.
//
// A straight segment is walkable when it stays in the closed free space
// (walls count as walkable), under one extra rule: where the boundary is
// pinched to a single point - two buildings touching at a corner - a path may
// not pass through that point from one side to the other. A pinch is not a
// doorway.

// Two ways to find which points see each other. Both use the same geometric
// tests, so any disagreement between them is a bug in one of the algorithms.
enum class visibility_algorithm : std::uint8_t {
  // Every pair tested against every wall: O(n^3).
  kPairwise,

  // Lee's rotational sweep (de Berg et al., Computational Geometry, ch. 15.2):
  // a ray turns around each point, keeping the walls it crosses in order of
  // distance, so each target needs only the nearest one tested. O(n^2 log n).
  kSweep
};

// How far a drawn path keeps off walls, fences and buildings, in meters.
// Nobody walks with a shoulder on the wall, and a path that grazes every
// corner looks it. With a clearance the free space shrinks by that much - the
// outline moves in, holes and barriers grow - so a gap narrower than twice it
// closes, and a point left outside (an entrance is *on* the outline) steps
// onto the free space.
//
// Layer 5 draws with it (see crossing.h). Everything measured against the
// geodesics - which areas need cells, what a crossing costs - uses no
// clearance: that is the true shortest walk, and the decisions are made on it.
constexpr auto kWallClearance = 0.25;

struct geodesic_options {
  double clearance_{0.0};
  visibility_algorithm algorithm_{visibility_algorithm::kSweep};
};

// The free space of one area, for testing many segments against it.
struct area_free_space {
  // `rings[0]` is the outline, `rings[1..]` are holes. Rings are closed
  // implicitly: the first point must not be repeated at the end. `obstacles`
  // are open polylines that block movement without bounding the area.
  area_free_space(std::vector<std::vector<geo::latlng>> const& rings,
                  std::vector<std::vector<geo::latlng>> const& obstacles = {},
                  double clearance = 0.0);
  ~area_free_space();
  area_free_space(area_free_space&&) noexcept;
  area_free_space& operator=(area_free_space&&) noexcept;

  // False if the outline is not a valid polygon (it crosses itself, say):
  // then nothing is walkable.
  bool valid() const;

  // Is the straight segment a -> b walkable? An end outside the free space
  // but within 1 m of it - on a wall, which the clearance leaves just outside
  // - first steps onto it.
  bool is_segment_inside(geo::latlng const& a, geo::latlng const& b) const;

  // The free space itself: the outline minus holes and barriers, cleaned and
  // snapped to the grid. Outer rings first in each piece, then its holes;
  // rings are not closed. For building something else on the same ground -
  // a skeleton, say.
  std::vector<std::vector<geo::latlng>> rings() const;

  struct impl;

private:
  std::unique_ptr<impl> impl_;
};

// Exact euclidean shortest paths ("geodesics") inside an area.
//
// The shortest path between two points in free space is a polyline that only
// ever bends at reflex corners, so running Dijkstra on the visibility graph
// over {connectors} u {reflex corners} yields the exact distance for every
// connector pair - no discretization, no approximation. It exists to
//   a) produce the per-pair distances the flat per-area cost averages over,
//   b) measure the error that flat cost introduces, per pair and end-to-end,
//   c) draw the pretty path across an area once a route is fixed (layer 5).
struct area_geodesics {
  using vertex_idx_t = std::uint32_t;

  static constexpr auto kNoVertex = std::numeric_limits<vertex_idx_t>::max();
  static constexpr auto kUnreachable = std::numeric_limits<float>::infinity();

  // `rings` and `obstacles` as for area_free_space.
  //
  // `connectors` are the points geodesics are computed between (entry/exit
  // nodes). They may sit on the boundary or anywhere inside; a connector
  // within 1 cm of a corner of the free space is that corner. A connector
  // outside the free space but within 1 m of it - an entrance on the outline,
  // which the clearance leaves just outside - steps onto it by the shortest
  // way: the step counts in its distances, and its paths start and end at the
  // connector itself. Connectors further out - outside the outline, deep in a
  // building - are unreachable.
  area_geodesics(std::vector<std::vector<geo::latlng>> const& rings,
                 std::vector<geo::latlng> const& connectors,
                 std::vector<std::vector<geo::latlng>> const& obstacles = {},
                 geodesic_options = {});

  std::size_t n_connectors() const noexcept { return n_connectors_; }

  // See area_free_space::valid. Everything is unreachable if not.
  bool valid() const noexcept { return valid_; }

  // Whether the visibility graph came from the sweep. It falls back to the
  // pairwise test if snapping to the grid left two walls crossing, which the
  // sweep's ordering of walls cannot cope with.
  bool swept() const noexcept { return swept_; }

  // In-area distance in meters between connector `i` and connector `j`,
  // `kUnreachable` if no path inside the area exists.
  float distance(std::size_t i, std::size_t j) const;

  // The geodesic itself, including both endpoints. Empty if unreachable.
  // Interior points are always corners of the free space.
  std::vector<geo::latlng> path(std::size_t i, std::size_t j) const;

  // Mean over all unordered connector pairs, ignoring unreachable ones. This
  // is the flat per-area crossing cost. NaN if no pair is reachable.
  float mean_pair_distance() const;

  // Largest |mean_pair_distance() - distance(i, j)| over all reachable pairs,
  // i.e. how far off the flat cost is at its worst in this area.
  float max_pair_error() const;

  bool is_connector_reachable(std::size_t i) const;

  // The visibility graph: pairs (a < b) of points that see each other, sorted,
  // and where those points are. For comparing the two algorithms, and for
  // drawing.
  std::vector<std::pair<vertex_idx_t, vertex_idx_t>> const& edges() const {
    return edges_;
  }
  std::vector<geo::latlng> const& vertices() const { return vertices_; }

  // The visibility predicate this is all built on, for a single segment (see
  // area_free_space, which is the better choice for many).
  static bool is_segment_inside(
      std::vector<std::vector<geo::latlng>> const& rings,
      geo::latlng const& a,
      geo::latlng const& b,
      std::vector<std::vector<geo::latlng>> const& obstacles = {});

private:
  std::size_t n_connectors_{0U};
  bool valid_{false};
  bool swept_{false};

  // Every point of the graph: connectors and corners of the free space.
  std::vector<geo::latlng> vertices_;

  // vertices_ index for each connector (connectors sharing a position map onto
  // the same point).
  std::vector<vertex_idx_t> connector_vertex_;

  // Each connector as given, and the length of its step onto the free space
  // (zero if it is on it).
  std::vector<geo::latlng> connector_pos_;
  std::vector<double> offset_;

  std::vector<std::pair<vertex_idx_t, vertex_idx_t>> edges_;

  // n_connectors_ x n_connectors_, row-major, meters.
  std::vector<float> dist_;

  // n_connectors_ x vertices_.size(), row-major: predecessor on the shortest
  // path from connector `i` to each point.
  std::vector<vertex_idx_t> pred_;
};

}  // namespace osr
