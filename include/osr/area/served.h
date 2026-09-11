#pragma once

#include <cstdint>
#include <optional>
#include <vector>

#include "geo/latlng.h"

#include "osr/area/walkable.h"
#include "osr/types.h"

namespace osr {

// Layer 1 of area routing: does an area need routing across it, or do the
// ways already mapped carry walkers well enough?
//
// Where footways are mapped across a plaza, routing follows them and the area
// gets no cells. It needs cells where the mapped ways - its own outline and
// whatever lies inside it - make crossing it a detour, or leave an entrance
// inside it that no mapped way reaches.

// A point where walkers enter or leave an area: a node of its outline or of a
// hole's ring that the street network meets, or a staircase, lift or footway
// stub ending inside it.
struct area_connector {
  std::int64_t node_{0};
  geo::latlng pos_{};
  bool on_ring_{true};  // on the outline or a hole's ring, not inside
  std::optional<std::size_t> outline_pos_{};  // index into rings_[0]
};

// An edge of a way mapped inside an area, between two OSM nodes.
struct network_edge {
  std::int64_t a_{0}, b_{0};
};

struct served_params {
  // The mapped ways may be at most this much longer than the shortest walk
  // across before the area counts as needing cells.
  double max_detour_{1.2};

  // A pair of outline connectors matters only if walking round the outline
  // is at least this much longer than crossing: otherwise the crossing gains
  // nothing worth modelling.
  double min_crossing_gain_{1.2};

  // Closer than this, two connectors are not a crossing.
  double min_pair_distance_{1.0};
};

// Walking distances between the connectors over the area's rings and the
// `mapped` ways, k x k in meters. Infinite where the mapped ways do not join
// two connectors - a stub ending inside an area cannot be reached without
// crossing it.
std::vector<double> network_distances(
    walkable_area const&,
    std::vector<network_edge> const& mapped,
    hash_map<std::int64_t, geo::latlng> const& node_pos,
    std::vector<area_connector> const&);

// Connectors inside the area that no mapped way reaches. Any one of them is
// reason enough to give the area cells.
std::size_t n_stranded(std::vector<area_connector> const&,
                       std::vector<double> const& network);

// Before any geometry: the shortest walk is never shorter than the straight
// line, so where every pair of ring connectors is within max_detour of the
// straight line over mapped ways, and nothing is stranded, the area is
// served - without building the costly geodesics (layer 2). Only a shortcut:
// false means "not known", not "needs cells".
bool served_without_geometry(std::vector<area_connector> const&,
                             std::vector<double> const& network,
                             served_params const& = {});

// Which pairs a crossing shortens enough to matter (see min_crossing_gain_),
// k x k. Pairs not both on the outline stay relevant. Also returns, for each
// outline pair it looked at, the detour of walking round the outline.
struct relevance {
  std::vector<bool> relevant_;
  std::vector<double> outline_detours_;
};
relevance relevant_pairs(walkable_area const&,
                         std::vector<area_connector> const&,
                         std::vector<float> const& geodesic,
                         served_params const& = {});

struct served_verdict {
  bool served_{false};
  double worst_detour_{0.0};  // mapped ways / shortest walk, relevant pairs
  std::size_t n_relevant_{0U};
  std::size_t n_stranded_{0U};
  std::vector<double> detours_{};  // per relevant ring pair
};

// With the geodesics (layer 2): served where every relevant pair of ring
// connectors is within max_detour of its shortest walk over mapped ways,
// and nothing is stranded. Only ring-to-ring pairs can judge the area: they
// are joined by its rings whatever else is mapped, while an interior stub is
// either reached by a mapped way or stranded.
served_verdict served_by_geodesics(std::vector<area_connector> const&,
                                   std::vector<double> const& network,
                                   std::vector<float> const& geodesic,
                                   std::vector<bool> const& relevant,
                                   served_params const& = {});

}  // namespace osr
