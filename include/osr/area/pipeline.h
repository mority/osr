#pragma once

#include <cstdint>
#include <filesystem>
#include <optional>
#include <vector>

#include "geo/latlng.h"

#include "osr/area/cells.h"
#include "osr/area/geodesic.h"
#include "osr/area/served.h"
#include "osr/area/walkable.h"
#include "osr/types.h"

namespace osr {

// Layers 0 to 3 of area routing (the layers are listed in
// osr/area/walkable.h) run over OSM data: collect_areas reads a file and
// finds its walkable areas with everything that matters about them (layer
// 0); prepare_area then decides, area by area, whether it needs cells and
// builds them (layers 1 to 3).

struct area_options {
  walkable_options walkable_{};
  bool strict_levels_{false};  // see area_levels::matches
  bool merge_{true};  // merge overlapping areas (see merge_overlapping)

  // Ways mapped inside an area count when deciding whether it is served.
  // Off: only its own outline does, to see area routing on its own.
  bool interior_ways_{true};

  // Areas with more ring and barrier vertices plus connectors than this are
  // declined: exact shortest paths (layer 2) are cubic in them.
  std::size_t max_vertices_{1000U};

  served_params served_{};
  double cells_threshold_{10.0};  // see area_cells_params::threshold_

  // Cells neighbour where their regions share a border; off: one neighbour
  // pair per cut, which always yields a tree.
  bool shared_borders_{true};

  cost_fit cell_fit_{cost_fit::kLeastSquares};  // see cost_fit
  std::size_t min_split_{4U};  // see area_cells_params::min_split_
  std::size_t direct_edges_{0U};  // see area_cells_params::direct_edges_
};

// Where a way stops inside an area rather than on its edge - the end of a
// staircase, a lift, a footway stub - or a lift or entrance node. Entries
// the boundary never sees.
struct loose_end {
  std::int64_t id_{0};
  geo::latlng pos_{};
  area_levels levels_{};
};

// What reading the OSM data yields: the walkable areas (merged), and for
// each the things inside it, plus what is known about the street network's
// nodes.
struct area_data {
  std::vector<walkable_area> areas_;
  merge_stats merged_{};

  // Indexed like areas_.
  std::vector<std::vector<network_edge>> interior_edges_;  // ways mapped inside
  std::vector<std::vector<std::vector<geo::latlng>>> barriers_;  // segments
  std::vector<std::vector<std::vector<geo::latlng>>> buildings_;  // open rings
  std::vector<std::vector<loose_end>> loose_ends_;

  // The areas each node is a ring node of.
  hash_map<std::int64_t, std::vector<std::uint32_t>> ring_areas_;

  // Nodes of the street network: their levels, positions (for those inside
  // an area), and which of them a way without access restrictions reaches.
  hash_map<std::int64_t, area_levels> node_levels_;
  hash_map<std::int64_t, geo::latlng> node_pos_;
  hash_set<std::int64_t> unrestricted_nodes_;

  // Counts, for reporting.
  int n_areas_with_level_{0};
  int n_ways_{0};
  int n_ways_with_level_{0};
  int n_restricted_ways_{0};
  int n_rejected_edges_{0};  // interior edges on another level or layer
};

// Reads `osm` (.osm.pbf; other formats if the program links their reader).
area_data collect_areas(std::filesystem::path const& osm, area_options const&);

enum class area_status : std::uint8_t {
  kTooFewConnectors,  // nothing to cross between
  kTooBig,  // over max_vertices_
  kNoCrossing,  // no two connectors reach each other inside
  kServed,  // the mapped ways suffice
  kNeedsCells,  // not served; cells not asked for
  kMeshed,  // cells built
  kUnreachable  // not served, but no cells could be built
};

char const* to_str(area_status);

struct connector_counts {
  int n_boundary_{0};  // ring nodes the network meets
  int n_restricted_{0};  // ...left out: only access=private/no ways there
  int n_interior_{0};  // stairs, lifts, stubs inside, and dissolved nodes
  int n_rejected_interior_{0};  // ...left out: on another level or layer
};

struct prepared_area {
  area_status status_{area_status::kTooFewConnectors};

  // The rings - outline, then clean holes (inner rings and buildings) - and
  // the barriers. Until the holes are cleaned (kTooFewConnectors) the rings
  // are the area's own.
  std::vector<std::vector<geo::latlng>> rings_;
  std::vector<std::vector<geo::latlng>> barriers_;
  std::size_t n_vertices_{0U};  // ring and barrier vertices

  std::vector<area_connector> connectors_;
  connector_counts counts_{};

  // Layer 1 without geometry.
  std::vector<double> network_;
  bool served_without_geometry_{false};

  // Layer 2, and layer 1 with it.
  std::optional<area_geodesics> geodesics_;
  std::vector<float> geodesic_distances_;  // k x k
  relevance relevance_{};
  served_verdict verdict_{};

  // Layer 3.
  std::optional<area_cells> cells_;
  std::vector<std::vector<area_cut_side>> regions_;
};

// Layers 1 to 3 for area `i`. Stops as soon as the area's fate is known -
// an area the geometry-free shortcut shows served needs no geodesics -
// unless `always_geodesics` (for statistics over every area). Cells are only
// built with `build_cells`.
prepared_area prepare_area(area_data const&,
                           std::size_t i,
                           area_options const&,
                           bool always_geodesics = false,
                           bool build_cells = true);

}  // namespace osr
