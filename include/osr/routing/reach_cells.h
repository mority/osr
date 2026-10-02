#pragma once

#include <cstdint>
#include <limits>
#include <optional>
#include <span>
#include <vector>

#include "geo/latlng.h"

#include "osr/location.h"
#include "osr/routing/parameters.h"
#include "osr/routing/profile.h"
#include "osr/types.h"

namespace osr {

struct ways;
struct lookup;

// Network Voronoi diagram: one forward search from all sources at once. Every
// node belongs to the source that settles it first (= cheapest source). The
// owned nodes are then rasterized and traced to one polygon set per source.

using source_idx_t = std::uint16_t;
constexpr auto const kNoSource = std::numeric_limits<source_idx_t>::max();

struct reach_node {
  node_idx_t node_;
  cost_t cost_;
  duration_t duration_;
  source_idx_t source_;
};

struct reach_cells_params {
  cost_t max_{7200U};
  double cell_size_m_{500.0};
  unsigned fill_steps_{4U};  // dilation rounds at the outer border
  unsigned smooth_iterations_{3U};  // Chaikin corner cutting, 0 = off
  duration_t band_size_{1800U};  // drive time bands, 0 = no bands

  // Meeting zones (nullopt = off, either needs one extra search per source):
  // - latest: all places where the latest arrival over all sources is at most
  //   `meet_margin_` later than at the place where it is earliest
  // - total: all places where the sum of all drive times is at most
  //   `total_margin_` more than at the place where it is smallest
  std::optional<duration_t> meet_margin_{};
  std::optional<duration_t> total_margin_{};
};

struct reach_grid {
  std::size_t idx(std::int32_t const x, std::int32_t const y) const {
    return static_cast<std::size_t>(y) * nx_ + static_cast<std::size_t>(x);
  }

  bool contains(std::int32_t const x, std::int32_t const y) const {
    return x >= 0 && y >= 0 && x < static_cast<std::int32_t>(nx_) &&
           y < static_cast<std::int32_t>(ny_);
  }

  source_idx_t owner(std::int32_t const x, std::int32_t const y) const {
    return contains(x, y) ? owner_[idx(x, y)] : kNoSource;
  }

  // Grid vertex (x, y) = south-west corner of cell (x, y).
  geo::latlng vertex_pos(std::int32_t const x, std::int32_t const y) const {
    return {origin_.lat_ + y * dlat_, origin_.lng_ + x * dlng_};
  }

  geo::latlng origin_{};  // south-west corner of the grid
  double dlat_{0.0};
  double dlng_{0.0};
  std::uint32_t nx_{0U};
  std::uint32_t ny_{0U};
  std::vector<source_idx_t> owner_;
  std::vector<cost_t> cost_;  // decides the owner
  std::vector<duration_t> duration_;  // drive time of the owner

  // Latest arrival / total drive time over all sources (min over the cell's
  // nodes, seconds), only set when a meeting zone is computed.
  std::vector<std::uint32_t> meet_latest_;
  std::vector<std::uint32_t> meet_total_;
};

struct grid_vertex {
  friend bool operator==(grid_vertex, grid_vertex) = default;
  std::int32_t x_, y_;
  bool pinned_{false};  // junction of 3+ areas: kept in place by smoothing
};

// Rings are not closed (last != first). Ring 0 is the outer ring
// (counter-clockwise), the others are holes (clockwise).
using grid_polygon = std::vector<std::vector<grid_vertex>>;

// GeoJSON conventions: rings are closed, outer ring counter-clockwise.
using reach_polygon = std::vector<std::vector<geo::latlng>>;

// Settled nodes with their cost and the source that reached them first.
std::vector<reach_node> compute_reach_owners(
    profile_parameters const&,
    ways const&,
    lookup const&,
    search_profile,
    std::vector<location> const& sources,
    cost_t max,
    double max_match_distance,
    bitvec<node_idx_t> const* blocked = nullptr);

enum class meet_criterion : std::uint8_t {
  kLatest,  // drive time of the last to arrive
  kTotal  // sum of all drive times
};

// A node every source reaches (times in seconds).
struct meet_node {
  std::uint32_t get(meet_criterion const c) const {
    return c == meet_criterion::kLatest ? latest_ : total_;
  }

  node_idx_t node_;
  std::uint32_t latest_;
  std::uint32_t total_;
};

struct meet_times {
  struct best {
    meet_node node_;
    std::vector<duration_t> per_source_;
  };

  std::vector<meet_node> nodes_;

  // Node where the latest arrival is earliest / the total is smallest.
  std::optional<best> latest_;
  std::optional<best> total_;
};

// Zone margins the per-source searches have to cover. Without any margin,
// every search runs up to `max`.
struct meet_bounds {
  std::optional<duration_t> latest_margin_{};
  std::optional<duration_t> total_margin_{};
};

// Searches run until the best value found so far + margin + this, in cost:
// cost is drive time plus penalties (U-turn, private gates).
constexpr auto const kMeetSearchSlack = cost_t{300U};

// One search per source, in parallel and in lockstep. With bounds, the
// searches stop once every place within a zone has been reached by all
// sources: `nodes_` then only contains the places reached by then (at least
// everything within the zones).
meet_times compute_meet_times(profile_parameters const&,
                              ways const&,
                              lookup const&,
                              search_profile,
                              std::vector<location> const& sources,
                              cost_t max,
                              double max_match_distance,
                              bitvec<node_idx_t> const* blocked = nullptr,
                              meet_bounds const& = {});

// Each cell gets the owner of its cheapest node, empty cells are filled from
// their neighbours `fill_steps` times, enclosed empty areas are filled
// completely (`close_holes`). `meet` (optional) sets
// `reach_grid::meet_latest_/meet_total_`, its nodes must be a subset of
// `nodes`.
reach_grid make_reach_grid(ways const&,
                           std::span<reach_node const> nodes,
                           std::span<meet_node const> meet,
                           double cell_size_m,
                           unsigned fill_steps);

void fill_gaps(reach_grid&, unsigned steps);

// Fills empty cells that are not connected to the grid border from the
// nearest owned cell.
void close_holes(reach_grid&);

// result[source] = polygons of all cells owned by `source`.
std::vector<std::vector<grid_polygon>> trace_reach_grid(reach_grid const&,
                                                        std::size_t n_sources);

// result[source * n_bands + band] = polygons of the cells owned by `source`
// with drive time in [band * band_size, (band + 1) * band_size). The last band
// is open-ended.
std::vector<std::vector<grid_polygon>> trace_reach_bands(reach_grid const&,
                                                         std::size_t n_sources,
                                                         duration_t band_size,
                                                         std::size_t n_bands);

// Polygons of the cells whose latest arrival / total drive time is <= limit.
std::vector<grid_polygon> trace_meet_zone(reach_grid const&,
                                          meet_criterion,
                                          std::uint32_t limit);

// Optionally smoothed (Chaikin); pinned vertices stay in place.
reach_polygon to_latlng(reach_grid const&,
                        grid_polygon const&,
                        unsigned smooth_iterations = 0U);

struct reach_band {
  source_idx_t source_;
  std::uint16_t
      band_;  // drive time in [band * band_size, (band+1) * band_size)
  reach_polygon polygon_;
};

struct reach_cells {
  // sources_[source] = polygons of the area reached first by `source`
  std::vector<std::vector<reach_polygon>> sources_;

  // the same areas, split into drive time bands (empty if band_size == 0)
  std::vector<reach_band> bands_;

  // only if requested and some place is reached by all sources
  struct meet {
    geo::latlng pos_;
    std::uint32_t latest_;  // latest arrival at pos_ [s]
    std::uint32_t total_;  // sum of all drive times to pos_ [s]
    std::vector<duration_t> per_source_;
    std::vector<reach_polygon> zone_;
  };
  std::optional<meet> meet_;  // best latest arrival
  std::optional<meet> total_meet_;  // best total drive time
};

reach_cells compute_reach_cells(profile_parameters const&,
                                ways const&,
                                lookup const&,
                                search_profile,
                                std::vector<location> const& sources,
                                reach_cells_params const&,
                                double max_match_distance,
                                bitvec<node_idx_t> const* blocked = nullptr);

}  // namespace osr
