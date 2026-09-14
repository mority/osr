#pragma once

#include <cstdint>
#include <limits>
#include <optional>
#include <utility>
#include <vector>

#include "geo/latlng.h"

namespace osr {

// Layer 3 of area routing (the layers are listed in osr/area/walkable.h): the
// subdivision.
//   in:  an area's connectors, their exact shortest distances (layer 2),
//        which pairs a crossing shortens enough to matter (layer 1), and its
//        rings with clean holes (layer 0)
//   out: area_cells - which cell each connector is in, what crossing a cell
//        costs, which cells neighbour - for routing (layer 4). Cell regions
//        (area_cut_side, region_rings) are for drawing only.

// The routing model for crossing one area: a handful of cells, each crossed
// at a flat cost.
//
// Every cell becomes a hub the router can pass through. A connector is joined
// to the hub of its cell, and two neighbouring hubs are joined to each other.
// Crossing cell `a` costs `c_a`, put half onto each edge touching its hub:
// connector -> hub is `c_a / 2`, hub `a` -> hub `b` is `(c_a + c_b) / 2`. A
// route through cells a, b, c therefore costs `c_a + c_b + c_c`, wherever in
// each cell it enters and leaves.
//
// That is all there is. Where two cells meet is not stored: a point on their
// border would only ever connect those two hubs, so it is the same thing as a
// direct hub edge, and the cell costs already fix what that edge costs. Which
// cells meet is a bit matrix. Hub positions are not stored either - they are a
// function of the members' positions and are derived when the data is loaded.
struct area_cells {
  using cell_idx_t = std::uint8_t;

  // For a connector that cannot reach any other connector inside the area -
  // typically stairs or a stub ending inside a building or behind a fence.
  // There is nothing for it to cross to, so it belongs to no cell.
  static constexpr auto kNoCell = std::numeric_limits<cell_idx_t>::max();

  // Each cut halves a cell, so the depth bounds the cell count at 2^depth.
  static constexpr auto kMaxDepth = std::size_t{6U};
  static constexpr auto kMaxCells = std::size_t{1U} << kMaxDepth;

  std::size_t n_cells() const noexcept { return cost_.size(); }
  std::size_t n_connectors() const noexcept { return connector_cell_.size(); }

  bool is_neighbour(cell_idx_t a, cell_idx_t b) const;

  // For rebuilding a stored subdivision: cost_ has to be set first.
  void set_neighbour(cell_idx_t a, cell_idx_t b);
  void set_link(cell_idx_t a, cell_idx_t b, float cost);

  // What the spoke from connector `i` to its hub costs, and the link between
  // two neighbouring hubs. From spoke_ and link_ where they are set, else from
  // the flat costs: half the cell's cost, half the sum of the two cells'.
  float spoke_cost(std::size_t i) const;
  float link_cost(cell_idx_t a, cell_idx_t b) const;

  // n_connectors() x n_connectors(), row-major: the modelled cost of crossing
  // from connector i to connector j - its spoke, the cheapest hub path, the
  // other spoke, or any cheaper way through the direct edges. Zero on the
  // diagonal; infinite where no path exists or a connector has no cell. This
  // is the model, whichever costs it has.
  std::vector<float> connector_distances() const;

  // n_cells() x n_cells(), row-major, for the flat costs only: the modelled
  // cost of entering the area in cell `a` and leaving it in cell `b` - the sum
  // of the costs of the cells on the cheapest cell path, both ends included.
  // The diagonal is each cell's own cost. Infinite where no cell path exists.
  std::vector<float> cell_distances() const;

  // The cell each connector belongs to, or kNoCell.
  std::vector<cell_idx_t> connector_cell_;

  // Flat crossing cost of each cell, in meters. With spoke_ set, only a
  // summary for drawing: twice the mean spoke of the cell's connectors.
  std::vector<float> cost_;

  // Per-connector costs (cost_fit::kSpokes); both empty for the flat model.
  // Each connector's spoke, and each link indexed like the neighbour bits.
  std::vector<float> spoke_;
  std::vector<float> link_;

  // A pair of connectors the cells price badly, kept exactly: the walk
  // between them at its true distance. Since that is a walk that exists, it
  // can only bring a modelled crossing down to the truth, never below it.
  struct direct_edge {
    std::uint32_t a_{0U}, b_{0U};
    float cost_{0.F};
  };
  std::vector<direct_edge> direct_;

  // The strict upper triangle (a < b) of the symmetric neighbour matrix,
  // row-major, packed into bits: n_cells() * (n_cells() - 1) / 2 of them.
  std::vector<std::uint64_t> neighbours_;
};

// One side of a cut, for drawing cells - build-time only; the router never
// needs it. A cut is a straight line: it keeps the points p with
//   key(p) = p.lng() * lng_scale_ * cx_ + p.lat() * cy_
// below threshold_ (below_ == true) or at or above it (below_ == false). A
// cell's region is the area polygon on the kept side of every cut on its
// path through the cut tree. In a concave area that can be several pieces,
// and the regions of separate walkable parts overlap.
struct area_cut_side {
  bool keeps(geo::latlng const& p) const {
    auto const key = p.lng() * lng_scale_ * cx_ + p.lat() * cy_;
    return below_ ? key < threshold_ : key >= threshold_;
  }

  double lng_scale_{0.0}, cx_{0.0}, cy_{0.0}, threshold_{0.0};
  bool below_{true};
};

// How the cell costs are fitted to the geodesics.
enum class cost_fit : std::uint8_t {
  // Least squares: small errors either way. The router picks the crossings
  // priced too low, so routes come out optimistic.
  kLeastSquares,

  // No crossing priced below its geodesic, with as little overpricing along
  // the routes taken as that allows: a linear program over the cell costs,
  // solved with HiGHS.
  kNoUnderpricing,

  // As kNoUnderpricing, but with a cost per connector - its spoke to the hub
  // - and per link between neighbouring hubs instead of one per cell, so a
  // short crossing does not pay for the longest one in its cell.
  kSpokes,

  // As kSpokes, but minimising first the worst overpricing of any pair in the
  // area, and only then, holding that, the overpricing overall.
  kSpokesWorstCase
};

struct area_cells_params {
  // A cell is split while some route through it misses the geodesic by more
  // than this many meters - so it bounds an error the router sees, not a
  // statistic of the cell.
  double threshold_{10.0};

  // Depth of the cut tree, at most area_cells::kMaxDepth.
  std::size_t max_depth_{area_cells::kMaxDepth};

  // The area's outline and holes (rings_[0] is the outline). With them, two
  // cells are neighbours wherever their regions share a border inside the
  // area. Without them, each cut makes just one pair of cells neighbours,
  // which always yields a tree - cells along the rest of a cut's line are
  // then not neighbours although they border each other.
  std::vector<std::vector<geo::latlng>> const* rings_{nullptr};

  cost_fit fit_{cost_fit::kLeastSquares};

  // Cells with fewer connectors than this are not split further. At 2, the
  // cut tree may go down to one connector per cell.
  std::size_t min_split_{4U};

  // How many pairs an area may keep exactly, as direct edges between two
  // connectors (see area_cells::direct_): the ones its cells price worst.
  // A spoke has to cover its connector's most distant partner, so where one
  // far connector exists every nearby pair inherits that; a handful of these
  // edges per area is what the tail costs to remove.
  std::size_t direct_edges_{0U};
};

// Subdivides an area and fits the cell costs.
//
// `dist` is the connectors' geodesic distance matrix, k x k and row-major, as
// computed by area_geodesics (kUnreachable where no path exists). It is the
// only thing about the area's geometry this reads: the cuts, the fit and the
// refinement are all arithmetic on it.
//
// `relevant` (k x k) marks the pairs a crossing actually shortens compared to
// walking round the edge. Where to cut is chosen on those pairs only.
//
// An area can fall apart into walkable parts that cannot reach each other: a
// building touching both sides, or a fence. Each part with two or more
// connectors is subdivided on its own; parts are never neighbours, so the
// neighbour relation is a forest. Connectors that reach nothing get kNoCell.
//
// Returns nullopt if no two connectors reach each other, or if the matrix
// contradicts itself (i reaches j and j reaches l, but i does not reach l),
// which no polygon can produce.
//
// If `cell_regions` is given, it receives for every cell the cut sides that
// bound its region (see area_cut_side), indexed like the cells.
std::optional<area_cells> build_area_cells(
    std::vector<geo::latlng> const& connectors,
    std::vector<float> const& dist,
    std::vector<bool> const& relevant,
    area_cells_params const& = {},
    std::vector<std::vector<area_cut_side>>* cell_regions = nullptr);

// Where each cell's hub sits: the centre of the smallest circle enclosing the
// cell's connectors. That keeps the longest spoke as short as possible, and
// the longest spoke is what an admissible spoke cost has to cover. Derived
// from the connectors' positions, never stored.
std::vector<geo::latlng> hub_positions(
    area_cells const&, std::vector<geo::latlng> const& connectors);

// Which cells border each other where one can walk: the pairs (a < b) whose
// regions share at least 1 m of border inside `rings` (outline first, then
// holes). A border running through a building, or across a notch outside a
// concave outline, does not count. `cell_regions` are the cells' paths
// through a cut tree, as build_area_cells reports them.
std::vector<std::pair<std::size_t, std::size_t>> cells_sharing_a_border(
    std::vector<std::vector<geo::latlng>> const& rings,
    std::vector<std::vector<area_cut_side>> const& cell_regions);

// Each cell's region as polygon rings, for drawing: every ring of the area -
// outline first, then holes - clipped to the kept side of each of the cell's
// cuts. The first ring is the piece of the outline; empty for a cell of which
// nothing of the outline is left. Rings are not closed.
std::vector<std::vector<std::vector<geo::latlng>>> region_rings(
    std::vector<std::vector<geo::latlng>> const& rings,
    std::vector<std::vector<area_cut_side>> const& cell_regions);

}  // namespace osr
