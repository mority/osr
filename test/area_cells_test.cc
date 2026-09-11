#include "gtest/gtest.h"

#include <cmath>
#include <cstdint>
#include <algorithm>
#include <numeric>
#include <random>
#include <vector>

#include "geo/constants.h"

#include "osr/area/cells.h"
#include "osr/area/geodesic.h"

using namespace osr;

namespace {

constexpr auto kLat0 = 49.0;
constexpr auto kLng0 = 8.0;

geo::latlng at(double const x, double const y) {
  return {kLat0 + y / geo::kApproxDistanceLatDegrees,
          kLng0 + x / (std::cos(kLat0 * geo::kPI / 180.0) *
                       geo::kApproxDistanceLatDegrees)};
}

// A convex, empty area: every geodesic is the straight line.
std::vector<float> straight_lines(std::vector<geo::latlng> const& pts) {
  auto const k = pts.size();
  auto d = std::vector<float>(k * k, 0.F);
  for (auto i = std::size_t{0U}; i != k; ++i) {
    for (auto j = std::size_t{0U}; j != k; ++j) {
      d[i * k + j] = static_cast<float>(geo::distance(pts[i], pts[j]));
    }
  }
  return d;
}

std::vector<bool> all_relevant(std::size_t const k) {
  return std::vector<bool>(k * k, true);
}

std::size_t n_neighbour_pairs(area_cells const& c) {
  auto n = std::size_t{0U};
  for (auto a = 0U; a != c.n_cells(); ++a) {
    for (auto b = a + 1U; b != c.n_cells(); ++b) {
      n += c.is_neighbour(static_cast<area_cells::cell_idx_t>(a),
                          static_cast<area_cells::cell_idx_t>(b))
               ? 1U
               : 0U;
    }
  }
  return n;
}

}  // namespace

TEST(area_cells, one_cell_when_within_threshold) {
  auto const pts =
      std::vector{at(0, 0), at(30, 0), at(30, 30), at(0, 30), at(15, 15)};
  auto const d = straight_lines(pts);

  auto const c = build_area_cells(pts, d, all_relevant(pts.size()),
                                  {.threshold_ = 1000.0});
  ASSERT_TRUE(c.has_value());
  ASSERT_EQ(1U, c->n_cells());
  EXPECT_TRUE(c->neighbours_.empty());

  // One cell is one unknown: least squares makes it the mean distance.
  auto sum = 0.0;
  auto n = 0;
  for (auto i = 0U; i != pts.size(); ++i) {
    for (auto j = i + 1U; j != pts.size(); ++j) {
      sum += d[i * pts.size() + j];
      ++n;
    }
  }
  EXPECT_NEAR(sum / n, c->cost_[0], 0.01);
}

TEST(area_cells, splits_two_distant_groups) {
  // Two 1m pairs 100m apart. A single flat cost cannot be both ~1m and ~100m,
  // so the cell has to split - and the best cut separates the groups.
  auto const pts = std::vector{at(0, 0), at(0, 1), at(100, 0), at(100, 1)};
  auto const d = straight_lines(pts);

  auto const c =
      build_area_cells(pts, d, all_relevant(pts.size()), {.threshold_ = 10.0});
  ASSERT_TRUE(c.has_value());
  ASSERT_EQ(2U, c->n_cells());
  EXPECT_EQ(c->connector_cell_[0], c->connector_cell_[1]);
  EXPECT_EQ(c->connector_cell_[2], c->connector_cell_[3]);
  EXPECT_NE(c->connector_cell_[0], c->connector_cell_[2]);
  EXPECT_TRUE(c->is_neighbour(0, 1));
  EXPECT_TRUE(c->is_neighbour(1, 0));

  // Going from one cell into the other pays both cells' costs.
  auto const cd = c->cell_distances();
  EXPECT_FLOAT_EQ(c->cost_[0] + c->cost_[1], cd[0 * 2 + 1]);
  EXPECT_FLOAT_EQ(c->cost_[0], cd[0]);
}

TEST(area_cells, cells_sharing_a_border_are_neighbours) {
  // Four clusters in the corners of a square area: two cuts make four cells.
  // The first cut borders two cells on each side. One pair per cut would give
  // three neighbour pairs; every two cells along a side share a border.
  auto const pts = std::vector{at(0, 0),    at(2, 0),    at(0, 2),
                               at(100, 0),  at(98, 0),   at(100, 2),
                               at(0, 100),  at(2, 100),  at(0, 98),
                               at(100, 100), at(98, 100), at(100, 98)};
  auto const square = std::vector<std::vector<geo::latlng>>{
      {at(-10, -10), at(110, -10), at(110, 110), at(-10, 110)}};
  auto const c = build_area_cells(pts, straight_lines(pts),
                                  all_relevant(pts.size()),
                                  {.threshold_ = 5.0, .rings_ = &square});
  ASSERT_TRUE(c.has_value());
  ASSERT_EQ(4U, c->n_cells());
  auto const cell = [&](std::size_t const i) { return c->connector_cell_[i]; };
  auto const sw = cell(0), se = cell(3), nw = cell(6), ne = cell(9);
  EXPECT_TRUE(c->is_neighbour(sw, se));
  EXPECT_TRUE(c->is_neighbour(nw, ne));
  EXPECT_TRUE(c->is_neighbour(sw, nw));
  EXPECT_TRUE(c->is_neighbour(se, ne));
}

namespace {

std::vector<geo::latlng> rect(double const x0,
                              double const y0,
                              double const x1,
                              double const y1) {
  return {at(x0, y0), at(x1, y0), at(x1, y1), at(x0, y1)};
}

// Cuts along x = x0 and y = y0, recorded as the cut tree records them.
area_cut_side cut_x(double const x0, bool const below) {
  return {.lng_scale_ = 1.0,
          .cx_ = 1.0,
          .cy_ = 0.0,
          .threshold_ = at(x0, 0).lng(),
          .below_ = below};
}

area_cut_side cut_y(double const y0, bool const below) {
  return {.lng_scale_ = 1.0,
          .cx_ = 0.0,
          .cy_ = 1.0,
          .threshold_ = at(0, y0).lat(),
          .below_ = below};
}

double area_m2(std::vector<geo::latlng> const& ring) {
  auto a = 0.0;
  auto const k = std::cos(kLat0 * geo::kPI / 180.0);
  for (auto i = std::size_t{0U}; i != ring.size(); ++i) {
    auto const& p = ring[i];
    auto const& q = ring[(i + 1U) % ring.size()];
    a += p.lng() * k * q.lat() - q.lng() * k * p.lat();
  }
  return std::abs(a) / 2.0 * geo::kApproxDistanceLatDegrees *
         geo::kApproxDistanceLatDegrees;
}

using pairs = std::vector<std::pair<std::size_t, std::size_t>>;

}  // namespace

TEST(area_cells, quadrants_neighbour_along_edges_not_at_the_corner) {
  // Cut at x = 50, then each half at y = 50.
  auto const quadrants = std::vector<std::vector<area_cut_side>>{
      {cut_x(50, true), cut_y(50, true)},
      {cut_x(50, true), cut_y(50, false)},
      {cut_x(50, false), cut_y(50, true)},
      {cut_x(50, false), cut_y(50, false)}};
  EXPECT_EQ((pairs{{0, 1}, {0, 2}, {1, 3}, {2, 3}}),
            cells_sharing_a_border({rect(0, 0, 100, 100)}, quadrants));
}

TEST(area_cells, a_building_on_the_border_keeps_cells_apart) {
  // Two halves, and a building over the line between them. With gaps of
  // 30 cm at either end nobody walks from one to the other there; with a
  // 10 m passage one does.
  auto const halves = std::vector<std::vector<area_cut_side>>{
      {cut_x(50, true)}, {cut_x(50, false)}};
  auto const outline = rect(0, 0, 100, 100);
  EXPECT_EQ((pairs{{0, 1}}), cells_sharing_a_border({outline}, halves));
  EXPECT_EQ((pairs{}),
            cells_sharing_a_border({outline, rect(45, 0.3, 55, 99.7)}, halves));
  EXPECT_EQ((pairs{{0, 1}}),
            cells_sharing_a_border({outline, rect(45, 10, 55, 99.7)}, halves));
}

TEST(area_cells, cells_across_a_notch_are_not_neighbours) {
  // A U-shaped area: base below y = 30, two arms above it. The line between
  // the arms runs through the gap between them, outside the area.
  auto const u = std::vector{at(0, 0),   at(100, 0),  at(100, 100), at(60, 100),
                             at(60, 30), at(40, 30),  at(40, 100),  at(0, 100)};
  auto const cells = std::vector<std::vector<area_cut_side>>{
      {cut_y(30, true)},                   // base
      {cut_y(30, false), cut_x(50, true)},  // left arm
      {cut_y(30, false), cut_x(50, false)}};  // right arm
  EXPECT_EQ((pairs{{0, 1}, {0, 2}}), cells_sharing_a_border({u}, cells));
}

TEST(area_cells, regions_cover_the_area_exactly) {
  // A square with a building in it, subdivided for real: the regions must
  // add up to the square without the building, no more, no less.
  auto const rings = std::vector<std::vector<geo::latlng>>{
      rect(0, 0, 200, 200), rect(80, 80, 120, 120)};
  auto rng = std::mt19937{7U};
  auto coord = std::uniform_real_distribution<double>{0.0, 200.0};
  auto pts = std::vector<geo::latlng>{};
  while (pts.size() != 30U) {
    auto const x = coord(rng);
    auto const y = coord(rng);
    if (x < 75.0 || x > 125.0 || y < 75.0 || y > 125.0) {
      pts.push_back(at(x, y));
    }
  }
  auto regions = std::vector<std::vector<area_cut_side>>{};
  auto const c = build_area_cells(pts, straight_lines(pts),
                                  all_relevant(pts.size()),
                                  {.threshold_ = 5.0, .rings_ = &rings},
                                  &regions);
  ASSERT_TRUE(c.has_value());
  ASSERT_GT(c->n_cells(), 2U);

  auto total = 0.0;
  for (auto const& cell : region_rings(rings, regions)) {
    for (auto i = std::size_t{0U}; i != cell.size(); ++i) {
      total += (i == 0U ? 1.0 : -1.0) * area_m2(cell[i]);
    }
  }
  EXPECT_NEAR(200.0 * 200.0 - 40.0 * 40.0, total, 1.0);
}

TEST(area_cells, with_a_cycle_the_cheaper_way_round_is_taken) {
  // Four cells in a ring, one of them expensive to cross: from cell 0 to
  // cell 2, round the cheap side.
  auto c = area_cells{};
  c.cost_ = {1.F, 100.F, 1.F, 1.F};
  c.connector_cell_ = {0, 1, 2, 3};
  // Neighbours 0-1, 1-2, 2-3, 0-3, in the packed upper triangle:
  // (0,1)=0 (0,2)=1 (0,3)=2 (1,2)=3 (1,3)=4 (2,3)=5.
  c.neighbours_ = {(1U << 0U) | (1U << 2U) | (1U << 3U) | (1U << 5U)};
  EXPECT_TRUE(c.is_neighbour(0, 3));
  EXPECT_FALSE(c.is_neighbour(0, 2));
  auto const d = c.cell_distances();
  EXPECT_FLOAT_EQ(3.F, d[0 * 4 + 2]);  // 0 -> 3 -> 2, not through 1
}

TEST(area_cells, contradictory_reachability_yields_nothing) {
  // 0 and 3 both reach 1 and 2, so inside a polygon they reach each other.
  auto const pts = std::vector{at(0, 0), at(10, 0), at(20, 0), at(30, 0)};
  auto d = straight_lines(pts);
  d[0 * 4 + 3] = d[3 * 4 + 0] = area_geodesics::kUnreachable;

  EXPECT_FALSE(
      build_area_cells(pts, d, all_relevant(pts.size())).has_value());
}

TEST(area_cells, stray_connector_gets_no_cell) {
  // Stairs ending inside a building: reaches none of the others.
  auto const pts =
      std::vector{at(0, 0), at(30, 0), at(30, 30), at(0, 30), at(15, 15)};
  auto d = straight_lines(pts);
  for (auto i = 0U; i != 4U; ++i) {
    d[i * 5 + 4] = d[4 * 5 + i] = area_geodesics::kUnreachable;
  }

  auto const c = build_area_cells(pts, d, all_relevant(pts.size()),
                                  {.threshold_ = 1e6});
  ASSERT_TRUE(c.has_value());
  ASSERT_EQ(1U, c->n_cells());
  EXPECT_EQ(area_cells::kNoCell, c->connector_cell_[4]);
  for (auto i = 0U; i != 4U; ++i) {
    EXPECT_EQ(0U, c->connector_cell_[i]);
  }
  EXPECT_EQ(1U, hub_positions(*c, pts).size());
}

TEST(area_cells, separate_parts_are_never_neighbours) {
  // Two groups a building cuts apart: each is meshed, none of their cells
  // neighbour the other's, and no cell path joins them.
  auto const pts = std::vector{at(0, 0),   at(0, 1),   at(40, 0), at(40, 1),
                               at(200, 0), at(200, 1), at(240, 0), at(240, 1)};
  auto d = straight_lines(pts);
  for (auto i = 0U; i != 4U; ++i) {
    for (auto j = 4U; j != 8U; ++j) {
      d[i * 8 + j] = d[j * 8 + i] = area_geodesics::kUnreachable;
    }
  }

  auto const c =
      build_area_cells(pts, d, all_relevant(pts.size()), {.threshold_ = 5.0});
  ASSERT_TRUE(c.has_value());
  auto const m = c->n_cells();
  auto const cd = c->cell_distances();
  for (auto i = 0U; i != 4U; ++i) {
    for (auto j = 4U; j != 8U; ++j) {
      auto const a = c->connector_cell_[i];
      auto const b = c->connector_cell_[j];
      ASSERT_LT(a, m);
      ASSERT_LT(b, m);
      EXPECT_NE(a, b);
      EXPECT_FALSE(c->is_neighbour(a, b));
      EXPECT_FALSE(std::isfinite(cd[a * m + b]));
    }
  }
}

TEST(area_cells, hub_is_centre_of_smallest_enclosing_circle) {
  // One cell: huge threshold. The square's centre is its circle's centre,
  // and the interior point must not pull it anywhere.
  auto const pts =
      std::vector{at(0, 0), at(40, 0), at(40, 40), at(0, 40), at(5, 5)};
  auto const c = build_area_cells(pts, straight_lines(pts),
                                  all_relevant(pts.size()),
                                  {.threshold_ = 1e6});
  ASSERT_TRUE(c.has_value());
  ASSERT_EQ(1U, c->n_cells());
  auto const hubs = hub_positions(*c, pts);
  ASSERT_EQ(1U, hubs.size());
  EXPECT_NEAR(0.0, geo::distance(at(20, 20), hubs[0]), 0.01);
}

TEST(area_cells, hub_of_collinear_connectors_is_midpoint_of_extremes) {
  auto const pts = std::vector{at(0, 0), at(10, 0), at(70, 0), at(30, 0)};
  auto const c = build_area_cells(pts, straight_lines(pts),
                                  all_relevant(pts.size()),
                                  {.threshold_ = 1e6});
  ASSERT_TRUE(c.has_value());
  ASSERT_EQ(1U, c->n_cells());
  EXPECT_NEAR(0.0, geo::distance(at(35, 0), hub_positions(*c, pts)[0]), 0.01);
}

TEST(area_cells, neighbours_form_a_tree_over_every_connector) {
  auto rng = std::mt19937{42U};
  auto coord = std::uniform_real_distribution<double>{0.0, 200.0};
  for (auto round = 0; round != 20; ++round) {
    auto pts = std::vector<geo::latlng>(40);
    for (auto& p : pts) {
      p = at(coord(rng), coord(rng));
    }
    auto const d = straight_lines(pts);

    auto regions = std::vector<std::vector<area_cut_side>>{};
    auto const c = build_area_cells(pts, d, all_relevant(pts.size()),
                                    {.threshold_ = 5.0}, &regions);
    ASSERT_TRUE(c.has_value());
    auto const m = c->n_cells();
    ASSERT_GE(m, 1U);
    ASSERT_LE(m, area_cells::kMaxCells);

    // Every connector lies in its own cell's region, and in no other.
    ASSERT_EQ(m, regions.size());
    for (auto i = std::size_t{0U}; i != pts.size(); ++i) {
      for (auto cell = std::size_t{0U}; cell != m; ++cell) {
        auto const inside = std::ranges::all_of(
            regions[cell],
            [&](area_cut_side const& s) { return s.keeps(pts[i]); });
        EXPECT_EQ(cell == c->connector_cell_[i], inside)
            << "connector " << i << " cell " << cell;
      }
    }

    // Every connector sits in a cell, and every cell has a connector.
    auto has_connector = std::vector<bool>(m, false);
    for (auto const cell : c->connector_cell_) {
      ASSERT_LT(cell, m);
      has_connector[cell] = true;
    }
    EXPECT_TRUE(std::ranges::all_of(has_connector, [](bool b) { return b; }));

    for (auto const cost : c->cost_) {
      EXPECT_GE(cost, 0.F);
    }

    // m - 1 neighbour pairs and connected: a spanning tree.
    EXPECT_EQ(m - 1U, n_neighbour_pairs(*c));
    auto const cd = c->cell_distances();
    for (auto a = std::size_t{0U}; a != m; ++a) {
      EXPECT_FLOAT_EQ(c->cost_[a], cd[a * m + a]);
      for (auto b = std::size_t{0U}; b != m; ++b) {
        EXPECT_TRUE(std::isfinite(cd[a * m + b]));
        EXPECT_FLOAT_EQ(cd[a * m + b], cd[b * m + a]);
      }
    }
  }
}
