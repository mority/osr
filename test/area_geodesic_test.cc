#include "gtest/gtest.h"

#include <cmath>
#include <algorithm>

#include "geo/constants.h"

#include "osr/area/geodesic.h"

using namespace osr;

namespace {

// Synthetic polygons are specified in meters and converted to coordinates
// here, so the expected distances below can be computed by hand.
constexpr auto kLat0 = 49.0;
constexpr auto kLng0 = 8.0;

geo::latlng at(double const x, double const y) {
  return {kLat0 + y / geo::kApproxDistanceLatDegrees,
          kLng0 + x / (std::cos(kLat0 * geo::kPI / 180.0) *
                       geo::kApproxDistanceLatDegrees)};
}

}  // namespace

TEST(area_geodesic, straight_across_convex_area) {
  auto const square = std::vector<std::vector<geo::latlng>>{
      {at(0, 0), at(100, 0), at(100, 100), at(0, 100)}};

  auto const g = area_geodesics{square, {at(0, 0), at(100, 100)}};

  // Nothing in the way: the geodesic is the diagonal, not the two walls.
  EXPECT_NEAR(std::sqrt(2.0) * 100.0, g.distance(0, 1), 0.01);
  EXPECT_EQ(2U, g.path(0, 1).size());
}

TEST(area_geodesic, bends_around_reflex_corner) {
  // L-shaped area: the missing quadrant is x > 40 && y > 40.
  auto const l_shape = std::vector<std::vector<geo::latlng>>{
      {at(0, 0), at(100, 0), at(100, 40), at(40, 40), at(40, 100), at(0, 100)}};

  auto const g = area_geodesics{l_shape, {at(90, 10), at(10, 90)}};

  // The straight line passes through (50, 50), which is outside the L, so the
  // path has to bend at the reflex corner (40, 40) - and only there.
  auto const expected = std::hypot(50.0, 30.0) + std::hypot(30.0, 50.0);
  EXPECT_NEAR(expected, g.distance(0, 1), 0.01);

  auto const path = g.path(0, 1);
  ASSERT_EQ(3U, path.size());
  EXPECT_NEAR(at(40, 40).lat(), path[1].lat(), 1e-7);
  EXPECT_NEAR(at(40, 40).lng(), path[1].lng(), 1e-7);
}

TEST(area_geodesic, routes_around_hole) {
  auto const with_hole = std::vector<std::vector<geo::latlng>>{
      {at(0, 0), at(100, 0), at(100, 100), at(0, 100)},
      {at(40, 40), at(60, 40), at(60, 60), at(40, 60)}};

  auto const g = area_geodesics{with_hole, {at(0, 50), at(100, 50)}};

  // Straight through would cross the hole; the geodesic hugs two of its
  // corners.
  auto const expected = 2.0 * std::hypot(40.0, 10.0) + 20.0;
  EXPECT_NEAR(expected, g.distance(0, 1), 0.01);
  EXPECT_EQ(4U, g.path(0, 1).size());
}

TEST(area_geodesic, symmetric_and_zero_on_diagonal) {
  auto const square = std::vector<std::vector<geo::latlng>>{
      {at(0, 0), at(100, 0), at(100, 100), at(0, 100)}};

  auto const g =
      area_geodesics{square, {at(0, 0), at(100, 0), at(100, 100), at(0, 100)}};

  for (auto i = 0U; i != g.n_connectors(); ++i) {
    EXPECT_FLOAT_EQ(0.F, g.distance(i, i));
    for (auto j = 0U; j != g.n_connectors(); ++j) {
      EXPECT_FLOAT_EQ(g.distance(i, j), g.distance(j, i));
    }
  }
}

TEST(area_geodesic, flat_cost_error_on_square) {
  auto const square = std::vector<std::vector<geo::latlng>>{
      {at(0, 0), at(100, 0), at(100, 100), at(0, 100)}};

  // Four corners: four side pairs at 100m, two diagonal pairs at ~141.42m.
  auto const g =
      area_geodesics{square, {at(0, 0), at(100, 0), at(100, 100), at(0, 100)}};

  auto const mean = (4.0 * 100.0 + 2.0 * std::sqrt(2.0) * 100.0) / 6.0;
  EXPECT_NEAR(mean, g.mean_pair_distance(), 0.01);

  // The flat cost is worst on the pairs furthest from the mean - here the
  // diagonals, which it underestimates.
  EXPECT_NEAR(std::sqrt(2.0) * 100.0 - mean, g.max_pair_error(), 0.01);
}

TEST(area_geodesic, chord_leaving_through_a_vertex_is_not_visible) {
  // The bottom boundary dents inward at (50, 10), so the interior lies above
  // the polyline (0,0)-(50,10)-(100,0) and the chord straight along y = 0 is
  // outside for its whole length.
  //
  // The chord shares an endpoint with each of the two edges it would have to
  // cross, so it never PROPERLY crosses the boundary anywhere - it leaves and
  // re-enters through vertices. Any visibility test built only on segment
  // crossings calls this visible; tg_geom_covers() does exactly that, while
  // separately reporting the chord's own midpoint as outside.
  auto const dented = std::vector<std::vector<geo::latlng>>{
      {at(0, 0), at(50, 10), at(100, 0), at(100, 100), at(0, 100)}};

  EXPECT_FALSE(area_geodesics::is_segment_inside(dented, at(0, 0), at(100, 0)));

  // So the way from corner to corner has to follow the dent, not cut under it.
  auto const g = area_geodesics{dented, {at(0, 0), at(100, 0)}};
  auto const expected = 2.0 * std::hypot(50.0, 10.0);
  EXPECT_NEAR(expected, g.distance(0, 1), 0.01);
  EXPECT_EQ(3U, g.path(0, 1).size());
}

TEST(area_geodesic, fence_across_the_area_has_to_be_walked_around) {
  auto const square = std::vector<std::vector<geo::latlng>>{
      {at(0, 0), at(100, 0), at(100, 100), at(0, 100)}};

  // A fence from the left wall to x = 70, leaving a 30m gap on the right.
  auto const fence =
      std::vector<std::vector<geo::latlng>>{{at(0, 50), at(70, 50)}};

  auto const from = at(10, 10);
  auto const to = at(10, 90);

  auto const open = area_geodesics{square, {from, to}};
  EXPECT_NEAR(80.0, open.distance(0, 1), 0.01);

  auto const blocked = area_geodesics{square, {from, to}, fence};
  EXPECT_FALSE(area_geodesics::is_segment_inside(square, from, to, fence));

  // Round the free end of the fence at (70, 50) and back. The barrier is
  // given a width so that it can separate its two sides at all, so the answer
  // sits just outside the zero-width ideal rather than on it.
  auto const expected = 2.0 * std::hypot(60.0, 40.0);
  EXPECT_NEAR(expected, blocked.distance(0, 1), 1.0);
  EXPECT_GT(blocked.distance(0, 1), expected);

  // ... and the turn happens at the fence's free end, not at its anchor.
  auto const path = blocked.path(0, 1);
  ASSERT_GE(path.size(), 3U);
  auto const near_free_end = std::ranges::any_of(
      path,
      [&](geo::latlng const& c) { return geo::distance(c, at(70, 50)) < 1.0; });
  auto const near_anchor = std::ranges::any_of(path, [&](geo::latlng const& c) {
    return geo::distance(c, at(0, 50)) < 1.0;
  });
  EXPECT_TRUE(near_free_end);
  EXPECT_FALSE(near_anchor);
}

TEST(area_geodesic, fence_does_not_change_what_counts_as_inside) {
  auto const square = std::vector<std::vector<geo::latlng>>{
      {at(0, 0), at(100, 0), at(100, 100), at(0, 100)}};

  // An obstacle is not an outline: a fence running past the area must not
  // make the ground beyond it count as part of the area.
  auto const fence =
      std::vector<std::vector<geo::latlng>>{{at(-50, 150), at(150, 150)}};

  auto const g = area_geodesics{square, {at(50, 50), at(50, 200)}, fence};
  EXPECT_EQ(area_geodesics::kUnreachable, g.distance(0, 1));
}

TEST(area_geodesic, connector_outside_stays_unreachable) {
  auto const square = std::vector<std::vector<geo::latlng>>{
      {at(0, 0), at(100, 0), at(100, 100), at(0, 100)}};

  auto const g = area_geodesics{square, {at(50, 50), at(200, 200)}};

  EXPECT_EQ(area_geodesics::kUnreachable, g.distance(0, 1));
  EXPECT_TRUE(g.path(0, 1).empty());
  EXPECT_FALSE(g.is_connector_reachable(1));
}
