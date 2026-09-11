#include "gtest/gtest.h"

#include <cmath>
#include <vector>

#include "geo/constants.h"

#include "osr/area/geodesic.h"
#include "osr/area/served.h"

using namespace osr;

// Layer 1 of area routing: does an area need cells, or do the ways already
// mapped serve it? Each test is a small version of a case from Berlin.

namespace {

constexpr auto kLat0 = 52.52;
constexpr auto kLng0 = 13.41;

geo::latlng at(double const x, double const y) {
  return {kLat0 + y / geo::kApproxDistanceLatDegrees,
          kLng0 + x / (std::cos(kLat0 * geo::kPI / 180.0) *
                       geo::kApproxDistanceLatDegrees)};
}

// An area from its outline corners, node ids 1, 2, ... in order.
walkable_area area_of(std::vector<geo::latlng> const& corners) {
  auto a = walkable_area{};
  auto& r = a.rings_.emplace_back();
  for (auto i = std::size_t{0U}; i != corners.size(); ++i) {
    r.ids_.push_back(static_cast<std::int64_t>(i + 1U));
    r.points_.push_back(corners[i]);
  }
  r.ids_.push_back(r.ids_.front());  // closed, like osmium's rings
  r.points_.push_back(r.points_.front());
  return a;
}

// The connector at outline corner `i` (0-based).
area_connector corner(walkable_area const& a, std::size_t const i) {
  return {.node_ = a.rings_[0].ids_[i],
          .pos_ = a.rings_[0].points_[i],
          .on_ring_ = true,
          .outline_pos_ = i};
}

std::vector<float> geodesics(walkable_area const& a,
                             std::vector<area_connector> const& connectors) {
  auto rings = std::vector<std::vector<geo::latlng>>{};
  for (auto const& r : a.rings_) {
    rings.emplace_back(r.points_.begin(), r.points_.end() - 1);
  }
  auto pts = std::vector<geo::latlng>{};
  for (auto const& c : connectors) {
    pts.push_back(c.pos_);
  }
  auto const g = area_geodesics{rings, pts};
  auto const k = connectors.size();
  auto d = std::vector<float>(k * k);
  for (auto i = std::size_t{0U}; i != k; ++i) {
    for (auto j = std::size_t{0U}; j != k; ++j) {
      d[i * k + j] = g.distance(i, j);
    }
  }
  return d;
}

auto const square = std::vector{at(0, 0), at(100, 0), at(100, 100), at(0, 100)};

}  // namespace

TEST(area_served, crossing_a_plaza_only_round_its_outline_needs_cells) {
  // Alexanderplatz: nothing mapped across. Opposite corners are 141 m apart,
  // 200 m round the outline.
  auto const a = area_of(square);
  auto const c = std::vector{corner(a, 0), corner(a, 2)};
  auto const network = network_distances(a, {}, {}, c);
  EXPECT_NEAR(200.0, network[0 * 2 + 1], 0.5);

  auto const geo = geodesics(a, c);
  auto const rel = relevant_pairs(a, c, geo);
  EXPECT_TRUE(rel.relevant_[0 * 2 + 1]);

  auto const v = served_by_geodesics(c, network, geo, rel.relevant_);
  EXPECT_FALSE(v.served_);
  EXPECT_NEAR(std::sqrt(2.0), v.worst_detour_, 0.01);
  EXPECT_FALSE(served_without_geometry(c, network));
}

TEST(area_served, a_mapped_footway_across_serves_it) {
  // Hansaplatz: footways mapped across the square. Here one, corner to
  // corner.
  auto const a = area_of(square);
  auto const c = std::vector{corner(a, 0), corner(a, 2)};
  auto const pos = hash_map<std::int64_t, geo::latlng>{{1, square[0]},
                                                       {3, square[2]}};
  auto const network = network_distances(a, {{.a_ = 1, .b_ = 3}}, pos, c);
  EXPECT_NEAR(100.0 * std::sqrt(2.0), network[0 * 2 + 1], 0.5);

  auto const geo = geodesics(a, c);
  auto const v =
      served_by_geodesics(c, network, geo, relevant_pairs(a, c, geo).relevant_);
  EXPECT_TRUE(v.served_);
  EXPECT_TRUE(served_without_geometry(c, network));
}

TEST(area_served, neighbouring_entrances_need_nothing_but_the_outline) {
  // Two entrances 10 m apart along the same edge: the outline is the
  // straight line, known without any geometry.
  auto a = area_of({at(0, 0), at(10, 0), at(100, 0), at(100, 100), at(0, 100)});
  auto const c = std::vector{corner(a, 0), corner(a, 1)};
  auto const network = network_distances(a, {}, {}, c);
  EXPECT_TRUE(served_without_geometry(c, network));

  // Nor does the pair count as one a crossing would shorten.
  auto const rel = relevant_pairs(a, c, geodesics(a, c));
  EXPECT_FALSE(rel.relevant_[0 * 2 + 1]);
}

TEST(area_served, a_staircase_ending_inside_forces_cells) {
  // Alexanderplatz: stairs from the station come up in the middle of the
  // plaza, and no mapped way leads there.
  auto const a = area_of(square);
  auto const stairs =
      area_connector{.node_ = 99, .pos_ = at(50, 50), .on_ring_ = false};
  auto const c = std::vector{corner(a, 0), corner(a, 1), stairs};
  auto const network = network_distances(a, {}, {}, c);
  EXPECT_EQ(1U, n_stranded(c, network));
  EXPECT_FALSE(served_without_geometry(c, network));

  auto const geo = geodesics(a, c);
  auto const v =
      served_by_geodesics(c, network, geo, relevant_pairs(a, c, geo).relevant_);
  EXPECT_FALSE(v.served_);  // the outline pair is fine; the stairs are not
  EXPECT_EQ(1U, v.n_stranded_);
}

TEST(area_served, a_staircase_a_mapped_way_leads_to_is_not_stranded) {
  auto const a = area_of(square);
  auto const stairs =
      area_connector{.node_ = 99, .pos_ = at(50, 50), .on_ring_ = false};
  auto const c = std::vector{corner(a, 0), stairs};
  auto const pos = hash_map<std::int64_t, geo::latlng>{{1, square[0]},
                                                       {99, at(50, 50)}};
  auto const network = network_distances(a, {{.a_ = 1, .b_ = 99}}, pos, c);
  EXPECT_EQ(0U, n_stranded(c, network));
}

TEST(area_served, the_shortcut_only_proves_served_never_the_opposite) {
  // An L-shaped area, crossed from one arm to the other. The straight line
  // cuts the missing quadrant, so the outline looks like a big detour next to
  // it - but the real shortest walk bends round the corner and is barely
  // shorter than the outline. Without geometry the answer is "not known";
  // with it, "served".
  auto const a = area_of({at(0, 0), at(100, 0), at(100, 40), at(40, 40),
                          at(40, 100), at(0, 100)});
  auto const c = std::vector{corner(a, 2), corner(a, 4)};  // (100,40), (40,100)
  auto const network = network_distances(a, {}, {}, c);
  EXPECT_NEAR(120.0, network[0 * 2 + 1], 0.5);  // via the inner corner

  EXPECT_FALSE(served_without_geometry(c, network));  // 120 > 1.2 * 84.9

  auto const geo = geodesics(a, c);
  EXPECT_NEAR(120.0, geo[0 * 2 + 1], 0.5);
  auto const v =
      served_by_geodesics(c, network, geo, relevant_pairs(a, c, geo).relevant_);
  EXPECT_TRUE(v.served_);
}
