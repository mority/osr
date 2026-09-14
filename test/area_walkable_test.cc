#include "gtest/gtest.h"

#include <cmath>
#include <algorithm>
#include <map>
#include <memory>
#include <string>

#include "geo/constants.h"

#include "osr/area/geodesic.h"
#include "osr/area/walkable.h"
#include "osr/util/levels.h"

using namespace osr;

// Layer 0 of area routing, one test per case where the OSM data turned out
// to mean something other than it seemed to. The names say which place in
// Berlin each one comes from.

namespace {

constexpr auto kLat0 = 52.52;
constexpr auto kLng0 = 13.41;

geo::latlng at(double const x, double const y) {
  return {kLat0 + y / geo::kApproxDistanceLatDegrees,
          kLng0 + x / (std::cos(kLat0 * geo::kPI / 180.0) *
                       geo::kApproxDistanceLatDegrees)};
}

tag_lookup tags(std::map<std::string, std::string, std::less<>> kv) {
  auto const m = std::make_shared<decltype(kv)>(std::move(kv));
  return [m](std::string_view const key) -> std::optional<std::string_view> {
    auto const it = m->find(key);
    return it == end(*m) ? std::nullopt
                         : std::optional<std::string_view>{it->second};
  };
}

std::vector<geo::latlng> rect(double const x0,
                              double const y0,
                              double const x1,
                              double const y1) {
  return {at(x0, y0), at(x1, y0), at(x1, y1), at(x0, y1)};
}

area_ring ring_of(std::vector<geo::latlng> const& pts, std::int64_t first_id) {
  auto r = area_ring{};
  for (auto const& p : pts) {
    r.ids_.push_back(first_id++);
    r.points_.push_back(p);
  }
  r.ids_.push_back(r.ids_.front());  // closed, like osmium's rings
  r.points_.push_back(r.points_.front());
  return r;
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

}  // namespace

// --- what is a walkable area ------------------------------------------------

TEST(area_walkable, bare_place_square_is_not_walkable) {
  // Hansaplatz: place=square and nothing else, with Altonaer Strasse
  // running through the polygon.
  auto const t = tags({{"place", "square"}, {"name", "Hansaplatz"}});
  EXPECT_FALSE(is_walkable_area(t, true));
  EXPECT_TRUE(is_walkable_area(t, true, {.place_square_ = true}));
}

TEST(area_walkable, place_square_with_a_surface_is_walkable) {
  EXPECT_TRUE(is_walkable_area(
      tags({{"place", "square"}, {"highway", "pedestrian"}, {"area", "yes"}}),
      true));
}

TEST(area_walkable, pedestrian_multipolygon_is_walkable) {
  // Alexanderplatz' actual surface, relation 131761.
  EXPECT_TRUE(is_walkable_area(tags({{"highway", "pedestrian"}}), false));
}

TEST(area_walkable, closed_footway_needs_area_yes) {
  EXPECT_FALSE(is_walkable_area(tags({{"highway", "footway"}}), true));
  EXPECT_TRUE(
      is_walkable_area(tags({{"highway", "footway"}, {"area", "yes"}}), true));
}

TEST(area_walkable, station_and_area_no_are_not_walkable) {
  EXPECT_FALSE(is_walkable_area(
      tags({{"public_transport", "station"}, {"highway", "pedestrian"}}),
      false));
  EXPECT_FALSE(is_walkable_area(
      tags({{"highway", "pedestrian"}, {"area", "no"}}), false));
}

TEST(area_walkable, platform_is_walkable) {
  // U Hansaplatz' platform, level -1.
  EXPECT_TRUE(is_walkable_area(
      tags({{"railway", "platform"}, {"area", "yes"}, {"level", "-1"}}),
      true));
}

TEST(area_walkable, linear_ways_are_highways_that_are_not_areas) {
  EXPECT_TRUE(is_linear_way(tags({{"highway", "footway"}})));
  EXPECT_FALSE(is_linear_way(tags({{"highway", "pedestrian"}, {"area", "yes"}})));
  EXPECT_FALSE(is_linear_way(tags({{"building", "yes"}})));
}

// --- levels -------------------------------------------------------------------

TEST(area_walkable, levels_are_parsed_like_osr_tags) {
  auto const bit = [](float const l) {
    return static_cast<level_bits_t>(1) << to_idx(level_t{l});
  };
  EXPECT_EQ(bit(0.F), parse_levels("0"));
  EXPECT_EQ(bit(-1.F), parse_levels("-1"));
  EXPECT_EQ(bit(0.F) | bit(1.F), parse_levels("0;1"));
  EXPECT_EQ(bit(1.5F), parse_levels("1.5"));

  auto const l = area_levels::of(tags({{"indoor:level", "-2"}, {"layer", "-1"}}));
  EXPECT_FALSE(l.any_);
  EXPECT_EQ(bit(-2.F), l.bits_);
  EXPECT_EQ(-1, l.layer_);
  EXPECT_TRUE(area_levels::of(tags({})).any_);
}

TEST(area_walkable, untagged_joins_anything_unless_strict) {
  auto const area = area_levels::of(tags({{"level", "-1"}}));
  auto const untagged = area_levels::of(tags({}));
  EXPECT_TRUE(area.matches(untagged, false));
  EXPECT_FALSE(area.matches(untagged, true));
  EXPECT_FALSE(area.matches(area_levels::of(tags({{"level", "0"}})), false));
}

// --- what blocks ----------------------------------------------------------------

TEST(area_walkable, building_on_another_layer_still_blocks) {
  // The Park Inn at Alexanderplatz: layer=1, standing over the station.
  auto const hotel = blocker_of(tags({{"building", "hotel"}, {"layer", "1"}}));
  EXPECT_EQ(blocker::kind::kBuilding, hotel.kind_);
  EXPECT_TRUE(hotel.stands_in(area_levels::of(tags({})), false));
}

TEST(area_walkable, fence_on_a_bridge_does_not_block_the_plaza_below) {
  auto const fence = blocker_of(tags({{"barrier", "fence"}, {"layer", "1"}}));
  EXPECT_EQ(blocker::kind::kBarrier, fence.kind_);
  EXPECT_FALSE(fence.stands_in(area_levels::of(tags({})), false));
  EXPECT_TRUE(
      fence.stands_in(area_levels::of(tags({{"layer", "1"}})), false));
}

TEST(area_walkable, upper_floor_on_columns_leaves_the_ground_open) {
  // Way 23723131 at Alexanderplatz: building:min_level=1 - one walks under
  // it, but it is there on the first floor.
  auto const upper = blocker_of(tags({{"building", "commercial"},
                                      {"building:min_level", "1"},
                                      {"layer", "1"}}));
  EXPECT_EQ(blocker::kind::kBuilding, upper.kind_);
  EXPECT_FALSE(upper.stands_in(area_levels::of(tags({})), false));
  EXPECT_FALSE(upper.stands_in(area_levels::of(tags({{"level", "0"}})), false));
  EXPECT_TRUE(upper.stands_in(area_levels::of(tags({{"level", "1"}})), false));

  auto const lifted = blocker_of(tags({{"building", "yes"}, {"min_height", "5"}}));
  EXPECT_FALSE(lifted.stands_in(area_levels::of(tags({})), false));
}

TEST(area_walkable, roof_and_building_no_do_not_block) {
  EXPECT_EQ(blocker::kind::kNone, blocker_of(tags({{"building", "roof"}})).kind_);
  EXPECT_EQ(blocker::kind::kNone, blocker_of(tags({{"building", "no"}})).kind_);
  EXPECT_EQ(blocker::kind::kNone,
            blocker_of(tags({{"building:part", "yes"}})).kind_);
}

TEST(area_walkable, building_on_another_level_does_not_block) {
  auto const shop = blocker_of(tags({{"building", "retail"}, {"level", "-1"}}));
  EXPECT_FALSE(shop.stands_in(area_levels::of(tags({{"level", "0"}})), false));
  EXPECT_TRUE(shop.stands_in(area_levels::of(tags({{"level", "-1"}})), false));
}

// --- clean holes --------------------------------------------------------------

TEST(area_walkable, building_cut_out_twice_stays_cut_out) {
  // Galeria at Alexanderplatz was a hole twice over, and under the even-odd
  // rule two holes cancel: the building became walkable again.
  auto const outer = rect(0, 0, 100, 100);
  auto const galeria = rect(30, 30, 70, 70);
  auto const holes = clean_holes(outer, {galeria, galeria});
  ASSERT_EQ(1U, holes.size());
  EXPECT_NEAR(1600.0, area_m2(holes[0]), 1.0);

  // The geodesics no longer depend on it: they union the holes themselves, so
  // the doubled hole is as solid as the clean one. (Under the old even-odd
  // test, the diagonal from corner to corner - touching walls without crossing
  // them - counted as open ground.)
  auto rings = std::vector<std::vector<geo::latlng>>{outer};
  rings.insert(end(rings), begin(holes), end(holes));
  EXPECT_FALSE(area_geodesics::is_segment_inside({outer, galeria, galeria},
                                                 at(30, 30), at(70, 70)));
  EXPECT_FALSE(area_geodesics::is_segment_inside(rings, at(30, 30), at(70, 70)));

  auto const clean = area_geodesics{rings, {at(50, 5), at(50, 95)}};
  EXPECT_GT(clean.distance(0, 1), 100.0);  // round it
}

TEST(area_walkable, overlapping_buildings_become_one_hole) {
  auto const holes =
      clean_holes(rect(0, 0, 100, 100), {rect(20, 20, 60, 60), rect(40, 40, 80, 80)});
  ASSERT_EQ(1U, holes.size());
  EXPECT_NEAR(1600.0 + 1600.0 - 400.0, area_m2(holes[0]), 1.0);
}

TEST(area_walkable, building_reaching_past_the_outline_is_clipped) {
  auto const holes = clean_holes(rect(0, 0, 100, 100), {rect(80, 40, 140, 60)});
  ASSERT_EQ(1U, holes.size());
  EXPECT_NEAR(20.0 * 20.0, area_m2(holes[0]), 1.0);
}

TEST(area_walkable, courtyard_inside_a_building_is_dropped) {
  // A building with a courtyard, as an outline with its own inner ring,
  // unions to a polygon with a hole; only its outer boundary is a hole of
  // the area - the courtyard cannot be reached from it.
  auto const holes = clean_holes(
      rect(0, 0, 100, 100), {rect(20, 20, 80, 80), rect(20, 20, 80, 30),
                             rect(20, 70, 80, 80), rect(20, 20, 30, 80),
                             rect(70, 20, 80, 80)});
  ASSERT_EQ(1U, holes.size());
  EXPECT_NEAR(3600.0, area_m2(holes[0]), 1.0);
}

// --- merging ------------------------------------------------------------------

TEST(area_walkable, overlapping_areas_on_one_level_merge) {
  // Alexanderplatz: the station entrance areas lie within the pedestrian
  // surface.
  auto areas = std::vector<walkable_area>{
      {.id_ = 1, .from_way_ = false, .name_ = "plaza",
       .rings_ = {ring_of(rect(0, 0, 100, 100), 100)}},
      {.id_ = 2, .from_way_ = true, .name_ = "entrance",
       .rings_ = {ring_of(rect(90, 40, 120, 60), 200)}}};
  auto const stats = merge_overlapping(areas);
  EXPECT_EQ(2, stats.n_members_);
  ASSERT_EQ(1U, areas.size());
  auto const& merged = areas.front();
  EXPECT_EQ(1, merged.id_);  // named after the larger member
  EXPECT_EQ((std::vector<std::string>{"relation/1", "way/2"}), merged.members_);

  // Corners that survive keep their node ids; where the outlines cross, the
  // new points get negative ids.
  auto const& ids = merged.rings_.front().ids_;
  EXPECT_NE(end(ids), std::ranges::find(ids, 100));  // plaza corner (0, 0)
  EXPECT_NE(end(ids), std::ranges::find(ids, 201));  // entrance corner (120, 40)
  EXPECT_TRUE(std::ranges::any_of(ids, [](auto const id) { return id < 0; }));

  // The entrance's corners inside the plaza are gone from the outline, but a
  // way meeting them there still enters the area.
  auto const dissolved_ids = [&]() {
    auto out = std::vector<std::int64_t>{};
    for (auto const& [id, pos] : merged.dissolved_) {
      out.push_back(id);
    }
    return out;
  }();
  EXPECT_NE(end(dissolved_ids), std::ranges::find(dissolved_ids, 200));
  EXPECT_NE(end(dissolved_ids), std::ranges::find(dissolved_ids, 203));
}

TEST(area_walkable, areas_on_other_levels_or_only_touching_stay_apart) {
  auto areas = std::vector<walkable_area>{
      {.id_ = 1, .rings_ = {ring_of(rect(0, 0, 100, 100), 100)}},
      {.id_ = 2,
       .levels_ = area_levels::of(tags({{"level", "-1"}})),
       .rings_ = {ring_of(rect(50, 50, 150, 150), 200)}},
      {.id_ = 3, .rings_ = {ring_of(rect(100, 0, 200, 40), 300)}}};
  auto const stats = merge_overlapping(areas);
  EXPECT_EQ(0, stats.n_members_);
  EXPECT_EQ(3U, areas.size());
}
