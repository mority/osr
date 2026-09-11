#include <algorithm>
#include <filesystem>

#include "gtest/gtest.h"

#include "osr/area/pipeline.h"

#include "xml_to_pbf.h"

using namespace osr;

// Layers 0 to 3 of area routing run over OSM data - test/area-pipeline.osm.

namespace {

std::filesystem::path pbf() {
  static auto const path = test::osm_to_pbf("test/area-pipeline.osm");
  return path;
}

std::size_t find(area_data const& d, std::int64_t const id) {
  auto const it = std::ranges::find_if(
      d.areas_, [&](walkable_area const& a) { return a.id_ == id; });
  EXPECT_NE(end(d.areas_), it) << "no area " << id;
  return static_cast<std::size_t>(std::distance(begin(d.areas_), it));
}

}  // namespace

TEST(area_pipeline, the_square_is_ignored_and_the_entrance_merged) {
  auto const d = collect_areas(pbf(), {});
  // Plaza A (with the entrance area merged into it) and plaza B; the bare
  // place=square is no walkable area.
  ASSERT_EQ(2U, d.areas_.size());
  EXPECT_EQ(2, d.merged_.n_members_);
  auto const& a = d.areas_[find(d, 2001)];
  EXPECT_EQ((std::vector<std::string>{"way/2001", "way/2105"}), a.members_);

  // With --place-square the square covers plaza A and is merged too.
  auto const with_square =
      collect_areas(pbf(), {.walkable_ = {.place_square_ = true}});
  EXPECT_EQ(2U, with_square.areas_.size());
  EXPECT_EQ(3, with_square.merged_.n_members_);
}

TEST(area_pipeline, plaza_a_needs_cells) {
  auto const d = collect_areas(pbf(), {});
  auto const i = find(d, 2001);

  // The building stands in it.
  ASSERT_EQ(1U, d.buildings_[i].size());

  auto const pa = prepare_area(d, i, {});
  // Entered at two corners, plus the stairs' end inside.
  EXPECT_EQ(2, pa.counts_.n_boundary_);
  EXPECT_EQ(1, pa.counts_.n_interior_);
  ASSERT_EQ(3U, pa.connectors_.size());
  EXPECT_TRUE(std::ranges::any_of(pa.connectors_, [](auto const& c) {
    return c.node_ == 12 && !c.on_ring_;
  }));

  // Outline and the building as a hole.
  EXPECT_EQ(2U, pa.rings_.size());

  // No way leads to the stairs' end, and the corners are a detour apart
  // round the outline: cells.
  EXPECT_FALSE(pa.served_without_geometry_);
  EXPECT_EQ(1U, pa.verdict_.n_stranded_);
  EXPECT_FALSE(pa.verdict_.served_);
  EXPECT_EQ(area_status::kMeshed, pa.status_);
  ASSERT_TRUE(pa.cells_.has_value());
  EXPECT_EQ(3U, pa.cells_->n_connectors());
  EXPECT_EQ(pa.cells_->n_cells(), pa.regions_.size());

  // Without cells asked for, it stops at the verdict.
  EXPECT_EQ(area_status::kNeedsCells,
            prepare_area(d, i, {}, false, false).status_);
}

TEST(area_pipeline, plaza_b_is_served_by_the_footway_across) {
  auto const d = collect_areas(pbf(), {});
  auto const i = find(d, 3001);
  EXPECT_FALSE(d.interior_edges_[i].empty());

  // Known without any geometry: the geodesics are never built.
  auto const pa = prepare_area(d, i, {});
  EXPECT_EQ(area_status::kServed, pa.status_);
  EXPECT_TRUE(pa.served_without_geometry_);
  EXPECT_FALSE(pa.geodesics_.has_value());

  // For statistics they can be built anyway; the verdict is the same.
  auto const full = prepare_area(d, i, {}, true);
  EXPECT_EQ(area_status::kServed, full.status_);
  EXPECT_TRUE(full.geodesics_.has_value());
  EXPECT_TRUE(full.verdict_.served_);

  // Ignoring the mapped footway, it needs cells.
  auto const alone = prepare_area(d, i, {.interior_ways_ = false});
  EXPECT_EQ(area_status::kMeshed, alone.status_);
}

TEST(area_pipeline, areas_over_the_vertex_cap_are_declined) {
  auto const d = collect_areas(pbf(), {});
  auto const pa = prepare_area(d, find(d, 2001), {.max_vertices_ = 5U});
  EXPECT_EQ(area_status::kTooBig, pa.status_);
  EXPECT_FALSE(pa.geodesics_.has_value());
}
