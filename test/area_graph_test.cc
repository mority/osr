#include <algorithm>
#include <filesystem>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "gtest/gtest.h"

#include "geo/constants.h"

#include "osr/area/area_graph.h"
#include "osr/extract/extract.h"
#include "osr/lookup.h"
#include "osr/routing/profiles/foot.h"
#include "osr/routing/route.h"
#include "osr/routing/tracking.h"
#include "osr/ways.h"

#include "xml_to_pbf.h"

namespace fs = std::filesystem;
using namespace osr;

// Layer 4 of area routing: the cells as a graph the router walks.

namespace {

tag_lookup tags(std::map<std::string, std::string, std::less<>> kv) {
  auto const m = std::make_shared<decltype(kv)>(std::move(kv));
  return [m](std::string_view const key) -> std::optional<std::string_view> {
    auto const it = m->find(key);
    return it == end(*m) ? std::nullopt
                         : std::optional<std::string_view>{it->second};
  };
}

// A lookup that knows OSM nodes 1 and 2 as routing nodes 5 and 6.
std::optional<node_idx_t> fake_lookup(std::int64_t const osm) {
  switch (osm) {
    case 1: return node_idx_t{5U};
    case 2: return node_idx_t{6U};
    default: return std::nullopt;
  }
}

distance_t edge(sharing_data const& s, node_idx_t const from,
                node_idx_t const to) {
  auto const it = s.additional_edges_.find(from);
  if (it != end(s.additional_edges_)) {
    for (auto const& e : it->second) {
      if (e.to_ == to) {
        return e.distance_;
      }
    }
  }
  return std::numeric_limits<distance_t>::max();
}

// Two cells, 10 m and 30 m to cross, neighbours; three connectors, the
// third on an OSM node that is no routing node.
area_crossing_input two_cells() {
  auto a = area_crossing_input{.osm_ = "way/7"};
  a.cells_.cost_ = {10.F, 30.F};
  a.cells_.connector_cell_ = {0, 1, 1};
  a.cells_.set_neighbour(0, 1);
  a.connector_nodes_ = {1, 2, 3};
  a.connector_pos_ = {{49.0, 8.0}, {49.001, 8.001}, {49.0015, 8.001}};
  return a;
}

}  // namespace

TEST(area_graph, untagged_areas_are_at_street_level) {
  EXPECT_EQ(level_t{0.F}, area_floor(area_levels::of(tags({}))));
  EXPECT_EQ(level_t{-1.F}, area_floor(area_levels::of(tags({{"level", "-1"}}))));
  // Stairs mapped as an area: open from any floor.
  EXPECT_EQ(kNoLevel, area_floor(area_levels::of(tags({{"level", "0;1"}}))));
}

TEST(area_graph, connectors_and_neighbouring_hubs_are_joined_at_half_cost) {
  auto const g = area_graph{100U, {two_cells()}, fake_lookup};
  auto const s = g.sharing();
  auto const hub0 = node_idx_t{100U};
  auto const hub1 = node_idx_t{101U};

  ASSERT_EQ(2U, g.n_hubs());
  EXPECT_EQ(100U, s.additional_node_offset_);
  EXPECT_EQ(2U, s.additional_node_coordinates_.size());
  EXPECT_TRUE(g.is_hub(hub0));
  EXPECT_FALSE(g.is_hub(node_idx_t{5U}));

  EXPECT_EQ(5U, edge(s, node_idx_t{5U}, hub0));   // 10 / 2
  EXPECT_EQ(5U, edge(s, hub0, node_idx_t{5U}));
  EXPECT_EQ(15U, edge(s, node_idx_t{6U}, hub1));  // 30 / 2
  EXPECT_EQ(20U, edge(s, hub0, hub1));             // (10 + 30) / 2
  EXPECT_EQ(20U, edge(s, hub1, hub0));

  // The third connector's node is no routing node: counted, left out.
  EXPECT_EQ(3U, g.n_connectors_);
  EXPECT_EQ(1U, g.n_without_routing_node_);

  EXPECT_EQ(0U, g.area_of(hub1));
  EXPECT_EQ(std::optional<std::size_t>{1U}, g.connector_of(0, node_idx_t{6U}));
  EXPECT_EQ(std::nullopt, g.connector_of(0, node_idx_t{7U}));

  // Every hub on its area's floor.
  ASSERT_NE(nullptr, s.additional_node_levels_);
  EXPECT_EQ(level_t{0.F}, (*s.additional_node_levels_)[0]);
}

TEST(area_graph, cells_of_separate_parts_are_not_linked) {
  // An area in two walkable parts, one cell each: no neighbour pair at all,
  // so the neighbour bits were never allocated. (This crashed the backend.)
  auto a = two_cells();
  a.cells_.neighbours_.clear();
  a.level_ = level_t{-1.F};
  auto const g = area_graph{100U, {a}, fake_lookup};
  auto const s = g.sharing();
  EXPECT_EQ(std::numeric_limits<distance_t>::max(),
            edge(s, node_idx_t{100U}, node_idx_t{101U}));
  EXPECT_EQ(5U, edge(s, node_idx_t{5U}, node_idx_t{100U}));
  EXPECT_EQ(level_t{-1.F}, g.level_of(node_idx_t{101U}));
}

TEST(area_graph, hubs_of_further_areas_follow_on) {
  auto stray = two_cells();
  stray.cells_.connector_cell_[2] = area_cells::kNoCell;  // reaches nothing
  auto const g = area_graph{100U, {two_cells(), stray}, fake_lookup};
  EXPECT_EQ(4U, g.n_hubs());
  EXPECT_EQ(1U, g.area_of(node_idx_t{103U}));
  EXPECT_EQ(3U + 2U, g.n_connectors_);  // the stray is no connector of a cell
}

// --- routing, end to end ------------------------------------------------------

namespace {

struct area_routing_test : public ::testing::Test {
  static void SetUpTestSuite() {
    dir_ = fs::temp_directory_path() / "osr-area-routing-test";
    auto ec = std::error_code{};
    fs::remove_all(dir_, ec);
    fs::create_directories(dir_, ec);
    extract(false, test::osm_to_pbf("test/area-routing.osm"), dir_, {});
    ways_ = std::make_unique<ways>(dir_, cista::mmap::protection::READ);
    lookup_ =
        std::make_unique<lookup>(*ways_, dir_, cista::mmap::protection::READ);
  }

  static void TearDownTestSuite() {
    lookup_.reset();
    ways_.reset();
    auto ec = std::error_code{};
    fs::remove_all(dir_, ec);
  }

  // The plaza as one cell between entrances A (OSM node 2) and B (5),
  // crossed at the straight-line distance.
  static area_graph plaza(level_t const level) {
    auto const& w = *ways_;
    auto const a_pos = w.get_node_pos(w.get_node_idx(osm_node_idx_t{2U}));
    auto const b_pos = w.get_node_pos(w.get_node_idx(osm_node_idx_t{5U}));
    auto a = area_crossing_input{.osm_ = "plaza", .level_ = level};
    a.cells_.cost_ = {static_cast<float>(
        geo::distance(a_pos.as_latlng(), b_pos.as_latlng()))};
    a.cells_.connector_cell_ = {0, 0};
    a.connector_nodes_ = {2, 5};
    a.connector_pos_ = {a_pos.as_latlng(), b_pos.as_latlng()};
    return area_graph{w.n_nodes(), {a}, [&](std::int64_t const osm) {
                        return w.find_node_idx(to_osm_node_idx(osm));
                      }};
  }

  static std::optional<path> walk(sharing_data const* sharing,
                                  routing_algorithm const algo) {
    auto const from = location{49.0, 7.999863, kNoLevel};  // D
    auto const to = location{49.0, 8.001507, kNoLevel};    // E
    return route(foot<false, elevator_tracking>::parameters{}, *ways_,
                 *lookup_, search_profile::kFoot, from, to, 3600,
                 direction::kForward, 25.0, nullptr, sharing, nullptr, algo);
  }

  static inline fs::path dir_{};
  static inline std::unique_ptr<ways> ways_{};
  static inline std::unique_ptr<lookup> lookup_{};
};

}  // namespace

TEST_F(area_routing_test, without_cells_the_walk_goes_round) {
  auto const p = walk(nullptr, routing_algorithm::kDijkstra);
  ASSERT_TRUE(p.has_value());
  EXPECT_GT(p->dist_, 390.0);
}

TEST_F(area_routing_test, the_route_crosses_through_the_hub) {
  auto const g = plaza(level_t{0.F});
  auto const s = g.sharing();
  for (auto const algo :
       {routing_algorithm::kDijkstra, routing_algorithm::kAStarBi}) {
    auto const p = walk(&s, algo);
    ASSERT_TRUE(p.has_value());
    EXPECT_LT(p->dist_, 150.0);  // 10 + 102 + 10, not 10 + 400 + 10
    EXPECT_TRUE(std::ranges::any_of(p->segments_, [&](auto const& seg) {
      return g.is_hub(seg.to_);
    }));
  }
}

TEST_F(area_routing_test, an_area_on_another_floor_is_not_entered) {
  // The ways are untagged, i.e. at street level; the plaza is a station hall
  // on level -1.
  auto const g = plaza(level_t{-1.F});
  auto const s = g.sharing();
  auto const p = walk(&s, routing_algorithm::kDijkstra);
  ASSERT_TRUE(p.has_value());
  EXPECT_GT(p->dist_, 390.0);
}
