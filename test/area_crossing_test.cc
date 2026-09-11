#include <cmath>
#include <algorithm>
#include <filesystem>
#include <memory>
#include <optional>
#include <vector>

#include "gtest/gtest.h"

#include "geo/constants.h"

#include "osr/area/crossing.h"
#include "osr/extract/extract.h"
#include "osr/lookup.h"
#include "osr/routing/profiles/foot.h"
#include "osr/routing/route.h"
#include "osr/routing/tracking.h"
#include "osr/ways.h"

#include "xml_to_pbf.h"

namespace fs = std::filesystem;
using namespace osr;

// Layer 5 of area routing: drawing a routed crossing as the shortest
// walkable line.

namespace {

constexpr auto kLat0 = 49.0;
constexpr auto kLng0 = 8.0;

geo::latlng at(double const x, double const y) {
  return {kLat0 + y / geo::kApproxDistanceLatDegrees,
          kLng0 + x / (std::cos(kLat0 * geo::kPI / 180.0) *
                       geo::kApproxDistanceLatDegrees)};
}

std::vector<geo::latlng> rect(double const x0, double const y0,
                              double const x1, double const y1) {
  return {at(x0, y0), at(x1, y0), at(x1, y1), at(x0, y1)};
}

// A 100 m square plaza with a building in the middle; entrances W and E
// halfway up its sides, OSM nodes 1 and 2, routing nodes 5 and 6.
std::optional<node_idx_t> fake_lookup(std::int64_t const osm) {
  switch (osm) {
    case 1: return node_idx_t{5U};
    case 2: return node_idx_t{6U};
    default: return std::nullopt;
  }
}

area_crossing_input plaza_cells(std::size_t const n_cells) {
  auto a = area_crossing_input{.osm_ = "way/1"};
  a.cells_.cost_ = n_cells == 1U ? std::vector{100.F} : std::vector{50.F, 50.F};
  a.cells_.connector_cell_ = {0, static_cast<area_cells::cell_idx_t>(n_cells - 1U)};
  if (n_cells == 2U) {
    a.cells_.set_neighbour(0, 1);
  }
  a.connector_nodes_ = {1, 2};
  a.connector_pos_ = {at(0, 50), at(100, 50)};
  return a;
}

area_geometry plaza_geometry() {
  return {.osm_ = "way/1",
          .rings_ = {rect(0, 0, 100, 100), rect(40, 40, 60, 60)},
          .barriers_ = {},
          .connectors_ = {at(0, 50), at(100, 50)}};
}

path::segment seg(std::uint32_t const from, std::uint32_t const to,
                  distance_t const dist) {
  return {.polyline_ = {},
          .from_ = node_idx_t{from},
          .to_ = node_idx_t{to},
          .cost_ = static_cast<cost_t>(dist),
          .dist_ = dist};
}

bool inside_building(geo::latlng const& p) {
  auto const x = (p.lng() - kLng0) * std::cos(kLat0 * geo::kPI / 180.0) *
                 geo::kApproxDistanceLatDegrees;
  auto const y = (p.lat() - kLat0) * geo::kApproxDistanceLatDegrees;
  return x > 40.01 && x < 59.99 && y > 40.01 && y < 59.99;
}

}  // namespace

TEST(area_crossing, a_crossing_is_drawn_round_the_building) {
  auto const g = area_graph{100U, {plaza_cells(1U)}, fake_lookup};
  auto const d = crossing_drawer{g, {plaza_geometry()}};

  // Walk to W, through the hub, out at E, walk on.
  auto p = path{};
  p.segments_ = {seg(4, 5, 10), seg(5, 100, 50), seg(100, 6, 50), seg(6, 7, 10)};

  auto const cs = d.crossings(p);
  ASSERT_EQ(1U, cs.size());
  auto const& c = cs.front();
  EXPECT_EQ(1U, c.first_segment_);
  EXPECT_EQ(2U, c.last_segment_);
  EXPECT_EQ(0U, c.area_);
  EXPECT_EQ(level_t{0.F}, c.level_);
  EXPECT_NEAR(100.0, c.model_distance_, 0.01);

  // Round a corner of the building: 2 * hypot(40, 10) + 20.
  EXPECT_NEAR(2.0 * std::hypot(40.0, 10.0) + 20.0, c.geodesic_distance_, 0.1);
  ASSERT_GE(c.geodesic_.size(), 3U);
  EXPECT_NEAR(0.0, geo::distance(at(0, 50), c.geodesic_.front()), 0.01);
  EXPECT_NEAR(0.0, geo::distance(at(100, 50), c.geodesic_.back()), 0.01);
  for (auto k = std::size_t{1U}; k != c.geodesic_.size(); ++k) {
    auto const& a = c.geodesic_[k - 1U];
    auto const& b = c.geodesic_[k];
    EXPECT_FALSE(inside_building({0.5 * (a.lat() + b.lat()),
                                  0.5 * (a.lng() + b.lng())}));
  }
}

TEST(area_crossing, the_drawn_path_walks_the_line) {
  auto const g = area_graph{100U, {plaza_cells(1U)}, fake_lookup};
  auto const d = crossing_drawer{g, {plaza_geometry()}};
  auto p = path{};
  p.segments_ = {seg(4, 5, 10), seg(5, 100, 50), seg(100, 6, 50), seg(6, 7, 10)};
  p.dist_ = 120.0;

  auto const drawn = d.drawn(p);
  ASSERT_EQ(3U, drawn.segments_.size());
  EXPECT_EQ(node_idx_t{4U}, drawn.segments_[0].from_);  // untouched
  EXPECT_EQ(node_idx_t{5U}, drawn.segments_[1].from_);
  EXPECT_EQ(node_idx_t{6U}, drawn.segments_[1].to_);
  EXPECT_EQ(d.crossings(p).front().geodesic_, drawn.segments_[1].polyline_);
  EXPECT_EQ(102U, drawn.segments_[1].dist_);
  EXPECT_EQ(cost_t{100U}, drawn.segments_[1].cost_);
  EXPECT_EQ(node_idx_t{7U}, drawn.segments_[2].to_);  // untouched
  EXPECT_DOUBLE_EQ(120.0, drawn.dist_);  // the router's totals stay
}

TEST(area_crossing, through_two_cells_is_one_crossing) {
  auto const g = area_graph{100U, {plaza_cells(2U)}, fake_lookup};
  auto const d = crossing_drawer{g, {plaza_geometry()}};
  auto p = path{};
  p.segments_ = {seg(5, 100, 25), seg(100, 101, 50), seg(101, 6, 25)};
  auto const cs = d.crossings(p);
  ASSERT_EQ(1U, cs.size());
  EXPECT_EQ(0U, cs[0].first_segment_);
  EXPECT_EQ(2U, cs[0].last_segment_);
  EXPECT_EQ(1U, d.drawn(p).segments_.size());
}

TEST(area_crossing, every_crossing_is_found_and_nothing_else) {
  auto const g = area_graph{100U, {plaza_cells(1U)}, fake_lookup};
  auto const d = crossing_drawer{g, {plaza_geometry()}};

  auto there_and_back = path{};
  there_and_back.segments_ = {seg(5, 100, 50), seg(100, 6, 50), seg(6, 8, 10),
                              seg(8, 6, 10), seg(6, 100, 50), seg(100, 5, 50)};
  EXPECT_EQ(2U, d.crossings(there_and_back).size());

  auto no_area = path{};
  no_area.segments_ = {seg(4, 5, 10), seg(5, 7, 10)};
  EXPECT_TRUE(d.crossings(no_area).empty());
  EXPECT_EQ(2U, d.drawn(no_area).segments_.size());
}

// --- end to end ---------------------------------------------------------------

namespace {

struct area_crossing_routing_test : public ::testing::Test {
  static void SetUpTestSuite() {
    dir_ = fs::temp_directory_path() / "osr-area-crossing-test";
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

  static geo::latlng pos(std::uint64_t const osm) {
    return ways_->get_node_pos(ways_->get_node_idx(osm_node_idx_t{osm}))
        .as_latlng();
  }

  static inline fs::path dir_{};
  static inline std::unique_ptr<ways> ways_{};
  static inline std::unique_ptr<lookup> lookup_{};
};

}  // namespace

TEST_F(area_crossing_routing_test, the_routed_crossing_is_drawn_from_a_to_b) {
  // test/area-routing.osm: the plaza lies between A (2), N1 (3), N2 (4) and
  // B (5), with one cell from entrance A to entrance B.
  auto const& w = *ways_;
  auto const a = pos(2U);
  auto const b = pos(5U);
  auto in = area_crossing_input{.osm_ = "plaza"};
  in.cells_.cost_ = {static_cast<float>(geo::distance(a, b))};
  in.cells_.connector_cell_ = {0, 0};
  in.connector_nodes_ = {2, 5};
  in.connector_pos_ = {a, b};
  auto const g = area_graph{w.n_nodes(), {in}, [&](std::int64_t const osm) {
                              return w.find_node_idx(to_osm_node_idx(osm));
                            }};
  auto const d = crossing_drawer{
      g, {{.osm_ = "plaza",
           // N1 and N2 lie on one way each, so they are no routing nodes:
           // their positions are the ones in the file.
           .rings_ = {{a, b, geo::latlng{49.00135, 8.00137},
                       geo::latlng{49.00135, 8.0}}},
           .barriers_ = {},
           .connectors_ = {a, b}}}};

  auto const s = g.sharing();
  auto const p = route(foot<false, elevator_tracking>::parameters{}, w,
                       *lookup_, search_profile::kFoot,
                       location{49.0, 7.999863, kNoLevel},
                       location{49.0, 8.001507, kNoLevel}, 3600,
                       direction::kForward, 25.0, nullptr, &s, nullptr);
  ASSERT_TRUE(p.has_value());

  auto const cs = d.crossings(*p);
  ASSERT_EQ(1U, cs.size());
  EXPECT_NEAR(geo::distance(a, b), cs[0].geodesic_distance_, 0.5);

  auto const drawn = d.drawn(*p);
  auto const it = std::ranges::find_if(drawn.segments_, [&](auto const& x) {
    return x.polyline_ == cs[0].geodesic_;
  });
  ASSERT_NE(end(drawn.segments_), it);
  EXPECT_NEAR(0.0, geo::distance(a, it->polyline_.front()), 0.01);
  EXPECT_NEAR(0.0, geo::distance(b, it->polyline_.back()), 0.01);
  EXPECT_TRUE(std::ranges::none_of(drawn.segments_, [&](auto const& x) {
    return g.is_hub(x.from_) || g.is_hub(x.to_);
  }));
}
