#include <cmath>
#include <cstdio>
#include <algorithm>
#include <random>
#include <vector>

#include "gtest/gtest.h"

#include "geo/constants.h"

#include "osr/area/geodesic.h"

using namespace osr;

// Layer 2's two visibility algorithms - every pair against every wall, and
// Lee's rotational sweep - share their geometric tests, so they must produce
// the same visibility graph edge for edge. And whatever points are added
// inside an area, the distances between its connectors must not change: a
// shortest path bends only at reflex corners. Both on random areas made to
// be awkward - buildings clipped against the outline, buildings touching at
// corners, fences, coordinates on a coarse grid.

namespace {

constexpr auto kLat0 = 52.5;
constexpr auto kLng0 = 13.4;

geo::latlng at(double const x, double const y) {
  return {kLat0 + y / geo::kApproxDistanceLatDegrees,
          kLng0 + x / (std::cos(kLat0 * geo::kPI / 180.0) *
                       geo::kApproxDistanceLatDegrees)};
}

struct scenario {
  std::vector<std::vector<geo::latlng>> rings_;
  std::vector<std::vector<geo::latlng>> fences_;
  std::vector<geo::latlng> connectors_;
  std::vector<geo::latlng> interior_;  // points added inside
};

scenario make_scenario(std::mt19937& rng, double const grid) {
  auto u = std::uniform_real_distribution<double>{0.0, 1.0};
  auto const snap = [&](double const v) {
    return grid > 0.0 ? std::round(v / grid) * grid : v;
  };
  auto const point = [&](double const x, double const y) {
    return at(snap(x), snap(y));
  };

  auto s = scenario{};

  // A star-shaped outline: never self-intersecting.
  auto const n = 6 + static_cast<int>(u(rng) * 10.0);
  auto outline = std::vector<std::pair<double, double>>{};
  for (auto i = 0; i != n; ++i) {
    auto const a = 2.0 * geo::kPI * i / n;
    auto const r = 60.0 + 40.0 * u(rng);
    outline.emplace_back(r * std::cos(a), r * std::sin(a));
  }
  auto& ring = s.rings_.emplace_back();
  for (auto const& [x, y] : outline) {
    ring.push_back(point(x, y));
  }

  // Buildings: axis-parallel boxes, some reaching past the outline (the
  // oracle cuts them out of it), some sharing a corner or a wall.
  auto const n_buildings = static_cast<int>(u(rng) * 5.0);
  for (auto i = 0; i != n_buildings; ++i) {
    auto const x = -80.0 + 160.0 * u(rng);
    auto const y = -80.0 + 160.0 * u(rng);
    auto const w = 5.0 + 25.0 * u(rng);
    auto const h = 5.0 + 25.0 * u(rng);
    s.rings_.push_back({point(x, y), point(x + w, y), point(x + w, y + h),
                        point(x, y + h)});
    if (u(rng) < 0.3) {  // a neighbour touching at the corner
      s.rings_.push_back({point(x + w, y + h), point(x + 2 * w, y + h),
                          point(x + 2 * w, y + 2 * h), point(x + w, y + 2 * h)});
    }
  }

  auto const n_fences = static_cast<int>(u(rng) * 3.0);
  for (auto i = 0; i != n_fences; ++i) {
    auto& f = s.fences_.emplace_back();
    auto const k = 2 + static_cast<int>(u(rng) * 3.0);
    for (auto j = 0; j != k; ++j) {
      f.push_back(point(-90.0 + 180.0 * u(rng), -90.0 + 180.0 * u(rng)));
    }
  }

  // Connectors on outline corners, plus a few inside. (Not through `ring`:
  // adding the buildings may have moved it.)
  for (auto i = 0; i < n; i += 1 + static_cast<int>(u(rng) * 3.0)) {
    s.connectors_.push_back(s.rings_.front()[static_cast<std::size_t>(i)]);
  }
  for (auto i = 0; i != 3; ++i) {
    s.connectors_.push_back(point(-50.0 + 100.0 * u(rng), -50.0 + 100.0 * u(rng)));
    s.interior_.push_back(point(-50.0 + 100.0 * u(rng), -50.0 + 100.0 * u(rng)));
  }
  return s;
}

struct tally {
  int n_areas_{0};
  int n_edges_{0};
  int n_edge_mismatches_{0};
  int n_distance_mismatches_{0};
  int n_pairs_{0};
  int worst_seed_{-1};
  std::string examples_;  // the first few disagreements, for the log
};

std::string meters(geo::latlng const& p) {
  auto const x = (p.lng() - kLng0) * std::cos(kLat0 * geo::kPI / 180.0) *
                 geo::kApproxDistanceLatDegrees;
  auto const y = (p.lat() - kLat0) * geo::kApproxDistanceLatDegrees;
  char buf[64];
  std::snprintf(buf, sizeof(buf), "(%.4f, %.4f)", x, y);
  return buf;
}

tally run(double const grid, int const n_cases,
          double const clearance = 0.0) {
  auto t = tally{};
  for (auto seed = 0; seed != n_cases; ++seed) {
    auto rng = std::mt19937{static_cast<unsigned>(seed)};
    auto const s = make_scenario(rng, grid);

    auto const pairwise = area_geodesics{
        s.rings_, s.connectors_, s.fences_,
        {clearance, visibility_algorithm::kPairwise}};
    auto const sweep =
        area_geodesics{s.rings_, s.connectors_, s.fences_,
                       {clearance, visibility_algorithm::kSweep}};
    if (!pairwise.valid()) {
      continue;
    }
    ++t.n_areas_;
    t.n_edges_ += static_cast<int>(pairwise.edges().size());

    auto a = pairwise.edges();
    auto b = sweep.edges();
    auto only = std::vector<std::pair<std::uint32_t, std::uint32_t>>{};
    std::ranges::set_symmetric_difference(a, b, std::back_inserter(only));
    if (!only.empty()) {
      t.n_edge_mismatches_ += static_cast<int>(only.size());
      t.worst_seed_ = seed;
      auto const& v = pairwise.vertices();
      for (auto const& [x, y] : only) {
        if (std::ranges::count(t.examples_, '\n') >= 12) {
          break;
        }
        auto const in_pairwise = std::ranges::binary_search(a, std::pair{x, y});
        t.examples_ += "  seed " + std::to_string(seed) + ": " + meters(v[x]) +
                       " - " + meters(v[y]) + " only " +
                       (in_pairwise ? "pairwise" : "sweep") + "\n";
      }
    }

    auto with = s.connectors_;
    with.insert(end(with), begin(s.interior_), end(s.interior_));
    auto const more = area_geodesics{s.rings_, with, s.fences_,
                                     {clearance}};
    for (auto i = std::size_t{0U}; i != s.connectors_.size(); ++i) {
      for (auto j = i + 1U; j != s.connectors_.size(); ++j) {
        auto const d1 = sweep.distance(i, j);
        auto const d2 = more.distance(i, j);
        ++t.n_pairs_;
        if (std::isinf(d1) != std::isinf(d2) ||
            (!std::isinf(d1) && std::abs(d1 - d2) > 1e-3F)) {
          ++t.n_distance_mismatches_;
          t.worst_seed_ = seed;
        }
      }
    }
  }
  return t;
}

}  // namespace

TEST(area_visibility, sweep_matches_pairwise_on_continuous_coordinates) {
  auto const t = run(0.0, 400);
  std::printf("[          ] %d areas, %d edges, %d pairs\n", t.n_areas_,
              t.n_edges_, t.n_pairs_);
  ASSERT_GT(t.n_areas_, 300);
  EXPECT_EQ(0, t.n_edge_mismatches_) << "seed " << t.worst_seed_ << "\n"
                                      << t.examples_;
  EXPECT_EQ(0, t.n_distance_mismatches_) << "seed " << t.worst_seed_;
}

TEST(area_visibility, sweep_matches_pairwise_on_a_coarse_grid) {
  // On a 1 m grid, walls in line, points on walls and corners touching are
  // everywhere.
  auto const t = run(1.0, 400);
  std::printf("[          ] %d areas, %d edges, %d pairs\n", t.n_areas_,
              t.n_edges_, t.n_pairs_);
  ASSERT_GT(t.n_areas_, 300);
  EXPECT_EQ(0, t.n_edge_mismatches_) << "seed " << t.worst_seed_ << "\n"
                                      << t.examples_;
  EXPECT_EQ(0, t.n_distance_mismatches_) << "seed " << t.worst_seed_;
}

TEST(area_visibility, sweep_matches_pairwise_with_wall_clearance) {
  // The free space shrunk by kWallClearance: mitred corners everywhere, and
  // every connector on the outline steps onto it.
  auto const t = run(0.0, 400, kWallClearance);
  std::printf("[          ] %d areas, %d edges, %d pairs\n", t.n_areas_,
              t.n_edges_, t.n_pairs_);
  ASSERT_GT(t.n_areas_, 300);
  EXPECT_EQ(0, t.n_edge_mismatches_) << "seed " << t.worst_seed_ << "\n"
                                      << t.examples_;
  EXPECT_EQ(0, t.n_distance_mismatches_) << "seed " << t.worst_seed_;
}

// For looking at one disagreement: OSR_SEED=<seed> OSR_GRID=<m> prints the
// area (in meters) and the edges only one algorithm has, as JSON.
TEST(area_visibility, DISABLED_dump_seed) {
  auto const* const seed_env = std::getenv("OSR_SEED");
  auto const* const grid_env = std::getenv("OSR_GRID");
  if (seed_env == nullptr) {
    GTEST_SKIP() << "set OSR_SEED";
  }
  auto rng = std::mt19937{static_cast<unsigned>(std::atoi(seed_env))};
  auto const s =
      make_scenario(rng, grid_env == nullptr ? 0.0 : std::atof(grid_env));
  auto const pairwise = area_geodesics{s.rings_, s.connectors_, s.fences_,
                                       {.algorithm_ = visibility_algorithm::kPairwise}};
  auto const sweep = area_geodesics{s.rings_, s.connectors_, s.fences_,
                                    {.algorithm_ = visibility_algorithm::kSweep}};
  auto const lines = [](std::vector<std::vector<geo::latlng>> const& v) {
    auto out = std::string{"["};
    for (auto const& r : v) {
      out += out.size() == 1U ? "[" : ",[";
      for (auto i = std::size_t{0U}; i != r.size(); ++i) {
        auto m = meters(r[i]);
        std::ranges::replace(m, '(', '[');
        std::ranges::replace(m, ')', ']');
        out += (i == 0U ? "" : ",") + m;
      }
      out += "]";
    }
    return out + "]";
  };
  auto only = std::vector<std::pair<std::uint32_t, std::uint32_t>>{};
  std::ranges::set_symmetric_difference(pairwise.edges(), sweep.edges(),
                                        std::back_inserter(only));
  auto edges = std::vector<std::vector<geo::latlng>>{};
  for (auto const& [a, b] : only) {
    edges.push_back({pairwise.vertices()[a], pairwise.vertices()[b]});
  }
  std::printf("{\"rings\":%s,\"fences\":%s,\"only_sweep\":%s,\"swept\":%d}\n",
              lines(s.rings_).c_str(), lines(s.fences_).c_str(),
              lines(edges).c_str(), sweep.swept() ? 1 : 0);
}

TEST(area_visibility, buildings_touching_at_a_corner_are_no_doorway) {
  // Two buildings meet at (50, 50), one to the lower left, one to the upper
  // right. Walking diagonally from the lower right corner area to the upper
  // left one passes exactly through the meeting point: not allowed.
  auto const rings = std::vector<std::vector<geo::latlng>>{
      {at(0, 0), at(100, 0), at(100, 100), at(0, 100)},
      {at(20, 20), at(50, 20), at(50, 50), at(20, 50)},
      {at(50, 50), at(80, 50), at(80, 80), at(50, 80)}};
  EXPECT_FALSE(area_geodesics::is_segment_inside(rings, at(20, 80), at(80, 20)));
  EXPECT_TRUE(area_geodesics::is_segment_inside(rings, at(10, 90), at(10, 10)));
}

TEST(area_visibility, walking_along_a_wall_is_fine) {
  auto const square = std::vector<std::vector<geo::latlng>>{
      {at(0, 0), at(100, 0), at(100, 100), at(0, 100)}};
  EXPECT_TRUE(area_geodesics::is_segment_inside(square, at(0, 0), at(100, 0)));
  EXPECT_TRUE(area_geodesics::is_segment_inside(square, at(10, 0), at(90, 0)));
}
