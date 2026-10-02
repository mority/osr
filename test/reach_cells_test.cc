#include "gtest/gtest.h"

#include <cmath>
#include <filesystem>
#include <memory>
#include <random>
#include <set>
#include <tuple>
#include <utility>

#include "utl/enumerate.h"

#include "osr/extract/extract.h"
#include "osr/lookup.h"
#include "osr/routing/parameters.h"
#include "osr/routing/reach_cells.h"
#include "osr/ways.h"

namespace fs = std::filesystem;
using namespace osr;

namespace {

reach_grid make_grid(std::vector<std::string> const& rows) {
  // rows[0] is the northernmost row, '.' = no owner, '0'..'9' = source
  auto g = reach_grid{};
  g.nx_ = static_cast<std::uint32_t>(rows.front().size());
  g.ny_ = static_cast<std::uint32_t>(rows.size());
  g.owner_.assign(g.nx_ * g.ny_, kNoSource);
  g.cost_.assign(g.nx_ * g.ny_, kInfeasible);
  g.duration_.assign(g.nx_ * g.ny_, kMaxDuration);
  for (auto y = 0U; y != g.ny_; ++y) {
    auto const& row = rows[g.ny_ - 1U - y];
    for (auto x = 0U; x != g.nx_; ++x) {
      if (row[x] != '.') {
        g.owner_[g.idx(static_cast<std::int32_t>(x),
                       static_cast<std::int32_t>(y))] =
            static_cast<source_idx_t>(row[x] - '0');
      }
    }
  }
  return g;
}

std::int64_t area2(std::vector<grid_vertex> const& ring) {
  auto a = std::int64_t{0};
  for (auto i = 0U; i != ring.size(); ++i) {
    auto const& p = ring[i];
    auto const& q = ring[(i + 1U) % ring.size()];
    a += static_cast<std::int64_t>(p.x_) * q.y_ -
         static_cast<std::int64_t>(q.x_) * p.y_;
  }
  return a;
}

// Polygon validity checks that hold for every traced grid:
// outer ring CCW, holes CW, axis-parallel edges, no ring revisits a vertex,
// and the enclosed area equals the number of owned cells.
void check_polygons(reach_grid const& g,
                    std::vector<std::vector<grid_polygon>> const& polygons) {
  for (auto source = 0U; source != polygons.size(); ++source) {
    auto cells = std::int64_t{0};
    for (auto const o : g.owner_) {
      cells += o == source ? 1 : 0;
    }

    auto area = std::int64_t{0};
    for (auto const& p : polygons[source]) {
      ASSERT_FALSE(p.empty());
      for (auto const [i, ring] : utl::enumerate(p)) {
        ASSERT_GE(ring.size(), 4U);
        auto const a = area2(ring);
        if (i == 0U) {
          EXPECT_GT(a, 0) << "outer ring must be counter-clockwise";
        } else {
          EXPECT_LT(a, 0) << "hole must be clockwise";
        }
        area += a;

        auto seen = std::set<std::pair<std::int32_t, std::int32_t>>{};
        for (auto j = 0U; j != ring.size(); ++j) {
          auto const& v = ring[j];
          auto const& next = ring[(j + 1U) % ring.size()];
          EXPECT_TRUE(v.x_ == next.x_ || v.y_ == next.y_);
          EXPECT_TRUE(seen.emplace(v.x_, v.y_).second)
              << "ring revisits vertex (" << v.x_ << ", " << v.y_ << ")";
        }
      }
    }
    EXPECT_EQ(cells * 2, area) << "source " << source;
  }
}

// Bounded meeting searches: same best places, same zones as the full search
// `full`. Returns how many fewer nodes the bounded searches reached.
std::size_t expect_bounded_equals_full(ways const& w,
                                       lookup const& l,
                                       std::vector<location> const& sources,
                                       cost_t const max,
                                       double const max_match_distance,
                                       meet_times const& full) {
  auto const params = get_parameters(search_profile::kCar);
  auto saved = std::size_t{0U};
  for (auto const [latest_margin, total_margin] :
       {std::pair{60U, 120U}, std::pair{0U, 300U}}) {
    auto const bounded = compute_meet_times(
        params, w, l, search_profile::kCar, sources, max, max_match_distance,
        nullptr,
        {.latest_margin_ =
             duration_t{static_cast<std::uint16_t>(latest_margin)},
         .total_margin_ =
             duration_t{static_cast<std::uint16_t>(total_margin)}});
    EXPECT_EQ(full.latest_.has_value(), bounded.latest_.has_value());
    if (!full.latest_.has_value() || !bounded.latest_.has_value()) {
      continue;
    }
    EXPECT_EQ(full.latest_->node_.latest_, bounded.latest_->node_.latest_);
    EXPECT_EQ(full.total_->node_.total_, bounded.total_->node_.total_);
    EXPECT_EQ(full.latest_->per_source_, bounded.latest_->per_source_);
    EXPECT_EQ(full.total_->per_source_, bounded.total_->per_source_);

    auto const zone = [&](meet_times const& m) {
      auto in = std::set<
          std::tuple<node_idx_t::value_t, std::uint32_t, std::uint32_t>>{};
      for (auto const& n : m.nodes_) {
        if (n.latest_ <= full.latest_->node_.latest_ + latest_margin ||
            n.total_ <= full.total_->node_.total_ + total_margin) {
          in.emplace(to_idx(n.node_), n.latest_, n.total_);
        }
      }
      return in;
    };
    EXPECT_EQ(zone(full), zone(bounded));
    EXPECT_LE(bounded.nodes_.size(), full.nodes_.size());
    saved += full.nodes_.size() -
             std::min(full.nodes_.size(), bounded.nodes_.size());
  }
  return saved;
}

}  // namespace

TEST(reach_cells, single_cell) {
  auto const g = make_grid({"0"});
  auto const p = trace_reach_grid(g, 1U);
  ASSERT_EQ(1U, p[0].size());
  ASSERT_EQ(1U, p[0][0].size());
  EXPECT_EQ((std::vector<grid_vertex>{{0, 0}, {1, 0}, {1, 1}, {0, 1}}),
            p[0][0][0]);
  check_polygons(g, p);
}

TEST(reach_cells, l_shape) {
  auto const g = make_grid({"0.",  //
                            "00"});
  auto const p = trace_reach_grid(g, 1U);
  ASSERT_EQ(1U, p[0].size());
  ASSERT_EQ(1U, p[0][0].size());
  EXPECT_EQ(6U, p[0][0][0].size());
  check_polygons(g, p);
}

TEST(reach_cells, hole) {
  auto const g = make_grid({"000",  //
                            "010",  //
                            "000"});
  auto const p = trace_reach_grid(g, 2U);
  ASSERT_EQ(1U, p[0].size());
  ASSERT_EQ(2U, p[0][0].size());
  EXPECT_EQ(4U, p[0][0][0].size());
  EXPECT_EQ(4U, p[0][0][1].size());
  ASSERT_EQ(1U, p[1].size());
  ASSERT_EQ(1U, p[1][0].size());
  check_polygons(g, p);
}

TEST(reach_cells, diagonal_cells_are_separate_polygons) {
  auto const g = make_grid({".0",  //
                            "0."});
  auto const p = trace_reach_grid(g, 1U);
  EXPECT_EQ(2U, p[0].size());
  check_polygons(g, p);
}

TEST(reach_cells, island_in_hole) {
  auto const g = make_grid({"00000",  //
                            "01110",  //
                            "01010",  //
                            "01110",  //
                            "00000"});
  auto const p = trace_reach_grid(g, 2U);
  ASSERT_EQ(2U, p[0].size());
  auto const outer_with_hole = p[0][0].size() == 2U ? 0U : 1U;
  EXPECT_EQ(2U, p[0][outer_with_hole].size());
  EXPECT_EQ(1U, p[0][1U - outer_with_hole].size());
  ASSERT_EQ(1U, p[1].size());
  EXPECT_EQ(2U, p[1][0].size());
  check_polygons(g, p);
}

TEST(reach_cells, random_grids) {
  auto rng = std::mt19937{42U};
  auto dist = std::uniform_int_distribution<int>{0, 3};
  for (auto round = 0U; round != 200U; ++round) {
    auto rows = std::vector<std::string>(12U, std::string(15U, '.'));
    for (auto& row : rows) {
      for (auto& c : row) {
        auto const o = dist(rng);
        c = o == 3 ? '.' : static_cast<char>('0' + o);
      }
    }
    auto const g = make_grid(rows);
    check_polygons(g, trace_reach_grid(g, 3U));
  }
}

TEST(reach_cells, fill_gaps) {
  auto g = make_grid({".....",  //
                      ".0.1.",  //
                      "....."});
  g.cost_[g.idx(1, 1)] = 10U;
  g.cost_[g.idx(3, 1)] = 20U;
  g.duration_[g.idx(1, 1)] = duration_t{7U};
  fill_gaps(g, 1U);
  EXPECT_EQ(duration_t{7U}, g.duration_[g.idx(2, 1)]);
  // the column between both goes to the cheaper source
  EXPECT_EQ(0U, g.owner(2, 1));
  EXPECT_EQ(0U, g.owner(0, 0));
  EXPECT_EQ(1U, g.owner(4, 2));
}

TEST(reach_cells, close_holes) {
  auto g = make_grid({"......",  //
                      ".000..",  //
                      ".0..0.",  //
                      ".0001.",  //
                      "......"});
  close_holes(g);
  // enclosed cells are filled, the open area on the right stays empty
  EXPECT_NE(kNoSource, g.owner(2, 2));
  EXPECT_NE(kNoSource, g.owner(3, 2));
  EXPECT_EQ(kNoSource, g.owner(4, 3));
  EXPECT_EQ(kNoSource, g.owner(0, 0));
  auto const polygons = trace_reach_grid(g, 2U);
  check_polygons(g, polygons);
  for (auto const& poly : polygons[0]) {
    EXPECT_EQ(1U, poly.size()) << "no holes left";
  }
}

// Smoothed borders between two sources must be identical in both polygons
// (reversed), so there are no gaps or overlaps. Only the outer border of the
// fully covered grid has no counterpart.
TEST(reach_cells, smoothing_is_watertight) {
  auto rng = std::mt19937{3U};
  auto dist = std::uniform_int_distribution<int>{0, 2};
  for (auto round = 0U; round != 100U; ++round) {
    auto rows = std::vector<std::string>(10U, std::string(12U, '0'));
    for (auto& row : rows) {
      for (auto& c : row) {
        c = static_cast<char>('0' + dist(rng));
      }
    }
    auto g = make_grid(rows);
    g.dlat_ = 1.0;
    g.dlng_ = 1.0;  // lat = y, lng = x

    using seg = std::pair<std::pair<double, double>, std::pair<double, double>>;
    auto segments = std::set<seg>{};
    auto all = std::vector<seg>{};
    for (auto const& polys : trace_reach_grid(g, 3U)) {
      for (auto const& poly : polys) {
        for (auto const& ring : to_latlng(g, poly, 2U)) {
          for (auto i = 1U; i < ring.size(); ++i) {
            auto const s = seg{{ring[i - 1U].lng_, ring[i - 1U].lat_},
                               {ring[i].lng_, ring[i].lat_}};
            EXPECT_TRUE(segments.insert(s).second);
            all.push_back(s);
          }
        }
      }
    }

    // The outer border of the fully covered grid has no counterpart. It stays
    // close to the grid frame (only its corners are rounded).
    auto const on_border = [&](std::pair<double, double> const& p) {
      auto const [x, y] = p;
      return std::min({x, y, g.nx_ - x, g.ny_ - y}) <= 0.5;
    };
    for (auto const& [a, b] : all) {
      if (!segments.contains(seg{b, a})) {
        EXPECT_TRUE(on_border(a) && on_border(b))
            << "unmatched inner segment (" << a.first << ", " << a.second
            << ") -> (" << b.first << ", " << b.second << ")";
      }
    }
  }
}

TEST(reach_cells, bands) {
  // source 0: durations 0..4 along the row, band size 2 -> bands 0,0,1,1,2
  // with 3 bands, durations >= 6 fall into the last band
  auto g = make_grid({"00000",  //
                      "11.11"});
  for (auto x = 0; x != 5; ++x) {
    g.duration_[g.idx(x, 1)] = duration_t{static_cast<std::uint16_t>(x)};
    g.duration_[g.idx(x, 0)] = duration_t{100U};
  }
  auto const p = trace_reach_bands(g, 2U, duration_t{2U}, 3U);
  ASSERT_EQ(6U, p.size());
  auto const cells = [](std::vector<grid_polygon> const& polys) {
    auto a = std::int64_t{0};
    for (auto const& poly : polys) {
      for (auto const& ring : poly) {
        a += area2(ring);
      }
    }
    return a / 2;
  };
  EXPECT_EQ(2, cells(p[0]));  // source 0, band 0
  EXPECT_EQ(2, cells(p[1]));
  EXPECT_EQ(1, cells(p[2]));
  EXPECT_EQ(0, cells(p[3]));  // source 1, band 0, 1
  EXPECT_EQ(0, cells(p[4]));
  EXPECT_EQ(4, cells(p[5]));  // source 1, last band (open-ended), 2 polygons
  EXPECT_EQ(2U, p[5].size());
}

TEST(reach_cells, owner_is_cheapest_source) {
  auto const dir = fs::temp_directory_path() / "osr_reach_cells_test";
  auto ec = std::error_code{};
  fs::remove_all(dir, ec);
  fs::create_directories(dir, ec);
  extract(false, "test/luisenplatz-darmstadt.osm.pbf", dir, {});
  auto const w = ways{dir, cista::mmap::protection::READ};
  auto const l = lookup{w, dir, cista::mmap::protection::READ};

  constexpr auto const kMax = cost_t{900U};
  constexpr auto const kMaxMatchDistance = 100.0;
  auto const params = get_parameters(search_profile::kCar);

  auto rng = std::mt19937{7U};
  auto node_dist =
      std::uniform_int_distribution<node_idx_t::value_t>{0U, w.n_nodes() - 1U};
  for (auto round = 0U; round != 5U; ++round) {
    auto sources = std::vector<location>{};
    for (auto i = 0U; i != 3U; ++i) {
      auto const n = node_idx_t{node_dist(rng)};
      sources.push_back({w.get_node_pos(n).as_latlng(), kNoLevel});
    }

    auto const multi = compute_reach_owners(params, w, l, search_profile::kCar,
                                            sources, kMax, kMaxMatchDistance);

    auto single = std::vector<hash_map<node_idx_t, cost_t>>{};
    auto reached = hash_set<node_idx_t>{};
    for (auto const& s : sources) {
      auto& costs = single.emplace_back();
      for (auto const& n :
           compute_reach_owners(params, w, l, search_profile::kCar, {s}, kMax,
                                kMaxMatchDistance)) {
        EXPECT_EQ(0U, n.source_);
        costs.emplace(n.node_, n.cost_);
        reached.emplace(n.node_);
      }
    }

    EXPECT_EQ(reached.size(), multi.size());
    for (auto const& n : multi) {
      ASSERT_LT(n.source_, sources.size());
      auto min = kInfeasible;
      for (auto const& costs : single) {
        if (auto const it = costs.find(n.node_); it != end(costs)) {
          min = std::min(min, it->second);
        }
      }
      EXPECT_EQ(min, n.cost_) << "node " << n.node_;
      auto const it = single[n.source_].find(n.node_);
      ASSERT_NE(it, end(single[n.source_])) << "node " << n.node_;
      EXPECT_EQ(n.cost_, it->second) << "node " << n.node_;
    }

    // Meeting times: latest arrival = max over the single source searches,
    // only for nodes reached by every source.
    auto const meet = compute_meet_times(params, w, l, search_profile::kCar,
                                         sources, kMax, kMaxMatchDistance);
    auto expected_reached_by_all = 0U;
    for (auto const& n : multi) {
      auto all = true;
      for (auto const& costs : single) {
        all = all && costs.contains(n.node_);
      }
      expected_reached_by_all += all ? 1U : 0U;
    }
    EXPECT_EQ(expected_reached_by_all, meet.nodes_.size());
    ASSERT_EQ(meet.latest_.has_value(), meet.total_.has_value());
    if (meet.latest_.has_value()) {
      for (auto const& m : meet.nodes_) {
        EXPECT_LE(meet.latest_->node_.latest_, m.latest_);
        EXPECT_LE(meet.total_->node_.total_, m.total_);
        EXPECT_LE(m.latest_, m.total_);
      }
      for (auto const* b : {&*meet.latest_, &*meet.total_}) {
        ASSERT_EQ(sources.size(), b->per_source_.size());
        auto latest = 0U;
        auto total = 0U;
        for (auto const d : b->per_source_) {
          latest = std::max(latest, static_cast<unsigned>(d.count()));
          total += d.count();
        }
        EXPECT_EQ(b->node_.latest_, latest);
        EXPECT_EQ(b->node_.total_, total);
      }
    }

    expect_bounded_equals_full(w, l, sources, kMax, kMaxMatchDistance, meet);

    auto const g = make_reach_grid(w, multi, meet.nodes_, 50.0, 2U);
    check_polygons(g, trace_reach_grid(g, sources.size()));
    if (meet.latest_.has_value()) {
      // the zones at the optimum contain the cell of the best node
      EXPECT_FALSE(trace_meet_zone(g, meet_criterion::kLatest,
                                   meet.latest_->node_.latest_)
                       .empty());
      EXPECT_FALSE(
          trace_meet_zone(g, meet_criterion::kTotal, meet.total_->node_.total_)
              .empty());
    }
  }
}

// On a larger extract with nearby sources, the bounded searches stop long
// before `max` and still find the same best places and zones.
TEST(reach_cells, bounded_meet_search) {
  auto const dir = fs::temp_directory_path() / "osr_reach_cells_bounded_test";
  auto ec = std::error_code{};
  fs::remove_all(dir, ec);
  fs::create_directories(dir, ec);
  extract(false, "test/london-northern-line.osm.pbf", dir, {});
  auto const w = ways{dir, cista::mmap::protection::READ};
  auto const l = lookup{w, dir, cista::mmap::protection::READ};

  constexpr auto const kMax = cost_t{3600U};
  constexpr auto const kMaxMatchDistance = 100.0;
  auto const params = get_parameters(search_profile::kCar);

  auto rng = std::mt19937{11U};
  auto offset = std::uniform_real_distribution<double>{-0.004, 0.004};
  auto saved = std::size_t{0U};
  for (auto round = 0U; round != 4U; ++round) {
    // around Kennington, sources ~ 300m apart
    auto sources = std::vector<location>{};
    for (auto i = 0U; i != 3U; ++i) {
      sources.push_back(
          {geo::latlng{51.488 + offset(rng), -0.111 + offset(rng)}, kNoLevel});
    }
    auto const full = compute_meet_times(params, w, l, search_profile::kCar,
                                         sources, kMax, kMaxMatchDistance);
    ASSERT_TRUE(full.latest_.has_value());
    saved += expect_bounded_equals_full(w, l, sources, kMax, kMaxMatchDistance,
                                        full);
  }
  EXPECT_GT(saved, 0U) << "bounded searches never stopped early";
}
