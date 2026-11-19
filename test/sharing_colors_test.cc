#include "gtest/gtest.h"

#include <chrono>
#include <algorithm>
#include <array>
#include <random>
#include <span>
#include <string>

#include "osr/lookup.h"
#include "osr/routing/parameters.h"
#include "osr/routing/profiles/car_sharing.h"
#include "osr/routing/route.h"
#include "osr/routing/sharing_data.h"
#include "osr/types.h"
#include "osr/ways.h"

#include "test_graph.h"

using namespace osr;
using test::graph;

// A colored car_sharing search (sharing_data::colors_) has to give, for every
// color c, the same rides as an uncolored search whose far side - end_allowed_
// forward, start_allowed_ backward - is just colors_[c]: one search instead of
// one per color, without one color's cheaper paths pruning another's.

namespace {

constexpr auto const kMaxMatchDistance = 100.0;
constexpr auto const kMax = std::chrono::seconds{1800};
constexpr auto const kNumDestinations = 120U;
constexpr auto const kNumColors = 3U;
constexpr auto const kNumStations = 24U;

constexpr auto const kCarAccessible = [](way_properties const& p) {
  return p.is_car_accessible();
};

node_idx_t random_node(ways const& w, std::mt19937& prng) {
  return node_idx_t{static_cast<std::uint32_t>(prng() % w.n_nodes())};
}

node_idx_t random_car_node(ways const& w, std::mt19937& prng) {
  while (true) {
    auto const n = random_node(w, prng);
    for (auto const way : w.r_->node_ways_[n]) {
      if (kCarAccessible(w.r_->way_properties_[way])) {
        return n;
      }
    }
  }
}

// Near side (where the search picks the vehicle up in search order): every
// street node. Far side: kNumColors overlapping colors, either longitude bands
// over the street nodes (flex areas) or additional nodes (location group
// stops), some of them with two colors.
struct colored_sharing {
  colored_sharing(ways const& w, std::mt19937& prng, bool const stations) {
    auto const n = static_cast<bitvec<node_idx_t>::size_type>(
        w.n_nodes() + (stations ? kNumStations : 0U));
    near_.resize(n);
    through_.resize(n);
    far_union_.resize(n);
    near_.zero_out();
    through_.one_out();
    far_union_.zero_out();
    for (auto i = 0U; i != w.n_nodes(); ++i) {
      near_.set(node_idx_t{i}, true);
    }
    for (auto& c : colors_) {
      c.resize(n);
      c.zero_out();
    }

    auto const color = [&](node_idx_t const x, unsigned const c) {
      colors_[c].set(x, true);
      far_union_.set(x, true);
    };

    if (stations) {
      for (auto s = 0U; s != kNumStations; ++s) {
        auto const street = random_car_node(w, prng);
        auto const a = node_idx_t{w.n_nodes() + s};
        coordinates_.push_back(w.get_node_pos(street).as_latlng());
        edges_[street].push_back({.to_ = a, .distance_ = 0U});
        edges_[a].push_back({.to_ = street, .distance_ = 0U});
        color(a, s % kNumColors);
        if (s % 4U == 0U) {
          color(a, (s + 1U) % kNumColors);
        }
      }
    } else {
      auto min_lng = 180.0;
      auto max_lng = -180.0;
      for (auto i = 0U; i != w.n_nodes(); ++i) {
        auto const lng = w.get_node_pos(node_idx_t{i}).lng();
        min_lng = std::min(min_lng, lng);
        max_lng = std::max(max_lng, lng);
      }
      // Bands [0, .4), (.3, .7), (.6, 1]: overlaps get two colors.
      for (auto i = 0U; i != w.n_nodes(); ++i) {
        auto const x = node_idx_t{i};
        auto const t =
            (w.get_node_pos(x).lng() - min_lng) / (max_lng - min_lng);
        if (t < 0.4) {
          color(x, 0U);
        }
        if (t > 0.3 && t < 0.7) {
          color(x, 1U);
        }
        if (t > 0.6) {
          color(x, 2U);
        }
      }
    }
    for (auto c = 0U; c != kNumColors; ++c) {
      color_ptrs_[c] = &colors_[c];
    }
  }

  sharing_data view(ways const& w,
                    direction const dir,
                    bitvec<node_idx_t> const* far,
                    bool const colored) const {
    auto s = sharing_data{
        .start_allowed_ = dir == direction::kForward ? &near_ : far,
        .end_allowed_ = dir == direction::kForward ? far : &near_,
        .through_allowed_ = &through_,
        .additional_node_offset_ = w.n_nodes(),
        .additional_node_coordinates_ = coordinates_,
        .additional_edges_ = edges_};
    if (colored) {
      s.colors_ = std::span{color_ptrs_};
    }
    return s;
  }

  bitvec<node_idx_t> near_, through_, far_union_;
  std::array<bitvec<node_idx_t>, kNumColors> colors_;
  std::array<bitvec<node_idx_t> const*, kNumColors> color_ptrs_{};
  std::vector<geo::latlng> coordinates_;
  hash_map<node_idx_t, std::vector<additional_edge>> edges_;
};

std::optional<std::pair<duration_t, duration_t>> ride_from_segments(
    path const& p) {
  auto pos = std::optional<std::pair<duration_t, duration_t>>{};
  auto at = duration_t{0U};
  for (auto const& s : p.segments_) {
    if (s.mode_ == mode::kCar) {
      if (!pos.has_value()) {
        pos = std::pair{at, at};
      }
      pos->second = at + s.duration_;
    }
    at += s.duration_;
  }
  return pos;
}

struct outcome {
  unsigned n_rides_{0U};  // colored (destination, color) rides
  unsigned n_slower_kept_{0U};  // ... slower than the uncolored best ride
};

outcome check(graph const& g,
              direction const dir,
              bool const stations,
              unsigned const seed) {
  auto const& w = *g.w_;
  auto const& l = *g.l_;
  auto const profile = search_profile::kCarSharing;
  auto const params =
      profile_parameters{car_sharing<track_node_tracking>::parameters{}};
  auto prng = std::mt19937{seed};
  auto const cs = colored_sharing{w, prng, stations};

  auto const from = location{w.get_node_pos(random_car_node(w, prng))};
  auto to = std::vector<location>{};
  for (auto i = 0U; i != kNumDestinations; ++i) {
    to.push_back(location{w.get_node_pos(random_node(w, prng))});
  }
  auto from_m = match_result{};
  l.match(params, from, false, dir, kMaxMatchDistance, nullptr, profile, {},
          from_m);
  auto to_m = match_result{};
  for (auto const& x : to) {
    l.match(params, x, true, dir, kMaxMatchDistance, nullptr, profile, {},
            to_m);
  }
  auto const from_match = from_m[match_idx_t{0U}];
  EXPECT_FALSE(from_match.empty());
  if (from_match.empty()) {
    return {};
  }

  auto const route = [&](sharing_data const& s) {
    return route_one_to_many(params, w, l, profile, from, to, from_match, to_m,
                             kMax, dir, nullptr, &s);
  };

  auto const colored_sharing = cs.view(w, dir, &cs.far_union_, true);
  auto const colored = route(colored_sharing);
  EXPECT_EQ(to.size() * kNumColors, colored->results().size());

  auto const union_sharing = cs.view(w, dir, &cs.far_union_, false);
  auto const uncolored = route(union_sharing);

  auto out = outcome{};
  for (auto c = 0U; c != kNumColors; ++c) {
    SCOPED_TRACE("color " + std::to_string(c));
    auto const single_sharing = cs.view(w, dir, &cs.colors_[c], false);
    auto const single = route(single_sharing);

    for (auto i = 0U; i != to.size(); ++i) {
      auto const k = i * kNumColors + c;
      auto const& r = colored->results()[k];
      auto const rc = colored->rental_durations(k);

      // Colored results are rides only.
      EXPECT_EQ(r.has_value(), rc.has_value()) << i;

      // Same ride as the search restricted to this color. If that search
      // walks (cheaper than any ride), the colored ride costs at least as much.
      auto const& s = single->results()[i];
      if (single->rental_durations(i).has_value()) {
        EXPECT_TRUE(r.has_value()) << i;
        if (r.has_value()) {
          EXPECT_EQ(s->cost_, r->cost_) << i;
          EXPECT_EQ(s->duration_, r->duration_) << i;
        }
      } else if (r.has_value()) {
        EXPECT_TRUE(s.has_value()) << i;
        if (s.has_value()) {
          EXPECT_LE(s->cost_, r->cost_) << i;
        }
      }
      if (!r.has_value()) {
        continue;
      }
      ++out.n_rides_;

      // The uncolored search over the union keeps only its cheapest path.
      auto const& u = uncolored->results()[i];
      EXPECT_TRUE(u.has_value()) << i;
      if (u.has_value()) {
        EXPECT_LE(u->cost_, r->cost_) << i;
        if (u->cost_ < r->cost_) {
          ++out.n_slower_kept_;
        }
      }

      // Reconstruction of a colored result: same duration, same ride.
      auto const p = colored->reconstruct(w, l, k, &colored_sharing);
      EXPECT_TRUE(p.has_value()) << i;
      if (p.has_value()) {
        EXPECT_EQ(r->duration_, p->duration_) << i;
        auto const ride = ride_from_segments(*p);
        EXPECT_TRUE(ride.has_value()) << i;
        if (ride.has_value()) {
          EXPECT_EQ(ride->first, rc->ride_start_) << i;
          EXPECT_EQ(ride->second, rc->ride_end_) << i;
        }
      }
    }
  }

  // The best uncolored ride is the best ride over all colors.
  for (auto i = 0U; i != to.size(); ++i) {
    if (!uncolored->rental_durations(i).has_value()) {
      continue;
    }
    auto best = std::optional<cost_t>{};
    for (auto c = 0U; c != kNumColors; ++c) {
      if (auto const& r = colored->results()[i * kNumColors + c];
          r.has_value()) {
        best = std::min(best.value_or(r->cost_), r->cost_);
      }
    }
    EXPECT_EQ(uncolored->results()[i]->cost_, best) << i;
  }
  return out;
}

void check_both_directions(bool const stations, unsigned const seed) {
  auto const g = graph{"test/miraustr.osm.pbf", "osr_sharing_colors_test"};
  for (auto const dir : {direction::kForward, direction::kBackward}) {
    SCOPED_TRACE(std::string{to_str(dir)});
    auto out = outcome{};
    for (auto i = 0U; i != 3U; ++i) {
      auto const o = check(g, dir, stations, seed + i);
      out.n_rides_ += o.n_rides_;
      out.n_slower_kept_ += o.n_slower_kept_;
    }
    // The case the colors exist for has to occur.
    EXPECT_GT(out.n_rides_, 20U);
    EXPECT_GT(out.n_slower_kept_, 5U);
  }
}

}  // namespace

TEST(sharing_colors, car_sharing_areas) { check_both_directions(false, 11U); }

TEST(sharing_colors, car_sharing_stations) { check_both_directions(true, 21U); }
