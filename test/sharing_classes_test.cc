#include "gtest/gtest.h"

#include <chrono>
#include <random>
#include <string>
#include <vector>

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

// sharing_data::classes_ splits the walk after the vehicle (search direction)
// into classes, so that one car_sharing search replaces one search per class
// that differ only in where the vehicle may be left (forward) or taken
// (backward). The classed search has to report, for every class and
// destination, exactly what the search for that class alone reports.

namespace {

constexpr auto const kMaxMatchDistance = 100.0;
constexpr auto const kMax = std::chrono::seconds{1800};
constexpr auto const kNumDestinations = 200U;
constexpr auto const kNumNearNodes = 30U;
constexpr auto const kNumClasses = 4U;
constexpr auto const kNodesPerClass = 15U;
constexpr auto const kSharedNodes = 5U;  // in classes 0 and 1

using profile_t = car_sharing<track_node_tracking>;

node_idx_t random_car_node(ways const& w, std::mt19937& prng) {
  while (true) {
    auto const n = node_idx_t{static_cast<std::uint32_t>(prng() % w.n_nodes())};
    for (auto const way : w.r_->node_ways_[n]) {
      if (w.r_->way_properties_[way].is_car_accessible()) {
        return n;
      }
    }
  }
}

// The near side (where the search's own vehicle part begins) is a set of
// street nodes; the far side has kNumClasses classes of street nodes or - with
// `stations` - of additional nodes, each tied to a street node.
struct setup {
  setup(ways const& w,
        std::mt19937& prng,
        direction const dir,
        bool const stations)
      : dir_{dir}, n_nodes_{w.n_nodes()} {
    auto const size = static_cast<bitvec<node_idx_t>::size_type>(
        w.n_nodes() +
        (stations ? kNumClasses * kNodesPerClass + kSharedNodes : 0U));
    near_.resize(size);
    far_union_.resize(size);
    through_.resize(size);
    near_.zero_out();
    far_union_.zero_out();
    through_.one_out();
    far_.resize(kNumClasses);
    for (auto& f : far_) {
      f.resize(size);
      f.zero_out();
    }

    for (auto i = 0U; i != kNumNearNodes; ++i) {
      near_.set(random_car_node(w, prng), true);
    }

    auto next_station = node_idx_t{w.n_nodes()};
    auto const add_far = [&](std::vector<unsigned> const& classes) {
      auto n = random_car_node(w, prng);
      if (stations) {
        auto const a = next_station++;
        coordinates_.push_back(w.get_node_pos(n).as_latlng());
        edges_[n].push_back({.to_ = a, .distance_ = 0U});
        edges_[a].push_back({.to_ = n, .distance_ = 0U});
        n = a;
      }
      far_union_.set(n, true);
      for (auto const c : classes) {
        far_[c].set(n, true);
        classes_[n] |= sharing_class_mask_t{1U} << c;
      }
    };
    for (auto c = 0U; c != kNumClasses; ++c) {
      for (auto i = 0U; i != kNodesPerClass; ++i) {
        add_far({c});
      }
    }
    for (auto i = 0U; i != kSharedNodes; ++i) {
      add_far({0U, 1U});
    }
  }

  // start / end allowed in travel order: the far side is where the vehicle is
  // left (forward search) or taken (backward search).
  sharing_data view(bitvec<node_idx_t> const* far) const {
    auto const forward = dir_ == direction::kForward;
    return {.start_allowed_ = forward ? &near_ : far,
            .end_allowed_ = forward ? far : &near_,
            .through_allowed_ = &through_,
            .additional_node_offset_ = n_nodes_,
            .additional_node_coordinates_ = coordinates_,
            .additional_edges_ = edges_};
  }

  sharing_data classed() const {
    auto s = view(&far_union_);
    s.classes_ = &classes_;
    s.n_classes_ = kNumClasses;
    return s;
  }

  sharing_data single(unsigned const c) const { return view(&far_[c]); }

  // What a caller keeps for reconstruction after the search: additional
  // nodes and classes, no start / end sets.
  sharing_data for_reconstruction() const {
    return {.additional_node_offset_ = n_nodes_,
            .additional_node_coordinates_ = coordinates_,
            .additional_edges_ = edges_,
            .classes_ = &classes_,
            .n_classes_ = kNumClasses};
  }

  direction dir_;
  node_idx_t::value_t n_nodes_;
  bitvec<node_idx_t> near_, far_union_, through_;
  std::vector<bitvec<node_idx_t>> far_;
  hash_map<node_idx_t, sharing_class_mask_t> classes_;
  std::vector<geo::latlng> coordinates_;
  hash_map<node_idx_t, std::vector<additional_edge>> edges_;
};

struct ride_position {
  duration_t ride_start_, ride_end_;
};

std::optional<ride_position> ride_from_segments(path const& p) {
  auto pos = std::optional<ride_position>{};
  auto at = duration_t{0U};
  for (auto const& s : p.segments_) {
    if (s.mode_ == mode::kCar) {
      if (!pos.has_value()) {
        pos = ride_position{.ride_start_ = at, .ride_end_ = at};
      }
      pos->ride_end_ = at + s.duration_;
    }
    at += s.duration_;
  }
  return pos;
}

struct outcome {
  unsigned n_rides_{0U}, n_class_differs_{0U};
};

outcome check(graph const& g,
              direction const dir,
              bool const stations,
              unsigned const seed) {
  auto const& w = *g.w_;
  auto const& l = *g.l_;
  auto const params = profile_parameters{profile_t::parameters{}};
  auto prng = std::mt19937{seed};
  auto const s = setup{w, prng, dir, stations};

  auto const from = location{w.get_node_pos(random_car_node(w, prng))};
  auto to = std::vector<location>{};
  for (auto i = 0U; i != kNumDestinations; ++i) {
    to.push_back(location{w.get_node_pos(random_car_node(w, prng))});
  }
  auto from_m = match_result{};
  l.match(params, from, false, dir, kMaxMatchDistance, nullptr,
          search_profile::kCarSharing, {}, from_m);
  auto to_m = match_result{};
  for (auto const& x : to) {
    l.match(params, x, true, dir, kMaxMatchDistance, nullptr,
            search_profile::kCarSharing, {}, to_m);
  }
  auto const from_match = from_m[match_idx_t{0U}];
  EXPECT_FALSE(from_match.empty());

  auto const route = [&](sharing_data const& sharing) {
    return route_one_to_many(params, w, l, search_profile::kCarSharing, from,
                             to, from_match, to_m, kMax, dir, nullptr,
                             &sharing);
  };

  auto const classed_sharing = s.classed();
  auto const classed = route(classed_sharing);
  EXPECT_EQ(kNumClasses * to.size(), classed->results().size());

  auto out = outcome{};
  auto const reconstruction_sharing = s.for_reconstruction();
  for (auto c = 0U; c != kNumClasses; ++c) {
    SCOPED_TRACE("class " + std::to_string(c));
    auto const single_sharing = s.single(c);
    auto const single = route(single_sharing);
    EXPECT_EQ(to.size(), single->results().size());
    for (auto i = 0U; i != to.size(); ++i) {
      auto const k = c * to.size() + i;
      auto const& expected = single->results()[i];
      auto const& actual = classed->results()[k];
      EXPECT_EQ(expected.has_value(), actual.has_value()) << i;
      if (!expected.has_value() || !actual.has_value()) {
        continue;
      }
      EXPECT_EQ(expected->cost_, actual->cost_) << i;
      EXPECT_EQ(expected->duration_, actual->duration_) << i;

      auto const expected_ride = single->rental_durations(i);
      auto const actual_ride = classed->rental_durations(k);
      EXPECT_EQ(expected_ride.has_value(), actual_ride.has_value()) << i;
      if (!expected_ride.has_value() || !actual_ride.has_value()) {
        continue;
      }
      ++out.n_rides_;
      EXPECT_EQ(expected_ride->ride_start_, actual_ride->ride_start_) << i;
      EXPECT_EQ(expected_ride->ride_end_, actual_ride->ride_end_) << i;

      // The classed path can be laid out again from additional nodes and
      // classes alone, and shows the ride where the search tracked it.
      auto const p = classed->reconstruct(w, l, k, &reconstruction_sharing);
      EXPECT_TRUE(p.has_value()) << i;
      if (!p.has_value()) {
        continue;
      }
      EXPECT_EQ(actual->duration_, p->duration_) << i;
      auto const ride = ride_from_segments(*p);
      EXPECT_TRUE(ride.has_value()) << i;
      if (ride.has_value()) {
        EXPECT_EQ(actual_ride->ride_start_, ride->ride_start_) << i;
        EXPECT_EQ(actual_ride->ride_end_, ride->ride_end_) << i;
      }

      if (c != 0U) {
        auto const& class_0 = classed->results()[i];
        if (!class_0.has_value() || class_0->cost_ != actual->cost_) {
          ++out.n_class_differs_;
        }
      }
    }
  }
  return out;
}

void check_both_directions(bool const stations) {
  auto const g = graph{"test/miraustr.osm.pbf", "osr_sharing_classes_test"};
  for (auto const dir : {direction::kForward, direction::kBackward}) {
    SCOPED_TRACE(std::string{to_str(dir)});
    auto out = outcome{};
    for (auto seed = 0U; seed != 4U; ++seed) {
      auto const o = check(g, dir, stations, 10U + seed);
      out.n_rides_ += o.n_rides_;
      out.n_class_differs_ += o.n_class_differs_;
    }
    // Rides have to occur, and the classes have to matter.
    EXPECT_GT(out.n_rides_, 50U);
    EXPECT_GT(out.n_class_differs_, 10U);
  }
}

}  // namespace

TEST(sharing_classes, street_nodes) { check_both_directions(false); }

TEST(sharing_classes, stations) { check_both_directions(true); }
