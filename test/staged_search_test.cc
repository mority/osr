#include "gtest/gtest.h"

#include <chrono>
#include <filesystem>
#include <iostream>
#include <memory>
#include <random>

#include "osr/extract/extract.h"
#include "osr/lookup.h"
#include "osr/routing/parameters.h"
#include "osr/routing/profiles/bike.h"
#include "osr/routing/profiles/bike_sharing.h"
#include "osr/routing/profiles/car.h"
#include "osr/routing/profiles/car_sharing.h"
#include "osr/routing/profiles/foot.h"
#include "osr/routing/route.h"
#include "osr/routing/sharing_data.h"
#include "osr/routing/staged_search.h"
#include "osr/ways.h"

namespace fs = std::filesystem;
using namespace osr;

namespace {

constexpr auto const kMaxMatchDistance = 100.0;
constexpr auto const kMax = cost_t{1800U};
constexpr auto const kNumDestinations = 200U;
constexpr auto const kNumStations = 40U;

using footp = foot<false>;
using bikep = bike<bike_costing::kSafe, kElevationNoCost>;
using carp = car;
using car_sharingp = car_sharing<track_node_tracking>;

struct staged_search_test : public ::testing::Test {
  static void SetUpTestSuite() {
    dir_ = fs::temp_directory_path() / "osr-staged-search-test";
    auto ec = std::error_code{};
    fs::remove_all(dir_, ec);
    fs::create_directories(dir_, ec);
    extract(false, "test/miraustr.osm.pbf", dir_, {});
    w_ = std::make_unique<ways>(dir_, cista::mmap::protection::READ);
    l_ = std::make_unique<lookup>(*w_, dir_, cista::mmap::protection::READ);
  }

  static void TearDownTestSuite() {
    l_.reset();
    w_.reset();
    auto ec = std::error_code{};
    fs::remove_all(dir_, ec);
  }

  static inline fs::path dir_{};
  static inline std::unique_ptr<ways> w_{};
  static inline std::unique_ptr<lookup> l_{};
};

// One rental product: `n_stations` stations at random nodes (additional node
// attached to a real node with a short edge in both directions), optionally a
// free-floating zone (random real nodes where picking up / dropping off is
// allowed directly). Zone nodes are restricted to nodes that are accessible
// for every mode so that the node feasibility check of the staged search does
// not make a difference (see the header comment of staged_search.h).
struct product {
  product(ways const& w,
          std::mt19937& prng,
          unsigned const n_stations,
          double const zone_start_fraction,
          double const zone_end_fraction)
      : offset_{w.n_nodes()} {
    auto const size =
        static_cast<bitvec<node_idx_t>::size_type>(w.n_nodes() + n_stations);
    start_allowed_.resize(size);
    end_allowed_.resize(size);
    through_allowed_.resize(size);
    through_allowed_.one_out();

    auto distr =
        std::uniform_int_distribution<std::uint32_t>{0, w.n_nodes() - 1};
    auto unit = std::uniform_real_distribution<double>{0.0, 1.0};

    for (auto i = 0U; i != n_stations; ++i) {
      auto const s = node_idx_t{distr(prng)};
      auto const a = node_idx_t{offset_ + i};
      coords_.push_back(w.get_node_pos(s).as_latlng());
      edges_[s].push_back(additional_edge{.to_ = a, .distance_ = 10U});
      edges_[a].push_back(additional_edge{.to_ = s, .distance_ = 10U});
      start_allowed_.set(a, true);
      end_allowed_.set(a, true);
    }

    if (zone_start_fraction > 0.0 || zone_end_fraction > 0.0) {
      for (auto n = node_idx_t{0U}; n != node_idx_t{w.n_nodes()}; ++n) {
        auto const p = w.r_->node_properties_[n];
        if (!p.is_walk_accessible() || !p.is_bike_accessible() ||
            !p.is_car_accessible()) {
          continue;
        }
        if (unit(prng) < zone_start_fraction) {
          start_allowed_.set(n, true);
        }
        if (unit(prng) < zone_end_fraction) {
          end_allowed_.set(n, true);
        }
      }
    }
  }

  sharing_data view() const {
    verify_additional_edge_count(edges_, offset_);
    return {.start_allowed_ = &start_allowed_,
            .end_allowed_ = &end_allowed_,
            .through_allowed_ = &through_allowed_,
            .additional_node_offset_ = offset_,
            .additional_node_coordinates_ = coords_,
            .additional_edges_ = edges_};
  }

  node_idx_t::value_t offset_;
  bitvec<node_idx_t> start_allowed_{};
  bitvec<node_idx_t> end_allowed_{};
  bitvec<node_idx_t> through_allowed_{};
  std::vector<geo::latlng> coords_{};
  hash_map<node_idx_t, std::vector<additional_edge>> edges_{};
};

struct query {
  query(ways const& w, unsigned const seed) {
    auto prng = std::mt19937{seed};
    auto distr =
        std::uniform_int_distribution<std::uint32_t>{0, w.n_nodes() - 1};
    from_ = location{w.get_node_pos(node_idx_t{distr(prng)})};
    for (auto i = 0U; i != kNumDestinations; ++i) {
      to_.push_back(location{w.get_node_pos(node_idx_t{distr(prng)})});
    }
  }

  // Both searches get the same matches (the joint profiles match like foot).
  template <typename P>
  void match(lookup const& l, direction const dir) {
    auto const pp = typename P::parameters{};
    from_m_ = match_result{};
    l.match<P>(pp, from_, false, dir, kMaxMatchDistance, nullptr, from_m_,
               std::nullopt);
    to_m_ = match_result{};
    for (auto const& x : to_) {
      l.match<P>(pp, x, true, dir, kMaxMatchDistance, nullptr, to_m_,
                 std::nullopt);
    }
  }

  match_view_t from_match() const { return from_m_[match_idx_t{0U}]; }

  location from_;
  std::vector<location> to_;
  match_result from_m_;
  match_result to_m_;
};

template <typename Sharing>
std::vector<std::optional<path>> run_joint(ways const& w,
                                           lookup const& l,
                                           search_profile const profile,
                                           query const& q,
                                           sharing_data const& sharing,
                                           direction const dir) {
  auto const params = profile_parameters{typename Sharing::parameters{}};
  return route_one_to_many(params, w, l, profile, q.from_, q.to_,
                           q.from_match(), q.to_m_, kMax, dir, nullptr,
                           &sharing)
      ->results();
}

template <typename Vehicle>
std::unique_ptr<staged_search> run_staged(
    ways const& w,
    query const& q,
    std::vector<sharing_data> const& products,
    direction const dir,
    bool const direct_pickup) {
  auto s = std::make_unique<staged_search>();
  auto const first = s->add_stage<footp>({}, nullptr, true);
  auto vehicles = std::vector<std::size_t>{};
  for (auto const& p : products) {
    vehicles.push_back(s->add_stage<Vehicle>({}, p.through_allowed_));
  }
  auto const last = s->add_stage<footp>({}, nullptr, true);
  for (auto i = 0U; i != products.size(); ++i) {
    add_rental_transitions(*s, first, vehicles[i], last, products[i], dir,
                           direct_pickup);
  }
  s->run(w, q.from_, q.to_, q.from_match(), q.to_m_, kMax, dir);
  return s;
}

cost_t segment_cost_sum(path const& p) {
  auto sum = cost_t{0U};
  for (auto const& s : p.segments_) {
    sum += s.cost_;
  }
  return sum;
}

void expect_contiguous(path const& p) {
  ASSERT_FALSE(p.segments_.empty());
  for (auto i = 1U; i < p.segments_.size(); ++i) {
    if (p.segments_[i - 1U].to_ != node_idx_t::invalid() &&
        p.segments_[i].from_ != node_idx_t::invalid()) {
      EXPECT_EQ(p.segments_[i - 1U].to_, p.segments_[i].from_) << i;
    }
    EXPECT_FALSE(p.segments_[i].polyline_.empty()) << i;
  }
}

// Compares per destination: same reachability, same cost, and a consistent
// reconstructed path.
void compare(std::vector<std::optional<path>> const& expected,
             staged_search const& s,
             lookup const& l,
             unsigned& n_found) {
  auto const& got = s.results();
  ASSERT_EQ(expected.size(), got.size());
  for (auto k = 0U; k != expected.size(); ++k) {
    ASSERT_EQ(expected[k].has_value(), got[k].has_value()) << k;
    if (!got[k].has_value()) {
      continue;
    }
    ++n_found;
    EXPECT_EQ(expected[k]->cost_, got[k]->cost_) << k;

    auto const p = s.reconstruct(l, k);
    ASSERT_TRUE(p.has_value()) << k;
    EXPECT_EQ(got[k]->cost_, p->cost_) << k;
    // Later start candidates can improve the search state after the result
    // was recorded, so the rendered segments may be cheaper, never dearer.
    EXPECT_LE(segment_cost_sum(*p), p->cost_) << k;
    expect_contiguous(*p);
  }
}

}  // namespace

template <typename Sharing, typename Vehicle>
void check_single_product(ways const& w,
                          lookup const& l,
                          search_profile const profile,
                          bool const direct_pickup,
                          double const zone_start_fraction,
                          double const zone_end_fraction) {
  for (auto const dir : {direction::kForward, direction::kBackward}) {
    for (auto seed = 1U; seed != 6U; ++seed) {
      auto prng = std::mt19937{seed};
      auto const prod = product{w, prng, kNumStations, zone_start_fraction,
                                zone_end_fraction};
      auto const sharings = std::vector<sharing_data>{prod.view()};
      auto q = query{w, seed};
      q.match<Sharing>(l, dir);
      ASSERT_FALSE(q.from_match().empty());

      auto const expected =
          run_joint<Sharing>(w, l, profile, q, sharings.front(), dir);
      auto const s = run_staged<Vehicle>(w, q, sharings, dir, direct_pickup);

      auto n_found = 0U;
      compare(expected, *s, l, n_found);
      EXPECT_GT(n_found, 10U) << to_str(dir) << " seed=" << seed;
    }
  }
}

TEST_F(staged_search_test, bike_stations_only_matches_bike_sharing) {
  check_single_product<bike_sharing, bikep>(
      *w_, *l_, search_profile::kBikeSharing, false, 0.0, 0.0);
}

TEST_F(staged_search_test, bike_with_dropoff_zone_matches_bike_sharing) {
  check_single_product<bike_sharing, bikep>(
      *w_, *l_, search_profile::kBikeSharing, false, 0.0, 0.3);
}

TEST_F(staged_search_test, car_zone_pickup_matches_car_sharing) {
  check_single_product<car_sharingp, carp>(
      *w_, *l_, search_profile::kCarSharing, true, 0.3, 0.0);
}

TEST_F(staged_search_test, car_zone_pickup_and_dropoff_matches_car_sharing) {
  check_single_product<car_sharingp, carp>(
      *w_, *l_, search_profile::kCarSharing, true, 0.3, 0.3);
}

// Several products in one staged search (shared foot stages) give the best
// product per destination, i.e. the minimum over the individual searches.
TEST_F(staged_search_test, multiple_products_give_minimum) {
  auto const& w = *w_;
  auto const& l = *l_;
  constexpr auto const kNumProducts = 4U;
  for (auto const dir : {direction::kForward, direction::kBackward}) {
    for (auto seed = 1U; seed != 4U; ++seed) {
      auto prng = std::mt19937{seed};
      auto products = std::vector<product>{};
      for (auto i = 0U; i != kNumProducts; ++i) {
        products.emplace_back(w, prng, kNumStations / 2U, 0.0, 0.0);
      }
      auto sharings = std::vector<sharing_data>{};
      for (auto const& p : products) {
        sharings.push_back(p.view());
      }
      auto q = query{w, seed};
      q.match<bike_sharing>(l, dir);
      ASSERT_FALSE(q.from_match().empty());

      auto t0 = std::chrono::steady_clock::now();
      auto expected = std::vector<std::optional<path>>(q.to_.size());
      for (auto const& sharing : sharings) {
        auto const r = run_joint<bike_sharing>(
            w, l, search_profile::kBikeSharing, q, sharing, dir);
        for (auto k = 0U; k != r.size(); ++k) {
          if (r[k].has_value() &&
              (!expected[k].has_value() || r[k]->cost_ < expected[k]->cost_)) {
            expected[k] = r[k];
          }
        }
      }
      auto t1 = std::chrono::steady_clock::now();
      auto const s = run_staged<bikep>(w, q, sharings, dir, false);
      auto t2 = std::chrono::steady_clock::now();

      auto n_found = 0U;
      compare(expected, *s, l, n_found);
      EXPECT_GT(n_found, 10U) << to_str(dir) << " seed=" << seed;

      using ms = std::chrono::duration<double, std::milli>;
      std::cout << to_str(dir) << " seed=" << seed << ": " << kNumProducts
                << " joint searches "
                << std::chrono::duration_cast<ms>(t1 - t0).count()
                << " ms, staged "
                << std::chrono::duration_cast<ms>(t2 - t1).count() << " ms\n";
    }
  }
}
