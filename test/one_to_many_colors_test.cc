#include "gtest/gtest.h"

#include <chrono>
#include <cmath>
#include <filesystem>
#include <memory>
#include <random>
#include <vector>

#include "osr/extract/extract.h"
#include "osr/lookup.h"
#include "osr/routing/cost_search_limit.h"
#include "osr/routing/parameters.h"
#include "osr/routing/profiles/bike.h"
#include "osr/routing/profiles/foot.h"
#include "osr/routing/route.h"
#include "osr/ways.h"

namespace fs = std::filesystem;
using namespace osr;

// Colors: once a member of a color is found within the limits, the other
// members are searched past them, each up to the cheapest found member's cost
// plus the search limit for the straight-line gap between them.

namespace {

constexpr auto const kMaxMatchDistance = 100.0;
constexpr auto const kMax = std::chrono::seconds{300};
constexpr auto const kReferenceMax = std::chrono::seconds{3600};
constexpr auto const kNumDestinations = 400U;

struct graph {
  graph() {
    dir_ = fs::temp_directory_path() / "osr_one_to_many_colors_test";
    auto ec = std::error_code{};
    fs::remove_all(dir_, ec);
    fs::create_directories(dir_, ec);
    extract(false, "test/miraustr.osm.pbf", dir_, {});
    w_ = std::make_unique<ways>(dir_, cista::mmap::protection::READ);
    l_ = std::make_unique<lookup>(*w_, dir_, cista::mmap::protection::READ);
  }

  fs::path dir_;
  std::unique_ptr<ways> w_;
  std::unique_ptr<lookup> l_;
};

template <typename P>
cost_t budget(typename P::parameters const& params,
              profile_parameters const& pp,
              location const& found,
              location const& member) {
  auto const gap =
      P::upper_bound_heuristic(params, geo::distance(found.pos_, member.pos_));
  return cost_search_limit(
      pp, clamp_duration(static_cast<std::uint64_t>(std::ceil(gap))));
}

template <typename P>
void check(graph const& g,
           search_profile const profile,
           direction const dir,
           unsigned const seed,
           bool const expect_over_budget) {
  auto const& w = *g.w_;
  auto const& l = *g.l_;
  auto const params = profile_parameters{typename P::parameters{}};
  auto const& pp = std::get<typename P::parameters>(params);

  auto prng = std::mt19937{seed};
  auto distr = std::uniform_int_distribution<std::uint32_t>{0, w.n_nodes() - 1};
  auto const from = location{w.get_node_pos(node_idx_t{distr(prng)})};
  auto to = std::vector<location>{};
  for (auto i = 0U; i != kNumDestinations; ++i) {
    to.push_back(location{w.get_node_pos(node_idx_t{distr(prng)})});
  }

  auto from_m = match_result{};
  l.match<P>(pp, from, false, dir, kMaxMatchDistance, nullptr, from_m,
             std::nullopt);
  auto to_m = match_result{};
  for (auto const& x : to) {
    l.match<P>(pp, x, true, dir, kMaxMatchDistance, nullptr, to_m,
               std::nullopt);
  }
  auto const from_match = from_m[match_idx_t{0U}];
  ASSERT_FALSE(from_match.empty());

  auto const run = [&](std::chrono::seconds const max,
                       std::vector<color_idx_t> const& colors,
                       bool const reconstruct) {
    return route_one_to_many(
        params, w, l, profile, from, to, from_match, to_m, max, dir, nullptr,
        nullptr, nullptr, [&](path const&) { return reconstruct; },
        std::nullopt, colors);
  };

  auto const today = run(kMax, {}, false)->results();
  auto const reference = run(kReferenceMax, {}, false)->results();

  // No colors, or only colors with a single member: today's results.
  auto singletons = std::vector<color_idx_t>{};
  for (auto i = 0U; i != to.size(); ++i) {
    singletons.push_back(color_idx_t{i});
  }
  for (auto const& colors :
       {std::vector<color_idx_t>(to.size(), color_idx_t::invalid()),
        singletons}) {
    auto const res = run(kMax, colors, false)->results();
    for (auto i = 0U; i != to.size(); ++i) {
      ASSERT_EQ(today[i].has_value(), res[i].has_value()) << i;
      if (today[i].has_value()) {
        EXPECT_EQ(today[i]->cost_, res[i]->cost_) << i;
        EXPECT_EQ(today[i]->duration_, res[i]->duration_) << i;
      }
    }
  }

  // Pair each target found today with one that isn't; colors of two targets
  // that aren't found stay without results.
  auto found = std::vector<std::size_t>{};
  auto missing = std::vector<std::size_t>{};
  for (auto i = 0U; i != to.size(); ++i) {
    (today[i].has_value() ? found : missing).push_back(i);
  }
  ASSERT_GT(found.size(), 5U);
  ASSERT_GT(missing.size(), 5U);

  auto colors = std::vector<color_idx_t>(to.size(), color_idx_t::invalid());
  auto partner = std::vector<std::size_t>(to.size(), to.size());
  auto next_color = 0U;
  auto m = std::size_t{0U};
  for (auto const f : found) {
    if (m == missing.size()) {
      break;
    }
    colors[f] = colors[missing[m]] = color_idx_t{next_color++};
    partner[missing[m]] = f;
    ++m;
  }
  for (; m + 1U < missing.size(); m += 2U) {
    colors[missing[m]] = colors[missing[m + 1U]] = color_idx_t{next_color++};
  }

  auto const immediate = run(kMax, colors, true)->results();
  auto const state = run(kMax, colors, false);
  auto const& res = state->results();

  auto n_extended = 0U, n_over_budget = 0U;
  for (auto i = 0U; i != to.size(); ++i) {
    ASSERT_EQ(res[i].has_value(), immediate[i].has_value()) << i;
    if (today[i].has_value()) {  // results within the limits don't change
      ASSERT_TRUE(res[i].has_value()) << i;
      EXPECT_EQ(today[i]->cost_, res[i]->cost_) << i;
      EXPECT_EQ(today[i]->duration_, res[i]->duration_) << i;
      continue;
    }
    if (partner[i] == to.size()) {  // no member of its color found
      EXPECT_FALSE(res[i].has_value()) << i;
      continue;
    }
    auto const limit =
        clamp_cost(static_cast<std::uint64_t>(today[partner[i]]->cost_) +
                   budget<P>(pp, params, to[partner[i]], to[i]));
    if (!reference[i].has_value() || reference[i]->cost_ > limit) {
      EXPECT_FALSE(res[i].has_value()) << i;
      n_over_budget += reference[i].has_value();
      continue;
    }
    ASSERT_TRUE(res[i].has_value()) << i;
    EXPECT_EQ(reference[i]->cost_, res[i]->cost_) << i;
    EXPECT_EQ(reference[i]->duration_, res[i]->duration_) << i;

    // The second phase is part of the stored search.
    auto const p = state->reconstruct(w, l, i, nullptr);
    ASSERT_TRUE(p.has_value()) << i;
    EXPECT_EQ(res[i]->cost_, p->cost_) << i;
    EXPECT_EQ(immediate[i]->cost_, p->cost_) << i;
    EXPECT_EQ(immediate[i]->segments_.size(), p->segments_.size()) << i;
    EXPECT_FALSE(p->segments_.empty()) << i;
    ++n_extended;
  }
  EXPECT_GT(n_extended, 5U);
  if (expect_over_budget) {
    EXPECT_GT(n_over_budget, 0U);
  }
}

}  // namespace

TEST(one_to_many_colors, foot_forward) {
  auto const g = graph{};
  check<foot<false, elevator_tracking>>(g, search_profile::kFoot,
                                        direction::kForward, 1U, false);
}

TEST(one_to_many_colors, foot_backward) {
  auto const g = graph{};
  check<foot<false, elevator_tracking>>(g, search_profile::kFoot,
                                        direction::kBackward, 2U, false);
}

TEST(one_to_many_colors, bike_forward) {
  auto const g = graph{};
  check<bike<bike_costing::kSafe, kElevationNoCost>>(
      g, search_profile::kBike, direction::kForward, 3U, true);
}
