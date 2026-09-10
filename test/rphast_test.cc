#include "gtest/gtest.h"

#include <chrono>
#include <filesystem>
#include <random>
#include <vector>

#include "cista/mmap.h"

#include "fmt/core.h"

#include "osr/extract/extract.h"
#include "osr/location.h"
#include "osr/lookup.h"
#include "osr/routing/algorithms.h"
#include "osr/routing/cch/cch.h"
#include "osr/routing/profiles/car.h"
#include "osr/routing/route.h"
#include "osr/types.h"
#include "osr/ways.h"

namespace fs = std::filesystem;
using namespace osr;

namespace {

constexpr auto const kMaxMatchDistance = 100.0;

void load(std::string_view raw_data, std::string_view data_dir) {
  if (fs::exists(raw_data)) {
    auto const p = fs::path{data_dir};
    auto ec = std::error_code{};
    fs::remove_all(p, ec);
    fs::create_directories(p, ec);
    osr::extract(false, raw_data, data_dir, fs::path{});
  }
}

// Both algorithms are handed the same matched candidates, so any difference is
// the routing and not the matching. The Dijkstra is the oracle: it is the
// one-to-many that production used before the hierarchy.
//
// The `kBackward` cases matter: `cch_supported` admits backward one to many
// queries, so they compare lazy RPHAST over the reversed graph against the
// Dijkstra rather than the Dijkstra against itself.
void run(ways const& w,
         lookup const& l,
         unsigned const n_sources,
         unsigned const n_targets,
         cost_t const max_cost,
         direction const dir) {
  auto prng = std::mt19937{42U};
  auto distr =
      std::uniform_int_distribution<std::uint32_t>{0U, w.n_nodes() - 1U};
  auto const pos = [&](node_idx_t const n) {
    return location{w.get_node_pos(n).as_latlng(), level_t{}};
  };

  auto n_src = 0U, n_cmp = 0U, n_agree = 0U, n_reachable = 0U,
       n_dur_differ = 0U;

  // Timed per arm: a speedup of one means the hierarchy is not being used at
  // all -- which is what a `cch_supported` guard rejecting the query looks
  // like, since `route` then answers both arms with the Dijkstra.
  auto reference_time = std::chrono::steady_clock::duration{};
  auto experiment_time = std::chrono::steady_clock::duration{};

  for (auto i = 0U; i != n_sources; ++i) {
    auto const from_loc = pos(node_idx_t{distr(prng)});
    auto to_locs = std::vector<location>{};
    to_locs.reserve(n_targets);
    for (auto k = 0U; k != n_targets; ++k) {
      to_locs.push_back(pos(node_idx_t{distr(prng)}));
    }

    auto from_m = match_result{};
    l.match<car>(car::parameters{}, from_loc, false, dir, kMaxMatchDistance,
                 nullptr, from_m);
    if (from_m.empty()) {
      continue;
    }
    auto const from_match = from_m[match_idx_t{0U}];
    if (from_match.empty()) {
      continue;
    }

    auto to_match = match_result{};
    for (auto const& x : to_locs) {
      l.match<car>(car::parameters{}, x, true, dir, kMaxMatchDistance, nullptr,
                   to_match);
    }
    ASSERT_EQ(to_locs.size(), to_match.size());

    ++n_src;

    auto const t0 = std::chrono::steady_clock::now();
    auto const reference =
        route(car::parameters{}, w, l, search_profile::kCar, from_loc, to_locs,
              from_match, to_match, max_cost, dir, nullptr, nullptr, nullptr,
              [](path const&) { return false; }, std::nullopt,
              routing_algorithm::kDijkstra);
    auto const t1 = std::chrono::steady_clock::now();
    auto const experiment =
        route(car::parameters{}, w, l, search_profile::kCar, from_loc, to_locs,
              from_match, to_match, max_cost, dir, nullptr, nullptr, nullptr,
              [](path const&) { return false; }, std::nullopt,
              routing_algorithm::kCCH);
    reference_time += t1 - t0;
    experiment_time += std::chrono::steady_clock::now() - t1;

    ASSERT_EQ(reference.size(), experiment.size());
    for (auto k = 0U; k != reference.size(); ++k) {
      ++n_cmp;
      auto const& a = reference[k];
      auto const& b = experiment[k];
      EXPECT_EQ(a.has_value(), b.has_value())
          << "source " << i << " target " << k << ": dijkstra "
          << (a.has_value() ? std::to_string(a->cost_) : "none") << " rphast "
          << (b.has_value() ? std::to_string(b->cost_) : "none");
      if (a.has_value() && b.has_value()) {
        ++n_reachable;
        EXPECT_EQ(a->cost_, b->cost_) << "source " << i << " target " << k;
        if (a->cost_ == b->cost_) {
          ++n_agree;
        }
        // The duration reported by rphast is derived from cost and is a bound,
        // not the exact figure the reference tracks separately. It must never
        // come in under the truth.
        EXPECT_GE(b->duration_.count(), a->duration_.count())
            << "source " << i << " target " << k << ": rphast duration "
            << b->duration_.count() << " under reference "
            << a->duration_.count();
        if (b->duration_.count() != a->duration_.count()) {
          ++n_dur_differ;
        }
      } else if (!a.has_value() && !b.has_value()) {
        ++n_agree;
      }
    }
  }

  fmt::println(
      "lazy rphast vs dijkstra ({}): {} sources, {} comparisons, {} reachable, "
      "{} agree, {} duration bounds above exact",
      dir == direction::kForward ? "fwd" : "bwd", n_src, n_cmp, n_reachable,
      n_agree, n_dur_differ);
  fmt::println("  speedup over dijkstra: {:.2f}",
               static_cast<double>(reference_time.count()) /
                   static_cast<double>(experiment_time.count()));
  EXPECT_GT(n_src, 0U);
  EXPECT_EQ(n_cmp, n_agree);
}

}  // namespace

TEST(rphast, monaco) {
  auto const raw_data = "test/monaco.osm.pbf";
  auto const data_dir = "test/monaco";

  if (!fs::exists(raw_data) && !fs::exists(data_dir)) {
    GTEST_SKIP() << raw_data << " not found";
  }
  load(raw_data, data_dir);
  if (!cch::exists(data_dir)) {
    GTEST_SKIP() << "no cch in " << data_dir;
  }
  auto const w = osr::ways{data_dir, cista::mmap::protection::READ};
  auto const l = osr::lookup{w, data_dir, cista::mmap::protection::READ};

  run(w, l, 200U, 25U, 2 * 3600U, direction::kForward);
}

TEST(rphast, hamburg) {
  auto const raw_data = "test/hamburg.osm.pbf";
  auto const data_dir = "test/hamburg";

  if (!fs::exists(raw_data) && !fs::exists(data_dir)) {
    GTEST_SKIP() << raw_data << " not found";
  }
  load(raw_data, data_dir);
  if (!cch::exists(data_dir)) {
    GTEST_SKIP() << "no cch in " << data_dir;
  }
  auto const w = osr::ways{data_dir, cista::mmap::protection::READ};
  auto const l = osr::lookup{w, data_dir, cista::mmap::protection::READ};

  run(w, l, 60U, 50U, 4 * 3600U, direction::kForward);
}

TEST(rphast, monaco_backward) {
  auto const data_dir = "test/monaco";
  if (!fs::exists(data_dir) || !cch::exists(data_dir)) {
    GTEST_SKIP() << "no cch in " << data_dir;
  }
  auto const w = osr::ways{data_dir, cista::mmap::protection::READ};
  auto const l = osr::lookup{w, data_dir, cista::mmap::protection::READ};

  run(w, l, 200U, 25U, 2 * 3600U, direction::kBackward);
}

TEST(rphast, hamburg_backward) {
  auto const data_dir = "test/hamburg";
  if (!fs::exists(data_dir) || !cch::exists(data_dir)) {
    GTEST_SKIP() << "no cch in " << data_dir;
  }
  auto const w = osr::ways{data_dir, cista::mmap::protection::READ};
  auto const l = osr::lookup{w, data_dir, cista::mmap::protection::READ};

  run(w, l, 60U, 50U, 4 * 3600U, direction::kBackward);
}

TEST(rphast, switzerland) {
  auto const data_dir = "test/scale/switzerland";

  if (!fs::exists(data_dir)) {
    GTEST_SKIP() << data_dir << " not found";
  }
  if (!cch::exists(data_dir)) {
    GTEST_SKIP() << "no cch in " << data_dir;
  }
  auto const w = osr::ways{data_dir, cista::mmap::protection::READ};
  auto const l = osr::lookup{w, data_dir, cista::mmap::protection::READ};

  run(w, l, 15U, 50U, 12 * 3600U, direction::kForward);
}
