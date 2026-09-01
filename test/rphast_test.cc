#include "gtest/gtest.h"

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
#include "osr/routing/cch/customize.h"
#include "osr/routing/cch/rphast.h"
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
// one-to-many that production uses today.
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

    auto const reference =
        route(car::parameters{}, w, l, search_profile::kCar, from_loc, to_locs,
              from_match, to_match, max_cost, dir, nullptr, nullptr, nullptr,
              [](path const&) { return false; }, std::nullopt,
              routing_algorithm::kDijkstra);
    auto const experiment =
        route(car::parameters{}, w, l, search_profile::kCar, from_loc, to_locs,
              from_match, to_match, max_cost, dir, nullptr, nullptr, nullptr,
              [](path const&) { return false; }, std::nullopt,
              routing_algorithm::kCCH);

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
      "rphast vs dijkstra ({}): {} sources, {} comparisons, {} reachable, "
      "{} agree, {} duration bounds above exact",
      dir == direction::kForward ? "fwd" : "bwd", n_src, n_cmp, n_reachable,
      n_agree, n_dur_differ);
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

// The packed sweep is opt-in and therefore never exercised through `route()`,
// so it is validated directly against the default sweep on identical input.
// The default sweep is congruent with the Dijkstra above, so agreement here
// carries that guarantee over to the packed path.
namespace {

void run_packed_equivalence(ways const& w,
                            lookup const& l,
                            cch const& c,
                            cch_metric const& m,
                            unsigned const n_sources,
                            unsigned const n_targets,
                            cost_t const max_cost,
                            bool const backward) {
  auto prng = std::mt19937{7U};
  auto distr =
      std::uniform_int_distribution<std::uint32_t>{0U, w.n_nodes() - 1U};
  auto const params = car::parameters{};

  auto targets = std::vector<node_idx_t>{};
  for (auto k = 0U; k != n_targets; ++k) {
    auto const n = node_idx_t{distr(prng)};
    if (c.contains(n)) {
      targets.push_back(n);
    }
  }
  ASSERT_FALSE(targets.empty());

  auto n_cmp = 0U, n_agree = 0U, n_finite = 0U;

  for (auto i = 0U; i != n_sources; ++i) {
    auto const src = node_idx_t{distr(prng)};
    auto const from_loc = location{w.get_node_pos(src).as_latlng(), level_t{}};
    auto from_m = match_result{};
    l.match<car>(params, from_loc, false, direction::kForward, kMaxMatchDistance,
                 nullptr, from_m);
    if (from_m.empty() || from_m[match_idx_t{0U}].empty()) {
      continue;
    }
    auto const fm = from_m[match_idx_t{0U}];

    auto const seed = [&](rphast<car>& rp) {
      rp.clear();
      for (auto const t : targets) {
        rp.add_target(c, t);
      }
      for (auto j = 0U; j != fm.size(); ++j) {
        auto const way = fm.way_[j];
        auto const left = fm.left(j);
        auto const right = fm.right(j);
        for (auto const* nc : {&left, &right}) {
          if (!nc->valid() || nc->cost_ >= max_cost) {
            continue;
          }
          auto const sc = car::way_cost(
              params, *w.r_, w.timezones_, way, w.r_->way_properties_[way],
              flip(direction::kForward, nc->way_dir_),
              static_cast<distance_t>(nc->dist_to_node_), std::nullopt,
              duration_t{0}, direction::kForward);
          if (sc.cost_ == kInfeasible || sc.cost_ >= max_cost) {
            continue;
          }
          car::resolve_start_node(
              *w.r_, way, nc->node_, level_t{}, direction::kForward,
              [&](auto const node) {
                rp.add_start(c, node.n_, make_port(node.way_, node.dir_),
                             sc.cost_);
              });
        }
      }
    };

    auto plain = std::make_unique<rphast<car>>();
    seed(*plain);
    if (plain->no_starts()) {
      continue;
    }
    plain->select(*w.r_, c, /* pack */ false, backward);
    plain->run(params, w, c, m, max_cost);

    auto packed = std::make_unique<rphast<car>>();
    seed(*packed);
    packed->select(*w.r_, c, /* pack */ true, backward);
    packed->materialize(params, w, m);
    packed->run(params, w, c, m, max_cost);
    ASSERT_TRUE(packed->packed());

    for (auto const t : targets) {
      for (auto p = port_t{0U}; p != kMaxPorts; ++p) {
        auto const a = plain->get(c, t, p);
        auto const b = packed->get(c, t, p);
        ++n_cmp;
        EXPECT_EQ(a, b) << "source " << i << " target " << to_idx(t) << " port "
                        << static_cast<unsigned>(p);
        if (a == b) {
          ++n_agree;
        }
        if (a != kInfeasible) {
          ++n_finite;
        }
      }
    }
  }

  fmt::println("packed vs default sweep ({}): {} comparisons, {} finite, "
               "{} agree",
               backward ? "bwd" : "fwd", n_cmp, n_finite, n_agree);
  EXPECT_GT(n_finite, 0U);
  EXPECT_EQ(n_cmp, n_agree);
}

}  // namespace

TEST(rphast, packed_equivalence_monaco) {
  auto const data_dir = "test/monaco";
  if (!fs::exists(data_dir) || !cch::exists(data_dir)) {
    GTEST_SKIP() << "no cch in " << data_dir;
  }
  auto const w = osr::ways{data_dir, cista::mmap::protection::READ};
  auto const l = osr::lookup{w, data_dir, cista::mmap::protection::READ};
  auto const c = cch::read(data_dir);
  auto const m = cch_metric::read(data_dir);

  run_packed_equivalence(w, l, *c, *m, 25U, 40U, 2 * 3600U, false);
  run_packed_equivalence(w, l, *c, *m, 25U, 40U, 2 * 3600U, true);
}

TEST(rphast, packed_equivalence_hamburg) {
  auto const data_dir = "test/hamburg";
  if (!fs::exists(data_dir) || !cch::exists(data_dir)) {
    GTEST_SKIP() << "no cch in " << data_dir;
  }
  auto const w = osr::ways{data_dir, cista::mmap::protection::READ};
  auto const l = osr::lookup{w, data_dir, cista::mmap::protection::READ};
  auto const c = cch::read(data_dir);
  auto const m = cch_metric::read(data_dir);

  run_packed_equivalence(w, l, *c, *m, 10U, 60U, 4 * 3600U, false);
  run_packed_equivalence(w, l, *c, *m, 10U, 60U, 4 * 3600U, true);
}
