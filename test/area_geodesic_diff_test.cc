#include "gtest/gtest.h"

#include <cmath>
#include <cstdio>
#include <algorithm>
#include <array>
#include <limits>
#include <random>
#include <vector>

#include "tg.h"

#include "geo/constants.h"
#include "geo/latlng.h"

#include "osr/area/geodesic.h"

using namespace osr;

namespace {

// Differential test for area_geodesics.
//
// src/area/geodesic.cc decides visibility with its own orientation predicates
// and a midpoint containment test. This builds the same visibility graph using
// tg's `covers` instead - a different library, a different algorithm, the same
// predicate Valhalla uses through GEOS - runs Dijkstra over it, and requires
// the two to produce the same connector-pair distances on random polygons with
// holes.
//
// A disagreement is a real bug in one of them, and the direction says which:
// a shorter reference distance means the reference accepted an edge that
// leaves the polygon, a longer one means the implementation refused an edge it
// should have allowed.
//
// Note what the reference may NOT be built from: tg_geom_covers(polygon,
// linestring) reports a segment as covered when it leaves through a vertex
// without properly crossing any edge - it says "covered" for a chord whose own
// midpoint it separately reports as outside. That is precisely the case this
// is meant to check, so the reference samples the segment and asks tg about
// each point instead, which is a predicate tg gets right.
//
// Scope: this tests the visibility and shortest-path logic only. It replicates
// the same projection the implementation uses, so a mistake in the projection
// would be invisible to it.

constexpr auto kLat0 = 48.86;
constexpr auto kLng0 = 2.34;

struct xy {
  double x_{0.0}, y_{0.0};
};

geo::latlng to_latlng(xy const p) {
  return {kLat0 + p.y_ / geo::kApproxDistanceLatDegrees,
          kLng0 + p.x_ / (std::cos(kLat0 * geo::kPI / 180.0) *
                          geo::kApproxDistanceLatDegrees)};
}

// Same local plane as src/area/geodesic.cc, so both sides see bit-identical
// coordinates and any difference is algorithmic.
struct plane {
  explicit plane(geo::latlng const& ref)
      : lat0_{ref.lat()},
        lng0_{ref.lng()},
        lng_scale_{std::cos(ref.lat() * geo::kPI / 180.0) *
                   geo::kApproxDistanceLatDegrees} {}

  xy operator()(geo::latlng const& p) const {
    return {(p.lng() - lng0_) * lng_scale_,
            (p.lat() - lat0_) * geo::kApproxDistanceLatDegrees};
  }

  double lat0_, lng0_, lng_scale_;
};

double dist(xy const a, xy const b) {
  return std::hypot(b.x_ - a.x_, b.y_ - a.y_);
}

// A star-shaped ring around `center`: sorted angles with jittered radii can
// never self-intersect, so every generated polygon is simple by construction.
std::vector<xy> star_ring(std::mt19937& rng,
                          xy const center,
                          double const r_min,
                          double const r_max,
                          std::size_t const n,
                          bool const clockwise,
                          double const grid) {
  auto radius = std::uniform_real_distribution<double>{r_min, r_max};
  auto out = std::vector<xy>{};
  out.reserve(n);
  for (auto i = std::size_t{0U}; i != n; ++i) {
    auto const a =
        2.0 * geo::kPI * static_cast<double>(i) / static_cast<double>(n);
    auto const r = radius(rng);
    auto p = xy{center.x_ + r * std::cos(a), center.y_ + r * std::sin(a)};
    if (grid > 0.0) {
      // Snapping to a grid manufactures collinear triples and exact vertex
      // hits, which is where these predicates are most likely to disagree.
      p = {std::round(p.x_ / grid) * grid, std::round(p.y_ / grid) * grid};
    }
    out.push_back(p);
  }
  if (clockwise) {
    std::reverse(begin(out), end(out));
  }
  return out;
}

struct scenario {
  std::vector<std::vector<geo::latlng>> rings_;
  std::vector<geo::latlng> connectors_;
};

std::optional<scenario> make_scenario(std::mt19937& rng, double const grid) {
  auto const outer_r = 100.0;
  auto n_outer = std::uniform_int_distribution<std::size_t>{6U, 18U};
  auto n_hole = std::uniform_int_distribution<std::size_t>{4U, 8U};
  auto n_holes = std::uniform_int_distribution<int>{0, 3};
  auto unit = std::uniform_real_distribution<double>{-1.0, 1.0};
  auto hole_r = std::uniform_real_distribution<double>{6.0, 16.0};

  auto rings_xy = std::vector<std::vector<xy>>{};
  rings_xy.push_back(star_ring(rng, xy{0.0, 0.0}, 0.65 * outer_r, outer_r,
                               n_outer(rng), false, grid));

  // Holes are kept well clear of the outer ring and of each other, so the
  // polygon stays valid without needing a general validity check.
  auto centers = std::vector<std::pair<xy, double>>{};
  auto const n = n_holes(rng);
  for (auto i = 0; i != n; ++i) {
    for (auto attempt = 0; attempt != 20; ++attempt) {
      auto const r = hole_r(rng);
      auto const c = xy{unit(rng) * 0.45 * outer_r, unit(rng) * 0.45 * outer_r};
      if (dist(c, xy{0.0, 0.0}) + r > 0.5 * outer_r) {
        continue;
      }
      auto const clear = std::ranges::all_of(centers, [&](auto const& o) {
        return dist(c, o.first) > r + o.second + 6.0;
      });
      if (!clear) {
        continue;
      }
      centers.emplace_back(c, r);
      rings_xy.push_back(
          star_ring(rng, c, 0.6 * r, r, n_hole(rng), true, grid));
      break;
    }
  }

  auto s = scenario{};
  for (auto const& r : rings_xy) {
    auto& out = s.rings_.emplace_back();
    for (auto const& p : r) {
      out.push_back(to_latlng(p));
    }
  }

  // Connectors are outer-ring vertices, which is what real entry points are.
  auto const& outer = s.rings_.front();
  auto pick = std::uniform_int_distribution<std::size_t>{0U, outer.size() - 1U};
  auto n_connectors = std::uniform_int_distribution<std::size_t>{2U, 6U};
  auto chosen = std::vector<std::size_t>{};
  auto const want = std::min(n_connectors(rng), outer.size());
  for (auto attempt = 0; attempt != 50 && chosen.size() != want; ++attempt) {
    auto const i = pick(rng);
    if (std::ranges::find(chosen, i) == end(chosen)) {
      chosen.push_back(i);
    }
  }
  if (chosen.size() < 2U) {
    return std::nullopt;
  }
  for (auto const i : chosen) {
    s.connectors_.push_back(outer[i]);
  }
  return s;
}

struct geom_deleter {
  void operator()(tg_geom* g) const { tg_geom_free(g); }
};
using geom_ptr = std::unique_ptr<tg_geom, geom_deleter>;

// The polygon as a tg geometry, in the same plane coordinates. tg's predicates
// are planar, so feeding projected metres straight in is exactly right.
geom_ptr make_tg_polygon(std::vector<std::vector<xy>> const& rings) {
  auto to_ring = [](std::vector<xy> const& r) {
    auto pts = std::vector<tg_point>{};
    pts.reserve(r.size() + 1U);
    for (auto const& p : r) {
      pts.push_back(tg_point{p.x_, p.y_});
    }
    pts.push_back(pts.front());
    return tg_ring_new(pts.data(), static_cast<int>(pts.size()));
  };

  auto* outer = to_ring(rings.front());
  auto holes = std::vector<tg_ring*>{};
  for (auto i = std::size_t{1U}; i != rings.size(); ++i) {
    holes.push_back(to_ring(rings[i]));
  }
  auto* poly = tg_poly_new(outer, holes.data(), static_cast<int>(holes.size()));
  auto geom = geom_ptr{tg_geom_new_polygon(poly)};
  tg_poly_free(poly);
  tg_ring_free(outer);
  for (auto* h : holes) {
    tg_ring_free(h);
  }
  return geom;
}

// Densely sampled point-in-polygon. Endpoints are skipped: they are polygon
// vertices, so they sit exactly on the boundary either way and say nothing.
//
// The residual risk is an excursion narrower than the sample step, which would
// slip through unseen; at 25cm against generated features no smaller than 6m
// there is a wide margin. tg_geom_covers on the whole line is ANDed in on top,
// since it catches proper crossings exactly even though it misses vertex
// exits.
bool tg_visible(tg_geom const* poly, xy const a, xy const b) {
  auto const len = dist(a, b);
  if (len < 1e-9) {
    return true;
  }

  auto const pts = std::array{tg_point{a.x_, a.y_}, tg_point{b.x_, b.y_}};
  auto* line = tg_line_new(pts.data(), 2);
  auto* geom = tg_geom_new_linestring(line);
  auto const covered = tg_geom_covers(poly, geom);
  tg_geom_free(geom);
  tg_line_free(line);
  if (!covered) {
    return false;
  }

  auto const n = std::max(32, static_cast<int>(len / 0.25));
  auto const dx = (b.x_ - a.x_) / len;
  auto const dy = (b.y_ - a.y_) / len;
  for (auto i = 1; i != n; ++i) {
    auto const s = len * static_cast<double>(i) / static_cast<double>(n);
    if (!tg_geom_intersects_xy(poly, a.x_ + dx * s, a.y_ + dy * s)) {
      return false;
    }
  }
  return true;
}

// Reference distances: tg for visibility, plain Dijkstra on top.
std::vector<double> reference_distances(
    std::vector<std::vector<geo::latlng>> const& rings,
    std::vector<geo::latlng> const& connectors) {
  auto const p = plane{rings.front().front()};

  auto pts = std::vector<xy>{};
  auto connector_idx = std::vector<std::size_t>{};
  auto rings_xy = std::vector<std::vector<xy>>{};
  auto ring_idx = std::vector<std::vector<std::size_t>>{};

  auto const add = [&](geo::latlng const& c) {
    auto const q = p(c);
    for (auto i = std::size_t{0U}; i != pts.size(); ++i) {
      if (dist(pts[i], q) <= 1e-2) {
        return i;
      }
    }
    pts.push_back(q);
    return pts.size() - 1U;
  };

  for (auto const& c : connectors) {
    connector_idx.push_back(add(c));
  }
  for (auto const& r : rings) {
    auto& idx = ring_idx.emplace_back();
    auto& xys = rings_xy.emplace_back();
    for (auto const& c : r) {
      idx.push_back(add(c));
      xys.push_back(p(c));
    }
  }

  auto const poly = make_tg_polygon(rings_xy);

  auto const n = pts.size();
  auto adj = std::vector<std::vector<std::pair<std::size_t, double>>>(n);
  auto is_ring_edge = std::vector<bool>(n * n, false);
  for (auto const& idx : ring_idx) {
    for (auto i = std::size_t{0U}; i != idx.size(); ++i) {
      auto const a = idx[i];
      auto const b = idx[(i + 1U) % idx.size()];
      if (a != b && !is_ring_edge[a * n + b]) {
        is_ring_edge[a * n + b] = is_ring_edge[b * n + a] = true;
        adj[a].emplace_back(b, dist(pts[a], pts[b]));
        adj[b].emplace_back(a, dist(pts[a], pts[b]));
      }
    }
  }

  for (auto i = std::size_t{0U}; i != n; ++i) {
    for (auto j = i + 1U; j != n; ++j) {
      if (is_ring_edge[i * n + j]) {
        continue;
      }
      if (tg_visible(poly.get(), pts[i], pts[j])) {
        auto const d = dist(pts[i], pts[j]);
        adj[i].emplace_back(j, d);
        adj[j].emplace_back(i, d);
      }
    }
  }

  auto const inf = std::numeric_limits<double>::infinity();
  auto out = std::vector<double>(connectors.size() * connectors.size(), inf);
  for (auto s = std::size_t{0U}; s != connectors.size(); ++s) {
    auto d = std::vector<double>(n, inf);
    auto done = std::vector<bool>(n, false);
    d[connector_idx[s]] = 0.0;
    for (;;) {
      auto u = n;
      auto best = inf;
      for (auto i = std::size_t{0U}; i != n; ++i) {
        if (!done[i] && d[i] < best) {
          best = d[i];
          u = i;
        }
      }
      if (u == n) {
        break;
      }
      done[u] = true;
      for (auto const& [v, w] : adj[u]) {
        d[v] = std::min(d[v], d[u] + w);
      }
    }
    for (auto t = std::size_t{0U}; t != connectors.size(); ++t) {
      out[s * connectors.size() + t] = d[connector_idx[t]];
    }
  }
  return out;
}

struct mismatch_report {
  int n_cases_{0};
  int n_pairs_{0};
  int n_too_long_{0};  // implementation refused a valid edge
  int n_too_short_{0};  // implementation allowed an invalid edge
  double worst_{0.0};
  unsigned worst_seed_{0U};
};

mismatch_report run(double const grid, int const n_cases) {
  auto report = mismatch_report{};
  for (auto seed = 0; seed != n_cases; ++seed) {
    auto rng =
        std::mt19937{static_cast<unsigned>(seed) + (grid > 0 ? 1U << 20U : 0U)};
    auto const s = make_scenario(rng, grid);
    if (!s.has_value()) {
      continue;
    }

    auto const g = area_geodesics{s->rings_, s->connectors_};
    auto const ref = reference_distances(s->rings_, s->connectors_);

    ++report.n_cases_;
    auto const k = s->connectors_.size();
    for (auto i = std::size_t{0U}; i != k; ++i) {
      for (auto j = i + 1U; j != k; ++j) {
        auto const mine = static_cast<double>(g.distance(i, j));
        auto const theirs = ref[i * k + j];
        if (std::isinf(mine) && std::isinf(theirs)) {
          continue;
        }
        ++report.n_pairs_;
        auto const tol = std::max(0.01, 1e-3 * std::max(mine, theirs));
        auto const diff = mine - theirs;
        if (std::abs(diff) <= tol) {
          continue;
        }
        if (diff > 0.0) {
          ++report.n_too_long_;
        } else {
          ++report.n_too_short_;
        }
        if (std::abs(diff) > report.worst_) {
          report.worst_ = std::abs(diff);
          report.worst_seed_ = static_cast<unsigned>(seed);
        }
      }
    }
  }
  return report;
}

}  // namespace

TEST(area_geodesic_diff, continuous_coordinates) {
  auto const r = run(0.0, 300);
  std::printf("[          ] %d polygons, %d connector pairs compared\n",
              r.n_cases_, r.n_pairs_);
  ASSERT_GT(r.n_pairs_, 1000);
  EXPECT_EQ(0, r.n_too_long_ + r.n_too_short_)
      << r.n_cases_ << " polygons, " << r.n_pairs_ << " pairs, "
      << r.n_too_long_ << " too long, " << r.n_too_short_ << " too short, "
      << "worst " << r.worst_ << "m at seed " << r.worst_seed_;
}

TEST(area_geodesic_diff, grid_snapped_coordinates) {
  // Snapped to 1m, so collinear triples and exact vertex hits are common.
  auto const r = run(1.0, 300);
  std::printf("[          ] %d polygons, %d connector pairs compared\n",
              r.n_cases_, r.n_pairs_);
  ASSERT_GT(r.n_pairs_, 1000);
  EXPECT_EQ(0, r.n_too_long_ + r.n_too_short_)
      << r.n_cases_ << " polygons, " << r.n_pairs_ << " pairs, "
      << r.n_too_long_ << " too long, " << r.n_too_short_ << " too short, "
      << "worst " << r.worst_ << "m at seed " << r.worst_seed_;
}
