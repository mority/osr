// Measures how well a flat per-area crossing cost can model real OSM areas.
//
// For every routable area in an OSM file this computes the exact geodesic
// between each pair of its connectors and reports the *spread* of those
// distances (longest crossing minus shortest). Spread, not size, is what a
// single scalar per area has to cover: a cell whose crossings are all about the
// same length is modelled exactly by one number however large it is, while a
// cell mixing a 2m hop with a 200m crossing cannot be, at any threshold.
//
// It also reports what the perimeter-following route costs against the
// geodesic, which is the size of the prize, and what merging near-coincident
// connectors does to the spread, which is the cheapest known fix for the
// degenerate pairs that would otherwise force subdivision down to threshold
// scale everywhere.

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <queue>
#include <string>
#include <string_view>
#include <vector>

#include "fmt/core.h"
#include "fmt/format.h"

#include "utl/enumerate.h"
#include "utl/parser/arg_parser.h"
#include "utl/parser/cstr.h"

// The library reads .osm.pbf; registering the XML reader here lets the
// harness take .osm files too.
#include "osmium/io/xml_input.hpp"
#include "osmium/osm/types.hpp"

#include "geo/box.h"
#include "geo/latlng.h"

#include "osr/area/cells.h"
#include "osr/area/geodesic.h"
#include "osr/area/pipeline.h"
#include "osr/area/served.h"
#include "osr/area/walkable.h"
#include "osr/types.h"

namespace {

// Layers 0 to 3 - what a walkable area is, whether it needs cells, its
// geodesics, its cells - live in the library (osr/area/pipeline.h ties them
// together) and are tested there. These are the harness's names for layer 0.
using level_set = osr::area_levels;
using ring = osr::area_ring;
using area_record = osr::walkable_area;

bool g_strict_levels = false;

// Overlapping areas on the same level and layer are merged before anything
// else looks at them, unless --no-merge.
bool g_merge = true;

// Whether a place=square without any surface tag is an area (--place-square).
bool g_place_square = false;

// Whether ways mapped inside an area count as a way across it when deciding
// if it is already served. On by default: where footways are mapped across an
// area, routing follows them and the area gets no cells. --no-interior-ways
// counts only the area's outline, to see area routing on its own.
bool g_interior_ways = true;

// Cells neighbour where their regions share a border, unless this is set
// (--tree-neighbours): then each cut joins just one pair, as it used to.
bool g_tree_neighbours = false;

// Areas above this are declined as too big (--max-vertices). The strategy
// experiments keep kMaxVertices whatever this is set to: they build the
// visibility graph up to 21 times per area, the cells model once.
std::size_t g_max_vertices = 400U;

// The visibility graph is O(n^3) in the ring vertices, so huge areas have to
// be declined rather than waited for. OSRM stops its mesher at 100 obstacle
// vertices and its query-time solver at 256; this is an offline measurement so
// it can afford more, but not unboundedly - Ile-de-France has a pedestrian
// area of 2739 nodes.
constexpr auto kMaxVertices = std::size_t{400U};

// star_fit is O(k^3) in connectors and is only a side experiment, so it gets
// its own, tighter limit.
constexpr auto kMaxStarConnectors = std::size_t{120U};

constexpr auto kGroupingDistance = 3.0;  // Valhalla's kEntranceGroupingMeters.

// Depth of the placeholder subdivision used to measure error composition.
constexpr auto kMaxSplitDepth = std::size_t{3U};

// The greedy build tests O(n^2) candidate diagonals, each with an O(n)
// visibility check, so it needs a tighter cap than the oracle itself.
constexpr auto kMaxTriangulationVertices = std::size_t{110U};

// Adaptive subdivision needs headroom beyond the fixed-depth experiments.
constexpr auto kAdaptiveMaxDepth = std::size_t{6U};

// The fit alternates between deriving routes and refitting costs. It has
// converged by four rounds - sixteen produces identical output to the decimal
// on Berlin - so this is not a knob worth turning.
constexpr auto kFitIterations = 4;

// How many portals a single cell border may carry. One portal is the
// assumption the whole flat-cost framework was built on; it is only neutral
// for routes already inside the same pair of cells.
std::size_t g_portals_per_border = 1U;

struct pair_stats {
  std::size_t n_pairs_{0U};
  double min_{0.0}, max_{0.0}, mean_{0.0};

  double spread() const { return max_ - min_; }
  // The centre that minimises the worst-case error is the midrange, not the
  // mean; the mean minimises squared error instead.
  double midrange() const { return 0.5 * (min_ + max_); }
  double max_error_midrange() const { return 0.5 * spread(); }
  double max_error_mean() const { return std::max(mean_ - min_, max_ - mean_); }
};

pair_stats compute_pair_stats(osr::area_geodesics const& g,
                              std::vector<std::size_t> const& idx,
                              std::vector<bool> const* mask = nullptr,
                              std::size_t const stride = 0U) {
  auto s = pair_stats{};
  auto sum = 0.0;
  for (auto i = std::size_t{0U}; i != idx.size(); ++i) {
    for (auto j = i + 1U; j != idx.size(); ++j) {
      auto const d = g.distance(idx[i], idx[j]);
      if (d == osr::area_geodesics::kUnreachable) {
        continue;
      }
      if (mask != nullptr && !(*mask)[idx[i] * stride + idx[j]]) {
        continue;
      }
      if (s.n_pairs_ == 0U) {
        s.min_ = s.max_ = d;
      } else {
        s.min_ = std::min<double>(s.min_, d);
        s.max_ = std::max<double>(s.max_, d);
      }
      sum += d;
      ++s.n_pairs_;
    }
  }
  if (s.n_pairs_ != 0U) {
    s.mean_ = sum / static_cast<double>(s.n_pairs_);
  }
  return s;
}

// Fits a star metric with one radius per connector: g(a,b) ~ r_a + r_b. The
// uniform flat cost is the special case r_a = c/2 for every connector, so this
// says what a second float per connector buys over a single float per area.
//
// Each triple (a,b,c) determines r_a exactly if the metric really is a star:
// r_a = (g_ab + g_ac - g_bc)/2. Averaging over all triples is the standard
// least-squares fit and needs nothing but the distance matrix.
std::vector<double> star_fit(osr::area_geodesics const& g,
                             std::vector<std::size_t> const& idx) {
  auto r = std::vector<double>(idx.size(), 0.0);
  auto const d = [&](std::size_t const a, std::size_t const b) {
    return static_cast<double>(g.distance(idx[a], idx[b]));
  };
  if (idx.size() < 2U) {
    return r;
  }
  if (idx.size() == 2U) {
    // A two-connector area is a star exactly: half the crossing each.
    r[0] = r[1] = 0.5 * d(0, 1);
    return r;
  }
  for (auto a = std::size_t{0U}; a != idx.size(); ++a) {
    auto sum = 0.0;
    auto n = std::size_t{0U};
    for (auto b = std::size_t{0U}; b != idx.size(); ++b) {
      if (b == a || d(a, b) == osr::area_geodesics::kUnreachable) {
        continue;
      }
      for (auto c = b + 1U; c != idx.size(); ++c) {
        if (c == a || d(a, c) == osr::area_geodesics::kUnreachable ||
            d(b, c) == osr::area_geodesics::kUnreachable) {
          continue;
        }
        sum += 0.5 * (d(a, b) + d(a, c) - d(b, c));
        ++n;
      }
    }
    r[a] = n == 0U ? 0.0 : std::max(0.0, sum / static_cast<double>(n));
  }

  // The triple average is a least-squares fit, but the uniform scalar it is
  // being compared against is the minimax-optimal midrange. Shifting every
  // radius by the same delta moves every pair by 2*delta, so the best such
  // shift is available in closed form and at least removes that bias.
  auto lo = std::numeric_limits<double>::max();
  auto hi = std::numeric_limits<double>::lowest();
  for (auto a = std::size_t{0U}; a != idx.size(); ++a) {
    for (auto b = a + 1U; b != idx.size(); ++b) {
      if (d(a, b) == osr::area_geodesics::kUnreachable) {
        continue;
      }
      auto const resid = r[a] + r[b] - d(a, b);
      lo = std::min(lo, resid);
      hi = std::max(hi, resid);
    }
  }
  if (lo <= hi) {
    auto const shift = -0.25 * (lo + hi);
    for (auto& x : r) {
      x = std::max(0.0, x + shift);
    }
  }
  return r;
}

double star_max_error(osr::area_geodesics const& g,
                      std::vector<std::size_t> const& idx,
                      std::vector<double> const& r) {
  auto worst = 0.0;
  for (auto a = std::size_t{0U}; a != idx.size(); ++a) {
    for (auto b = a + 1U; b != idx.size(); ++b) {
      auto const d = g.distance(idx[a], idx[b]);
      if (d == osr::area_geodesics::kUnreachable) {
        continue;
      }
      worst = std::max(worst, std::abs(r[a] + r[b] - static_cast<double>(d)));
    }
  }
  return worst;
}

// Greedy merge of connectors within kGroupingDistance of each other, keeping
// one representative per group.
std::vector<std::size_t> group_connectors(
    std::vector<geo::latlng> const& connectors) {
  auto reps = std::vector<std::size_t>{};
  for (auto i = std::size_t{0U}; i != connectors.size(); ++i) {
    auto const merged = std::ranges::any_of(reps, [&](std::size_t const r) {
      return geo::distance(connectors[i], connectors[r]) <= kGroupingDistance;
    });
    if (!merged) {
      reps.push_back(i);
    }
  }
  return reps;
}

// ---------------------------------------------------------------------------
// End-to-end error: how per-cell flat costs compose along a route.
//
// A threshold T bounds the error of ONE cell. A route crossing k cells pays k
// flat costs, so the worst case is k*T - unless the errors cancel, which they
// will not if the cuts bias the same way. Nothing measured so far says which
// happens, and it decides whether T can be set per cell at all.
//
// The subdivision here is a placeholder: recursive bisection of the connector
// set on its widest axis, with a portal placed between the two halves. A
// better cut lowers the constant; it does not change how error accumulates
// with the number of cells crossed, which is what this is for.
// ---------------------------------------------------------------------------

struct cell_tree {
  std::vector<std::size_t> connectors_;  // indices into the combined point list
  // Two candidates: the point between the halves, and the nearest polygon
  // vertex to it as a fallback for when that point falls outside a concave
  // area. Which one is used is decided once the geodesics are known.
  std::size_t portal_{std::numeric_limits<std::size_t>::max()};
  // Used instead of portal_ by the geodesic cut, which needs no new points:
  // the cut's endpoints are connectors already, shared by both children.
  std::vector<std::size_t> portals_;
  // Whether this cut is materialised. The adaptive loop below refines by
  // flipping these, so the whole tree can be built once - portals and all -
  // and the geodesics computed once, while the decomposition still varies.
  bool expanded_{false};
  std::unique_ptr<cell_tree> lhs_, rhs_;

  bool leaf() const { return lhs_ == nullptr; }
};

// Splits on the widest axis at the median, which is the least clever cut that
// is not arbitrary.
void split(cell_tree& c,
           std::vector<geo::latlng> const& pos,
           std::vector<std::vector<geo::latlng>> const& rings,
           std::vector<geo::latlng>& portals,
           std::size_t const depth) {
  if (depth == 0U || c.connectors_.size() < 4U) {
    return;
  }

  auto box = geo::box{};
  for (auto const i : c.connectors_) {
    box.extend(pos[i]);
  }
  auto const span_lat = box.max_.lat() - box.min_.lat();
  auto const span_lng = (box.max_.lng() - box.min_.lng()) *
                        std::cos(box.min_.lat() * geo::kPI / 180.0);
  auto const by_lat = span_lat > span_lng;

  auto sorted = c.connectors_;
  std::ranges::sort(sorted, [&](std::size_t const a, std::size_t const b) {
    return by_lat ? pos[a].lat() < pos[b].lat() : pos[a].lng() < pos[b].lng();
  });
  auto const mid = sorted.size() / 2U;

  c.lhs_ = std::make_unique<cell_tree>();
  c.rhs_ = std::make_unique<cell_tree>();
  c.lhs_->connectors_.assign(begin(sorted), begin(sorted) + mid);
  c.rhs_->connectors_.assign(begin(sorted) + mid, end(sorted));

  // The portal sits between the two halves. It belongs to both, which is what
  // makes them neighbours - no adjacency is stored anywhere.
  auto const centroid = [&](std::vector<std::size_t> const& v) {
    auto lat = 0.0;
    auto lng = 0.0;
    for (auto const i : v) {
      lat += pos[i].lat();
      lng += pos[i].lng();
    }
    return geo::latlng{lat / static_cast<double>(v.size()),
                       lng / static_cast<double>(v.size())};
  };
  auto const a = centroid(c.lhs_->connectors_);
  auto const b = centroid(c.rhs_->connectors_);
  auto const between =
      geo::latlng{0.5 * (a.lat() + b.lat()), 0.5 * (a.lng() + b.lng())};
  auto nearest_vertex = between;
  auto best = std::numeric_limits<double>::max();
  for (auto const& r : rings) {
    for (auto const& v : r) {
      if (auto const d = geo::distance(v, between); d < best) {
        best = d;
        nearest_vertex = v;
      }
    }
  }

  c.portal_ = pos.size() + portals.size();
  portals.push_back(between);
  portals.push_back(nearest_vertex);

  split(*c.lhs_, pos, rings, portals, depth - 1U);
  split(*c.rhs_, pos, rings, portals, depth - 1U);
}

// The spread of a set of connectors over the pairs the area model is actually
// responsible for.
double spread_of(osr::area_geodesics const& g,
                 std::vector<std::size_t> const& idx,
                 std::vector<bool> const& relevant,
                 std::size_t const n_c) {
  auto lo = std::numeric_limits<double>::max();
  auto hi = 0.0;
  auto any = false;
  for (auto i = std::size_t{0U}; i != idx.size(); ++i) {
    for (auto j = i + 1U; j != idx.size(); ++j) {
      if (idx[i] >= n_c || idx[j] >= n_c || !relevant[idx[i] * n_c + idx[j]]) {
        continue;
      }
      auto const d = g.distance(idx[i], idx[j]);
      if (d == osr::area_geodesics::kUnreachable) {
        continue;
      }
      lo = std::min<double>(lo, d);
      hi = std::max<double>(hi, d);
      any = true;
    }
  }
  return any ? hi - lo : 0.0;
}

// Cuts geometrically - so cells stay contiguous and a route crosses them in
// spatial order - but chooses WHERE to cut by the spread criterion, taking the
// line that minimises the worse of the two children. Bisecting the widest axis
// at the median ignores the distances entirely; free clustering optimises them
// but abandons contiguity. This is the middle.
void split_smart(cell_tree& c,
                 std::vector<geo::latlng> const& pos,
                 std::vector<std::vector<geo::latlng>> const& rings,
                 std::vector<geo::latlng>& portals,
                 std::size_t const depth,
                 osr::area_geodesics const& g,
                 std::vector<bool> const& relevant,
                 std::size_t const n_c,
                 double const threshold) {
  if (depth == 0U || c.connectors_.size() < 4U) {
    return;
  }
  // The actual design: stop as soon as this cell is within tolerance, rather
  // than cutting to a fixed depth. A cell whose spread is already under 2T has
  // a flat cost accurate to T and nothing to gain from splitting.
  if (threshold > 0.0 &&
      spread_of(g, c.connectors_, relevant, n_c) <= 2.0 * threshold) {
    return;
  }

  auto const lat0 = pos[c.connectors_.front()].lat();
  auto const lng_scale = std::cos(lat0 * geo::kPI / 180.0);
  auto best_score = std::numeric_limits<double>::max();
  auto best_l = std::vector<std::size_t>{};
  auto best_r = std::vector<std::size_t>{};

  // Exhaustive where it is affordable: every split position along each
  // direction, rather than five fixed fractions. The scan is O(k^3) per
  // direction, so large connector sets keep the coarse grid.
  auto const exhaustive = c.connectors_.size() <= 40U;
  auto const kDirections = exhaustive ? 16U : 8U;
  for (auto d = 0U; d != kDirections; ++d) {
    auto const theta =
        geo::kPI * static_cast<double>(d) / static_cast<double>(kDirections);
    auto const cx = std::cos(theta);
    auto const cy = std::sin(theta);
    auto sorted = c.connectors_;
    std::ranges::sort(sorted, [&](std::size_t const a, std::size_t const b) {
      auto const ka = pos[a].lng() * lng_scale * cx + pos[a].lat() * cy;
      auto const kb = pos[b].lng() * lng_scale * cx + pos[b].lat() * cy;
      return ka < kb;
    });

    auto positions = std::vector<std::size_t>{};
    if (exhaustive) {
      for (auto m = std::size_t{1U}; m != sorted.size(); ++m) {
        positions.push_back(m);
      }
    } else {
      for (auto const f : {0.25, 0.375, 0.5, 0.625, 0.75}) {
        positions.push_back(static_cast<std::size_t>(
            f * static_cast<double>(sorted.size()) + 0.5));
      }
    }

    for (auto const m : positions) {
      if (m < 1U || m >= sorted.size()) {
        continue;
      }
      auto l = std::vector<std::size_t>{begin(sorted), begin(sorted) + m};
      auto r = std::vector<std::size_t>{begin(sorted) + m, end(sorted)};
      auto const score = std::max(spread_of(g, l, relevant, n_c),
                                  spread_of(g, r, relevant, n_c));
      if (score < best_score) {
        best_score = score;
        best_l = std::move(l);
        best_r = std::move(r);
      }
    }
  }

  if (best_l.empty() || best_r.empty()) {
    return;
  }

  c.lhs_ = std::make_unique<cell_tree>();
  c.rhs_ = std::make_unique<cell_tree>();
  c.lhs_->connectors_ = std::move(best_l);
  c.rhs_->connectors_ = std::move(best_r);

  auto const centroid = [&](std::vector<std::size_t> const& v) {
    auto lat = 0.0;
    auto lng = 0.0;
    for (auto const i : v) {
      lat += pos[i].lat();
      lng += pos[i].lng();
    }
    return geo::latlng{lat / static_cast<double>(v.size()),
                       lng / static_cast<double>(v.size())};
  };
  auto const a = centroid(c.lhs_->connectors_);
  auto const b = centroid(c.rhs_->connectors_);
  auto const between =
      geo::latlng{0.5 * (a.lat() + b.lat()), 0.5 * (a.lng() + b.lng())};
  auto nearest_vertex = between;
  auto best = std::numeric_limits<double>::max();
  for (auto const& r : rings) {
    for (auto const& v : r) {
      if (auto const dd = geo::distance(v, between); dd < best) {
        best = dd;
        nearest_vertex = v;
      }
    }
  }
  c.portal_ = pos.size() + portals.size();
  portals.push_back(between);
  portals.push_back(nearest_vertex);

  split_smart(*c.lhs_, pos, rings, portals, depth - 1U, g, relevant, n_c,
              threshold);
  split_smart(*c.rhs_, pos, rings, portals, depth - 1U, g, relevant, n_c,
              threshold);
}

// Cuts along the geodesic between two boundary connectors.
//
// The shortest path inside a polygon bends only at polygon vertices, so the
// cut is made entirely of points that already exist - and its two endpoints
// are connectors, so both children inherit them and no portal has to be
// invented. It also cannot cross a barrier, being a legal path, so neither
// half is severed.
//
// The pair is the median-length one, so the cut follows the grain of the area
// rather than an axis of its bounding box.
void split_geodesic(cell_tree& c,
                    std::vector<geo::latlng> const& pos,
                    std::vector<char> const& on_ring,
                    std::size_t const depth,
                    osr::area_geodesics const& g,
                    std::vector<bool> const& relevant,
                    std::size_t const n_c) {
  if (depth == 0U || c.connectors_.size() < 5U) {
    return;
  }

  auto pairs = std::vector<std::tuple<double, std::size_t, std::size_t>>{};
  for (auto const i : c.connectors_) {
    for (auto const j : c.connectors_) {
      if (i >= j || i >= n_c || j >= n_c || on_ring[i] == 0 ||
          on_ring[j] == 0 || !relevant[i * n_c + j]) {
        continue;
      }
      auto const d = g.distance(i, j);
      if (d != osr::area_geodesics::kUnreachable && d > 1.0) {
        pairs.emplace_back(d, i, j);
      }
    }
  }
  if (pairs.empty()) {
    return;
  }
  std::ranges::sort(pairs);
  auto const [len, a, b] = pairs[pairs.size() / 2U];

  auto const path = g.path(a, b);
  if (path.size() < 2U) {
    return;
  }

  // Which side of the cut each connector falls on, by the nearest segment.
  auto const side_of = [&](geo::latlng const& p) {
    auto best = std::numeric_limits<double>::max();
    auto sign = 0;
    for (auto i = std::size_t{0U}; i + 1U < path.size(); ++i) {
      auto const closest = geo::closest_on_segment(p, path[i], path[i + 1U]);
      if (auto const d = geo::distance(p, closest); d < best) {
        best = d;
        auto const cross =
            (path[i + 1U].lng() - path[i].lng()) * (p.lat() - path[i].lat()) -
            (path[i + 1U].lat() - path[i].lat()) * (p.lng() - path[i].lng());
        sign = cross > 0.0 ? 1 : -1;
      }
    }
    return sign;
  };

  auto lhs = std::vector<std::size_t>{};
  auto rhs = std::vector<std::size_t>{};
  for (auto const i : c.connectors_) {
    if (i == a || i == b) {
      continue;
    }
    (side_of(pos[i]) >= 0 ? lhs : rhs).push_back(i);
  }
  if (lhs.empty() || rhs.empty()) {
    return;
  }

  // The cut's endpoints lie on the border, so they belong to BOTH children -
  // that is what makes them the portal. Dropping them instead starves the
  // recursion, since each child would shed two connectors per cut.
  lhs.push_back(a);
  lhs.push_back(b);
  rhs.push_back(a);
  rhs.push_back(b);

  c.lhs_ = std::make_unique<cell_tree>();
  c.rhs_ = std::make_unique<cell_tree>();
  c.lhs_->connectors_ = std::move(lhs);
  c.rhs_->connectors_ = std::move(rhs);

  split_geodesic(*c.lhs_, pos, on_ring, depth - 1U, g, relevant, n_c);
  split_geodesic(*c.rhs_, pos, on_ring, depth - 1U, g, relevant, n_c);
}

// Cuts ACROSS the geodesic between two boundary connectors, at its midpoint.
//
// Cutting along that path fails because the path is the area's busiest
// corridor, so it separates the connectors that most need each other. Cutting
// perpendicular to it separates the corridor's two ENDS instead, and the
// portal lands on the corridor itself - the one place a route crossing the
// border was always going to pass through.
//
// The portal is a new point (a geodesic midpoint is rarely a polygon vertex),
// but it is guaranteed to lie inside the area, being a point on a legal path.
void split_orthogonal(cell_tree& c,
                      std::vector<geo::latlng> const& pos,
                      std::vector<char> const& on_ring,
                      std::vector<geo::latlng>& portals,
                      std::size_t const depth,
                      osr::area_geodesics const& g,
                      std::vector<bool> const& relevant,
                      std::size_t const n_c,
                      double const threshold,
                      std::size_t const max_portals) {
  if (depth == 0U || c.connectors_.size() < 4U) {
    return;
  }
  if (threshold > 0.0 &&
      spread_of(g, c.connectors_, relevant, n_c) <= 2.0 * threshold) {
    return;
  }

  auto pairs = std::vector<std::tuple<double, std::size_t, std::size_t>>{};
  for (auto const i : c.connectors_) {
    for (auto const j : c.connectors_) {
      if (i >= j || i >= n_c || j >= n_c || on_ring[i] == 0 ||
          on_ring[j] == 0 || !relevant[i * n_c + j]) {
        continue;
      }
      auto const d = g.distance(i, j);
      if (d != osr::area_geodesics::kUnreachable && d > 1.0) {
        pairs.emplace_back(d, i, j);
      }
    }
  }
  if (pairs.empty()) {
    return;
  }
  std::ranges::sort(pairs);
  auto const [len, a, b] = pairs[pairs.size() / 2U];

  auto const path = g.path(a, b);
  if (path.size() < 2U) {
    return;
  }

  // Walk to half the arc length: that point, and the heading there.
  auto total = 0.0;
  for (auto i = std::size_t{0U}; i + 1U < path.size(); ++i) {
    total += geo::distance(path[i], path[i + 1U]);
  }
  auto walked = 0.0;
  auto mid = path.front();
  auto seg = std::size_t{0U};
  for (auto i = std::size_t{0U}; i + 1U < path.size(); ++i) {
    auto const step = geo::distance(path[i], path[i + 1U]);
    if (walked + step >= 0.5 * total || i + 2U == path.size()) {
      auto const f = step < 1e-9 ? 0.0 : (0.5 * total - walked) / step;
      mid = {path[i].lat() + f * (path[i + 1U].lat() - path[i].lat()),
             path[i].lng() + f * (path[i + 1U].lng() - path[i].lng())};
      seg = i;
      break;
    }
    walked += step;
  }

  auto const lng_scale = std::cos(mid.lat() * geo::kPI / 180.0);
  auto const dx = (path[seg + 1U].lng() - path[seg].lng()) * lng_scale;
  auto const dy = path[seg + 1U].lat() - path[seg].lat();
  auto const norm = std::hypot(dx, dy);
  if (norm < 1e-12) {
    return;
  }

  // Which end of the corridor each connector is on.
  auto lhs = std::vector<std::size_t>{};
  auto rhs = std::vector<std::size_t>{};
  for (auto const i : c.connectors_) {
    auto const px = (pos[i].lng() - mid.lng()) * lng_scale;
    auto const py = pos[i].lat() - mid.lat();
    ((px * dx + py * dy) >= 0.0 ? rhs : lhs).push_back(i);
  }
  if (lhs.empty() || rhs.empty()) {
    return;
  }

  c.lhs_ = std::make_unique<cell_tree>();
  c.rhs_ = std::make_unique<cell_tree>();
  c.lhs_->connectors_ = std::move(lhs);
  c.rhs_->connectors_ = std::move(rhs);
  // Portals where the traffic actually is: every other shortest path in this
  // cell is walked, and wherever one crosses the cut line that crossing point
  // becomes a portal candidate. A single portal in the middle of the corridor
  // serves the corridor and taxes everything else; these serve what is there.
  auto const signed_off = [&](geo::latlng const& p) {
    return (p.lng() - mid.lng()) * lng_scale * dx + (p.lat() - mid.lat()) * dy;
  };
  auto crossings = std::vector<std::pair<geo::latlng, int>>{};
  auto const note = [&](geo::latlng const& p) {
    for (auto& [q, cnt] : crossings) {
      if (geo::distance(p, q) < 5.0) {
        ++cnt;
        return;
      }
    }
    crossings.emplace_back(p, 1);
  };
  note(mid);
  for (auto const& [plen, pi, pj] : pairs) {
    auto const pp = g.path(pi, pj);
    for (auto k = std::size_t{0U}; k + 1U < pp.size(); ++k) {
      auto const s0 = signed_off(pp[k]);
      auto const s1 = signed_off(pp[k + 1U]);
      if ((s0 > 0.0) == (s1 > 0.0)) {
        continue;
      }
      auto const f = std::abs(s1 - s0) < 1e-12 ? 0.0 : s0 / (s0 - s1);
      note({pp[k].lat() + f * (pp[k + 1U].lat() - pp[k].lat()),
            pp[k].lng() + f * (pp[k + 1U].lng() - pp[k].lng())});
      break;
    }
  }
  std::ranges::sort(crossings, [](auto const& x, auto const& y) {
    return x.second > y.second;
  });
  crossings.resize(std::min(crossings.size(), max_portals));

  for (auto const& [p, cnt] : crossings) {
    c.portals_.push_back(pos.size() + portals.size());
    portals.push_back(p);
  }

  split_orthogonal(*c.lhs_, pos, on_ring, portals, depth - 1U, g, relevant, n_c,
                   threshold, max_portals);
  split_orthogonal(*c.rhs_, pos, on_ring, portals, depth - 1U, g, relevant, n_c,
                   threshold, max_portals);
}

// A portal that fell outside the polygon is swapped for the boundary vertex
// beside it.
void resolve_portals(cell_tree& c, osr::area_geodesics const& g) {
  if (c.leaf() || !c.portals_.empty()) {
    return;
  }
  if (!g.is_connector_reachable(c.portal_)) {
    ++c.portal_;
  }
  resolve_portals(*c.lhs_, g);
  resolve_portals(*c.rhs_, g);
}

// The cells a route actually sees at this depth, each with the points that
// enter it (its own connectors plus the portals of every cut at or below it).
// Like collect_cells, but the decomposition is whatever the expanded_ flags
// say rather than a uniform depth.
void reset_expanded(cell_tree& c) {
  c.expanded_ = false;
  if (!c.leaf()) {
    reset_expanded(*c.lhs_);
    reset_expanded(*c.rhs_);
  }
}

std::pair<std::size_t, std::size_t> collect_frontier(
    cell_tree& c,
    std::vector<geo::latlng> const& pts,
    std::vector<std::vector<std::size_t>>& out,
    std::vector<cell_tree*>& owner) {
  if (c.leaf() || !c.expanded_) {
    out.push_back(c.connectors_);
    owner.push_back(&c);
    return {out.size() - 1U, out.size()};
  }

  auto const [l0, l1] = collect_frontier(*c.lhs_, pts, out, owner);
  auto const [r0, r1] = collect_frontier(*c.rhs_, pts, out, owner);
  auto const nearest_to = [&](std::size_t const from, std::size_t const to,
                              std::size_t const target) {
    auto best = from;
    auto best_d = std::numeric_limits<double>::max();
    for (auto i = from; i != to; ++i) {
      if (out[i].empty()) {
        continue;
      }
      auto lat = 0.0;
      auto lng = 0.0;
      for (auto const p : out[i]) {
        lat += pts[p].lat();
        lng += pts[p].lng();
      }
      auto const nn = static_cast<double>(out[i].size());
      auto const dd = geo::distance(pts[target], {lat / nn, lng / nn});
      if (dd < best_d) {
        best_d = dd;
        best = i;
      }
    }
    return best;
  };
  if (!c.portals_.empty()) {
    for (auto const p : c.portals_) {
      out[nearest_to(l0, l1, p)].push_back(p);
      out[nearest_to(r0, r1, p)].push_back(p);
    }
  } else if (c.portal_ != std::numeric_limits<std::size_t>::max()) {
    out[nearest_to(l0, l1, c.portal_)].push_back(c.portal_);
    out[nearest_to(r0, r1, c.portal_)].push_back(c.portal_);
  }
  return {l0, r1};
}

std::pair<std::size_t, std::size_t> collect_cells(
    cell_tree const& c,
    std::size_t const depth,
    std::vector<geo::latlng> const& pts,
    std::vector<std::vector<std::size_t>>& out) {
  if (c.leaf() || depth == 0U) {
    out.push_back(c.connectors_);
    return {out.size() - 1U, out.size()};
  }

  auto const [l0, l1] = collect_cells(*c.lhs_, depth - 1U, pts, out);
  auto const [r0, r1] = collect_cells(*c.rhs_, depth - 1U, pts, out);

  // The border this cut created is shared by exactly two cells - the one on
  // each side that the portal actually lies against. Handing it to every cell
  // in the subtree instead would make every pair of cells neighbours, and no
  // route would ever cross more than two.
  auto const nearest_to = [&](std::size_t const from, std::size_t const to,
                              std::size_t const target) {
    auto best = from;
    auto best_d = std::numeric_limits<double>::max();
    for (auto i = from; i != to; ++i) {
      if (out[i].empty()) {
        continue;
      }
      auto lat = 0.0;
      auto lng = 0.0;
      for (auto const p : out[i]) {
        lat += pts[p].lat();
        lng += pts[p].lng();
      }
      auto const n = static_cast<double>(out[i].size());
      auto const d = geo::distance(pts[target], {lat / n, lng / n});
      if (d < best_d) {
        best_d = d;
        best = i;
      }
    }
    return best;
  };
  auto const nearest = [&](std::size_t const from, std::size_t const to) {
    return nearest_to(from, to, c.portal_);
  };
  if (!c.portals_.empty()) {
    for (auto const p : c.portals_) {
      out[nearest_to(l0, l1, p)].push_back(p);
      out[nearest_to(r0, r1, p)].push_back(p);
    }
  } else if (c.portal_ == std::numeric_limits<std::size_t>::max()) {
    // Shared points already sit in both children's own connector lists.
  } else {
    out[nearest(l0, l1)].push_back(c.portal_);
    out[nearest(r0, r1)].push_back(c.portal_);
  }
  return {l0, r1};
}

void check_portals(cell_tree const& c, osr::area_geodesics const& g, bool& ok) {
  if (c.leaf()) {
    return;
  }
  // A cut carries EITHER a list of portals (split_orthogonal) or a single one
  // (split_smart) or none at all (split_geodesic shares its endpoints through
  // the children's own connector lists). Reading portal_ unconditionally means
  // reading its SIZE_MAX sentinel for the first and last of those.
  auto const check = [&](std::size_t const p) {
    ok = ok && p < g.n_connectors() && g.is_connector_reachable(p);
  };
  if (!c.portals_.empty()) {
    for (auto const p : c.portals_) {
      check(p);
    }
  } else if (c.portal_ != std::numeric_limits<std::size_t>::max()) {
    check(c.portal_);
  }
  check_portals(*c.lhs_, g, ok);
  check_portals(*c.rhs_, g, ok);
}

// Groups connectors so that the pairwise distances inside a group are as
// alike as possible - which is exactly what a single flat cost per group can
// represent. Note this is NOT "put nearby connectors together": four
// connectors all ~100m apart are a perfect group, while a 2m pair beside a
// 200m pair is a bad one however tight its bounding box.
//
// Greedy agglomerative: repeatedly merge whichever two groups produce the
// smallest resulting spread. Groups are free to be any subset, so this is an
// upper bound on what a geometric cut can achieve - if it barely beats
// bisection, contiguity costs nothing.
std::vector<std::vector<std::size_t>> cluster_by_spread(
    osr::area_geodesics const& g,
    std::size_t const k,
    std::size_t const target,
    std::vector<bool> const& relevant) {
  auto groups = std::vector<std::vector<std::size_t>>{};
  for (auto i = std::size_t{0U}; i != k; ++i) {
    groups.push_back({i});
  }
  if (k <= target) {
    return groups;
  }

  auto const inf = std::numeric_limits<double>::infinity();
  auto const n = groups.size();
  // Per-group extremes, and the extremes between each pair of groups. A merge
  // only has to combine six numbers, so the whole run stays O(m^3).
  auto lo = std::vector<double>(n, inf);
  auto hi = std::vector<double>(n, 0.0);
  auto cross_lo = std::vector<double>(n * n, inf);
  auto cross_hi = std::vector<double>(n * n, 0.0);
  auto alive = std::vector<char>(n, 1);
  // Reachability and relevance are tracked apart. A pair the perimeter already
  // serves must not drag the spread down - those pairs are why an untreated
  // spread is ~0 at the bottom end for almost every area - but the two groups
  // are still joinable.
  auto linked = std::vector<char>(n * n, 0);
  for (auto i = std::size_t{0U}; i != n; ++i) {
    for (auto j = i + 1U; j != n; ++j) {
      auto const d = g.distance(i, j);
      if (d == osr::area_geodesics::kUnreachable) {
        continue;
      }
      linked[i * n + j] = linked[j * n + i] = 1;
      if (!relevant[i * n + j]) {
        continue;
      }
      cross_lo[i * n + j] = cross_lo[j * n + i] = d;
      cross_hi[i * n + j] = cross_hi[j * n + i] = d;
    }
  }

  auto n_alive = n;
  while (n_alive > target) {
    auto best_i = n;
    auto best_j = n;
    auto best_spread = inf;
    for (auto i = std::size_t{0U}; i != n; ++i) {
      if (alive[i] == 0) {
        continue;
      }
      for (auto j = i + 1U; j != n; ++j) {
        if (alive[j] == 0 || linked[i * n + j] == 0) {
          continue;
        }
        auto const l = std::min({lo[i], lo[j], cross_lo[i * n + j]});
        auto const h = std::max({hi[i], hi[j], cross_hi[i * n + j]});
        if (auto const spread = h - (l == inf ? h : l); spread < best_spread) {
          best_spread = spread;
          best_i = i;
          best_j = j;
        }
      }
    }
    if (best_i == n) {
      break;
    }

    lo[best_i] =
        std::min({lo[best_i], lo[best_j], cross_lo[best_i * n + best_j]});
    hi[best_i] =
        std::max({hi[best_i], hi[best_j], cross_hi[best_i * n + best_j]});
    for (auto t = std::size_t{0U}; t != n; ++t) {
      if (alive[t] == 0 || t == best_i || t == best_j) {
        continue;
      }
      auto const l =
          std::min(cross_lo[best_i * n + t], cross_lo[best_j * n + t]);
      auto const h =
          std::max(cross_hi[best_i * n + t], cross_hi[best_j * n + t]);
      cross_lo[best_i * n + t] = cross_lo[t * n + best_i] = l;
      cross_hi[best_i * n + t] = cross_hi[t * n + best_i] = h;
      linked[best_i * n + t] = linked[t * n + best_i] =
          (linked[best_i * n + t] != 0 || linked[best_j * n + t] != 0) ? 1 : 0;
    }
    groups[best_i].insert(end(groups[best_i]), begin(groups[best_j]),
                          end(groups[best_j]));
    groups[best_j].clear();
    alive[best_j] = 0;
    --n_alive;
  }

  auto out = std::vector<std::vector<std::size_t>>{};
  for (auto& gp : groups) {
    if (!gp.empty()) {
      out.push_back(std::move(gp));
    }
  }
  return out;
}

struct composition_sample {
  int n_cells_{0};
  double error_{0.0};
};

// Builds the flat-cost model at a given depth and compares it, pair by pair,
// against the exact geodesic - recording how many cells each modelled route
// crossed, so the error can be read as a function of that.
struct composition_result {
  std::vector<composition_sample> samples_;
  std::vector<double> cell_errors_;  // the per-cell error a threshold bounds
  // How often each cell carried a route whose modelled cost missed by more
  // than the threshold - which is what says where to subdivide next.
  std::vector<int> cell_violations_;
};

// Solves the small dense normal equations of a least-squares fit.
bool solve_symmetric(std::vector<double>& a,
                     std::vector<double>& b,
                     std::size_t const n) {
  for (auto i = std::size_t{0U}; i != n; ++i) {
    auto piv = i;
    for (auto r = i + 1U; r != n; ++r) {
      if (std::abs(a[r * n + i]) > std::abs(a[piv * n + i])) {
        piv = r;
      }
    }
    if (std::abs(a[piv * n + i]) < 1e-9) {
      return false;
    }
    if (piv != i) {
      for (auto c = std::size_t{0U}; c != n; ++c) {
        std::swap(a[i * n + c], a[piv * n + c]);
      }
      std::swap(b[i], b[piv]);
    }
    for (auto r = i + 1U; r != n; ++r) {
      auto const f = a[r * n + i] / a[i * n + i];
      for (auto c = i; c != n; ++c) {
        a[r * n + c] -= f * a[i * n + c];
      }
      b[r] -= f * b[i];
    }
  }
  for (auto ii = n; ii-- > 0U;) {
    auto x = b[ii];
    for (auto c = ii + 1U; c != n; ++c) {
      x -= a[ii * n + c] * b[c];
    }
    b[ii] = x / a[ii * n + ii];
  }
  return true;
}

composition_result evaluate_composition(
    osr::area_geodesics const& g,
    std::size_t const n_connectors,
    std::size_t const n_points,
    std::vector<std::vector<std::size_t>> const& cells,
    bool const global_fit = false,
    double const violation_threshold = 0.0) {
  auto result = composition_result{};
  result.cell_violations_.assign(cells.size(), 0);
  auto const inf = std::numeric_limits<double>::infinity();

  // One hub per cell, appended after the points.
  auto const n = n_points + cells.size();
  auto adj = std::vector<std::vector<std::pair<std::size_t, double>>>(n);
  for (auto const [ci, pts] : utl::enumerate(cells)) {
    // Flat cost: the midrange over the cell's own pairs, which is the centre
    // that minimises worst-case error.
    auto lo = inf;
    auto hi = 0.0;
    for (auto i = std::size_t{0U}; i != pts.size(); ++i) {
      for (auto j = i + 1U; j != pts.size(); ++j) {
        // A portal is a token, not a place: it exists to make two cells
        // neighbours, and has no distance to anything. Only real connectors
        // seed the cost, and the global fit adjusts from there.
        if (pts[i] >= g.n_connectors() || pts[j] >= g.n_connectors()) {
          continue;
        }
        auto const d = g.distance(pts[i], pts[j]);
        if (d == osr::area_geodesics::kUnreachable) {
          return {};
        }
        lo = std::min<double>(lo, d);
        hi = std::max<double>(hi, d);
      }
    }
    // A cell may have no connector pair of its own to seed from - one real
    // connector plus portals is enough to be a cell, and portals are tokens
    // with no distances. Its hub still has to exist and be connected, or the
    // cell drops out of the graph entirely and nothing can route through it.
    // The fit sets its cost from the routes that pass through it.
    auto const c = lo == inf ? 0.0 : 0.5 * (lo + hi);
    result.cell_errors_.push_back(0.5 * (hi - lo));
    auto const hub = n_points + ci;
    for (auto const p : pts) {
      adj[p].emplace_back(hub, 0.5 * c);
      adj[hub].emplace_back(p, 0.5 * c);
    }
  }

  // The cell costs are free parameters of the model, not something the
  // geometry dictates: the topology says which cells a route crosses, and the
  // costs are whatever best reproduces the true distances. Fitting each cell
  // in isolation from its own point set - which is what the midrange above
  // does - is a local guess, and it is the only reason a portal's POSITION
  // ever affected the result. Refit them together instead.
  if (global_fit) {
    auto const m = cells.size();
    auto cost = std::vector<double>(m, 0.0);
    for (auto const [ci, pts] : utl::enumerate(cells)) {
      auto const hub = n_points + ci;
      cost[ci] = adj[hub].empty() ? 0.0 : 2.0 * adj[hub].front().second;
    }

    for (auto iter = 0; iter != kFitIterations; ++iter) {
      // Which cells each pair's modelled route crosses, under current costs.
      auto rows = std::vector<std::pair<std::vector<std::size_t>, double>>{};
      for (auto src = std::size_t{0U}; src != n_connectors; ++src) {
        auto d = std::vector<double>(n, inf);
        auto pred = std::vector<std::size_t>(n, n);
        auto q =
            std::priority_queue<std::pair<double, std::size_t>,
                                std::vector<std::pair<double, std::size_t>>,
                                std::greater<>>{};
        d[src] = 0.0;
        q.emplace(0.0, src);
        while (!q.empty()) {
          auto const [cst, u] = q.top();
          q.pop();
          if (cst > d[u]) {
            continue;
          }
          for (auto const& [v, w] : adj[u]) {
            if (auto const nx = cst + w; nx < d[v]) {
              d[v] = nx;
              pred[v] = u;
              q.emplace(nx, v);
            }
          }
        }
        for (auto dst = src + 1U; dst != n_connectors; ++dst) {
          auto const truth = g.distance(src, dst);
          if (!std::isfinite(d[dst]) ||
              truth == osr::area_geodesics::kUnreachable || truth < 1.0) {
            continue;
          }
          auto used = std::vector<std::size_t>{};
          for (auto v = dst; v != src && v != n; v = pred[v]) {
            if (v >= n_points) {
              used.push_back(v - n_points);
            }
          }
          if (!used.empty()) {
            rows.emplace_back(std::move(used), static_cast<double>(truth));
          }
        }
      }
      if (rows.size() < m) {
        break;
      }

      auto ata = std::vector<double>(m * m, 0.0);
      auto atb = std::vector<double>(m, 0.0);
      for (auto const& [used, truth] : rows) {
        for (auto const i : used) {
          atb[i] += truth;
          for (auto const j : used) {
            ata[i * m + j] += 1.0;
          }
        }
      }
      for (auto i = std::size_t{0U}; i != m; ++i) {
        ata[i * m + i] += 1e-3;  // keep it solvable when a cell is unused
      }
      if (!solve_symmetric(ata, atb, m)) {
        break;
      }
      for (auto i = std::size_t{0U}; i != m; ++i) {
        cost[i] = std::max(0.0, atb[i]);
      }

      for (auto const [ci, pts] : utl::enumerate(cells)) {
        auto const hub = n_points + ci;
        adj[hub].clear();
        for (auto const p : pts) {
          std::erase_if(adj[p], [&](auto const& e) { return e.first == hub; });
          adj[p].emplace_back(hub, 0.5 * cost[ci]);
          adj[hub].emplace_back(p, 0.5 * cost[ci]);
        }
      }
    }
  }

  using entry = std::tuple<double, std::size_t, int>;
  for (auto src = std::size_t{0U}; src != n_connectors; ++src) {
    auto d = std::vector<double>(n, inf);
    auto hops = std::vector<int>(n, 0);
    auto pred = std::vector<std::size_t>(n, n);
    auto q = std::priority_queue<entry, std::vector<entry>, std::greater<>>{};
    d[src] = 0.0;
    q.emplace(0.0, src, 0);
    while (!q.empty()) {
      auto const [cost, u, h] = q.top();
      q.pop();
      if (cost > d[u]) {
        continue;
      }
      for (auto const& [v, w] : adj[u]) {
        if (auto const next = cost + w; next < d[v]) {
          d[v] = next;
          pred[v] = u;
          hops[v] = h + (v >= n_points ? 1 : 0);
          q.emplace(next, v, hops[v]);
        }
      }
    }
    for (auto dst = src + 1U; dst != n_connectors; ++dst) {
      auto const truth = g.distance(src, dst);
      if (!std::isfinite(d[dst]) ||
          truth == osr::area_geodesics::kUnreachable || truth < 1.0) {
        continue;
      }
      auto const err = std::abs(d[dst] - static_cast<double>(truth));
      result.samples_.push_back({hops[dst], err});
      if (violation_threshold > 0.0 && err > violation_threshold) {
        for (auto v = dst; v != src && v != n; v = pred[v]) {
          if (v >= n_points) {
            ++result.cell_violations_[v - n_points];
          }
        }
      }
    }
  }
  return result;
}

// ---------------------------------------------------------------------------
// Triangulation as the routing structure.
//
// Attractive because it adds no vertices: every diagonal joins two existing
// polygon vertices, so a portal is a real edge between two real OSM nodes
// rather than the synthetic point the cut-based cells need. The cost is that
// paths then follow triangle edges instead of the geodesic, and that is what
// this measures.
//
// Built greedily - shortest non-crossing diagonal first - rather than by ear
// clipping. Greedy is a well-known minimum-weight-triangulation heuristic, it
// needs no hole bridging because visibility already accounts for holes, and it
// reuses the oracle's own predicate so barriers block diagonals for free.
// ---------------------------------------------------------------------------

// Do open segments a-b and c-d cross at an interior point of both?
bool segments_cross(geo::latlng const& a,
                    geo::latlng const& b,
                    geo::latlng const& c,
                    geo::latlng const& d) {
  auto const side = [](geo::latlng const& p, geo::latlng const& q,
                       geo::latlng const& r) {
    auto const v = (q.lng() - p.lng()) * (r.lat() - p.lat()) -
                   (q.lat() - p.lat()) * (r.lng() - p.lng());
    return v > 1e-12 ? 1 : (v < -1e-12 ? -1 : 0);
  };
  auto const d1 = side(a, b, c);
  auto const d2 = side(a, b, d);
  auto const d3 = side(c, d, a);
  auto const d4 = side(c, d, b);
  return d1 * d2 < 0 && d3 * d4 < 0;
}

struct triangulation {
  std::vector<geo::latlng> pts_;
  std::vector<std::pair<std::size_t, std::size_t>> edges_;
};

triangulation greedy_triangulation(
    std::vector<std::vector<geo::latlng>> const& rings,
    std::vector<std::vector<geo::latlng>> const& barriers) {
  auto t = triangulation{};
  auto ring_span = std::vector<std::pair<std::size_t, std::size_t>>{};
  for (auto const& r : rings) {
    auto const from = t.pts_.size();
    t.pts_.insert(end(t.pts_), begin(r), end(r));
    ring_span.emplace_back(from, t.pts_.size());
  }

  auto is_ring_edge = osr::hash_set<std::uint64_t>{};
  auto const key = [](std::size_t const a, std::size_t const b) {
    return (static_cast<std::uint64_t>(std::min(a, b)) << 32U) |
           static_cast<std::uint64_t>(std::max(a, b));
  };
  for (auto const& [from, to] : ring_span) {
    auto const n = to - from;
    for (auto i = std::size_t{0U}; i != n; ++i) {
      auto const a = from + i;
      auto const b = from + (i + 1U) % n;
      t.edges_.emplace_back(a, b);
      is_ring_edge.insert(key(a, b));
    }
  }

  // Every diagonal a triangulation could use, shortest first.
  auto candidates = std::vector<std::tuple<double, std::size_t, std::size_t>>{};
  for (auto i = std::size_t{0U}; i != t.pts_.size(); ++i) {
    for (auto j = i + 1U; j != t.pts_.size(); ++j) {
      if (!is_ring_edge.contains(key(i, j))) {
        candidates.emplace_back(geo::distance(t.pts_[i], t.pts_[j]), i, j);
      }
    }
  }
  std::ranges::sort(candidates);

  // n + 3h - 3 diagonals complete a triangulation of a polygon with h holes.
  auto const h = rings.size() - 1U;
  auto const needed = t.pts_.size() + 3U * h - 3U;
  auto diagonals = std::vector<std::pair<std::size_t, std::size_t>>{};
  for (auto const& [len, i, j] : candidates) {
    if (diagonals.size() >= needed) {
      break;
    }
    auto crosses = false;
    for (auto const& [a, b] : diagonals) {
      if (a != i && a != j && b != i && b != j &&
          segments_cross(t.pts_[i], t.pts_[j], t.pts_[a], t.pts_[b])) {
        crosses = true;
        break;
      }
    }
    if (crosses) {
      continue;
    }
    if (!osr::area_geodesics::is_segment_inside(rings, t.pts_[i], t.pts_[j],
                                                barriers)) {
      continue;
    }
    diagonals.emplace_back(i, j);
  }

  t.edges_.insert(end(t.edges_), begin(diagonals), end(diagonals));
  return t;
}

struct percentiles {
  double p50_{0.0}, p90_{0.0}, p95_{0.0}, max_{0.0};
};

percentiles quantiles_of(std::vector<double> v) {
  if (v.empty()) {
    return {};
  }
  std::ranges::sort(v);
  auto const at = [&](double const q) {
    auto const i =
        static_cast<std::size_t>(q * static_cast<double>(v.size() - 1U) + 0.5);
    return v[std::min(i, v.size() - 1U)];
  };
  return {at(0.50), at(0.90), at(0.95), v.back()};
}

}  // namespace

// The chosen model as osr::build_area_cells builds it - no portals, one fit.
// Route errors come out in the pair order evaluate_composition reports them,
// so the two implementations can be compared entry by entry.
struct library_run {
  std::size_t n_cells_{0U};
  std::vector<double> errors_;
};

std::optional<library_run> run_library(
    std::vector<geo::latlng> const& connectors,
    std::vector<float> const& dist,
    std::vector<bool> const& relevant,
    double const threshold) {
  auto const cells = osr::build_area_cells(connectors, dist, relevant,
                                           {.threshold_ = threshold});
  if (!cells.has_value()) {
    return std::nullopt;
  }
  auto const k = connectors.size();
  auto const m = cells->n_cells();
  auto const cd = cells->cell_distances();
  auto run = library_run{.n_cells_ = m, .errors_ = {}};
  for (auto src = std::size_t{0U}; src != k; ++src) {
    for (auto dst = src + 1U; dst != k; ++dst) {
      auto const truth = dist[src * k + dst];
      auto const a = cells->connector_cell_[src];
      auto const b = cells->connector_cell_[dst];
      if (truth == osr::area_geodesics::kUnreachable || truth < 1.0 ||
          a == osr::area_cells::kNoCell || b == osr::area_cells::kNoCell ||
          !std::isfinite(cd[a * m + b])) {
        continue;
      }
      auto const model = cd[a * m + b];
      run.errors_.push_back(
          std::abs(static_cast<double>(model) - static_cast<double>(truth)));
    }
  }
  return run;
}

// Every area the harness looks at, as GeoJSON for osr-backend's area layer:
// its outline with a status saying what the funnel decided, and for meshed
// areas the cells - connectors coloured by cell, hubs, spokes, and the
// neighbour links between hubs.
struct cells_dump {
  static std::string escape(std::string_view const s) {
    auto out = std::string{};
    for (auto const ch : s) {
      if (ch == '"' || ch == '\\') {
        out.push_back('\\');
        out.push_back(ch);
      } else if (static_cast<unsigned char>(ch) < 0x20U) {
        out += fmt::format("\\u{:04x}", static_cast<unsigned>(ch));
      } else {
        out.push_back(ch);
      }
    }
    return out;
  }

  static std::string coord(geo::latlng const& p) {
    return fmt::format("[{:.7f},{:.7f}]", p.lng(), p.lat());
  }

  static std::string coords(std::vector<geo::latlng> const& line,
                            bool const close) {
    auto out = std::string{"["};
    for (auto const& p : line) {
      if (out.size() != 1U) {
        out.push_back(',');
      }
      out += coord(p);
    }
    if (close && !line.empty() && line.front() != line.back()) {
      out += "," + coord(line.front());
    }
    return out + "]";
  }

  void feature(std::string const& geometry, std::string const& properties) {
    features_.push_back(fmt::format(
        R"({{"type":"Feature","geometry":{},"properties":{{{}}}}})", geometry,
        properties));
  }

  void area(area_record const& a,
            std::vector<std::vector<geo::latlng>> const& rings,
            std::vector<std::vector<geo::latlng>> const& barriers,
            std::vector<geo::latlng> const& connectors,
            std::vector<osmium::object_id_type> const& connector_ids,
            std::vector<char> const& on_ring,
            std::string_view const status,
            osr::area_cells const* cells = nullptr,
            std::vector<float> const* dist = nullptr,
            std::vector<std::vector<osr::area_cut_side>> const* regions =
                nullptr) {
    // Every feature carries its area's id and levels, so a whole area can be
    // selected - and filtered by level - from any of its features. An area
    // without a level tag gets no "levels" at all.
    // "area" tells apart the pieces of a merge, which share their "osm".
    auto osm = fmt::format(R"("osm":"{}/{}","area":{})",
                           a.from_way_ ? "way" : "relation", a.id_,
                           n_areas_++);
    if (!a.levels_.any_) {
      auto levels = std::string{};
      for (auto b = 0U; b != 64U; ++b) {
        if (((a.levels_.bits_ >> b) & 1U) != 0U) {
          levels += fmt::format("{}{}", levels.empty() ? "" : ",",
                                osr::level_t{static_cast<std::uint8_t>(b)}
                                    .to_float());
        }
      }
      osm += fmt::format(R"(,"levels":[{}])", levels);
    }

    auto polygon = std::string{};
    for (auto const& r : rings) {
      if (r.size() >= 3U) {
        polygon += (polygon.empty() ? "" : ",") + coords(r, true);
      }
    }
    auto extra = std::string{};
    if (!a.members_.empty()) {
      auto members = std::string{};
      for (auto const& m : a.members_) {
        members += (members.empty() ? "" : " ") + m;
      }
      extra += fmt::format(R"(,"members":"{}")", members);
    }
    if (cells != nullptr) {
      auto const k = connectors.size();
      auto const m = cells->n_cells();
      auto const cd = cells->cell_distances();
      auto errors = std::vector<double>{};
      for (auto i = std::size_t{0U}; i != k; ++i) {
        for (auto j = i + 1U; j != k; ++j) {
          auto const truth = (*dist)[i * k + j];
          auto const a = cells->connector_cell_[i];
          auto const b = cells->connector_cell_[j];
          if (truth != osr::area_geodesics::kUnreachable && truth >= 1.0 &&
              a != osr::area_cells::kNoCell && b != osr::area_cells::kNoCell) {
            errors.push_back(std::abs(static_cast<double>(cd[a * m + b]) -
                                      static_cast<double>(truth)));
          }
        }
      }
      extra += fmt::format(
          R"(,"cells":{},"unreachable_connectors":{})", m,
          std::ranges::count(cells->connector_cell_, osr::area_cells::kNoCell));
      all_errors_.insert(end(all_errors_), begin(errors), end(errors));
      cells_per_area_.push_back(static_cast<double>(m));
      if (!errors.empty()) {
        auto const q = quantiles_of(errors);
        extra += fmt::format(R"(,"err_p50":{:.1f},"err_max":{:.1f})", q.p50_,
                             q.max_);
      }
    }
    feature(fmt::format(R"({{"type":"Polygon","coordinates":[{}]}})", polygon),
            fmt::format(R"("kind":"area",{},"name":"{}","status":"{}",)"
                        R"("connectors":{}{})",
                        osm, escape(a.name_), status, connectors.size(),
                        extra));

    // Cell regions (osr::region_rings), for drawing.
    if (cells != nullptr && regions != nullptr) {
      auto const polygons = osr::region_rings(rings, *regions);
      for (auto c = std::size_t{0U}; c != polygons.size(); ++c) {
        auto region = std::string{};
        for (auto const& r : polygons[c]) {
          region += (region.empty() ? "" : ",") + coords(r, true);
        }
        if (!region.empty()) {
          feature(
              fmt::format(R"({{"type":"Polygon","coordinates":[{}]}})", region),
              fmt::format(R"("kind":"cell",{},"cell":{},"cost":{:.1f})", osm,
                          c, cells->cost_[c]));
        }
      }
    }

    for (auto const& b : barriers) {
      feature(fmt::format(R"({{"type":"LineString","coordinates":{}}})",
                          coords(b, false)),
              fmt::format(R"("kind":"barrier",{})", osm));
    }

    auto hubs = std::vector<geo::latlng>{};
    if (cells != nullptr) {
      hubs = osr::hub_positions(*cells, connectors);
      auto members = std::vector<int>(cells->n_cells(), 0);
      for (auto const c : cells->connector_cell_) {
        if (c != osr::area_cells::kNoCell) {
          ++members[c];
        }
      }
      for (auto c = std::size_t{0U}; c != hubs.size(); ++c) {
        feature(fmt::format(R"({{"type":"Point","coordinates":{}}})",
                            coord(hubs[c])),
                fmt::format(R"("kind":"hub",{},"cell":{},"cost":{:.1f},)"
                            R"("members":{})",
                            osm, c, cells->cost_[c], members[c]));
        for (auto d = c + 1U; d != hubs.size(); ++d) {
          if (cells->is_neighbour(static_cast<osr::area_cells::cell_idx_t>(c),
                                  static_cast<osr::area_cells::cell_idx_t>(d))) {
            feature(fmt::format(R"({{"type":"LineString","coordinates":{}}})",
                                coords({hubs[c], hubs[d]}, false)),
                    fmt::format(R"("kind":"neighbour",{},"from":{},"to":{},)"
                                R"("cost":{:.1f})",
                                osm, c, d,
                                0.5 * (cells->cost_[c] + cells->cost_[d])));
          }
        }
      }
    }

    for (auto i = std::size_t{0U}; i != connectors.size(); ++i) {
      auto cell = std::string{};
      if (cells != nullptr && cells->connector_cell_[i] == osr::area_cells::kNoCell) {
        cell = R"(,"unreachable":true)";
      } else if (cells != nullptr) {
        auto const c = cells->connector_cell_[i];
        cell = fmt::format(R"(,"cell":{})", c);
        feature(fmt::format(R"({{"type":"LineString","coordinates":{}}})",
                            coords({connectors[i], hubs[c]}, false)),
                fmt::format(R"("kind":"spoke",{},"cell":{},"cost":{:.1f})",
                            osm, c, 0.5 * cells->cost_[c]));
      }
      feature(fmt::format(R"({{"type":"Point","coordinates":{}}})",
                          coord(connectors[i])),
              fmt::format(R"("kind":"connector",{},"node":{},"interior":{}{})",
                          osm, connector_ids[i], on_ring[i] == 0, cell));
    }
  }

  void write(std::filesystem::path const& p) const {
    auto out = std::ofstream{p};
    out << R"({"type":"FeatureCollection","features":[)";
    for (auto i = std::size_t{0U}; i != features_.size(); ++i) {
      out << (i == 0U ? "" : ",\n") << features_[i];
    }
    out << "]}\n";
    fmt::print("wrote {} features to {}\n", features_.size(), p.string());
    if (!all_errors_.empty()) {
      auto const e = quantiles_of(all_errors_);
      auto const c = quantiles_of(cells_per_area_);
      fmt::print(
          "  the model as written ({}): cells/area p50 {:.0f} p90 {:.0f} | "
          "route err p50 {:.1f} p90 {:.1f} p95 {:.1f}\n",
          g_tree_neighbours ? "one neighbour pair per cut"
                            : "cells neighbour where they share a border",
          c.p50_, c.p90_, e.p50_, e.p90_, e.p95_);
    }
  }

  std::vector<std::string> features_;
  std::size_t n_areas_{0U};
  std::vector<double> all_errors_;
  std::vector<double> cells_per_area_;
};

int main(int argc, char** argv) {
  if (argc < 2) {
    fmt::print(stderr,
               "usage: {} [--strict-levels] [--cells-out <geojson>] "
               "[--cells-t <meters>] [--max-vertices <n>] "
               "[--no-interior-ways] [--place-square] [--no-merge] "
               "[--tree-neighbours] <osm-file> [more-osm-files...]\n",
               argv[0]);
    return 1;
  }

  // Options apply to the files after them.
  auto cells_out = std::optional<std::filesystem::path>{};
  auto cells_t = 10.0;
  auto dump = cells_dump{};
  for (auto arg = 1; arg != argc; ++arg) {
    auto const opt = std::string_view{argv[arg]};
    if (opt == "--strict-levels") {
      g_strict_levels = true;
      continue;
    }
    if (opt == "--cells-out" && arg + 1 != argc) {
      cells_out = argv[++arg];
      continue;
    }
    if (opt == "--cells-t" && arg + 1 != argc) {
      cells_t = std::stod(argv[++arg]);
      continue;
    }
    if (opt == "--no-interior-ways") {
      g_interior_ways = false;
      continue;
    }
    if (opt == "--place-square") {
      g_place_square = true;
      continue;
    }
    if (opt == "--no-merge") {
      g_merge = false;
      continue;
    }
    if (opt == "--tree-neighbours") {
      g_tree_neighbours = true;
      continue;
    }
    if (opt == "--max-vertices" && arg + 1 != argc) {
      g_max_vertices = std::stoul(argv[++arg]);
      continue;
    }
    auto const path = std::string{argv[arg]};
    fmt::print("\n===== {} =====\n",
               std::filesystem::path{path}.filename().string());

    // Layer 0 (osr/area/pipeline.h): the walkable areas and what is in them.
    auto const options =
        osr::area_options{.walkable_ = {.place_square_ = g_place_square},
                          .strict_levels_ = g_strict_levels,
                          .merge_ = g_merge,
                          .interior_ways_ = g_interior_ways,
                          .max_vertices_ = g_max_vertices,
                          .cells_threshold_ = cells_t,
                          .shared_borders_ = !g_tree_neighbours};
    auto const data = osr::collect_areas(path, options);
    if (g_merge) {
      auto const& m = data.merged_;
      fmt::print(
          "merged {} overlapping areas (same level and layer) into {}, "
          "{} groups left unmerged because the union failed\n",
          m.n_members_, m.n_merged_, m.n_failed_);
    }

    auto n_no_crossing = 0;
    auto n_reported = 0;
    auto n_boundary_connectors = 0;
    auto n_interior_connectors = 0;
    auto n_too_big = 0;
    auto big_vertex_counts = std::vector<double>{};
    auto network_ratios = std::vector<double>{};
    auto n_stranded_connectors = 0;
    auto n_no_interior_ways = 0;
    auto n_cheap_served = 0;
    auto n_no_ways_but_served = 0;
    auto n_rejected_interior = 0;
    auto n_portal_outside = 0;
    auto n_no_inner_rings = 0;
    auto tri_ratio = std::vector<double>{};
    auto tri_err = std::vector<double>{};
    auto n_triangulation_edges = std::vector<double>{};
    auto n_triangulation_verts = std::vector<double>{};
    auto n_simple_polygon = 0;
    auto by_cells = std::vector<std::vector<double>>(9);
    auto cell_err_by_depth =
        std::vector<std::vector<double>>(kMaxSplitDepth + 1U);
    auto route_err_by_depth =
        std::vector<std::vector<double>>(kMaxSplitDepth + 1U);
    auto cells_crossed_by_depth =
        std::vector<std::vector<double>>(kMaxSplitDepth + 1U);
    auto gf_route_err = std::map<std::size_t, std::vector<double>>{};
    auto og_route_err = std::map<std::size_t, std::vector<double>>{};
    auto or_cell_err = std::map<std::size_t, std::vector<double>>{};
    auto or_route_err = std::map<std::size_t, std::vector<double>>{};
    auto oa_cells = std::map<double, std::vector<double>>{};
    auto oa_route_err = std::map<double, std::vector<double>>{};
    auto rf_cells = std::map<double, std::vector<double>>{};
    auto rf_route_err = std::map<double, std::vector<double>>{};
    auto lib_cells = std::map<double, std::vector<double>>{};
    auto lib_route_err = std::map<double, std::vector<double>>{};
    auto lib_compared = 0;
    auto lib_mismatches = 0;
    auto lib_extra = 0;
    auto ad_cells = std::map<double, std::vector<double>>{};
    auto ad_route_err = std::map<double, std::vector<double>>{};
    auto ad_route_err_global = std::map<double, std::vector<double>>{};
    auto gd_cell_err = std::map<std::size_t, std::vector<double>>{};
    auto gd_route_err = std::map<std::size_t, std::vector<double>>{};
    auto sm_cell_err = std::map<std::size_t, std::vector<double>>{};
    auto sm_route_err = std::map<std::size_t, std::vector<double>>{};
    auto cl_cell_err = std::map<unsigned, std::vector<double>>{};
    auto cl_route_err = std::map<unsigned, std::vector<double>>{};
    auto n_areas_with_barrier = 0;
    auto n_barrier_edges = 0;
    auto n_areas_with_building = 0;
    auto n_buildings_inside = 0;
    auto n_restricted_connectors = 0;

    auto n_areas_with_stranded = 0;
    auto worst_network = std::vector<double>{};
    auto needs_mesh = std::vector<char>{};
    auto spreads = std::vector<double>{};
    auto spreads_grouped = std::vector<double>{};
    auto detours = std::vector<double>{};
    auto vertex_counts = std::vector<double>{};
    auto connector_counts = std::vector<double>{};
    auto mins = std::vector<double>{};
    auto maxes = std::vector<double>{};
    auto star_errors = std::vector<double>{};
    auto spreads_relevant = std::vector<double>{};
    auto star_errors_relevant = std::vector<double>{};
    struct worst_entry {
      double spread_;
      double longest_;
      osmium::object_id_type id_;
      bool from_way_;
      std::size_t n_vertices_;
      std::size_t n_connectors_;
      std::string name_;
    };
    auto worst = std::vector<worst_entry>{};

    for (auto const [area_i, a] : utl::enumerate(data.areas_)) {
      // Layers 1 to 3 (osr/area/pipeline.h). The geodesics are built for
      // every area, served or not, for the statistics below.
      auto const pa = osr::prepare_area(data, area_i, options, true,
                                        cells_out.has_value());
      auto const& rings = pa.rings_;
      auto const& barriers = pa.barriers_;
      auto connectors = std::vector<geo::latlng>{};
      auto connector_ids = std::vector<osmium::object_id_type>{};
      auto on_ring = std::vector<char>{};
      for (auto const& c : pa.connectors_) {
        connectors.push_back(c.pos_);
        connector_ids.push_back(c.node_);
        on_ring.push_back(c.on_ring_ ? 1 : 0);
      }
      n_boundary_connectors += pa.counts_.n_boundary_;
      n_restricted_connectors += pa.counts_.n_restricted_;
      n_interior_connectors += pa.counts_.n_interior_;
      n_rejected_interior += pa.counts_.n_rejected_interior_;

      if (pa.status_ == osr::area_status::kTooFewConnectors) {
        ++n_no_crossing;
        dump.area(a, rings, barriers, connectors, connector_ids, on_ring,
                  "too_few_connectors");
        continue;
      }
      // Epstein & Sack's counting/sampling DP is defined over subpolygons
      // P(i,j) of a SIMPLE polygon, so it does not survive holes - and every
      // barrier and building is a hole.
      auto const& buildings = data.buildings_[area_i];
      n_no_inner_rings += a.rings_.size() == 1U ? 1 : 0;
      n_simple_polygon +=
          (a.rings_.size() == 1U && barriers.empty() && buildings.empty()) ? 1
                                                                          : 0;
      if (!barriers.empty()) {
        ++n_areas_with_barrier;
        n_barrier_edges += static_cast<int>(barriers.size());
      }
      if (!buildings.empty()) {
        ++n_areas_with_building;
        n_buildings_inside += static_cast<int>(buildings.size());
      }

      auto const n_vertices = pa.n_vertices_;
      if (pa.status_ == osr::area_status::kTooBig) {
        ++n_too_big;
        big_vertex_counts.push_back(static_cast<double>(n_vertices));
        dump.area(a, rings, barriers, connectors, connector_ids, on_ring,
                  "too_big");
        continue;
      }

      // Stage 1 (the binary): are any ways mapped into this area at all?
      auto const has_interior_ways = !data.interior_edges_[area_i].empty();
      n_no_interior_ways += has_interior_ways ? 0 : 1;

      // The shortcut that needs no geometry (layer 1). Counted only: the
      // geodesics are built for every area anyway.
      auto const cheap_served = pa.served_without_geometry_;
      n_cheap_served += cheap_served ? 1 : 0;
      // An area with no ways mapped into it can still be adequately served by
      // its own ring - two adjacent entrances, or a thin shape. Sending it
      // straight to subdivision on the binary alone would be wasted work.
      n_no_ways_but_served += (!has_interior_ways && cheap_served) ? 1 : 0;

      auto const& g = *pa.geodesics_;

      auto all = std::vector<std::size_t>(connectors.size());
      for (auto i = std::size_t{0U}; i != all.size(); ++i) {
        all[i] = i;
      }
      auto const s = compute_pair_stats(g, all);
      if (pa.status_ == osr::area_status::kNoCrossing) {
        ++n_no_crossing;
        dump.area(a, rings, barriers, connectors, connector_ids, on_ring,
                  "no_crossing");
        continue;
      }

      auto const grouped = group_connectors(connectors);
      auto const sg = compute_pair_stats(g, grouped);

      ++n_reported;
      spreads.push_back(s.spread());
      if (sg.n_pairs_ != 0U) {
        spreads_grouped.push_back(sg.spread());
      }
      vertex_counts.push_back(static_cast<double>(n_vertices));
      connector_counts.push_back(static_cast<double>(connectors.size()));
      mins.push_back(s.min_);
      maxes.push_back(s.max_);
      star_errors.push_back(all.size() > kMaxStarConnectors
                                ? 0.0
                                : star_max_error(g, all, star_fit(g, all)));

      auto const n_c = connectors.size();
      auto const& geo_dist = pa.geodesic_distances_;

      // Only pairs a crossing actually shortens have to be modelled, and
      // those are the ones the model is judged on (layer 1).
      detours.insert(end(detours), begin(pa.relevance_.outline_detours_),
                     end(pa.relevance_.outline_detours_));
      auto const& relevant = pa.relevance_.relevant_;

      // Layer 1's verdict: served where every relevant pair of ring
      // connectors is within max_detour_ of its shortest walk over mapped
      // ways, and no interior connector is stranded.
      auto const& verdict = pa.verdict_;
      for (auto const d : verdict.detours_) {
        network_ratios.push_back(std::min(d, 1000.0));
      }
      n_stranded_connectors += static_cast<int>(verdict.n_stranded_);
      n_areas_with_stranded += verdict.n_stranded_ != 0U ? 1 : 0;
      if (verdict.n_relevant_ != 0U) {
        worst_network.push_back(std::min(verdict.worst_detour_, 1000.0));
      }
      // Stays index-aligned with spreads_relevant below.
      needs_mesh.push_back(verdict.served_ ? 0 : 1);

      switch (pa.status_) {
        case osr::area_status::kServed:
          dump.area(a, rings, barriers, connectors, connector_ids, on_ring,
                    "served");
          break;
        case osr::area_status::kMeshed:
          dump.area(a, rings, barriers, connectors, connector_ids, on_ring,
                    "meshed", &*pa.cells_, &geo_dist, &pa.regions_);
          break;
        case osr::area_status::kUnreachable:
          dump.area(a, rings, barriers, connectors, connector_ids, on_ring,
                    "unreachable_pair");
          break;
        default: break;
      }

      auto const sr = compute_pair_stats(g, all, &relevant, n_c);
      if (sr.n_pairs_ != 0U) {
        spreads_relevant.push_back(sr.spread());
        auto rel_idx = std::vector<std::size_t>{};
        for (auto i = std::size_t{0U}; i != n_c; ++i) {
          auto const used = [&] {
            for (auto j = std::size_t{0U}; j != n_c; ++j) {
              if (j != i && relevant[i * n_c + j] &&
                  g.distance(i, j) != osr::area_geodesics::kUnreachable) {
                return true;
              }
            }
            return false;
          }();
          if (used) {
            rel_idx.push_back(i);
          }
        }
        star_errors_relevant.push_back(
            rel_idx.size() > kMaxStarConnectors
                ? 0.0
                : star_max_error(g, rel_idx, star_fit(g, rel_idx)));
      } else {
        spreads_relevant.push_back(0.0);
        star_errors_relevant.push_back(0.0);
      }

      // Triangulation: how much does routing on triangle edges cost against
      // the geodesic? Capped by vertex count because the greedy build tests
      // O(n^2) candidate diagonals.
      if (needs_mesh.back() == 1 && n_vertices <= kMaxTriangulationVertices) {
        auto const tri = greedy_triangulation(rings, barriers);
        auto idx_of = osr::hash_map<std::uint64_t, std::size_t>{};
        auto const pos_key = [](geo::latlng const& c) {
          return (static_cast<std::uint64_t>(
                      static_cast<std::int64_t>(c.lat() * 1e7))
                  << 32U) ^
                 static_cast<std::uint64_t>(
                     static_cast<std::int64_t>(c.lng() * 1e7));
        };
        for (auto i = std::size_t{0U}; i != tri.pts_.size(); ++i) {
          idx_of.emplace(pos_key(tri.pts_[i]), i);
        }

        auto pts = tri.pts_;
        auto adj = std::vector<std::vector<std::pair<std::size_t, double>>>{};
        auto const node_of = [&](geo::latlng const& c) {
          if (auto const it = idx_of.find(pos_key(c)); it != end(idx_of)) {
            return it->second;
          }
          pts.push_back(c);
          idx_of.emplace(pos_key(c), pts.size() - 1U);
          return pts.size() - 1U;
        };
        auto conn_node = std::vector<std::size_t>{};
        for (auto const& c : connectors) {
          conn_node.push_back(node_of(c));
        }

        adj.resize(pts.size());
        auto const add = [&](std::size_t const u, std::size_t const v) {
          auto const d = geo::distance(pts[u], pts[v]);
          adj[u].emplace_back(v, d);
          adj[v].emplace_back(u, d);
        };
        for (auto const& [u, v] : tri.edges_) {
          add(u, v);
        }
        // A connector inside a triangle is not a polygon vertex, so it hangs
        // off the three nearest corners it can see - no new vertices in the
        // polygon, three edges each.
        for (auto i = std::size_t{0U}; i != connectors.size(); ++i) {
          if (conn_node[i] < tri.pts_.size()) {
            continue;
          }
          auto near = std::vector<std::pair<double, std::size_t>>{};
          for (auto v = std::size_t{0U}; v != tri.pts_.size(); ++v) {
            near.emplace_back(geo::distance(connectors[i], tri.pts_[v]), v);
          }
          std::ranges::partial_sort(
              near, begin(near) + std::min<std::size_t>(12U, near.size()));
          auto joined = 0;
          for (auto const& [d, v] : near) {
            if (joined == 3) {
              break;
            }
            if (osr::area_geodesics::is_segment_inside(rings, connectors[i],
                                                       tri.pts_[v], barriers)) {
              add(conn_node[i], v);
              ++joined;
            }
          }
        }

        n_triangulation_edges.push_back(static_cast<double>(tri.edges_.size()));
        n_triangulation_verts.push_back(static_cast<double>(tri.pts_.size()));

        auto const inf = std::numeric_limits<double>::infinity();
        using entry = std::pair<double, std::size_t>;
        for (auto i = std::size_t{0U}; i != connectors.size(); ++i) {
          auto dd = std::vector<double>(pts.size(), inf);
          auto q =
              std::priority_queue<entry, std::vector<entry>, std::greater<>>{};
          dd[conn_node[i]] = 0.0;
          q.emplace(0.0, conn_node[i]);
          while (!q.empty()) {
            auto const [cost, u] = q.top();
            q.pop();
            if (cost > dd[u]) {
              continue;
            }
            for (auto const& [v, w] : adj[u]) {
              if (auto const nx = cost + w; nx < dd[v]) {
                dd[v] = nx;
                q.emplace(nx, v);
              }
            }
          }
          for (auto j = i + 1U; j != connectors.size(); ++j) {
            auto const truth = g.distance(i, j);
            if (truth == osr::area_geodesics::kUnreachable || truth < 1.0 ||
                !relevant[i * n_c + j] || !std::isfinite(dd[conn_node[j]])) {
              continue;
            }
            tri_ratio.push_back(dd[conn_node[j]] / static_cast<double>(truth));
            tri_err.push_back(
                std::abs(dd[conn_node[j]] - static_cast<double>(truth)));
          }
        }
      }

      // End-to-end error: subdivide with the placeholder cut and see how the
      // error behaves as a route crosses more cells.
      if (needs_mesh.back() == 1 && connectors.size() >= 4U &&
          n_vertices + connectors.size() + 16U <= kMaxVertices) {
        auto portals = std::vector<geo::latlng>{};
        auto root = cell_tree{};
        root.connectors_.resize(connectors.size());
        for (auto i = std::size_t{0U}; i != connectors.size(); ++i) {
          root.connectors_[i] = i;
        }
        split(root, connectors, rings, portals, kMaxSplitDepth);

        if (!portals.empty()) {
          auto all_pts = connectors;
          all_pts.insert(end(all_pts), begin(portals), end(portals));
          auto const g2 = osr::area_geodesics{rings, all_pts, barriers};

          // Same connector set, grouped by spread instead of cut by position.
          // Portals go on a spanning tree over the groups, so the two schemes
          // are compared with the same number of cells and the same kind of
          // border structure.
          for (auto const m : {2U, 4U, 8U}) {
            if (connectors.size() < m * 2U) {
              continue;
            }
            auto groups = cluster_by_spread(g, connectors.size(), m, relevant);
            if (groups.size() < 2U) {
              continue;
            }

            auto centroid = [&](std::vector<std::size_t> const& v) {
              auto lat = 0.0;
              auto lng = 0.0;
              for (auto const i : v) {
                lat += connectors[i].lat();
                lng += connectors[i].lng();
              }
              auto const cnt = static_cast<double>(v.size());
              return geo::latlng{lat / cnt, lng / cnt};
            };

            auto cpts = connectors;
            auto cells2 = groups;
            // Prim over group centroids: one portal per tree edge.
            auto in_tree = std::vector<char>(groups.size(), 0);
            in_tree[0] = 1;
            for (auto added = std::size_t{1U}; added != groups.size();
                 ++added) {
              auto best_a = groups.size();
              auto best_b = groups.size();
              auto best_d = std::numeric_limits<double>::max();
              for (auto x = std::size_t{0U}; x != groups.size(); ++x) {
                if (in_tree[x] == 0) {
                  continue;
                }
                for (auto y = std::size_t{0U}; y != groups.size(); ++y) {
                  if (in_tree[y] != 0) {
                    continue;
                  }
                  auto const d =
                      geo::distance(centroid(groups[x]), centroid(groups[y]));
                  if (d < best_d) {
                    best_d = d;
                    best_a = x;
                    best_b = y;
                  }
                }
              }
              if (best_a == groups.size()) {
                break;
              }
              auto const ca = centroid(groups[best_a]);
              auto const cb = centroid(groups[best_b]);
              cells2[best_a].push_back(cpts.size());
              cells2[best_b].push_back(cpts.size());
              cpts.push_back(
                  {0.5 * (ca.lat() + cb.lat()), 0.5 * (ca.lng() + cb.lng())});
              in_tree[best_b] = 1;
            }

            if (cpts.size() == connectors.size() ||
                n_vertices + cpts.size() > kMaxVertices) {
              continue;
            }
            auto const g3 = osr::area_geodesics{rings, cpts, barriers};
            auto ok = true;
            for (auto i = connectors.size(); i != cpts.size() && ok; ++i) {
              ok = g3.is_connector_reachable(i);
            }
            if (!ok) {
              continue;
            }
            auto const r3 = evaluate_composition(g3, connectors.size(),
                                                 cpts.size(), cells2);
            for (auto const& smp : r3.samples_) {
              cl_route_err[m].push_back(smp.error_);
            }
            for (auto const e : r3.cell_errors_) {
              cl_cell_err[m].push_back(e);
            }
          }

          // Orthogonal cut, both at fixed depth (comparable to the others) and
          // adaptive (comparable to the design), at each portal budget.
          for (auto const np : {1U, 2U, 4U}) {
            g_portals_per_border = np;
            for (auto const mode : {0, 1}) {
              for (auto const t :
                   (mode == 0 ? std::vector<double>{0.0}
                              : std::vector<double>{10.0, 20.0, 40.0})) {
                auto portals_o = std::vector<geo::latlng>{};
                auto root_o = cell_tree{};
                root_o.connectors_ = root.connectors_;
                split_orthogonal(root_o, connectors, on_ring, portals_o,
                                 mode == 0 ? kMaxSplitDepth : kAdaptiveMaxDepth,
                                 g, relevant, n_c, t, g_portals_per_border);
                if (portals_o.empty()) {
                  continue;
                }
                auto pts_o = connectors;
                pts_o.insert(end(pts_o), begin(portals_o), end(portals_o));
                if (n_vertices + pts_o.size() > kMaxVertices) {
                  continue;
                }
                auto const go = osr::area_geodesics{rings, pts_o, barriers};
                auto ok_o = true;
                check_portals(root_o, go, ok_o);
                if (!ok_o) {
                  continue;
                }
                if (mode == 0) {
                  for (auto depth = std::size_t{1U}; depth <= kMaxSplitDepth;
                       ++depth) {
                    auto cells_o = std::vector<std::vector<std::size_t>>{};
                    collect_cells(root_o, depth, pts_o, cells_o);
                    auto const ro = evaluate_composition(go, connectors.size(),
                                                         pts_o.size(), cells_o);
                    for (auto const& smp : ro.samples_) {
                      or_route_err[depth + 10U * np].push_back(smp.error_);
                    }
                    for (auto const e : ro.cell_errors_) {
                      or_cell_err[depth + 10U * np].push_back(e);
                    }
                    for (auto const& smp :
                         evaluate_composition(go, connectors.size(),
                                              pts_o.size(), cells_o, true)
                             .samples_) {
                      og_route_err[depth + 10U * np].push_back(smp.error_);
                    }
                  }
                } else {
                  auto cells_o = std::vector<std::vector<std::size_t>>{};
                  collect_cells(root_o, kAdaptiveMaxDepth, pts_o, cells_o);
                  auto const ro = evaluate_composition(go, connectors.size(),
                                                       pts_o.size(), cells_o);
                  if (ro.samples_.empty()) {
                    continue;
                  }
                  oa_cells[t + 1000.0 * np].push_back(
                      static_cast<double>(cells_o.size()));
                  for (auto const& smp : ro.samples_) {
                    oa_route_err[t + 1000.0 * np].push_back(smp.error_);
                  }
                }
              }
            }
          }
          g_portals_per_border = 1U;

          // Adaptive subdivision driven by the GLOBAL FIT's own residuals.
          //
          // The tree and its portals are built once to full depth and the
          // geodesics computed once; refinement then only flips expanded_
          // flags. A cell is split when it carries a route the fitted model
          // gets wrong by more than T - so T bounds what the router
          // experiences, rather than a per-cell statistic that the fit has
          // already been shown not to depend on.
          {
            // Portals as tokens: no positions, so no second oracle run, no
            // fallback for one landing outside the polygon, and nothing to
            // check. `g` - built from the connectors alone - is reused.
            auto portals_f = std::vector<geo::latlng>{};
            auto root_f = cell_tree{};
            root_f.connectors_ = root.connectors_;
            split_smart(root_f, connectors, rings, portals_f, kAdaptiveMaxDepth,
                        g, relevant, n_c, 0.0);
            auto const n_tokens = connectors.size() + portals_f.size();
            if (!portals_f.empty()) {
              // Portals keep their positions - collect_cells needs them to
              // decide which leaf on each side a border belongs to, which is
              // the topology itself. They are NOT added to the oracle: their
              // distances are never read, so `g` (connectors only) is reused
              // and the second area_geodesics build per strategy disappears.
              auto const& gf = g;
              auto pts_f = connectors;
              pts_f.insert(end(pts_f), begin(portals_f), end(portals_f));
              {
                auto dist_m = std::vector<float>(n_c * n_c);
                for (auto i = std::size_t{0U}; i != n_c; ++i) {
                  for (auto j = std::size_t{0U}; j != n_c; ++j) {
                    dist_m[i * n_c + j] = g.distance(i, j);
                  }
                }
                for (auto const t : {10.0, 20.0, 40.0}) {
                  reset_expanded(root_f);
                  auto cells_f = std::vector<std::vector<std::size_t>>{};
                  auto owner = std::vector<cell_tree*>{};
                  auto res = composition_result{};
                  for (auto it = std::size_t{0U}; it <= kAdaptiveMaxDepth;
                       ++it) {
                    cells_f.clear();
                    owner.clear();
                    collect_frontier(root_f, pts_f, cells_f, owner);
                    res = evaluate_composition(gf, connectors.size(),
                                               pts_f.size(), cells_f, true, t);
                    auto refined = false;
                    if (res.cell_violations_.size() != owner.size()) {
                      break;  // unreachable pair: no usable fit here
                    }
                    for (auto ci = std::size_t{0U}; ci != owner.size(); ++ci) {
                      if (res.cell_violations_[ci] > 0 && !owner[ci]->leaf()) {
                        owner[ci]->expanded_ = true;
                        refined = true;
                      }
                    }
                    if (!refined) {
                      break;
                    }
                  }

                  // The library must reproduce this area exactly: same cells,
                  // same error for every pair (up to float storage of costs).
                  auto const lib = run_library(connectors, dist_m, relevant, t);
                  auto const old_ok = !res.samples_.empty();
                  auto const lib_ok = lib.has_value() && !lib->errors_.empty();
                  auto max_diff = 0.0;
                  // Where the old code gives up on an unreachable pair, the
                  // library meshes the reachable parts: nothing to compare.
                  auto same = !old_ok || lib_ok;
                  lib_extra += (!old_ok && lib_ok) ? 1 : 0;
                  if (old_ok && lib_ok) {
                    same = lib->n_cells_ == cells_f.size() &&
                           lib->errors_.size() == res.samples_.size();
                    for (auto i = std::size_t{0U};
                         same && i != res.samples_.size(); ++i) {
                      max_diff = std::max(
                          max_diff,
                          std::abs(lib->errors_[i] - res.samples_[i].error_));
                    }
                    same = same && max_diff <= 1e-2;
                  }
                  ++lib_compared;
                  if (!same) {
                    if (++lib_mismatches <= 10) {
                      fmt::print(stderr,
                                 "library mismatch: area {} T={} cells {} vs "
                                 "{}, pairs {} vs {}, max diff {}\n",
                                 a.id_, t, cells_f.size(),
                                 lib ? lib->n_cells_ : 0U, res.samples_.size(),
                                 lib ? lib->errors_.size() : 0U, max_diff);
                    }
                  }
                  if (old_ok && lib_ok) {
                    lib_cells[t].push_back(static_cast<double>(lib->n_cells_));
                    for (auto const e : lib->errors_) {
                      lib_route_err[t].push_back(e);
                    }
                  }

                  if (res.samples_.empty()) {
                    continue;
                  }
                  rf_cells[t].push_back(static_cast<double>(cells_f.size()));
                  for (auto const& smp : res.samples_) {
                    rf_route_err[t].push_back(smp.error_);
                  }
                }
              }
            }
          }

          // Adaptive subdivision: cut until every cell is within T, then stop.
          for (auto const t : {10.0, 20.0, 40.0}) {
            auto portals_t = std::vector<geo::latlng>{};
            auto root_t = cell_tree{};
            root_t.connectors_ = root.connectors_;
            split_smart(root_t, connectors, rings, portals_t, kAdaptiveMaxDepth,
                        g, relevant, n_c, t);
            auto pts_t = connectors;
            pts_t.insert(end(pts_t), begin(portals_t), end(portals_t));
            if (n_vertices + pts_t.size() > kMaxVertices) {
              continue;
            }
            auto const gt = portals_t.empty()
                                ? g
                                : osr::area_geodesics{rings, pts_t, barriers};
            if (!portals_t.empty()) {
              resolve_portals(root_t, gt);
              auto ok_t = true;
              check_portals(root_t, gt, ok_t);
              if (!ok_t) {
                continue;
              }
            }
            auto cells_t = std::vector<std::vector<std::size_t>>{};
            collect_cells(root_t, kAdaptiveMaxDepth, pts_t, cells_t);
            auto const rt = evaluate_composition(gt, connectors.size(),
                                                 pts_t.size(), cells_t);
            if (rt.samples_.empty()) {
              continue;
            }
            ad_cells[t].push_back(static_cast<double>(cells_t.size()));
            for (auto const& smp : rt.samples_) {
              ad_route_err[t].push_back(smp.error_);
            }
            for (auto const& smp :
                 evaluate_composition(gt, connectors.size(), pts_t.size(),
                                      cells_t, true)
                     .samples_) {
              ad_route_err_global[t].push_back(smp.error_);
            }
          }

          // Geodesic cut: no extra points at all, so no second oracle run.
          {
            auto root3 = cell_tree{};
            root3.connectors_ = root.connectors_;
            split_geodesic(root3, connectors, on_ring, kMaxSplitDepth, g,
                           relevant, n_c);
            if (root3.lhs_ != nullptr) {
              for (auto depth = std::size_t{1U}; depth <= kMaxSplitDepth;
                   ++depth) {
                auto cells3 = std::vector<std::vector<std::size_t>>{};
                collect_cells(root3, depth, connectors, cells3);
                auto const r5 = evaluate_composition(g, connectors.size(),
                                                     connectors.size(), cells3);
                for (auto const& smp : r5.samples_) {
                  gd_route_err[depth].push_back(smp.error_);
                }
                for (auto const e : r5.cell_errors_) {
                  gd_cell_err[depth].push_back(e);
                }
              }
            }
          }

          // Same machinery, cut by the spread criterion instead of by median.
          {
            auto portals2 = std::vector<geo::latlng>{};
            auto root2 = cell_tree{};
            root2.connectors_ = root.connectors_;
            split_smart(root2, connectors, rings, portals2, kMaxSplitDepth, g,
                        relevant, n_c, 0.0);
            if (!portals2.empty() &&
                n_vertices + connectors.size() + portals2.size() <=
                    kMaxVertices) {
              auto pts2 = connectors;
              pts2.insert(end(pts2), begin(portals2), end(portals2));
              auto const g4 = osr::area_geodesics{rings, pts2, barriers};
              resolve_portals(root2, g4);
              auto ok2 = true;
              check_portals(root2, g4, ok2);
              if (ok2) {
                for (auto depth = std::size_t{1U}; depth <= kMaxSplitDepth;
                     ++depth) {
                  auto cells2 = std::vector<std::vector<std::size_t>>{};
                  collect_cells(root2, depth, pts2, cells2);
                  auto const r4 = evaluate_composition(g4, connectors.size(),
                                                       pts2.size(), cells2);
                  for (auto const& smp : r4.samples_) {
                    sm_route_err[depth].push_back(smp.error_);
                  }
                  for (auto const e : r4.cell_errors_) {
                    sm_cell_err[depth].push_back(e);
                  }
                  for (auto const& smp :
                       evaluate_composition(g4, connectors.size(), pts2.size(),
                                            cells2, true)
                           .samples_) {
                    gf_route_err[depth].push_back(smp.error_);
                  }
                }
              }
            }
          }

          resolve_portals(root, g2);
          auto usable = true;
          check_portals(root, g2, usable);
          if (!usable) {
            ++n_portal_outside;
          } else {
            for (auto depth = std::size_t{1U}; depth <= kMaxSplitDepth;
                 ++depth) {
              auto cells = std::vector<std::vector<std::size_t>>{};
              collect_cells(root, depth, all_pts, cells);
              auto const r = evaluate_composition(g2, connectors.size(),
                                                  all_pts.size(), cells);
              for (auto const& smp : r.samples_) {
                if (smp.n_cells_ >= 1 && smp.n_cells_ <= 8) {
                  by_cells[static_cast<std::size_t>(smp.n_cells_)].push_back(
                      smp.error_);
                }
                route_err_by_depth[depth].push_back(smp.error_);
                cells_crossed_by_depth[depth].push_back(
                    static_cast<double>(smp.n_cells_));
              }
              for (auto const e : r.cell_errors_) {
                cell_err_by_depth[depth].push_back(e);
              }
            }
          }
        }
      }

      worst.push_back(worst_entry{s.spread(), s.max_, a.id_, a.from_way_,
                                  n_vertices, connectors.size(), a.name_});
    }

    if (cells_out.has_value()) {
      dump.write(*cells_out);
    }

    fmt::print(
        "areas with a crossing: {}   (skipped, <2 usable connectors: {})\n",
        n_reported, n_no_crossing);
    fmt::print("  level rule: {}\n",
               g_strict_levels ? "strict (a levelled area needs levelled "
                                 "connectors)"
                               : "permissive (untagged joins anything)");
    fmt::print("  level tags: {}/{} areas, {}/{} linear ways\n",
               data.n_areas_with_level_, data.areas_.size(),
               data.n_ways_with_level_, data.n_ways_);
    if (n_too_big != 0) {
      auto const p = quantiles_of(big_vertex_counts);
      fmt::print(
          "  declined, over {} vertices: {} areas "
          "(p50 {:.0f}, max {:.0f} vertices)\n",
          g_max_vertices, n_too_big, p.p50_, p.max_);
    }
    fmt::print(
        "  connectors: {} on the boundary, {} interior "
        "(stairs/lifts/stubs ending inside)\n",
        n_boundary_connectors, n_interior_connectors);
    if (n_reported == 0) {
      continue;
    }

    auto const q = [](char const* label, std::vector<double> const& v,
                      char const* unit) {
      auto const p = quantiles_of(v);
      fmt::print(
          "  {:<26} p50 {:7.1f}  p90 {:7.1f}  p95 {:7.1f}  max {:8.1f} {}\n",
          label, p.p50_, p.p90_, p.p95_, p.max_, unit);
    };
    auto spreads_half = std::vector<double>{};
    for (auto const sp : spreads) {
      spreads_half.push_back(0.5 * sp);
    }
    q("ring vertices", vertex_counts, "");
    q("connectors", connector_counts, "");
    q("shortest crossing", mins, "m");
    q("longest crossing", maxes, "m");
    q("spread (max-min pair)", spreads, "m");
    q("spread, 3m grouping", spreads_grouped, "m");
    q("perimeter / geodesic", detours, "x");
    q("existing ways / geodesic", network_ratios, "x");
    q("  worst per area", worst_network, "x");

    for (auto const t : {1.1, 1.2, 1.5, 2.0}) {
      auto const served = std::ranges::count_if(
          worst_network, [&](double const r) { return r <= t; });
      fmt::print(
          "  already served within {:.1f}x by existing ways: "
          "{:4}/{:<4} areas\n",
          t, served, worst_network.size());
    }
    fmt::print(
        "  simple polygons (no holes at all): {}; no inner rings but "
        "barriers/buildings added: {}\n",
        n_simple_polygon, n_no_inner_rings - n_simple_polygon);
    fmt::print("  barriers inside: {} areas, {} fence/wall/hedge segments\n",
               n_areas_with_barrier, n_barrier_edges);
    fmt::print(
        "  buildings overlapping (not mapped as holes): {} areas, "
        "{} buildings\n",
        n_areas_with_building, n_buildings_inside);
    fmt::print(
        "  access-restricted connectors dropped: {} (from {} restricted linear "
        "ways)\n",
        n_restricted_connectors, data.n_restricted_ways_);
    fmt::print(
        "  rejected by layer (bridge over / tunnel under, no shared "
        "node): {} interior edges, {} interior connectors\n",
        data.n_rejected_edges_, n_rejected_interior);
    fmt::print("  no ways mapped into the area at all: {} areas\n",
               n_no_interior_ways);
    fmt::print(
        "  provably served by the straight-line bound (no visibility graph "
        "needed): {} areas ({} of them with no interior ways, served by the "
        "ring alone)\n",
        n_cheap_served, n_no_ways_but_served);
    fmt::print(
        "  interior connectors unreachable over existing ways: "
        "{} in {} areas\n",
        n_stranded_connectors, n_areas_with_stranded);
    auto const n_needs = std::ranges::count(needs_mesh, char{1});
    fmt::print("  -> {} areas need meshing, {} do not\n", n_needs,
               n_reported - n_needs);
    q("max err, uniform scalar", spreads_half, "m");
    q("max err, per-connector r", star_errors, "m");
    auto spreads_rel_half = std::vector<double>{};
    for (auto const sp : spreads_relevant) {
      spreads_rel_half.push_back(0.5 * sp);
    }
    q("  ... relevant pairs only", spreads_rel_half, "m");
    q("  ... rel. pairs + per-conn", star_errors_relevant, "m");

    for (auto const t : {5.0, 10.0, 20.0, 40.0}) {
      auto const ok = std::ranges::count_if(
          spreads, [&](double const s) { return s <= 2.0 * t; });
      auto const ok_g = std::ranges::count_if(
          spreads_grouped, [&](double const s) { return s <= 2.0 * t; });
      auto const ok_star = std::ranges::count_if(
          star_errors, [&](double const e) { return e <= t; });
      auto const ok_rel = std::ranges::count_if(
          spreads_relevant, [&](double const s) { return s <= 2.0 * t; });
      auto const ok_rel_star = std::ranges::count_if(
          star_errors_relevant, [&](double const e) { return e <= t; });
      auto ok_needed = 0;
      for (auto i = std::size_t{0U}; i != spreads_relevant.size(); ++i) {
        if (needs_mesh[i] == 1 && spreads_relevant[i] <= 2.0 * t) {
          ++ok_needed;
        }
      }
      fmt::print(
          "  T={:>4.0f}m: uniform {:4}/{:<4}  per-conn r {:4}"
          "  relevant-only {:4}  both {:4}  |  of the {} needing mesh:"
          " {:4}\n",
          t, ok, n_reported, ok_star, ok_rel, ok_rel_star, n_needs, ok_needed);
    }

    if (!tri_ratio.empty()) {
      auto const tr = quantiles_of(tri_ratio);
      auto const te = quantiles_of(tri_err);
      auto const tv = quantiles_of(n_triangulation_verts);
      auto const te2 = quantiles_of(n_triangulation_edges);
      fmt::print("  triangulation ({} areas, <= {} vertices):\n",
                 n_triangulation_verts.size(), kMaxTriangulationVertices);
      fmt::print(
          "    path/geodesic  p50 {:5.2f}x p90 {:5.2f}x p95 {:5.2f}x "
          "max {:6.2f}x\n",
          tr.p50_, tr.p90_, tr.p95_, tr.max_);
      fmt::print(
          "    error          p50 {:6.1f} p90 {:6.1f} p95 {:6.1f} "
          "max {:7.1f} m\n",
          te.p50_, te.p90_, te.p95_, te.max_);
      fmt::print(
          "    vertices       p50 {:5.0f} | edges p50 {:5.0f} p90 "
          "{:5.0f}\n",
          tv.p50_, te2.p50_, te2.p90_);
    }

    fmt::print("  per-cell error vs route error, by subdivision depth:\n");
    for (auto d = std::size_t{1U}; d <= kMaxSplitDepth; ++d) {
      if (route_err_by_depth[d].empty()) {
        continue;
      }
      auto const ce = quantiles_of(cell_err_by_depth[d]);
      auto const re = quantiles_of(route_err_by_depth[d]);
      auto const cc = quantiles_of(cells_crossed_by_depth[d]);
      fmt::print(
          "    depth {}: cell err p50 {:6.1f} p90 {:6.1f} | route err "
          "p50 {:6.1f} p90 {:6.1f} | amplification p50 {:4.2f}x "
          "p90 {:4.2f}x | cells crossed p50 {:.0f} p90 {:.0f}\n",
          d, ce.p50_, ce.p90_, re.p50_, re.p90_,
          ce.p50_ > 0.0 ? re.p50_ / ce.p50_ : 0.0,
          ce.p90_ > 0.0 ? re.p90_ / ce.p90_ : 0.0, cc.p50_, cc.p90_);
    }

    fmt::print("  orthogonal cut at geodesic midpoint:\n");
    for (auto const& [d, errs] : or_cell_err) {
      if (errs.empty()) {
        continue;
      }
      auto const ce = quantiles_of(errs);
      auto const re = quantiles_of(or_route_err[d]);
      auto const gq = quantiles_of(og_route_err[d]);
      fmt::print(
          "    {} portal(s) depth {}: cell err p50 {:6.1f} | route err local "
          "p50 {:6.1f} p90 {:6.1f} | GLOBAL FIT p50 {:6.1f} p90 {:6.1f}\n",
          d / 10U, d % 10U, ce.p50_, re.p50_, re.p90_, gq.p50_, gq.p90_);
    }
    for (auto const& [t, cells] : oa_cells) {
      if (cells.empty()) {
        continue;
      }
      auto const cq = quantiles_of(cells);
      auto const rq = quantiles_of(oa_route_err[t]);
      fmt::print(
          "    adaptive T={:4.0f}m: cells/area p50 {:4.0f} p90 {:4.0f} "
          "| route err p50 {:6.1f} p90 {:6.1f}\n",
          t, cq.p50_, cq.p90_, rq.p50_, rq.p90_);
    }

    fmt::print(
        "  adaptive on GLOBAL-FIT residuals (split cells the fitted "
        "model gets wrong by > T):\n");
    for (auto const& [t, cells] : rf_cells) {
      if (cells.empty()) {
        continue;
      }
      auto const cq = quantiles_of(cells);
      auto const rq = quantiles_of(rf_route_err[t]);
      fmt::print(
          "    T={:4.0f}m: cells/area p50 {:4.0f} p90 {:4.0f} max {:4.0f}"
          " | route err p50 {:6.1f} p90 {:6.1f} p95 {:6.1f}\n",
          t, cq.p50_, cq.p90_, cq.max_, rq.p50_, rq.p90_, rq.p95_);
    }

    fmt::print(
        "  same model via osr::build_area_cells (no portals, one fit) - "
        "{} area/T runs, {} mismatches, {} meshed only by the library:\n",
        lib_compared, lib_mismatches, lib_extra);
    for (auto const& [t, cells] : lib_cells) {
      if (cells.empty()) {
        continue;
      }
      auto const cq = quantiles_of(cells);
      auto const rq = quantiles_of(lib_route_err[t]);
      fmt::print(
          "    T={:4.0f}m: cells/area p50 {:4.0f} p90 {:4.0f} max {:4.0f}"
          " | route err p50 {:6.1f} p90 {:6.1f} p95 {:6.1f}\n",
          t, cq.p50_, cq.p90_, cq.max_, rq.p50_, rq.p90_, rq.p95_);
    }

    fmt::print("  adaptive subdivision (cut until every cell is within T):\n");
    for (auto const& [t, cells] : ad_cells) {
      if (cells.empty()) {
        continue;
      }
      auto const cq = quantiles_of(cells);
      auto const rq = quantiles_of(ad_route_err[t]);
      auto const gq = quantiles_of(ad_route_err_global[t]);
      fmt::print(
          "    T={:4.0f}m: cells/area p50 {:4.0f} p90 {:4.0f} max {:4.0f}"
          " | route err local p50 {:6.1f} p90 {:6.1f}"
          " | GLOBAL FIT p50 {:6.1f} p90 {:6.1f}\n",
          t, cq.p50_, cq.p90_, cq.max_, rq.p50_, rq.p90_, gq.p50_, gq.p90_);
    }

    fmt::print("  geodesic cut (no new points):\n");
    for (auto const& [d, errs] : gd_cell_err) {
      if (errs.empty()) {
        continue;
      }
      auto const ce = quantiles_of(errs);
      auto const re = quantiles_of(gd_route_err[d]);
      fmt::print(
          "    depth {}: cell err p50 {:6.1f} p90 {:6.1f} | route err "
          "p50 {:6.1f} p90 {:6.1f} | amplification p50 {:4.2f}x\n",
          d, ce.p50_, ce.p90_, re.p50_, re.p90_,
          ce.p50_ > 0.0 ? re.p50_ / ce.p50_ : 0.0);
    }
    fmt::print("  spread-guided geometric cut:\n");
    for (auto const& [d, errs] : sm_cell_err) {
      if (errs.empty()) {
        continue;
      }
      auto const ce = quantiles_of(errs);
      auto const re = quantiles_of(sm_route_err[d]);
      auto const gq = quantiles_of(gf_route_err[d]);
      fmt::print(
          "    depth {}: cell err p50 {:6.1f} | route err local p50 {:6.1f} "
          "p90 {:6.1f} | GLOBAL FIT p50 {:6.1f} p90 {:6.1f}\n",
          d, ce.p50_, re.p50_, re.p90_, gq.p50_, gq.p90_);
    }
    fmt::print("  spread-clustered connectors (free grouping):\n");
    for (auto const& [m, errs] : cl_cell_err) {
      if (errs.empty()) {
        continue;
      }
      auto const ce = quantiles_of(errs);
      auto const re = quantiles_of(cl_route_err[m]);
      fmt::print(
          "    {} cells: cell err p50 {:6.1f} p90 {:6.1f} | route err "
          "p50 {:6.1f} p90 {:6.1f} | amplification p50 {:4.2f}x\n",
          m, ce.p50_, ce.p90_, re.p50_, re.p90_,
          ce.p50_ > 0.0 ? re.p50_ / ce.p50_ : 0.0);
    }
    fmt::print(
        "  end-to-end error by cells crossed"
        " ({} areas skipped: portal outside the polygon):\n",
        n_portal_outside);
    for (auto k = std::size_t{1U}; k != by_cells.size(); ++k) {
      if (by_cells[k].empty()) {
        continue;
      }
      auto const q = quantiles_of(by_cells[k]);
      fmt::print(
          "    {} cell{}  n={:7}   p50 {:7.1f}  p90 {:7.1f}  "
          "p95 {:7.1f}  max {:8.1f} m\n",
          k, k == 1U ? " " : "s", by_cells[k].size(), q.p50_, q.p90_, q.p95_,
          q.max_);
    }

    std::ranges::sort(worst, [](auto const& a, auto const& b) {
      return a.spread_ > b.spread_;
    });
    fmt::print("  worst spreads:\n");
    for (auto i = std::size_t{0U}; i != std::min<std::size_t>(8U, worst.size());
         ++i) {
      auto const& w = worst[i];
      fmt::print(
          "    spread {:8.1f}m  longest {:8.1f}m  {}/{:<11} "
          "{:4} verts {:4} conn  {}\n",
          w.spread_, w.longest_, w.from_way_ ? "way" : "rel", w.id_,
          w.n_vertices_, w.n_connectors_, w.name_);
    }
  }

  return 0;
}
