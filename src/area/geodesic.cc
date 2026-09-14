#include "osr/area/geodesic.h"

#include <cmath>
#include <cstdint>
#include <algorithm>
#include <exception>
#include <map>
#include <optional>
#include <queue>
#include <set>
#include <utility>

#include "boost/geometry.hpp"

#include "geo/constants.h"

namespace osr {

namespace bg = boost::geometry;

namespace {

// The free space is snapped to this grid, in meters. Far below anything that
// matters for walking, far above double rounding: everything after the
// snapping is exact integer arithmetic.
constexpr auto kGrid = 1e-4;

// Half-width given to a barrier. A barrier drawn as a zero-width polyline
// cannot block anything: a path always has the option of pivoting on one of
// its vertices. The value only has to exceed the width of nothing; it must
// stay well under the narrowest passage worth keeping.
constexpr auto kObstacleHalfWidth = 0.15;

// A point this close to a corner of the free space is that corner - the same
// node, rounded differently.
constexpr auto kSnap = 1e-2;

// Everything blocked grows by this much before it is cut out of the outline.
// A building clipped against the outline has its edge within nanometres of the
// outline's; Boost cuts that cleanly, but snapping two near-identical edges to
// the grid can make them cross. Grown by a millimetre, they overlap by a
// millimetre instead and the cut stays clean after snapping.
constexpr auto kOverlap = 1e-3;

// A connector outside the free space - on the outline, which the wall
// clearance leaves just outside - steps onto it if it is at most this far.
constexpr auto kMaxStep = 1.0;

struct xy {
  double x_{0.0}, y_{0.0};
};

// Local equirectangular projection. Areas are at most a few hundred meters
// across, so straight lines stay straight and the distortion is orders of
// magnitude below the error thresholds any of this is measured against.
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

  geo::latlng inverse(xy const p) const {
    return {lat0_ + p.y_ / geo::kApproxDistanceLatDegrees,
            lng0_ + p.x_ / lng_scale_};
  }

  double lat0_, lng0_, lng_scale_;
};

// --- exact predicates on the grid -------------------------------------------

struct pt {
  std::int64_t x_{0}, y_{0};

  friend bool operator==(pt const&, pt const&) = default;
  friend auto operator<=>(pt const&, pt const&) = default;
};

pt operator-(pt const a, pt const b) { return {a.x_ - b.x_, a.y_ - b.y_}; }

using wide = __int128;

int sign(wide const v) { return v > 0 ? 1 : (v < 0 ? -1 : 0); }

wide cross(pt const a, pt const b) {
  return static_cast<wide>(a.x_) * b.y_ - static_cast<wide>(a.y_) * b.x_;
}

wide dot(pt const a, pt const b) {
  return static_cast<wide>(a.x_) * b.x_ + static_cast<wide>(a.y_) * b.y_;
}

// > 0: c lies left of the line a -> b.
int orient(pt const a, pt const b, pt const c) { return sign(cross(b - a, c - a)); }

// The interiors of the two segments cross at a single point.
bool properly_crosses(pt const p1, pt const p2, pt const q1, pt const q2) {
  return orient(p1, p2, q1) * orient(p1, p2, q2) < 0 &&
         orient(q1, q2, p1) * orient(q1, q2, p2) < 0;
}

// For v on the line a-b: strictly between a and b.
bool strictly_between(pt const a, pt const b, pt const v) {
  return dot(v - a, b - a) > 0 && dot(v - b, a - b) > 0;
}

// Free space seen from one point: the closed sector swept counter-clockwise
// from direction `from_` to direction `to_`. At a corner of a ring that has
// the free space on its left, `from_` points to the next corner and `to_` to
// the previous one.
struct cone {
  bool contains(pt const d) const {
    auto const turn = sign(cross(from_, to_));
    if (turn > 0) {  // less than a half turn
      return sign(cross(from_, d)) >= 0 && sign(cross(d, to_)) >= 0;
    }
    if (turn < 0) {  // more than a half turn: not strictly in the rest
      return !(sign(cross(to_, d)) > 0 && sign(cross(d, from_)) > 0);
    }
    if (dot(from_, to_) < 0) {  // a straight wall: the half plane left of it
      return sign(cross(from_, d)) >= 0;
    }
    return true;  // a full turn
  }

  // More than a half turn: a shortest path may bend here.
  bool reflex() const { return sign(cross(from_, to_)) < 0; }

  pt from_, to_;
};

// Where a point sits relative to the free space: outside, strictly inside, or
// on the boundary - then the directions into free space are the cones.
struct locus {
  bool inside_{false};
  std::vector<cone> cones_;  // empty: strictly inside

  bool free_towards(pt const d) const {
    return cones_.empty() || std::ranges::any_of(cones_, [&](cone const& c) {
             return c.contains(d);
           });
  }

  // Passing straight through: in along d1 reversed, out along d2. Both have
  // to be in the same cone - where the boundary is pinched to a point, going
  // from one cone to another is going through the pinch.
  bool passable(pt const d1, pt const d2) const {
    return cones_.empty() || std::ranges::any_of(cones_, [&](cone const& c) {
             return c.contains(d1) && c.contains(d2);
           });
  }
};

// --- building the free space -------------------------------------------------

using bpoint = bg::model::d2::point_xy<double>;
using bpolygon = bg::model::polygon<bpoint>;
using bmulti = bg::model::multi_polygon<bpolygon>;
using bline = bg::model::linestring<bpoint>;

// Outline minus holes and barriers, as Boost computes it. Nothing if the
// outline is not a valid polygon.
std::optional<bmulti> free_polygon(
    plane const& p,
    std::vector<std::vector<geo::latlng>> const& rings,
    std::vector<std::vector<geo::latlng>> const& obstacles,
    double const clearance) {
  auto const polygon_of = [&](std::vector<geo::latlng> const& r) {
    auto poly = bpolygon{};
    for (auto const& c : r) {
      auto const q = p(c);
      poly.outer().push_back(bpoint{q.x_, q.y_});
    }
    bg::correct(poly);
    return poly;
  };

  try {
    auto const outline = polygon_of(rings.front());
    if (!bg::is_valid(outline)) {
      return std::nullopt;
    }

    auto blocked = bmulti{};
    auto const block = [&](auto const& geometry) {
      auto next = bmulti{};
      bg::union_(blocked, geometry, next);
      blocked = std::move(next);
    };

    for (auto i = std::size_t{1U}; i < rings.size(); ++i) {
      if (rings[i].size() < 3U) {
        continue;
      }
      if (auto const hole = polygon_of(rings[i]); bg::is_valid(hole)) {
        block(hole);
      }
    }

    // A barrier becomes a band of kObstacleHalfWidth each side, lengthened by
    // as much at both ends so that a fence ending at a wall closes the gap.
    for (auto const& o : obstacles) {
      auto line = bline{};
      for (auto const& c : o) {
        auto const q = p(c);
        if (line.empty() || std::hypot(q.x_ - line.back().x(),
                                       q.y_ - line.back().y()) > 1e-6) {
          line.push_back(bpoint{q.x_, q.y_});
        }
      }
      if (line.size() < 2U) {
        continue;
      }
      auto const lengthen = [](bpoint& end, bpoint const& towards) {
        auto const dx = end.x() - towards.x();
        auto const dy = end.y() - towards.y();
        auto const len = std::hypot(dx, dy);
        end = bpoint{end.x() + dx / len * kObstacleHalfWidth,
                     end.y() + dy / len * kObstacleHalfWidth};
      };
      lengthen(line.front(), line[1]);
      lengthen(line.back(), line[line.size() - 2U]);

      auto band = bmulti{};
      bg::buffer(
          line, band,
          bg::strategy::buffer::distance_symmetric<double>{kObstacleHalfWidth},
          bg::strategy::buffer::side_straight{},
          bg::strategy::buffer::join_miter{},
          bg::strategy::buffer::end_flat{},
          bg::strategy::buffer::point_square{});
      block(band);
    }

    if (!blocked.empty()) {
      auto grown = bmulti{};
      bg::buffer(blocked, grown,
                 bg::strategy::buffer::distance_symmetric<double>{kOverlap},
                 bg::strategy::buffer::side_straight{},
                 bg::strategy::buffer::join_miter{},
                 bg::strategy::buffer::end_flat{},
                 bg::strategy::buffer::point_square{});
      blocked = std::move(grown);
    }

    auto free = bmulti{};
    bg::difference(outline, blocked, free);
    if (clearance > 0.0 && !free.empty()) {
      // Keep off the walls: the free space shrinks by the clearance - the
      // outline moves in, holes and barriers grow - so every wall is that far
      // from where one walks.
      auto shrunk = bmulti{};
      bg::buffer(free, shrunk,
                 bg::strategy::buffer::distance_symmetric<double>{-clearance},
                 bg::strategy::buffer::side_straight{},
                 bg::strategy::buffer::join_miter{},
                 bg::strategy::buffer::end_flat{},
                 bg::strategy::buffer::point_square{});
      free = std::move(shrunk);
    }
    return free;
  } catch (std::exception const&) {
    return std::nullopt;
  }
}

// Drops repeated points and spikes - a point where the ring turns straight
// back - which snapping can produce.
void clean_ring(std::vector<pt>& r) {
  auto changed = true;
  while (changed && r.size() >= 3U) {
    changed = false;
    for (auto i = std::size_t{0U}; i < r.size() && r.size() >= 3U;) {
      auto const n = r.size();
      auto const& a = r[(i + n - 1U) % n];
      auto const& b = r[i];
      auto const& c = r[(i + 1U) % n];
      if (b == c || (orient(a, b, c) == 0 && dot(b - a, c - b) <= 0)) {
        r.erase(begin(r) + static_cast<std::ptrdiff_t>(i));
        changed = true;
      } else {
        ++i;
      }
    }
  }
}

wide twice_area(std::vector<pt> const& r) {
  auto a = wide{0};
  for (auto i = std::size_t{0U}; i != r.size(); ++i) {
    a += cross(r[i], r[(i + 1U) % r.size()]);
  }
  return a;
}

// Is `q` in the closed region the rings bound? Exact.
bool covers(std::vector<std::vector<pt>> const& rings, pt const q) {
  auto in = false;
  for (auto const& r : rings) {
    for (auto i = std::size_t{0U}; i != r.size(); ++i) {
      auto const a = r[i];
      auto const b = r[(i + 1U) % r.size()];
      if (q == a || (orient(a, b, q) == 0 && strictly_between(a, b, q))) {
        return true;
      }
      if ((a.y_ > q.y_) != (b.y_ > q.y_)) {
        auto const s = orient(a, b, q);
        if (b.y_ > a.y_ ? s > 0 : s < 0) {
          in = !in;
        }
      }
    }
  }
  return in;
}

// The shortest step from `q` (in meters) onto the region, if one of at most
// kMaxStep exists: the grid point it reaches, and its length.
std::optional<std::pair<pt, double>> step_onto(
    std::vector<std::vector<pt>> const& rings, xy const q) {
  auto best = kMaxStep;
  auto found = false;
  auto target = xy{};
  for (auto const& r : rings) {
    for (auto i = std::size_t{0U}; i != r.size(); ++i) {
      auto const a = xy{static_cast<double>(r[i].x_) * kGrid,
                        static_cast<double>(r[i].y_) * kGrid};
      auto const& nb = r[(i + 1U) % r.size()];
      auto const b = xy{static_cast<double>(nb.x_) * kGrid,
                        static_cast<double>(nb.y_) * kGrid};
      auto const dx = b.x_ - a.x_;
      auto const dy = b.y_ - a.y_;
      auto const len2 = dx * dx + dy * dy;
      auto const t =
          len2 == 0.0
              ? 0.0
              : std::clamp(((q.x_ - a.x_) * dx + (q.y_ - a.y_) * dy) / len2,
                           0.0, 1.0);
      auto const c = xy{a.x_ + t * dx, a.y_ + t * dy};
      if (auto const d = std::hypot(c.x_ - q.x_, c.y_ - q.y_); d <= best) {
        best = d;
        target = c;
        found = true;
      }
    }
  }
  if (!found) {
    return std::nullopt;
  }
  // Onto the boundary, and on a hair further where rounding to the grid would
  // leave it just outside.
  auto const len = std::hypot(target.x_ - q.x_, target.y_ - q.y_);
  auto const ux = len == 0.0 ? 0.0 : (target.x_ - q.x_) / len;
  auto const uy = len == 0.0 ? 0.0 : (target.y_ - q.y_) / len;
  for (auto const further : {0.0, 1.0, 3.0, 10.0, 30.0}) {
    auto const g = pt{std::llround((target.x_ + ux * further * kGrid) / kGrid),
                      std::llround((target.y_ + uy * further * kGrid) / kGrid)};
    if (covers(rings, g)) {
      return std::pair{g, std::hypot(static_cast<double>(g.x_) * kGrid - q.x_,
                                     static_cast<double>(g.y_) * kGrid - q.y_)};
    }
  }
  return std::nullopt;
}

}  // namespace

// The free space on the grid: rings with the free space on their left (outer
// rings counter-clockwise, holes clockwise), and extra points - connectors -
// that are either corners of it or lie strictly inside or outside.
struct area_free_space::impl {
  impl(plane const& p,
       std::vector<std::vector<geo::latlng>> const& rings,
       std::vector<std::vector<geo::latlng>> const& obstacles,
       std::vector<geo::latlng> const& extra,
       double const clearance)
      : plane_{p}, clearance_{clearance} {
    if (rings.empty() || rings.front().size() < 3U) {
      return;
    }
    auto const poly = free_polygon(p, rings, obstacles, clearance);
    if (!poly.has_value()) {
      return;
    }
    valid_ = true;

    auto point_rings = std::vector<std::vector<pt>>{};
    auto const add_ring = [&](auto const& boost_ring, bool const outer) {
      auto r = std::vector<pt>{};
      for (auto const& q : boost_ring) {
        r.push_back(snap(xy{q.x(), q.y()}));
      }
      clean_ring(r);
      if (r.size() < 3U) {
        return;
      }
      auto const a = twice_area(r);
      if (a == 0) {
        return;
      }
      if ((a > 0) != outer) {
        std::reverse(begin(r), end(r));
      }
      point_rings.push_back(std::move(r));
    };
    for (auto const& q : *poly) {
      add_ring(q.outer(), true);
      for (auto const& inner : q.inners()) {
        add_ring(inner, false);
      }
    }

    for (auto const& r : point_rings) {
      for (auto const& q : r) {
        index_of(q);
      }
    }
    auto const n_boundary = pts_.size();
    for (auto const& c : extra) {
      auto q = snap(p(c));
      auto offset = 0.0;
      if (clearance > 0.0 && !covers(point_rings, q)) {
        if (auto const s = step_onto(point_rings, p(c)); s.has_value()) {
          q = s->first;
          offset = s->second;
        }
      }
      extra_offset_.push_back(offset);
      auto best = std::size_t{0U};
      auto best_d = kSnap;
      for (auto i = std::size_t{0U}; i != n_boundary; ++i) {
        auto const d = std::hypot(static_cast<double>(pts_[i].x_ - q.x_),
                                  static_cast<double>(pts_[i].y_ - q.y_)) *
                       kGrid;
        if (d <= best_d) {
          best = i;
          best_d = d;
        }
      }
      extra_.push_back(best_d < kSnap || (n_boundary != 0U && best_d == kSnap &&
                                          pts_[best] == q)
                           ? static_cast<std::uint32_t>(best)
                           : index_of(q));
    }

    point_rings_ = point_rings;

    // Node the rings: a point lying on an edge becomes a corner of it, so
    // that rings touch, and connectors sit on walls, only ever at corners.
    for (auto& r : point_rings) {
      auto noded = std::vector<std::uint32_t>{};
      for (auto i = std::size_t{0U}; i != r.size(); ++i) {
        auto const a = r[i];
        auto const b = r[(i + 1U) % r.size()];
        noded.push_back(index_of(a));
        auto on = std::vector<std::pair<wide, std::uint32_t>>{};
        for (auto v = std::uint32_t{0U}; v != pts_.size(); ++v) {
          if (orient(a, b, pts_[v]) == 0 && strictly_between(a, b, pts_[v])) {
            on.emplace_back(dot(pts_[v] - a, b - a), v);
          }
        }
        std::ranges::sort(on);
        for (auto const& [_, v] : on) {
          noded.push_back(v);
        }
      }
      rings_.push_back(std::move(noded));
    }

    cones_.resize(pts_.size());
    incident_.resize(pts_.size());
    auto leaving = std::vector<std::vector<pt>>(pts_.size());
    auto walls = std::vector<std::vector<pt>>(pts_.size());
    for (auto const& r : rings_) {
      for (auto i = std::size_t{0U}; i != r.size(); ++i) {
        auto const prev = r[(i + r.size() - 1U) % r.size()];
        auto const v = r[i];
        auto const next = r[(i + 1U) % r.size()];
        leaving[v].push_back(pts_[next] - pts_[v]);
        walls[v].push_back(pts_[next] - pts_[v]);
        walls[v].push_back(pts_[prev] - pts_[v]);
        incident_[v].push_back(static_cast<std::uint32_t>(edges_.size()));
        incident_[next].push_back(static_cast<std::uint32_t>(edges_.size()));
        edges_.emplace_back(v, next);
      }
    }

    // The free directions at a corner run counter-clockwise from each wall
    // leaving it to the next wall of any ring. Where rings touch at a corner,
    // that splits the free space around it into its separate wedges; each
    // ring's own sector would reach across the other ring and let a path
    // through the pinch.
    auto const turn_class = [](pt const from, pt const d) {
      auto const s = sign(cross(from, d));
      return s > 0 ? 0 : (s == 0 ? 1 : 2);  // under, exactly, over a half turn
    };
    for (auto v = std::size_t{0U}; v != pts_.size(); ++v) {
      for (auto const from : leaving[v]) {
        auto to = std::optional<pt>{};
        for (auto const d : walls[v]) {
          if (sign(cross(from, d)) == 0 && dot(from, d) > 0) {
            continue;  // this wall itself
          }
          if (!to.has_value() || turn_class(from, d) < turn_class(from, *to) ||
              (turn_class(from, d) == turn_class(from, *to) &&
               sign(cross(d, *to)) > 0)) {
            to = d;
          }
        }
        cones_[v].push_back(cone{.from_ = from, .to_ = to.value_or(from)});
      }
    }

    inside_.resize(pts_.size());
    for (auto v = std::size_t{0U}; v != pts_.size(); ++v) {
      inside_[v] = !cones_[v].empty() || strictly_inside(pts_[v]);
    }

    // Snapping can, in principle, make two walls that nearly touched cross.
    for (auto i = std::size_t{0U}; i != edges_.size() && crossing_free_; ++i) {
      for (auto j = i + 1U; j != edges_.size(); ++j) {
        auto const [a, b] = edges_[i];
        auto const [c, d] = edges_[j];
        if (properly_crosses(pts_[a], pts_[b], pts_[c], pts_[d])) {
          crossing_free_ = false;
          break;
        }
      }
    }
  }

  pt snap(xy const q) const {
    return {std::llround(q.x_ / kGrid), std::llround(q.y_ / kGrid)};
  }

  std::uint32_t index_of(pt const q) {
    auto const [it, fresh] =
        index_.try_emplace(q, static_cast<std::uint32_t>(pts_.size()));
    if (fresh) {
      pts_.push_back(q);
    }
    return it->second;
  }

  // Even-odd over all rings. Exact; for points on no edge.
  bool strictly_inside(pt const q) const {
    auto in = false;
    for (auto const& [u, w] : edges_) {
      auto const a = pts_[u];
      auto const b = pts_[w];
      if ((a.y_ > q.y_) != (b.y_ > q.y_)) {
        auto const s = orient(a, b, q);
        if (b.y_ > a.y_ ? s > 0 : s < 0) {
          in = !in;
        }
      }
    }
    return in;
  }

  locus locate(std::uint32_t const v) const {
    return {.inside_ = inside_[v] != 0, .cones_ = cones_[v]};
  }

  // Any point, whether a corner, on a wall, or neither.
  locus locate(pt const q) const {
    if (auto const it = index_.find(q); it != end(index_)) {
      return locate(it->second);
    }
    for (auto const& [u, w] : edges_) {
      if (orient(pts_[u], pts_[w], q) == 0 &&
          strictly_between(pts_[u], pts_[w], q)) {
        return {.inside_ = true,
                .cones_ = {cone{.from_ = pts_[w] - pts_[u],
                                .to_ = pts_[u] - pts_[w]}}};
      }
    }
    return {.inside_ = strictly_inside(q), .cones_ = {}};
  }

  // The segment a -> b stays in the free space: it leaves a and enters b
  // through free space, crosses no wall, and passes through every corner on
  // its way without going through a pinch.
  bool visible(pt const a,
               locus const& la,
               pt const b,
               locus const& lb) const {
    if (!la.inside_ || !lb.inside_) {
      return false;
    }
    if (a == b) {
      return true;
    }
    if (!la.free_towards(b - a) || !lb.free_towards(a - b)) {
      return false;
    }
    for (auto const& [u, w] : edges_) {
      if (properly_crosses(a, b, pts_[u], pts_[w])) {
        return false;
      }
    }
    for (auto v = std::size_t{0U}; v != pts_.size(); ++v) {
      auto const q = pts_[v];
      if (cones_[v].empty() || q == a || q == b || orient(a, b, q) != 0 ||
          !strictly_between(a, b, q)) {
        continue;
      }
      if (!locate(static_cast<std::uint32_t>(v)).passable(a - q, b - q)) {
        return false;
      }
    }
    return true;
  }

  bool visible(std::uint32_t const a, std::uint32_t const b) const {
    return visible(pts_[a], locate(a), pts_[b], locate(b));
  }

  // Which points `p` sees, by turning a ray around it once. The walls the ray
  // currently crosses are kept ordered by distance; a point on the ray is
  // visible if the nearest wall is beyond it. Points in line with `p` are
  // taken nearest first, each seen only if the one before it is and the way
  // on from it is free.
  std::vector<char> visible_from(std::uint32_t const p) const {
    auto out = std::vector<char>(pts_.size(), 0);
    if (!inside_[p]) {
      return out;
    }
    auto const o = pts_[p];
    auto const at_p = locate(p);

    auto order = std::vector<std::uint32_t>{};
    for (auto v = std::uint32_t{0U}; v != pts_.size(); ++v) {
      if (v != p && inside_[v]) {
        order.push_back(v);
      }
    }
    auto const upper = [](pt const d) {
      return d.y_ > 0 || (d.y_ == 0 && d.x_ > 0);
    };
    std::ranges::sort(order, [&](std::uint32_t const a, std::uint32_t const b) {
      auto const da = pts_[a] - o;
      auto const db = pts_[b] - o;
      if (upper(da) != upper(db)) {
        return upper(da);
      }
      if (auto const s = sign(cross(da, db)); s != 0) {
        return s > 0;
      }
      return dot(da, da) < dot(db, db);
    });

    // Of two walls the ray meets, which does it meet first? Walls do not
    // cross, so the answer is the same for every ray that meets both.
    auto const nearer = [&](std::uint32_t const e1, std::uint32_t const e2) {
      if (e1 == e2) {
        return false;
      }
      auto const a1 = pts_[edges_[e1].first];
      auto const b1 = pts_[edges_[e1].second];
      auto const a2 = pts_[edges_[e2].first];
      auto const b2 = pts_[edges_[e2].second];
      auto const s1 = orient(a1, b1, a2);
      auto const s2 = orient(a1, b1, b2);
      if (s1 == 0 && s2 == 0) {
        return e1 < e2;  // in line: the ray meets both only at a shared end
      }
      if (s1 * s2 >= 0) {  // e2 on one side of e1's line
        return (s1 != 0 ? s1 : s2) != orient(a1, b1, o);
      }
      auto const t1 = orient(a2, b2, a1);  // then e1 is on one side of e2's
      auto const t2 = orient(a2, b2, b1);
      return (t1 != 0 ? t1 : t2) == orient(a2, b2, o);
    };
    auto walls = std::set<std::uint32_t, decltype(nearer)>{nearer};

    // The ray starts along +x.
    for (auto e = std::uint32_t{0U}; e != edges_.size(); ++e) {
      auto const [u, w] = edges_[e];
      if (u == p || w == p) {
        continue;
      }
      auto const du = pts_[u] - o;
      auto const dw = pts_[w] - o;
      auto const su = sign(du.y_);
      auto const sw = sign(dw.y_);
      auto const on_ray = [](pt const d) { return d.y_ == 0 && d.x_ > 0; };
      if (su * sw < 0) {
        auto const lo = su < 0 ? du : dw;
        auto const hi = su < 0 ? dw : du;
        if (sign(cross(hi - lo, pt{} - lo)) > 0) {
          walls.insert(e);
        }
      } else if ((on_ray(du) && sw < 0) || (on_ray(dw) && su < 0)) {
        walls.insert(e);
      }
    }

    auto chain = std::vector<std::uint32_t>{};  // points on the current ray
    auto chain_visible = false;
    for (auto const w : order) {
      auto const d = pts_[w] - o;
      if (!chain.empty()) {
        auto const first = pts_[chain.front()] - o;
        if (sign(cross(first, d)) != 0 || dot(first, d) <= 0) {
          chain.clear();
        }
      }

      auto const at_w = locate(w);
      auto vis = at_p.free_towards(d) && at_w.free_towards(o - pts_[w]);
      if (vis && !chain.empty()) {
        auto const prev = chain.back();
        vis = chain_visible && locate(prev).passable(o - pts_[prev],
                                                     pts_[w] - pts_[prev]);
      }
      if (vis) {
        for (auto const e : walls) {
          auto const [u, x] = edges_[e];
          auto const touches = [&](std::uint32_t const c) {
            return u == c || x == c;
          };
          if (touches(w) || std::ranges::any_of(chain, touches)) {
            continue;  // meets the ray at w or before it, at a corner
          }
          vis = !properly_crosses(o, pts_[w], pts_[u], pts_[x]);
          break;
        }
      }
      out[w] = vis ? 1 : 0;

      // Turn the ray past w: first drop the walls it leaves behind, then take
      // on the ones ahead. Not the other way round - a wall behind and a wall
      // ahead that meet at w are never both crossed by the same ray, so there
      // is no nearer one between them, and letting them share the set for a
      // moment corrupts its order.
      for (auto const behind : {true, false}) {
        for (auto const e : incident_[w]) {
          auto const [u, x] = edges_[e];
          auto const other = u == w ? x : u;
          if (other == p) {
            continue;
          }
          auto const s = orient(o, pts_[w], pts_[other]);
          if (behind && s < 0) {
            walls.erase(e);
          } else if (!behind && s > 0) {
            walls.insert(e);
          }
        }
      }
      chain.push_back(w);
      chain_visible = vis;
    }
    return out;
  }

  plane plane_;
  double clearance_{0.0};
  bool valid_{false};
  bool crossing_free_{true};
  std::vector<std::vector<pt>> point_rings_;  // before noding
  std::vector<double> extra_offset_;  // step of each extra point onto it
  std::vector<pt> pts_;
  std::map<pt, std::uint32_t> index_;
  std::vector<std::uint32_t> extra_;  // point of each extra point
  std::vector<std::vector<std::uint32_t>> rings_;
  std::vector<std::pair<std::uint32_t, std::uint32_t>> edges_;
  std::vector<std::vector<cone>> cones_;
  std::vector<std::vector<std::uint32_t>> incident_;
  std::vector<char> inside_;
};

area_free_space::area_free_space(
    std::vector<std::vector<geo::latlng>> const& rings,
    std::vector<std::vector<geo::latlng>> const& obstacles,
    double const clearance)
    : impl_{rings.empty() || rings.front().empty()
                ? nullptr
                : std::make_unique<impl>(plane{rings.front().front()}, rings,
                                         obstacles, std::vector<geo::latlng>{},
                                         clearance)} {}

area_free_space::~area_free_space() = default;
area_free_space::area_free_space(area_free_space&&) noexcept = default;
area_free_space& area_free_space::operator=(area_free_space&&) noexcept =
    default;

bool area_free_space::valid() const { return impl_ != nullptr && impl_->valid_; }

std::vector<std::vector<geo::latlng>> area_free_space::rings() const {
  auto out = std::vector<std::vector<geo::latlng>>{};
  if (!valid()) {
    return out;
  }
  for (auto const& r : impl_->point_rings_) {
    auto& ring = out.emplace_back();
    for (auto const& q : r) {
      ring.push_back(
          impl_->plane_.inverse(xy{static_cast<double>(q.x_) * kGrid,
                                   static_cast<double>(q.y_) * kGrid}));
    }
  }
  return out;
}

bool area_free_space::is_segment_inside(geo::latlng const& a,
                                        geo::latlng const& b) const {
  if (!valid()) {
    return false;
  }
  auto const snapped = [&](geo::latlng const& c) {
    auto q = impl_->snap(impl_->plane_(c));
    if (impl_->clearance_ > 0.0 && !covers(impl_->point_rings_, q)) {
      if (auto const s = step_onto(impl_->point_rings_, impl_->plane_(c));
          s.has_value()) {
        q = s->first;
      }
    }
    for (auto const& v : impl_->pts_) {
      if (std::hypot(static_cast<double>(v.x_ - q.x_),
                     static_cast<double>(v.y_ - q.y_)) *
              kGrid <
          kSnap) {
        q = v;
        break;
      }
    }
    return q;
  };
  auto const pa = snapped(a);
  auto const pb = snapped(b);
  return impl_->visible(pa, impl_->locate(pa), pb, impl_->locate(pb));
}

bool area_geodesics::is_segment_inside(
    std::vector<std::vector<geo::latlng>> const& rings,
    geo::latlng const& a,
    geo::latlng const& b,
    std::vector<std::vector<geo::latlng>> const& obstacles) {
  return area_free_space{rings, obstacles}.is_segment_inside(a, b);
}

area_geodesics::area_geodesics(
    std::vector<std::vector<geo::latlng>> const& rings,
    std::vector<geo::latlng> const& connectors,
    std::vector<std::vector<geo::latlng>> const& obstacles,
    geodesic_options const options)
    : n_connectors_{connectors.size()}, connector_pos_{connectors} {
  dist_.assign(n_connectors_ * n_connectors_, kUnreachable);
  connector_vertex_.assign(n_connectors_, kNoVertex);
  if (rings.empty() || rings.front().size() < 3U || connectors.empty()) {
    return;
  }

  auto const p = plane{rings.front().front()};
  auto const fs =
      area_free_space::impl{p, rings, obstacles, connectors, options.clearance_};
  if (!fs.valid_) {
    return;
  }
  valid_ = true;
  offset_ = fs.extra_offset_;

  auto const n = fs.pts_.size();
  vertices_.resize(n);
  for (auto v = std::size_t{0U}; v != n; ++v) {
    vertices_[v] = p.inverse(xy{static_cast<double>(fs.pts_[v].x_) * kGrid,
                                static_cast<double>(fs.pts_[v].y_) * kGrid});
  }
  auto is_connector = std::vector<char>(n, 0);
  for (auto i = std::size_t{0U}; i != n_connectors_; ++i) {
    auto const v = fs.extra_[i];
    connector_vertex_[i] = v;
    if (is_connector[v] == 0) {
      if (offset_[i] == 0.0) {
        vertices_[v] = connectors[i];  // the node itself, not its grid point
      }
      is_connector[v] = 1;
    }
  }

  // A shortest path bends only at reflex corners, so only those and the
  // connectors are needed.
  auto keep = std::vector<char>(n, 0);
  for (auto v = std::size_t{0U}; v != n; ++v) {
    keep[v] = (is_connector[v] != 0 ||
               std::ranges::any_of(fs.cones_[v],
                                   [](cone const& c) { return c.reflex(); }))
                  ? 1
                  : 0;
  }

  swept_ =
      options.algorithm_ == visibility_algorithm::kSweep && fs.crossing_free_;
  for (auto a = vertex_idx_t{0U}; a != n; ++a) {
    if (keep[a] == 0) {
      continue;
    }
    if (swept_) {
      auto const seen = fs.visible_from(a);
      for (auto b = static_cast<vertex_idx_t>(a + 1U); b < n; ++b) {
        if (keep[b] != 0 && seen[b] != 0) {
          edges_.emplace_back(a, b);
        }
      }
    } else {
      for (auto b = static_cast<vertex_idx_t>(a + 1U); b < n; ++b) {
        if (keep[b] != 0 && fs.visible(a, b)) {
          edges_.emplace_back(a, b);
        }
      }
    }
  }

  auto adj = std::vector<std::vector<std::pair<vertex_idx_t, double>>>(n);
  for (auto const& [a, b] : edges_) {
    auto const d = std::hypot(static_cast<double>(fs.pts_[a].x_ - fs.pts_[b].x_),
                              static_cast<double>(fs.pts_[a].y_ - fs.pts_[b].y_)) *
                   kGrid;
    adj[a].emplace_back(b, d);
    adj[b].emplace_back(a, d);
  }

  pred_.assign(n_connectors_ * n, kNoVertex);
  auto d = std::vector<double>(n);
  using queue_entry = std::pair<double, vertex_idx_t>;
  for (auto i = std::size_t{0U}; i != n_connectors_; ++i) {
    auto const src = connector_vertex_[i];
    if (fs.inside_[src] == 0) {
      continue;
    }
    std::ranges::fill(d, std::numeric_limits<double>::infinity());
    auto* pred = &pred_[i * n];
    d[src] = 0.0;
    auto q = std::priority_queue<queue_entry, std::vector<queue_entry>,
                                 std::greater<>>{};
    q.emplace(0.0, src);
    while (!q.empty()) {
      auto const [cost, u] = q.top();
      q.pop();
      if (cost > d[u]) {
        continue;
      }
      for (auto const& [v, w] : adj[u]) {
        if (auto const next = cost + w; next < d[v]) {
          d[v] = next;
          pred[v] = u;
          q.emplace(next, v);
        }
      }
    }
    for (auto j = std::size_t{0U}; j != n_connectors_; ++j) {
      auto const dst = connector_vertex_[j];
      if (std::isfinite(d[dst])) {
        dist_[i * n_connectors_ + j] = static_cast<float>(
            i == j ? 0.0 : d[dst] + offset_[i] + offset_[j]);
      }
    }
  }
}

float area_geodesics::distance(std::size_t const i, std::size_t const j) const {
  return dist_[i * n_connectors_ + j];
}

bool area_geodesics::is_connector_reachable(std::size_t const i) const {
  if (i >= n_connectors_) {
    return false;
  }
  for (auto j = std::size_t{0U}; j != n_connectors_; ++j) {
    if (i != j && distance(i, j) != kUnreachable) {
      return true;
    }
  }
  return false;
}

std::vector<geo::latlng> area_geodesics::path(std::size_t const i,
                                              std::size_t const j) const {
  if (distance(i, j) == kUnreachable) {
    return {};
  }

  auto const n = vertices_.size();
  auto const src = connector_vertex_[i];
  auto out = std::vector<geo::latlng>{};
  for (auto v = connector_vertex_[j]; v != kNoVertex;
       v = (v == src ? kNoVertex : pred_[i * n + v])) {
    out.push_back(vertices_[v]);
  }
  std::reverse(begin(out), end(out));
  // A connector that stepped onto the free space starts and ends its paths
  // itself.
  if (offset_[i] > 0.0) {
    out.insert(begin(out), connector_pos_[i]);
  }
  if (offset_[j] > 0.0) {
    out.push_back(connector_pos_[j]);
  }
  return out;
}

float area_geodesics::mean_pair_distance() const {
  auto sum = 0.0;
  auto n = std::size_t{0U};
  for (auto i = std::size_t{0U}; i != n_connectors_; ++i) {
    for (auto j = i + 1U; j != n_connectors_; ++j) {
      if (auto const d = distance(i, j); d != kUnreachable) {
        sum += d;
        ++n;
      }
    }
  }
  return n == 0U ? std::numeric_limits<float>::quiet_NaN()
                 : static_cast<float>(sum / static_cast<double>(n));
}

float area_geodesics::max_pair_error() const {
  auto const mean = mean_pair_distance();
  if (std::isnan(mean)) {
    return std::numeric_limits<float>::quiet_NaN();
  }
  auto worst = 0.F;
  for (auto i = std::size_t{0U}; i != n_connectors_; ++i) {
    for (auto j = i + 1U; j != n_connectors_; ++j) {
      if (auto const d = distance(i, j); d != kUnreachable) {
        worst = std::max(worst, std::abs(mean - d));
      }
    }
  }
  return worst;
}

}  // namespace osr
