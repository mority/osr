#include "osr/area/geodesic.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <queue>
#include <utility>

#include "geo/constants.h"

#include "utl/enumerate.h"

#include "osr/types.h"

namespace osr {

namespace {

// Everything below works in a local metric plane. Tolerances are therefore
// lengths in meters, not coordinate deltas or areas, which keeps them
// meaningful independent of where on the globe the area sits.
constexpr auto kEps = 1e-6;  // 1 um - well below OSM coordinate resolution.
constexpr auto kMergeDistance = 1e-2;  // 1 cm - same node, different rounding.

// Half-width given to a barrier when it is turned into a hole.
//
// A barrier drawn as a zero-width polyline cannot actually block anything: a
// path always has the option of pivoting on one of its vertices, so a fence
// anchored to a wall gets stepped around exactly at the anchor. Giving it a
// width separates the two sides, which is why ppr buffers its obstacles too.
// The value only has to exceed the width of nothing; it must stay well under
// the narrowest passage worth keeping.
constexpr auto kObstacleHalfWidth = 0.15;

// How far beside a wall free space is looked for (see free_beside). Well
// under the narrowest passage worth keeping, well over coordinate noise.
constexpr auto kFreeBeside = 0.05;  // 5 cm

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

double signed_area(std::vector<xy> const& r) {
  auto a = 0.0;
  for (auto i = std::size_t{0U}; i != r.size(); ++i) {
    auto const& p = r[i];
    auto const& q = r[(i + 1U) % r.size()];
    a += p.x_ * q.y_ - q.x_ * p.y_;
  }
  return 0.5 * a;
}

double dist(xy const a, xy const b) {
  return std::hypot(b.x_ - a.x_, b.y_ - a.y_);
}

// Signed distance of `c` from the infinite line a->b, in meters. Normalizing
// the cross product by the segment length turns the orientation predicate into
// a length, so a single length tolerance covers it.
double side(xy const a, xy const b, xy const c) {
  auto const cross =
      (b.x_ - a.x_) * (c.y_ - a.y_) - (b.y_ - a.y_) * (c.x_ - a.x_);
  auto const len = dist(a, b);
  return len < kEps ? 0.0 : cross / len;
}

int sgn(double const v) { return v > kEps ? 1 : (v < -kEps ? -1 : 0); }

bool on_segment(xy const a, xy const b, xy const p) {
  if (sgn(side(a, b, p)) != 0) {
    return false;
  }
  auto const dx = b.x_ - a.x_;
  auto const dy = b.y_ - a.y_;
  auto const len = std::hypot(dx, dy);
  if (len < kEps) {
    return dist(a, p) <= kEps;
  }
  auto const s = ((p.x_ - a.x_) * dx + (p.y_ - a.y_) * dy) / len;
  return s >= -kEps && s <= len + kEps;
}

// Interiors intersect. Touching at endpoints or running along a shared line
// does not count - those cases are resolved by the midpoint tests instead.
bool properly_crosses(xy const p1, xy const p2, xy const q1, xy const q2) {
  auto const d1 = sgn(side(p1, p2, q1));
  auto const d2 = sgn(side(p1, p2, q2));
  auto const d3 = sgn(side(q1, q2, p1));
  auto const d4 = sgn(side(q1, q2, p2));
  return d1 * d2 < 0 && d3 * d4 < 0;
}

// A barrier polyline as a set of closed rings: one box per segment, extended
// at both ends so consecutive boxes overlap at the joints, plus a square at
// each interior joint so a sharp turn cannot leave a gap to slip through.
std::vector<std::vector<xy>> thicken(std::vector<xy> const& line,
                                     double const h) {
  auto out = std::vector<std::vector<xy>>{};
  for (auto i = std::size_t{0U}; i + 1U < line.size(); ++i) {
    auto const a = line[i];
    auto const b = line[i + 1U];
    auto const len = dist(a, b);
    if (len < kEps) {
      continue;
    }
    auto const ux = (b.x_ - a.x_) / len;
    auto const uy = (b.y_ - a.y_) / len;
    auto const ax = a.x_ - ux * h;
    auto const ay = a.y_ - uy * h;
    auto const bx = b.x_ + ux * h;
    auto const by = b.y_ + uy * h;
    out.push_back({xy{ax - uy * h, ay + ux * h}, xy{bx - uy * h, by + ux * h},
                   xy{bx + uy * h, by - ux * h}, xy{ax + uy * h, ay - ux * h}});
  }
  for (auto i = std::size_t{1U}; i + 1U < line.size(); ++i) {
    auto const& v = line[i];
    out.push_back({xy{v.x_ - h, v.y_ - h}, xy{v.x_ + h, v.y_ - h},
                   xy{v.x_ + h, v.y_ + h}, xy{v.x_ - h, v.y_ + h}});
  }
  return out;
}

struct polygon {
  // Every blocking edge. Obstacles have been turned into rings by this point,
  // so there is only one kind.
  template <typename Fn>
  void for_each_edge(Fn&& f) const {
    for (auto const& r : rings_) {
      for (auto i = std::size_t{0U}; i != r.size(); ++i) {
        f(r[i], r[(i + 1U) % r.size()]);
      }
    }
  }

  template <typename Fn>
  void for_each_vertex(Fn&& f) const {
    for (auto const& r : rings_) {
      for (auto const& v : r) {
        f(v);
      }
    }
  }

  bool on_boundary(xy const p) const {
    for (auto const& r : rings_) {
      for (auto i = std::size_t{0U}; i != r.size(); ++i) {
        auto const& a = r[i];
        auto const& b = r[(i + 1U) % r.size()];
        if (on_segment(a, b, p)) {
          return true;
        }
      }
    }
    return false;
  }

  // Even-odd ray cast over all rings at once: with holes wound as separate
  // rings, an odd crossing count means inside the outer ring and outside every
  // hole, whatever the ring orientations are.
  bool strictly_inside(xy const p) const {
    auto inside = false;
    for (auto const& r : rings_) {
      for (auto i = std::size_t{0U}; i != r.size(); ++i) {
        auto const& a = r[i];
        auto const& b = r[(i + 1U) % r.size()];
        if ((a.y_ > p.y_) != (b.y_ > p.y_)) {
          auto const x = a.x_ + (p.y_ - a.y_) / (b.y_ - a.y_) * (b.x_ - a.x_);
          if (p.x_ < x) {
            inside = !inside;
          }
        }
      }
    }
    return inside;
  }

  // Walls count as walkable, so a point on the boundary is inside.
  bool contains(xy const p) const {
    return on_boundary(p) || strictly_inside(p);
  }

  // Is there free space right beside `p`, on either side of the line a->b?
  // A wall is walkable along its length only where one can stand next to
  // it. Where a hole lies against the outline, their shared edge has the
  // hole on one side and the outside on the other: walking along it is
  // walking through the building. That happens wherever an area's outline
  // runs through a building - under an arcade, say - and the building is
  // clipped to it.
  bool free_beside(xy const a, xy const b, xy const p) const {
    auto const len = dist(a, b);
    if (len < kEps) {
      return strictly_inside(p);
    }
    auto const nx = -(b.y_ - a.y_) / len * kFreeBeside;
    auto const ny = (b.x_ - a.x_) / len * kFreeBeside;
    return strictly_inside(xy{p.x_ + nx, p.y_ + ny}) ||
           strictly_inside(xy{p.x_ - nx, p.y_ - ny});
  }

  // Is the straight segment a->b fully contained in the polygon?
  //
  // Two conditions, and both are needed: no boundary edge may cross it (that
  // would take it out through a wall), and it may not slip through a vertex
  // into the outside - crossing a concave notch or squeezing between two holes
  // that touch at a point never properly crosses an edge. The second case is
  // caught by splitting the segment at every vertex it passes through and
  // testing each piece's midpoint.
  bool is_visible(xy const a, xy const b) const {
    auto const len = dist(a, b);
    if (len < kEps) {
      return true;
    }

    auto blocked = false;
    for_each_edge([&](xy const p, xy const q) {
      blocked = blocked || properly_crosses(a, b, p, q);
    });
    if (blocked) {
      return false;
    }

    auto const dx = (b.x_ - a.x_) / len;
    auto const dy = (b.y_ - a.y_) / len;

    // Vertices the segment passes straight through.
    auto pierced = std::vector<xy>{};
    auto touch = std::vector<double>{0.0, len};
    for_each_vertex([&](xy const v) {
      if (sgn(side(a, b, v)) != 0) {
        return;
      }
      auto const s = (v.x_ - a.x_) * dx + (v.y_ - a.y_) * dy;
      if (s > kEps && s < len - kEps) {
        touch.push_back(s);
        pierced.push_back(v);
      }
    });
    std::sort(begin(touch), end(touch));

    // Passing through a vertex is only allowed if the blocking edges meeting
    // there all lie to one side - that is what going round the free end of a
    // fence looks like. If they lie on both sides the path is going through
    // the barrier, not round it: a fence anchored to a wall blocks movement
    // along the wall, and a boundary pinched to a point is not a doorway.
    for (auto const& v : pierced) {
      auto left = false;
      auto right = false;
      for_each_edge([&](xy const p, xy const q) {
        auto const at_p = dist(p, v) <= kMergeDistance;
        auto const at_q = dist(q, v) <= kMergeDistance;
        if (!at_p && !at_q) {
          return;
        }
        auto const far = at_p ? q : p;
        auto const sd = sgn(side(a, b, far));
        left = left || sd > 0;
        right = right || sd < 0;
      });
      if (left && right) {
        return false;
      }
    }

    for (auto i = std::size_t{1U}; i != touch.size(); ++i) {
      auto const gap = touch[i] - touch[i - 1U];
      if (gap <= kEps) {
        continue;
      }
      auto const s = touch[i - 1U] + gap / 2.0;
      auto const mid = xy{a.x_ + dx * s, a.y_ + dy * s};
      if (!strictly_inside(mid) &&
          !(on_boundary(mid) && free_beside(a, b, mid))) {
        return false;
      }
    }
    return true;
  }

  // Ring 0 bounds the area; the rest are holes, whether they came from inner
  // rings or from thickened barriers.
  std::vector<std::vector<xy>> rings_;
};

polygon make_polygon(plane const& p,
                     std::vector<std::vector<geo::latlng>> const& rings,
                     std::vector<std::vector<geo::latlng>> const& obstacles) {
  auto poly = polygon{};
  for (auto const& r : rings) {
    if (r.size() < 3U) {
      continue;
    }
    auto& ring = poly.rings_.emplace_back();
    for (auto const& c : r) {
      auto const q = p(c);
      if (ring.empty() || dist(ring.back(), q) > kMergeDistance) {
        ring.push_back(q);
      }
    }
    while (ring.size() > 1U &&
           dist(ring.front(), ring.back()) <= kMergeDistance) {
      ring.pop_back();
    }
    if (ring.size() < 3U) {
      poly.rings_.pop_back();
    }
  }
  for (auto const& o : obstacles) {
    auto line = std::vector<xy>{};
    for (auto const& c : o) {
      auto const q = p(c);
      if (line.empty() || dist(line.back(), q) > kMergeDistance) {
        line.push_back(q);
      }
    }
    for (auto& box : thicken(line, kObstacleHalfWidth)) {
      poly.rings_.push_back(std::move(box));
    }
  }
  return poly;
}

}  // namespace

bool area_geodesics::is_segment_inside(
    std::vector<std::vector<geo::latlng>> const& rings,
    geo::latlng const& a,
    geo::latlng const& b,
    std::vector<std::vector<geo::latlng>> const& obstacles) {
  if (rings.empty() || rings.front().size() < 3U) {
    return false;
  }
  auto const p = plane{rings.front().front()};
  return make_polygon(p, rings, obstacles).is_visible(p(a), p(b));
}

area_geodesics::area_geodesics(
    std::vector<std::vector<geo::latlng>> const& rings,
    std::vector<geo::latlng> const& connectors,
    std::vector<std::vector<geo::latlng>> const& obstacles)
    : n_connectors_{connectors.size()} {
  if (rings.empty() || rings.front().size() < 3U || connectors.empty()) {
    dist_.assign(n_connectors_ * n_connectors_, kUnreachable);
    connector_vertex_.assign(n_connectors_, kNoVertex);
    return;
  }

  auto const p = plane{rings.front().front()};

  auto pts = std::vector<xy>{};
  auto const add_vertex = [&](geo::latlng const& c) -> vertex_idx_t {
    auto const q = p(c);
    for (auto i = std::size_t{0U}; i != pts.size(); ++i) {
      if (dist(pts[i], q) <= kMergeDistance) {
        return static_cast<vertex_idx_t>(i);
      }
    }
    pts.push_back(q);
    vertices_.push_back(c);
    return static_cast<vertex_idx_t>(pts.size() - 1U);
  };

  connector_vertex_.reserve(n_connectors_);
  for (auto const& c : connectors) {
    connector_vertex_.push_back(add_vertex(c));
  }

  auto ring_vertices = std::vector<std::vector<vertex_idx_t>>{};
  auto poly = polygon{};
  for (auto const& r : rings) {
    if (r.size() < 3U) {
      continue;
    }
    auto& rv = ring_vertices.emplace_back();
    for (auto const& c : r) {
      auto const v = add_vertex(c);
      if (rv.empty() || rv.back() != v) {
        rv.push_back(v);
      }
    }
    while (rv.size() > 1U && rv.front() == rv.back()) {
      rv.pop_back();
    }
    if (rv.size() < 3U) {
      ring_vertices.pop_back();
      continue;
    }
    auto& ring = poly.rings_.emplace_back();
    ring.reserve(rv.size());
    for (auto const v : rv) {
      ring.push_back(pts[v]);
    }
  }

  // Barriers become thin holes. Their corners join the vertex set so a path
  // can bend round the free end of a fence, but they get no unconditional
  // edges the way the outline does: an edge along a barrier box that runs
  // outside the area would otherwise offer a way around it.
  // Indices, not a count: add_vertex deduplicates, so an obstacle corner that
  // coincides with an existing point does not append one.
  auto obstacle_vertices = std::vector<vertex_idx_t>{};
  auto obstacle_rings = std::vector<std::vector<xy>>{};
  for (auto const& o : obstacles) {
    auto line = std::vector<xy>{};
    for (auto const& c : o) {
      auto const q = p(c);
      if (line.empty() || dist(line.back(), q) > kMergeDistance) {
        line.push_back(q);
      }
    }
    for (auto& box : thicken(line, kObstacleHalfWidth)) {
      for (auto const& v : box) {
        obstacle_vertices.push_back(add_vertex(p.inverse(v)));
      }
      obstacle_rings.push_back(std::move(box));
    }
  }

  // A wall is walkable along its length - unless a barrier crosses it. A fence
  // anchored to a wall does exactly that, and without this the outline itself
  // becomes the way round the fence.
  auto const crosses_obstacle = [&](xy const u, xy const v) {
    for (auto const& box : obstacle_rings) {
      for (auto i = std::size_t{0U}; i != box.size(); ++i) {
        if (properly_crosses(u, v, box[i], box[(i + 1U) % box.size()])) {
          return true;
        }
      }
    }
    return false;
  };

  // A shortest path inside a polygon bends only where the boundary sticks INTO
  // the free space - at a reflex vertex. Convex vertices can therefore be left
  // out of the graph entirely: any path that appeared to turn at one is either
  // straight through it or beaten by the chord across it. Since a simple
  // polygon always has more convex vertices than reflex ones, this is the
  // difference between O(n^3) over every vertex and O(n^3) over a fraction of
  // them.
  //
  // Reflexness is measured against the free space, not the ring: a hole's
  // corner sticks into the free space even though it is convex for the hole.
  // Normalising the outer ring to counter-clockwise and holes to clockwise
  // puts the free space on the left of every ring, so one rule covers both.
  if (pts.size() > std::numeric_limits<vertex_idx_t>::max()) {
    std::fprintf(stderr, "area_geodesics: %zu points exceeds vertex_idx_t\n",
                 pts.size());
    std::abort();
  }
  auto keep = std::vector<char>(pts.size(), 0);
  for (auto const v : connector_vertex_) {
    if (v != kNoVertex) {
      keep[v] = 1;
    }
  }
  for (auto const [ring_i, rv] : utl::enumerate(ring_vertices)) {
    auto const& ring = poly.rings_[ring_i];
    auto const area = signed_area(ring);
    auto const s = (ring_i == 0U ? (area >= 0.0) : (area <= 0.0)) ? 1.0 : -1.0;
    for (auto i = std::size_t{0U}; i != rv.size(); ++i) {
      auto const& u = ring[(i + rv.size() - 1U) % rv.size()];
      auto const& v = ring[i];
      auto const& w = ring[(i + 1U) % rv.size()];
      auto const cross =
          (v.x_ - u.x_) * (w.y_ - v.y_) - (v.y_ - u.y_) * (w.x_ - v.x_);
      if (s * cross < 0.0) {
        keep[rv[i]] = 1;
      }
    }
  }
  // Barrier corners are hole corners, so they are reflex by construction and
  // are what a path bends at to get round a fence.
  for (auto const v : obstacle_vertices) {
    keep[v] = 1;
  }

  auto const n = pts.size();
  auto adj = std::vector<std::vector<std::pair<vertex_idx_t, float>>>(n);
  auto const add_edge = [&](vertex_idx_t const a, vertex_idx_t const b) {
    auto const d = static_cast<float>(dist(pts[a], pts[b]));
    adj[a].emplace_back(b, d);
    adj[b].emplace_back(a, d);
  };

  // Boundary edges are walkable by definition - adding them up front means the
  // visibility test never has to decide whether a segment lying exactly on a
  // wall is inside, and guarantees the boundary stays connected even where the
  // geometry is degenerate.
  auto is_ring_edge = hash_set<std::uint64_t>{};
  auto const pair_key = [](vertex_idx_t const a, vertex_idx_t const b) {
    return (static_cast<std::uint64_t>(std::min(a, b)) << 32U) |
           static_cast<std::uint64_t>(std::max(a, b));
  };
  for (auto const& rv : ring_vertices) {
    for (auto i = std::size_t{0U}; i != rv.size(); ++i) {
      auto const a = rv[i];
      auto const b = rv[(i + 1U) % rv.size()];
      // Recorded as a ring edge either way, so the visibility test never
      // gets to reconsider it; walkable only with free space beside it.
      if (a != b && keep[a] != 0 && keep[b] != 0 &&
          is_ring_edge.insert(pair_key(a, b)).second &&
          !crosses_obstacle(pts[a], pts[b]) &&
          poly.free_beside(pts[a], pts[b],
                           xy{0.5 * (pts[a].x_ + pts[b].x_),
                              0.5 * (pts[a].y_ + pts[b].y_)})) {
        add_edge(a, b);
      }
    }
  }

  for (auto& box : obstacle_rings) {
    poly.rings_.push_back(std::move(box));
  }

  for (auto i = vertex_idx_t{0U}; i != n; ++i) {
    if (keep[i] == 0) {
      continue;
    }
    for (auto j = static_cast<vertex_idx_t>(i + 1U); j != n; ++j) {
      if (keep[j] == 0 || is_ring_edge.contains(pair_key(i, j))) {
        continue;
      }
      if (poly.is_visible(pts[i], pts[j])) {
        add_edge(i, j);
      }
    }
  }

  // Connectors outside the polygon have no edges at all and stay unreachable
  // rather than silently snapping onto the nearest wall.
  dist_.assign(n_connectors_ * n_connectors_, kUnreachable);
  pred_.assign(n_connectors_ * n, kNoVertex);

  auto d = std::vector<float>(n);
  using queue_entry = std::pair<float, vertex_idx_t>;
  auto q = std::priority_queue<queue_entry, std::vector<queue_entry>,
                               std::greater<>>{};
  for (auto i = std::size_t{0U}; i != n_connectors_; ++i) {
    auto const src = connector_vertex_[i];
    if (src == kNoVertex) {
      continue;
    }

    std::fill(begin(d), end(d), kUnreachable);
    auto* pred = &pred_[i * n];
    d[src] = 0.F;
    q = {};
    q.emplace(0.F, src);
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
      if (dst != kNoVertex) {
        dist_[i * n_connectors_ + j] = d[dst];
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
