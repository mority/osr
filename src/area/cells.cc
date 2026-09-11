#include "osr/area/cells.h"

#include <cmath>
#include <algorithm>
#include <limits>
#include <numeric>
#include <optional>
#include <set>
#include <utility>

#include "geo/box.h"

#include "geo/constants.h"

#include "utl/verify.h"

#include "osr/area/geodesic.h"

namespace osr {

namespace {

constexpr auto kNone = std::numeric_limits<std::size_t>::max();
constexpr auto kInf = std::numeric_limits<double>::infinity();

// The fit alternates between deriving each pair's cell path and refitting the
// costs, until the paths stop changing. As long as the neighbour relation is a
// tree - and the cut tree only ever produces trees - there is exactly one path
// between two cells whatever the costs, so the second round confirms the
// first. The bound only matters for a relation with cycles.
constexpr auto kMaxFitRounds = 4;

// Keeps the normal equations solvable when a cell carries no pair.
constexpr auto kRidge = 1e-3;

// Closer than this, two connectors are not a crossing worth modelling.
constexpr auto kMinPairDistance = 1.0;

std::size_t bit_index(std::size_t const m, std::size_t a, std::size_t b) {
  if (a > b) {
    std::swap(a, b);
  }
  return a * (2U * m - a - 1U) / 2U + (b - a - 1U);
}

std::size_t n_words(std::size_t const m) {
  return m < 2U ? 0U : (m * (m - 1U) / 2U + 63U) / 64U;
}

// Bounds-safe: cells that are never neighbours - separate walkable parts, one
// cell each - leave the bits unallocated.
bool test_bit(std::vector<std::uint64_t> const& bits, std::size_t const i) {
  return i / 64U < bits.size() && ((bits[i / 64U] >> (i % 64U)) & 1U) != 0U;
}

void set_bit(std::vector<std::uint64_t>& bits, std::size_t const i) {
  bits[i / 64U] |= std::uint64_t{1U} << (i % 64U);
}

// Cheapest cell paths from `src` when a path costs the sum of the costs of the
// cells on it. There are at most kMaxCells cells, so a plain O(m^2) scan
// beats a heap.
void cell_paths(std::size_t const m,
                std::vector<std::uint64_t> const& neighbours,
                std::vector<double> const& cost,
                std::size_t const src,
                double* dist,
                std::size_t* pred) {
  auto done = std::vector<char>(m, 0);
  std::fill(dist, dist + m, kInf);
  std::fill(pred, pred + m, kNone);
  dist[src] = cost[src];
  for (auto round = std::size_t{0U}; round != m; ++round) {
    auto u = kNone;
    for (auto v = std::size_t{0U}; v != m; ++v) {
      if (done[v] == 0 && dist[v] != kInf && (u == kNone || dist[v] < dist[u])) {
        u = v;
      }
    }
    if (u == kNone) {
      break;
    }
    done[u] = 1;
    for (auto v = std::size_t{0U}; v != m; ++v) {
      if (v != u && done[v] == 0 &&
          test_bit(neighbours, bit_index(m, u, v)) &&
          dist[u] + cost[v] < dist[v]) {
        dist[v] = dist[u] + cost[v];
        pred[v] = u;
      }
    }
  }
}

// Gaussian elimination with partial pivoting on the small dense normal
// equations. Leaves the solution in `b`.
bool solve(std::vector<double>& a, std::vector<double>& b, std::size_t const n) {
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

// Shared borders between cells. Geometry is done in the (lng, lat) plane,
// the plane in which every cut is a straight line.
struct ll {
  double x_{0.0}, y_{0.0};  // lng, lat
};

// Borders shorter than this are two cells touching at a corner.
constexpr auto kMinSharedBorder = 1.0;  // m

double side_value(area_cut_side const& s, ll const p) {
  return p.x_ * s.lng_scale_ * s.cx_ + p.y_ * s.cy_ - s.threshold_;
}

// Sutherland-Hodgman against the kept side of `s`.
std::vector<ll> clip(std::vector<ll> const& poly, area_cut_side const& s) {
  auto const keep = [&](double const f) { return s.below_ ? f < 0.0 : f >= 0.0; };
  auto out = std::vector<ll>{};
  for (auto i = std::size_t{0U}; i != poly.size(); ++i) {
    auto const p = poly[i];
    auto const q = poly[(i + 1U) % poly.size()];
    auto const fp = side_value(s, p);
    auto const fq = side_value(s, q);
    if (keep(fp)) {
      out.push_back(p);
    }
    if (keep(fp) != keep(fq)) {
      auto const t = fp / (fp - fq);
      out.push_back({p.x_ + t * (q.x_ - p.x_), p.y_ + t * (q.y_ - p.y_)});
    }
  }
  return out;
}

// The line of a cut, as base point + s * direction.
struct cut_line {
  explicit cut_line(area_cut_side const& s)
      : a_{s.lng_scale_ * s.cx_}, b_{s.cy_}, t_{s.threshold_} {
    auto const n2 = a_ * a_ + b_ * b_;
    p0_ = {a_ * t_ / n2, b_ * t_ / n2};
    d_ = {-b_, a_};
  }

  ll at(double const s) const {
    return {p0_.x_ + s * d_.x_, p0_.y_ + s * d_.y_};
  }

  // Where the line runs through the convex polygon `poly`, as a parameter
  // interval.
  std::optional<std::pair<double, double>> through(
      std::vector<ll> const& poly) const {
    auto lo = kInf;
    auto hi = -kInf;
    auto const f = [&](ll const p) { return a_ * p.x_ + b_ * p.y_ - t_; };
    auto const param = [&](ll const p) {
      return ((p.x_ - p0_.x_) * d_.x_ + (p.y_ - p0_.y_) * d_.y_) /
             (d_.x_ * d_.x_ + d_.y_ * d_.y_);
    };
    for (auto i = std::size_t{0U}; i != poly.size(); ++i) {
      auto const p = poly[i];
      auto const q = poly[(i + 1U) % poly.size()];
      auto const fp = f(p);
      auto const fq = f(q);
      if ((fp < 0.0) == (fq < 0.0) && fp != 0.0 && fq != 0.0) {
        continue;
      }
      auto const u = fp == fq ? 0.0 : fp / (fp - fq);
      auto const s = param({p.x_ + u * (q.x_ - p.x_), p.y_ + u * (q.y_ - p.y_)});
      lo = std::min(lo, s);
      hi = std::max(hi, s);
    }
    return lo <= hi ? std::optional{std::pair{lo, hi}} : std::nullopt;
  }

  double a_, b_, t_;
  ll p0_, d_;
};

// Even-odd over all rings: inside the outline and outside every hole.
bool inside(std::vector<std::vector<ll>> const& rings, ll const p) {
  auto in = false;
  for (auto const& r : rings) {
    for (auto i = std::size_t{0U}; i != r.size(); ++i) {
      auto const& a = r[i];
      auto const& b = r[(i + 1U) % r.size()];
      if ((a.y_ > p.y_) != (b.y_ > p.y_) &&
          p.x_ < a.x_ + (p.y_ - a.y_) / (b.y_ - a.y_) * (b.x_ - a.x_)) {
        in = !in;
      }
    }
  }
  return in;
}

// How much of the segment r0 -> r1 lies inside the area, in meters.
double inside_length(std::vector<std::vector<ll>> const& rings,
                     ll const r0,
                     ll const r1) {
  auto ts = std::vector<double>{0.0, 1.0};
  auto const dx = r1.x_ - r0.x_;
  auto const dy = r1.y_ - r0.y_;
  for (auto const& r : rings) {
    for (auto i = std::size_t{0U}; i != r.size(); ++i) {
      auto const& a = r[i];
      auto const& b = r[(i + 1U) % r.size()];
      auto const ex = b.x_ - a.x_;
      auto const ey = b.y_ - a.y_;
      auto const den = dx * ey - dy * ex;
      if (den == 0.0) {
        continue;
      }
      auto const t = ((a.x_ - r0.x_) * ey - (a.y_ - r0.y_) * ex) / den;
      auto const u = ((a.x_ - r0.x_) * dy - (a.y_ - r0.y_) * dx) / den;
      if (t > 0.0 && t < 1.0 && u >= 0.0 && u <= 1.0) {
        ts.push_back(t);
      }
    }
  }
  std::ranges::sort(ts);
  auto const lat = 0.5 * (r0.y_ + r1.y_);
  auto const len =
      std::hypot(dx * std::cos(lat * geo::kPI / 180.0), dy) *
      geo::kApproxDistanceLatDegrees;
  auto total = 0.0;
  for (auto i = std::size_t{1U}; i != ts.size(); ++i) {
    auto const m = 0.5 * (ts[i - 1U] + ts[i]);
    if (inside(rings, {r0.x_ + m * dx, r0.y_ + m * dy})) {
      total += (ts[i] - ts[i - 1U]) * len;
    }
  }
  return total;
}

std::vector<std::vector<ll>> to_ll(
    std::vector<std::vector<geo::latlng>> const& rings) {
  auto out = std::vector<std::vector<ll>>{};
  for (auto const& r : rings) {
    auto& o = out.emplace_back();
    for (auto const& p : r) {
      o.push_back({p.lng(), p.lat()});
    }
    if (o.size() > 1U && o.front().x_ == o.back().x_ &&
        o.front().y_ == o.back().y_) {
      o.pop_back();
    }
    if (o.size() < 3U) {
      out.pop_back();
    }
  }
  return out;
}

// The shared-border rule on rings already in the (lng, lat) plane (see
// cells_sharing_a_border). Two cells are separated by the first cut on their
// paths where they take opposite sides. Each cell's region with that one cut
// dropped is convex (before intersecting with the area) and crosses the
// cut's line; where the two crossings overlap, and the overlap lies inside
// the area, the cells border each other.
std::vector<std::pair<std::size_t, std::size_t>> shared_border_pairs(
    std::vector<std::vector<ll>> const& rings,
    std::vector<std::vector<area_cut_side>> const& regions) {
  auto out = std::vector<std::pair<std::size_t, std::size_t>>{};
  if (rings.empty()) {
    return out;
  }

  auto box = geo::box{};
  for (auto const& p : rings.front()) {
    box.extend(geo::latlng{p.y_, p.x_});
  }
  auto const margin = 1e-4;
  auto const bounds =
      std::vector<ll>{{box.min_.lng() - margin, box.min_.lat() - margin},
                      {box.max_.lng() + margin, box.min_.lat() - margin},
                      {box.max_.lng() + margin, box.max_.lat() + margin},
                      {box.min_.lng() - margin, box.max_.lat() + margin}};

  auto const m = regions.size();
  for (auto a = std::size_t{0U}; a != m; ++a) {
    for (auto b = a + 1U; b != m; ++b) {
      auto const& ra = regions[a];
      auto const& rb = regions[b];
      auto k = std::size_t{0U};
      while (k < ra.size() && k < rb.size() && ra[k].below_ == rb[k].below_) {
        ++k;
      }
      if (k == ra.size() || k == rb.size()) {
        continue;  // not separated by a cut: cannot happen for two cells
      }

      auto const convex = [&](std::vector<area_cut_side> const& r) {
        auto poly = bounds;
        for (auto j = std::size_t{0U}; j != r.size() && poly.size() >= 3U;
             ++j) {
          if (j != k) {
            poly = clip(poly, r[j]);
          }
        }
        return poly;
      };
      auto const line = cut_line{ra[k]};
      auto const ia = line.through(convex(ra));
      auto const ib = line.through(convex(rb));
      if (!ia.has_value() || !ib.has_value()) {
        continue;
      }
      auto const lo = std::max(ia->first, ib->first);
      auto const hi = std::min(ia->second, ib->second);
      if (hi > lo &&
          inside_length(rings, line.at(lo), line.at(hi)) >= kMinSharedBorder) {
        out.emplace_back(a, b);
      }
    }
  }
  return out;
}

// A node of the cut tree. The tree is built once, to full depth; refinement
// only decides which of its cuts are used.
struct cut_node {
  bool is_leaf() const { return lhs_ == kNone; }

  std::vector<std::size_t> connectors_;

  // Midpoint between the two halves' centroids. Only used while building, to
  // decide which cell on each side of this cut the other one neighbours.
  geo::latlng mid_{};

  // The cut as a line, only kept for drawing: lhs_ has the connectors below
  // threshold_, rhs_ those at or above it (see area_cut_side).
  double lng_scale_{0.0}, cx_{0.0}, cy_{0.0}, threshold_{0.0};

  std::size_t lhs_{kNone}, rhs_{kNone};
  bool expanded_{false};
};

// The cells the tree currently resolves to.
struct frontier {
  // Cut-tree node each cell is.
  std::vector<std::size_t> owner_;

  // Positions used to pick neighbours: the cell's connectors, plus the
  // midpoints of the cuts already assigned to it.
  std::vector<std::vector<geo::latlng>> points_;

  std::vector<std::pair<std::size_t, std::size_t>> neighbours_;
};

struct fit {
  std::vector<double> cost_;
  std::vector<int> violations_;
};

struct builder {
  float distance(std::size_t const i, std::size_t const j) const {
    return dist_[i * k_ + j];
  }

  // Spread of the pairwise distances within a candidate cell, over the pairs
  // the crossing matters for.
  double spread(std::vector<std::size_t> const& idx) const {
    auto lo = std::numeric_limits<double>::max();
    auto hi = 0.0;
    auto any = false;
    for (auto i = std::size_t{0U}; i != idx.size(); ++i) {
      for (auto j = i + 1U; j != idx.size(); ++j) {
        if (!relevant_[idx[i] * k_ + idx[j]]) {
          continue;
        }
        auto const d = distance(idx[i], idx[j]);
        if (d == area_geodesics::kUnreachable) {
          continue;
        }
        lo = std::min<double>(lo, d);
        hi = std::max<double>(hi, d);
        any = true;
      }
    }
    return any ? hi - lo : 0.0;
  }

  geo::latlng centroid(std::vector<std::size_t> const& idx) const {
    auto lat = 0.0;
    auto lng = 0.0;
    for (auto const i : idx) {
      lat += pos_[i].lat();
      lng += pos_[i].lng();
    }
    return {lat / static_cast<double>(idx.size()),
            lng / static_cast<double>(idx.size())};
  }

  // Cuts geometrically, so cells stay contiguous and a route passes through
  // them in spatial order, but chooses where by the pairwise distances: of all
  // straight cuts tried, the one whose worse half has the smaller spread.
  void split(std::size_t const n, std::size_t const depth) {
    // Copied: nodes_ grows below.
    auto const conns = nodes_[n].connectors_;
    if (depth == 0U || conns.size() < 4U) {
      return;
    }

    auto const lat0 = pos_[conns.front()].lat();
    auto const lng_scale = std::cos(lat0 * geo::kPI / 180.0);
    auto best_score = std::numeric_limits<double>::max();
    auto best_l = std::vector<std::size_t>{};
    auto best_r = std::vector<std::size_t>{};
    auto best_cx = 0.0;
    auto best_cy = 0.0;
    auto best_threshold = 0.0;

    // Every split position along each direction where that is affordable -
    // the scan is O(k^3) per direction - and five fixed fractions otherwise.
    auto const exhaustive = conns.size() <= 40U;
    auto const n_directions = exhaustive ? 16U : 8U;
    for (auto d = 0U; d != n_directions; ++d) {
      auto const theta =
          geo::kPI * static_cast<double>(d) / static_cast<double>(n_directions);
      auto const cx = std::cos(theta);
      auto const cy = std::sin(theta);
      auto sorted = conns;
      std::ranges::sort(sorted, [&](std::size_t const a, std::size_t const b) {
        auto const ka = pos_[a].lng() * lng_scale * cx + pos_[a].lat() * cy;
        auto const kb = pos_[b].lng() * lng_scale * cx + pos_[b].lat() * cy;
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
        auto const score = std::max(spread(l), spread(r));
        if (score < best_score) {
          auto const key = [&](std::size_t const a) {
            return pos_[a].lng() * lng_scale * cx + pos_[a].lat() * cy;
          };
          best_score = score;
          best_l = std::move(l);
          best_r = std::move(r);
          best_cx = cx;
          best_cy = cy;
          best_threshold = 0.5 * (key(sorted[m - 1U]) + key(sorted[m]));
        }
      }
    }

    if (best_l.empty() || best_r.empty()) {
      return;
    }

    auto const a = centroid(best_l);
    auto const b = centroid(best_r);
    auto const lhs = nodes_.size();
    nodes_.push_back(cut_node{.connectors_ = std::move(best_l)});
    auto const rhs = nodes_.size();
    nodes_.push_back(cut_node{.connectors_ = std::move(best_r)});
    nodes_[n].lhs_ = lhs;
    nodes_[n].rhs_ = rhs;
    nodes_[n].mid_ =
        geo::latlng{0.5 * (a.lat() + b.lat()), 0.5 * (a.lng() + b.lng())};
    nodes_[n].lng_scale_ = lng_scale;
    nodes_[n].cx_ = best_cx;
    nodes_[n].cy_ = best_cy;
    nodes_[n].threshold_ = best_threshold;

    split(lhs, depth - 1U);
    split(rhs, depth - 1U);
  }

  // Collects the cells below `n` and the neighbour pairs between them.
  //
  // A cut makes exactly one pair of cells neighbours: the cell on each side
  // that its midpoint lies closest to. Making every cell on one side a
  // neighbour of every cell on the other would put any two cells next to each
  // other, and no route would ever pass through more than two. One pair per
  // cut is also why the relation is always a tree.
  std::pair<std::size_t, std::size_t> collect(std::size_t const n,
                                              frontier& f) const {
    auto const& node = nodes_[n];
    if (node.is_leaf() || !node.expanded_) {
      f.owner_.push_back(n);
      auto& pts = f.points_.emplace_back();
      for (auto const c : node.connectors_) {
        pts.push_back(pos_[c]);
      }
      return {f.owner_.size() - 1U, f.owner_.size()};
    }

    auto const [l0, l1] = collect(node.lhs_, f);
    auto const [r0, r1] = collect(node.rhs_, f);
    auto const nearest = [&](std::size_t const from, std::size_t const to) {
      auto best = from;
      auto best_d = std::numeric_limits<double>::max();
      for (auto i = from; i != to; ++i) {
        auto const& pts = f.points_[i];
        if (pts.empty()) {
          continue;
        }
        auto lat = 0.0;
        auto lng = 0.0;
        for (auto const& p : pts) {
          lat += p.lat();
          lng += p.lng();
        }
        auto const nn = static_cast<double>(pts.size());
        auto const d = geo::distance(node.mid_, {lat / nn, lng / nn});
        if (d < best_d) {
          best_d = d;
          best = i;
        }
      }
      return best;
    };
    auto const a = nearest(l0, l1);
    auto const b = nearest(r0, r1);
    f.points_[a].push_back(node.mid_);
    f.points_[b].push_back(node.mid_);
    f.neighbours_.emplace_back(a, b);
    return {l0, r1};
  }

  // The cut sides bounding each cell below `n`, in the same order collect()
  // lists the cells.
  void collect_regions(std::size_t const n,
                       std::vector<area_cut_side>& path,
                       std::vector<std::vector<area_cut_side>>& out) const {
    auto const& node = nodes_[n];
    if (node.is_leaf() || !node.expanded_) {
      out.push_back(path);
      return;
    }
    auto const side = [&](bool const below) {
      return area_cut_side{.lng_scale_ = node.lng_scale_,
                           .cx_ = node.cx_,
                           .cy_ = node.cy_,
                           .threshold_ = node.threshold_,
                           .below_ = below};
    };
    path.push_back(side(true));
    collect_regions(node.lhs_, path, out);
    path.back() = side(false);
    collect_regions(node.rhs_, path, out);
    path.pop_back();
  }

  // Makes every two cells of `f` neighbours whose regions share a border of
  // at least kMinSharedBorder inside the area. The pairs collect() found are
  // kept: they alone connect all cells, whatever the geometry.
  //
  // Two cells are separated by the first cut on their paths where they take
  // opposite sides. Each cell's region with that one cut dropped is convex
  // (before intersecting with the area), and crosses the cut's line; where
  // the two crossings overlap, and the overlap lies inside the area, the
  // cells border each other.
  void add_shared_borders(frontier& f,
                          std::vector<std::vector<ll>> const& rings) const {
    if (rings.empty()) {
      return;
    }
    auto regions = std::vector<std::vector<area_cut_side>>{};
    auto path = std::vector<area_cut_side>{};
    collect_regions(0U, path, regions);

    auto pairs = std::set<std::pair<std::size_t, std::size_t>>{};
    for (auto const& [a, b] : f.neighbours_) {
      pairs.emplace(std::min(a, b), std::max(a, b));
    }
    for (auto const& p : shared_border_pairs(rings, regions)) {
      pairs.insert(p);
    }
    f.neighbours_.assign(begin(pairs), end(pairs));
  }

  // Fits the cell costs of `f` against the geodesics by least squares, and
  // counts per cell the routes through it that miss by more than `threshold`.
  //
  // The costs are free parameters: the cells only fix which costs a pair's
  // route adds up, and the costs are whatever reproduces the true distances
  // best. Fitting each cell alone from its own pairs would be a guess.
  fit fit_costs(frontier const& f, double const threshold) const {
    auto const m = f.owner_.size();

    auto cell_of = std::vector<std::size_t>(k_, kNone);
    for (auto ci = std::size_t{0U}; ci != m; ++ci) {
      for (auto const c : nodes_[f.owner_[ci]].connectors_) {
        cell_of[c] = ci;
      }
    }

    auto neighbours = std::vector<std::uint64_t>(n_words(m), 0U);
    for (auto const& [a, b] : f.neighbours_) {
      set_bit(neighbours, bit_index(m, a, b));
    }

    // Seed: the midrange of each cell's own pairs, the flat cost that
    // minimises the worst error inside it. Only kept if there is too little
    // to fit. A cell with a single connector has no pair and starts at zero.
    auto result = fit{.cost_ = std::vector<double>(m, 0.0),
                      .violations_ = std::vector<int>(m, 0)};
    for (auto ci = std::size_t{0U}; ci != m; ++ci) {
      auto const& conns = nodes_[f.owner_[ci]].connectors_;
      auto lo = kInf;
      auto hi = 0.0;
      for (auto i = std::size_t{0U}; i != conns.size(); ++i) {
        for (auto j = i + 1U; j != conns.size(); ++j) {
          auto const d = distance(conns[i], conns[j]);
          lo = std::min<double>(lo, d);
          hi = std::max<double>(hi, d);
        }
      }
      result.cost_[ci] = lo == kInf ? 0.0 : 0.5 * (lo + hi);
    }

    auto dist = std::vector<double>(m * m);
    auto pred = std::vector<std::size_t>(m * m);
    auto const derive_paths = [&]() {
      for (auto src = std::size_t{0U}; src != m; ++src) {
        cell_paths(m, neighbours, result.cost_, src, &dist[src * m],
                   &pred[src * m]);
      }
    };
    auto const for_each_on_path = [&](std::size_t const a, std::size_t const b,
                                      auto&& fn) {
      for (auto v = b; v != kNone; v = (v == a ? kNone : pred[a * m + v])) {
        fn(v);
      }
    };
    // Every pair the model is measured on, in a fixed order.
    auto const for_each_pair = [&](auto&& fn) {
      for (auto src = std::size_t{0U}; src != k_; ++src) {
        for (auto dst = src + 1U; dst != k_; ++dst) {
          auto const truth = distance(src, dst);
          auto const a = cell_of[src];
          auto const b = cell_of[dst];
          if (truth == area_geodesics::kUnreachable ||
              truth < kMinPairDistance || dist[a * m + b] == kInf) {
            continue;
          }
          fn(a, b, static_cast<double>(truth));
        }
      }
    };

    auto prev_pred = std::vector<std::size_t>{};
    auto used = std::vector<std::size_t>{};
    for (auto round = 0; round != kMaxFitRounds; ++round) {
      derive_paths();
      if (pred == prev_pred) {
        break;  // same paths, same rows, same solution
      }
      prev_pred = pred;

      auto ata = std::vector<double>(m * m, 0.0);
      auto atb = std::vector<double>(m, 0.0);
      auto n_rows = std::size_t{0U};
      for_each_pair([&](std::size_t const a, std::size_t const b,
                        double const truth) {
        ++n_rows;
        used.clear();
        for_each_on_path(a, b, [&](std::size_t const c) { used.push_back(c); });
        for (auto const i : used) {
          atb[i] += truth;
          for (auto const j : used) {
            ata[i * m + j] += 1.0;
          }
        }
      });
      if (n_rows < m) {
        break;
      }
      for (auto i = std::size_t{0U}; i != m; ++i) {
        ata[i * m + i] += kRidge;
      }
      if (!solve(ata, atb, m)) {
        break;
      }
      for (auto i = std::size_t{0U}; i != m; ++i) {
        result.cost_[i] = std::max(0.0, atb[i]);
      }
    }

    derive_paths();
    for_each_pair(
        [&](std::size_t const a, std::size_t const b, double const truth) {
          if (std::abs(dist[a * m + b] - truth) > threshold) {
            for_each_on_path(a, b,
                             [&](std::size_t const c) { ++result.violations_[c]; });
          }
        });
    return result;
  }

  std::vector<geo::latlng> const& pos_;
  std::vector<float> const& dist_;
  std::vector<bool> const& relevant_;
  std::size_t k_;
  std::vector<cut_node> nodes_;
};

// Smallest enclosing circle, in a local metric plane.
struct xy {
  double x_{0.0}, y_{0.0};
};

struct circle {
  bool covers(xy const p) const {
    return std::hypot(p.x_ - c_.x_, p.y_ - c_.y_) <= r_ + 1e-6;
  }

  xy c_;
  double r_{0.0};
};

circle circle_through(xy const a, xy const b) {
  return {.c_ = {0.5 * (a.x_ + b.x_), 0.5 * (a.y_ + b.y_)},
          .r_ = 0.5 * std::hypot(b.x_ - a.x_, b.y_ - a.y_)};
}

circle circle_through(xy const a, xy const b, xy const c) {
  auto const bx = b.x_ - a.x_;
  auto const by = b.y_ - a.y_;
  auto const cx = c.x_ - a.x_;
  auto const cy = c.y_ - a.y_;
  auto const d = 2.0 * (bx * cy - by * cx);
  if (std::abs(d) < 1e-12) {
    // Collinear: the two points furthest apart span the other one.
    auto best = circle_through(a, b);
    for (auto const& x : {circle_through(a, c), circle_through(b, c)}) {
      if (x.r_ > best.r_) {
        best = x;
      }
    }
    return best;
  }
  auto const b2 = bx * bx + by * by;
  auto const c2 = cx * cx + cy * cy;
  auto const ux = (cy * b2 - by * c2) / d;
  auto const uy = (bx * c2 - cx * b2) / d;
  return {.c_ = {a.x_ + ux, a.y_ + uy}, .r_ = std::hypot(ux, uy)};
}

// Welzl's algorithm in its iterative form. Without shuffling the worst case
// is cubic, which a cell's few dozen connectors can afford, and the result
// stays deterministic.
circle smallest_enclosing_circle(std::vector<xy> const& p) {
  auto c = circle{.c_ = p.front(), .r_ = 0.0};
  for (auto i = std::size_t{1U}; i < p.size(); ++i) {
    if (c.covers(p[i])) {
      continue;
    }
    c = circle{.c_ = p[i], .r_ = 0.0};
    for (auto j = std::size_t{0U}; j != i; ++j) {
      if (c.covers(p[j])) {
        continue;
      }
      c = circle_through(p[i], p[j]);
      for (auto k = std::size_t{0U}; k != j; ++k) {
        if (!c.covers(p[k])) {
          c = circle_through(p[i], p[j], p[k]);
        }
      }
    }
  }
  return c;
}

}  // namespace

std::vector<geo::latlng> hub_positions(
    area_cells const& cells, std::vector<geo::latlng> const& connectors) {
  utl::verify(connectors.size() == cells.n_connectors(),
              "hub positions: {} connectors for {} cell assignments",
              connectors.size(), cells.n_connectors());

  auto members = std::vector<std::vector<geo::latlng>>(cells.n_cells());
  for (auto i = std::size_t{0U}; i != connectors.size(); ++i) {
    if (cells.connector_cell_[i] != area_cells::kNoCell) {
      members[cells.connector_cell_[i]].push_back(connectors[i]);
    }
  }

  auto hubs = std::vector<geo::latlng>{};
  hubs.reserve(members.size());
  for (auto const& m : members) {
    utl::verify(!m.empty(), "hub positions: cell without connectors");
    auto const lat0 = m.front().lat();
    auto const lng0 = m.front().lng();
    auto const lng_scale =
        std::cos(lat0 * geo::kPI / 180.0) * geo::kApproxDistanceLatDegrees;
    auto pts = std::vector<xy>{};
    pts.reserve(m.size());
    for (auto const& p : m) {
      pts.push_back({(p.lng() - lng0) * lng_scale,
                     (p.lat() - lat0) * geo::kApproxDistanceLatDegrees});
    }
    auto const c = smallest_enclosing_circle(pts).c_;
    hubs.emplace_back(lat0 + c.y_ / geo::kApproxDistanceLatDegrees,
                      lng0 + c.x_ / lng_scale);
  }
  return hubs;
}

std::vector<std::pair<std::size_t, std::size_t>> cells_sharing_a_border(
    std::vector<std::vector<geo::latlng>> const& rings,
    std::vector<std::vector<area_cut_side>> const& cell_regions) {
  return shared_border_pairs(to_ll(rings), cell_regions);
}

std::vector<std::vector<std::vector<geo::latlng>>> region_rings(
    std::vector<std::vector<geo::latlng>> const& rings,
    std::vector<std::vector<area_cut_side>> const& cell_regions) {
  auto const plane = to_ll(rings);
  auto out = std::vector<std::vector<std::vector<geo::latlng>>>{};
  for (auto const& sides : cell_regions) {
    auto& cell = out.emplace_back();
    for (auto const& r : plane) {
      auto clipped = r;
      for (auto const& side : sides) {
        clipped = clip(clipped, side);
      }
      if (clipped.size() >= 3U) {
        auto& ring = cell.emplace_back();
        for (auto const& p : clipped) {
          ring.emplace_back(p.y_, p.x_);
        }
      } else if (&r == &plane.front()) {
        break;  // nothing of the outline is left
      }
    }
  }
  return out;
}

bool area_cells::is_neighbour(cell_idx_t const a, cell_idx_t const b) const {
  return a != b && test_bit(neighbours_, bit_index(n_cells(), a, b));
}

void area_cells::set_neighbour(cell_idx_t const a, cell_idx_t const b) {
  utl::verify(a != b && a < n_cells() && b < n_cells(),
              "area cells: no neighbour pair {}-{} among {} cells", a, b,
              n_cells());
  neighbours_.resize(std::max(neighbours_.size(), n_words(n_cells())), 0U);
  set_bit(neighbours_, bit_index(n_cells(), a, b));
}

std::vector<float> area_cells::cell_distances() const {
  auto const m = n_cells();
  auto const cost = std::vector<double>(begin(cost_), end(cost_));
  auto dist = std::vector<double>(m * m);
  auto pred = std::vector<std::size_t>(m * m);
  for (auto src = std::size_t{0U}; src != m; ++src) {
    cell_paths(m, neighbours_, cost, src, &dist[src * m], &pred[src * m]);
  }
  return {begin(dist), end(dist)};
}

namespace {

// One walkable part: every pair of `connectors` reachable.
area_cells build_part(
    std::vector<geo::latlng> const& connectors,
    std::vector<float> const& dist,
    std::vector<bool> const& relevant,
    area_cells_params const& params,
    std::vector<std::vector<area_cut_side>>* cell_regions) {
  auto const k = connectors.size();
  auto b = builder{.pos_ = connectors,
                   .dist_ = dist,
                   .relevant_ = relevant,
                   .k_ = k,
                   .nodes_ = {}};
  auto& root = b.nodes_.emplace_back();
  root.connectors_.resize(k);
  for (auto i = std::size_t{0U}; i != k; ++i) {
    root.connectors_[i] = i;
  }
  b.split(0U, params.max_depth_);
  auto const rings = params.rings_ != nullptr ? to_ll(*params.rings_)
                                              : std::vector<std::vector<ll>>{};

  // Start from a single cell and split every cell that carries a route the
  // fitted model gets wrong by more than the threshold, until none does or
  // those cells cannot be split further.
  auto f = frontier{};
  auto result = fit{};
  for (auto round = std::size_t{0U}; round <= params.max_depth_; ++round) {
    f = frontier{};
    b.collect(0U, f);
    b.add_shared_borders(f, rings);
    result = b.fit_costs(f, params.threshold_);

    // Out of rounds, and this does happen: the fit is global, so a cell that
    // was within the threshold can start missing after a later refit. Keep
    // the frontier just evaluated, and leave the tree matching it.
    if (round == params.max_depth_) {
      break;
    }

    auto refined = false;
    for (auto ci = std::size_t{0U}; ci != f.owner_.size(); ++ci) {
      auto& node = b.nodes_[f.owner_[ci]];
      if (result.violations_[ci] > 0 && !node.is_leaf()) {
        node.expanded_ = true;
        refined = true;
      }
    }
    if (!refined) {
      break;
    }
  }

  auto const m = f.owner_.size();
  auto cells = area_cells{};
  cells.connector_cell_.resize(k);
  for (auto ci = std::size_t{0U}; ci != m; ++ci) {
    for (auto const c : b.nodes_[f.owner_[ci]].connectors_) {
      cells.connector_cell_[c] = static_cast<area_cells::cell_idx_t>(ci);
    }
  }
  cells.cost_.assign(begin(result.cost_), end(result.cost_));
  cells.neighbours_.assign(n_words(m), 0U);
  for (auto const& [x, y] : f.neighbours_) {
    set_bit(cells.neighbours_, bit_index(m, x, y));
  }
  if (cell_regions != nullptr) {
    auto path = std::vector<area_cut_side>{};
    auto const before = cell_regions->size();
    b.collect_regions(0U, path, *cell_regions);
    utl::verify(cell_regions->size() - before == m,
                "area cells: {} regions for {} cells",
                cell_regions->size() - before, m);
  }
  return cells;
}

}  // namespace

std::optional<area_cells> build_area_cells(
    std::vector<geo::latlng> const& connectors,
    std::vector<float> const& dist,
    std::vector<bool> const& relevant,
    area_cells_params const& params,
    std::vector<std::vector<area_cut_side>>* cell_regions) {
  auto const k = connectors.size();
  if (cell_regions != nullptr) {
    cell_regions->clear();
  }
  utl::verify(dist.size() == k * k, "area cells: {} distances for {} connectors",
              dist.size(), k);
  utl::verify(relevant.size() == k * k,
              "area cells: {} relevance flags for {} connectors",
              relevant.size(), k);
  utl::verify(params.max_depth_ <= area_cells::kMaxDepth,
              "area cells: depth {} exceeds {}", params.max_depth_,
              area_cells::kMaxDepth);

  // Walkable parts: groups of connectors that reach each other. Ordered by
  // their first connector, and each keeps its connectors' order, so an area
  // that is one part is built from exactly the input it was given.
  auto root = std::vector<std::size_t>(k);
  std::iota(begin(root), end(root), std::size_t{0U});
  auto const find = [&](std::size_t x) {
    while (root[x] != x) {
      x = root[x] = root[root[x]];
    }
    return x;
  };
  for (auto i = std::size_t{0U}; i != k; ++i) {
    for (auto j = i + 1U; j != k; ++j) {
      if (dist[i * k + j] != area_geodesics::kUnreachable) {
        root[find(i)] = find(j);
      }
    }
  }
  auto part_of_root = std::vector<std::size_t>(k, kNone);
  auto parts = std::vector<std::vector<std::size_t>>{};
  for (auto i = std::size_t{0U}; i != k; ++i) {
    auto& p = part_of_root[find(i)];
    if (p == kNone) {
      p = parts.size();
      parts.emplace_back();
    }
    parts[p].push_back(i);
  }

  auto cells = area_cells{};
  cells.connector_cell_.assign(k, area_cells::kNoCell);
  auto neighbours = std::vector<std::pair<std::size_t, std::size_t>>{};
  for (auto const& part : parts) {
    auto const n = part.size();
    if (n < 2U) {
      continue;  // a stray: nothing to cross to
    }

    auto pos = std::vector<geo::latlng>(n);
    auto part_dist = std::vector<float>(n * n);
    auto part_relevant = std::vector<bool>(n * n);
    for (auto a = std::size_t{0U}; a != n; ++a) {
      pos[a] = connectors[part[a]];
      for (auto b = std::size_t{0U}; b != n; ++b) {
        part_dist[a * n + b] = dist[part[a] * k + part[b]];
        part_relevant[a * n + b] = relevant[part[a] * k + part[b]];
        if (a != b && part_dist[a * n + b] == area_geodesics::kUnreachable) {
          return std::nullopt;
        }
      }
    }

    auto const sub =
        build_part(pos, part_dist, part_relevant, params, cell_regions);
    auto const offset = cells.cost_.size();
    utl::verify(offset + sub.n_cells() < area_cells::kNoCell,
                "area cells: more than {} cells", area_cells::kNoCell - 1U);
    for (auto a = std::size_t{0U}; a != n; ++a) {
      cells.connector_cell_[part[a]] = static_cast<area_cells::cell_idx_t>(
          offset + sub.connector_cell_[a]);
    }
    cells.cost_.insert(end(cells.cost_), begin(sub.cost_), end(sub.cost_));
    for (auto x = std::size_t{0U}; x != sub.n_cells(); ++x) {
      for (auto y = x + 1U; y != sub.n_cells(); ++y) {
        if (sub.is_neighbour(static_cast<area_cells::cell_idx_t>(x),
                             static_cast<area_cells::cell_idx_t>(y))) {
          neighbours.emplace_back(offset + x, offset + y);
        }
      }
    }
  }

  if (cells.cost_.empty()) {
    return std::nullopt;
  }
  cells.neighbours_.assign(n_words(cells.n_cells()), 0U);
  for (auto const& [x, y] : neighbours) {
    set_bit(cells.neighbours_, bit_index(cells.n_cells(), x, y));
  }
  return cells;
}

}  // namespace osr
