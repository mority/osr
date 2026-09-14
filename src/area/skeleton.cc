#include "osr/area/skeleton.h"

#include <cmath>
#include <cstdint>
#include <algorithm>
#include <limits>
#include <map>
#include <queue>
#include <utility>

#include "boost/polygon/point_data.hpp"
#include "boost/polygon/segment_data.hpp"
#include "boost/polygon/voronoi.hpp"

#include "geo/constants.h"

#include "osr/area/geodesic.h"

namespace osr {

namespace bp = boost::polygon;

namespace {

// The Voronoi diagram is built on integer coordinates: millimetres, so an
// area spanning kilometres stays far inside a 32-bit integer.
constexpr auto kUnit = 1e3;

// How many skeleton vertices a point tries before giving up on reaching it.
constexpr auto kNearest = std::size_t{8U};

constexpr auto kNoVertex = std::numeric_limits<std::size_t>::max();
constexpr auto kInf = std::numeric_limits<double>::infinity();

}  // namespace

struct area_skeleton::impl {
  impl(std::vector<std::vector<geo::latlng>> const& rings,
       std::vector<std::vector<geo::latlng>> const& obstacles)
      : free_{rings, obstacles} {
    if (!free_.valid()) {
      return;
    }
    auto const free_rings = free_.rings();
    if (free_rings.empty() || free_rings.front().size() < 3U) {
      return;
    }

    lat0_ = free_rings.front().front().lat();
    lng0_ = free_rings.front().front().lng();
    scale_ = std::cos(lat0_ * geo::kPI / 180.0) * geo::kApproxDistanceLatDegrees;

    auto segments = std::vector<bp::segment_data<std::int32_t>>{};
    for (auto const& r : free_rings) {
      auto& ring = rings_.emplace_back();
      for (auto const& p : r) {
        ring.push_back(to_meters(p));
      }
      for (auto i = std::size_t{0U}; i != ring.size(); ++i) {
        auto const& a = ring[i];
        auto const& b = ring[(i + 1U) % ring.size()];
        if (unit(a.first) != unit(b.first) || unit(a.second) != unit(b.second)) {
          segments.push_back({{unit(a.first), unit(a.second)},
                              {unit(b.first), unit(b.second)}});
        }
      }
    }
    if (segments.empty()) {
      return;
    }

    auto vd = bp::voronoi_diagram<double>{};
    auto const points = std::vector<bp::point_data<std::int32_t>>{};
    bp::construct_voronoi(points.begin(), points.end(), segments.begin(),
                          segments.end(), &vd);

    auto index = std::map<std::pair<std::int32_t, std::int32_t>, std::size_t>{};
    auto const node = [&](double const x, double const y) {
      auto const [it, fresh] =
          index.try_emplace(std::pair{unit(x), unit(y)}, pos_.size());
      if (fresh) {
        pos_.push_back(std::pair{x, y});
        adj_.emplace_back();
      }
      return it->second;
    };

    for (auto const& e : vd.edges()) {
      // Primary edges only: a secondary edge runs from a wall to one of its
      // own corners and is no part of the medial axis. Each edge is stored
      // both ways round, so only one of the pair is taken.
      if (!e.is_primary() || !e.is_finite() || &e > e.twin()) {
        continue;
      }
      auto const x0 = e.vertex0()->x() / kUnit;
      auto const y0 = e.vertex0()->y() / kUnit;
      auto const x1 = e.vertex1()->x() / kUnit;
      auto const y1 = e.vertex1()->y() / kUnit;
      if (!inside(x0, y0) || !inside(x1, y1) ||
          !inside(0.5 * (x0 + x1), 0.5 * (y0 + y1))) {
        continue;
      }
      auto const a = node(x0, y0);
      auto const b = node(x1, y1);
      if (a == b) {
        continue;
      }
      auto const d = std::hypot(x1 - x0, y1 - y0);
      adj_[a].emplace_back(b, d);
      adj_[b].emplace_back(a, d);
      ++n_edges_;
    }
  }

  std::int32_t unit(double const v) const {
    return static_cast<std::int32_t>(std::llround(v * kUnit));
  }

  std::pair<double, double> to_meters(geo::latlng const& p) const {
    return {(p.lng() - lng0_) * scale_,
            (p.lat() - lat0_) * geo::kApproxDistanceLatDegrees};
  }

  geo::latlng to_latlng(std::pair<double, double> const& p) const {
    return {lat0_ + p.second / geo::kApproxDistanceLatDegrees,
            lng0_ + p.first / scale_};
  }

  // Even-odd over the free space's own rings.
  bool inside(double const x, double const y) const {
    auto in = false;
    for (auto const& r : rings_) {
      for (auto i = std::size_t{0U}; i != r.size(); ++i) {
        auto const& [ax, ay] = r[i];
        auto const& [bx, by] = r[(i + 1U) % r.size()];
        if ((ay > y) != (by > y) && x < ax + (y - ay) / (by - ay) * (bx - ax)) {
          in = !in;
        }
      }
    }
    return in;
  }

  // The nearest vertex of the skeleton this point can walk to, and how far.
  std::pair<std::size_t, double> attach(geo::latlng const& p) const {
    auto order = std::vector<std::pair<double, std::size_t>>{};
    for (auto v = std::size_t{0U}; v != pos_.size(); ++v) {
      order.emplace_back(geo::distance(p, to_latlng(pos_[v])), v);
    }
    auto const n = std::min(kNearest, order.size());
    std::ranges::partial_sort(order, begin(order) + static_cast<std::ptrdiff_t>(n));
    for (auto i = std::size_t{0U}; i != n; ++i) {
      if (free_.is_segment_inside(p, to_latlng(pos_[order[i].second]))) {
        return {order[i].second, order[i].first};
      }
    }
    return {kNoVertex, kInf};
  }

  // Cheapest walks along the skeleton from `src`, with predecessors.
  void walk(std::size_t const src,
            std::vector<double>& d,
            std::vector<std::size_t>& pred) const {
    d.assign(pos_.size(), kInf);
    pred.assign(pos_.size(), kNoVertex);
    d[src] = 0.0;
    auto q = std::priority_queue<std::pair<double, std::size_t>,
                                 std::vector<std::pair<double, std::size_t>>,
                                 std::greater<>>{};
    q.emplace(0.0, src);
    while (!q.empty()) {
      auto const [cost, u] = q.top();
      q.pop();
      if (cost > d[u]) {
        continue;
      }
      for (auto const& [v, w] : adj_[u]) {
        if (cost + w < d[v]) {
          d[v] = cost + w;
          pred[v] = u;
          q.emplace(d[v], v);
        }
      }
    }
  }

  area_free_space free_;
  double lat0_{0.0}, lng0_{0.0}, scale_{1.0};
  std::vector<std::vector<std::pair<double, double>>> rings_;  // meters
  std::vector<std::pair<double, double>> pos_;  // meters
  std::vector<std::vector<std::pair<std::size_t, double>>> adj_;
  std::size_t n_edges_{0U};
};

area_skeleton::area_skeleton(
    std::vector<std::vector<geo::latlng>> const& rings,
    std::vector<std::vector<geo::latlng>> const& obstacles)
    : impl_{std::make_unique<impl>(rings, obstacles)} {}

area_skeleton::~area_skeleton() = default;
area_skeleton::area_skeleton(area_skeleton&&) noexcept = default;
area_skeleton& area_skeleton::operator=(area_skeleton&&) noexcept = default;

bool area_skeleton::empty() const { return impl_->pos_.empty(); }

std::size_t area_skeleton::n_vertices() const { return impl_->pos_.size(); }

std::size_t area_skeleton::n_edges() const { return impl_->n_edges_; }

std::vector<std::pair<geo::latlng, geo::latlng>> area_skeleton::edges() const {
  auto out = std::vector<std::pair<geo::latlng, geo::latlng>>{};
  for (auto u = std::size_t{0U}; u != impl_->pos_.size(); ++u) {
    for (auto const& [v, w] : impl_->adj_[u]) {
      if (u < v) {
        out.emplace_back(impl_->to_latlng(impl_->pos_[u]),
                         impl_->to_latlng(impl_->pos_[v]));
      }
    }
  }
  return out;
}

std::vector<geo::latlng> area_skeleton::path(geo::latlng const& a,
                                             geo::latlng const& b) const {
  if (empty()) {
    return {};
  }
  auto const [from, from_cost] = impl_->attach(a);
  auto const [to, to_cost] = impl_->attach(b);
  if (from == kNoVertex || to == kNoVertex) {
    return {};
  }
  auto d = std::vector<double>{};
  auto pred = std::vector<std::size_t>{};
  impl_->walk(from, d, pred);
  if (d[to] == kInf) {
    return {};
  }

  auto out = std::vector<geo::latlng>{};
  for (auto v = to; v != kNoVertex; v = (v == from ? kNoVertex : pred[v])) {
    out.push_back(impl_->to_latlng(impl_->pos_[v]));
  }
  std::reverse(begin(out), end(out));
  out.insert(begin(out), a);
  out.push_back(b);
  return out;
}

double area_skeleton::distance(geo::latlng const& a,
                               geo::latlng const& b) const {
  if (empty()) {
    return kInf;
  }
  auto const [from, from_cost] = impl_->attach(a);
  auto const [to, to_cost] = impl_->attach(b);
  if (from == kNoVertex || to == kNoVertex) {
    return kInf;
  }
  auto d = std::vector<double>{};
  auto pred = std::vector<std::size_t>{};
  impl_->walk(from, d, pred);
  return d[to] == kInf ? kInf : from_cost + d[to] + to_cost;
}

}  // namespace osr
