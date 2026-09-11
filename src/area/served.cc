#include "osr/area/served.h"

#include <cmath>
#include <algorithm>
#include <limits>
#include <queue>
#include <utility>

namespace osr {

namespace {

constexpr auto kInf = std::numeric_limits<double>::infinity();

}  // namespace

std::vector<double> network_distances(
    walkable_area const& a,
    std::vector<network_edge> const& mapped,
    hash_map<std::int64_t, geo::latlng> const& node_pos,
    std::vector<area_connector> const& connectors) {
  auto idx_of = hash_map<std::int64_t, std::size_t>{};
  auto pos = std::vector<geo::latlng>{};
  auto adj = std::vector<std::vector<std::pair<std::size_t, double>>>{};

  auto const add_node = [&](std::int64_t const id, geo::latlng const& p) {
    if (auto const it = idx_of.find(id); it != end(idx_of)) {
      return it->second;
    }
    idx_of[id] = pos.size();
    pos.push_back(p);
    adj.emplace_back();
    return pos.size() - 1U;
  };
  auto const add_edge = [&](std::size_t const u, std::size_t const v) {
    if (u == v) {
      return;
    }
    auto const d = geo::distance(pos[u], pos[v]);
    adj[u].emplace_back(v, d);
    adj[v].emplace_back(u, d);
  };

  for (auto const& r : a.rings_) {
    for (auto i = std::size_t{0U}; i != r.ids_.size(); ++i) {
      auto const j = (i + 1U) % r.ids_.size();
      add_edge(add_node(r.ids_[i], r.points_[i]),
               add_node(r.ids_[j], r.points_[j]));
    }
  }
  for (auto const& e : mapped) {
    auto const pa = node_pos.find(e.a_);
    auto const pb = node_pos.find(e.b_);
    if (pa != end(node_pos) && pb != end(node_pos)) {
      add_edge(add_node(e.a_, pa->second), add_node(e.b_, pb->second));
    }
  }

  auto const k = connectors.size();
  auto out = std::vector<double>(k * k, kInf);
  auto d = std::vector<double>{};
  using entry = std::pair<double, std::size_t>;
  for (auto i = std::size_t{0U}; i != k; ++i) {
    auto const it = idx_of.find(connectors[i].node_);
    if (it == end(idx_of)) {
      continue;
    }
    d.assign(pos.size(), kInf);
    d[it->second] = 0.0;
    auto q = std::priority_queue<entry, std::vector<entry>, std::greater<>>{};
    q.emplace(0.0, it->second);
    while (!q.empty()) {
      auto const [cost, u] = q.top();
      q.pop();
      if (cost > d[u]) {
        continue;
      }
      for (auto const& [v, w] : adj[u]) {
        if (auto const next = cost + w; next < d[v]) {
          d[v] = next;
          q.emplace(next, v);
        }
      }
    }
    for (auto j = std::size_t{0U}; j != k; ++j) {
      if (auto const jt = idx_of.find(connectors[j].node_);
          jt != end(idx_of)) {
        out[i * k + j] = d[jt->second];
      }
    }
  }
  return out;
}

std::size_t n_stranded(std::vector<area_connector> const& connectors,
                       std::vector<double> const& network) {
  auto const k = connectors.size();
  auto stranded = std::size_t{0U};
  for (auto i = std::size_t{0U}; i != k; ++i) {
    if (connectors[i].on_ring_) {
      continue;
    }
    auto reachable = false;
    for (auto j = std::size_t{0U}; j != k && !reachable; ++j) {
      reachable = i != j && std::isfinite(network[i * k + j]);
    }
    stranded += reachable ? 0U : 1U;
  }
  return stranded;
}

bool served_without_geometry(std::vector<area_connector> const& connectors,
                             std::vector<double> const& network,
                             served_params const& params) {
  auto const k = connectors.size();
  for (auto i = std::size_t{0U}; i != k; ++i) {
    for (auto j = i + 1U; j != k; ++j) {
      if (!connectors[i].on_ring_ || !connectors[j].on_ring_) {
        continue;
      }
      auto const straight =
          geo::distance(connectors[i].pos_, connectors[j].pos_);
      if (straight < params.min_pair_distance_) {
        continue;
      }
      if (network[i * k + j] > params.max_detour_ * straight) {
        return false;
      }
    }
  }
  return n_stranded(connectors, network) == 0U;
}

relevance relevant_pairs(walkable_area const& a,
                         std::vector<area_connector> const& connectors,
                         std::vector<float> const& geodesic,
                         served_params const& params) {
  // Distance along the outline, both ways round.
  auto const& outline = a.rings_.front().points_;
  auto cum = std::vector<double>(outline.size() + 1U, 0.0);
  for (auto i = std::size_t{0U}; i != outline.size(); ++i) {
    cum[i + 1U] = cum[i] + geo::distance(outline[i],
                                         outline[(i + 1U) % outline.size()]);
  }
  auto const round_outline = [&](std::size_t const x, std::size_t const y) {
    auto const along = std::abs(cum[x] - cum[y]);
    return std::min(along, cum.back() - along);
  };

  auto const k = connectors.size();
  auto r = relevance{.relevant_ = std::vector<bool>(k * k, true),
                     .outline_detours_ = {}};
  for (auto i = std::size_t{0U}; i != k; ++i) {
    for (auto j = i + 1U; j != k; ++j) {
      auto const d = static_cast<double>(geodesic[i * k + j]);
      auto const& pi = connectors[i].outline_pos_;
      auto const& pj = connectors[j].outline_pos_;
      if (!std::isfinite(d) || d < params.min_pair_distance_ ||
          !pi.has_value() || !pj.has_value()) {
        continue;
      }
      auto const detour = round_outline(*pi, *pj) / d;
      r.outline_detours_.push_back(detour);
      if (detour < params.min_crossing_gain_) {
        r.relevant_[i * k + j] = r.relevant_[j * k + i] = false;
      }
    }
  }
  return r;
}

served_verdict served_by_geodesics(std::vector<area_connector> const& connectors,
                                   std::vector<double> const& network,
                                   std::vector<float> const& geodesic,
                                   std::vector<bool> const& relevant,
                                   served_params const& params) {
  auto const k = connectors.size();
  auto v = served_verdict{};
  for (auto i = std::size_t{0U}; i != k; ++i) {
    for (auto j = i + 1U; j != k; ++j) {
      auto const d = static_cast<double>(geodesic[i * k + j]);
      if (!relevant[i * k + j] || !std::isfinite(d) ||
          d < params.min_pair_distance_ || !connectors[i].on_ring_ ||
          !connectors[j].on_ring_) {
        continue;
      }
      auto const detour = network[i * k + j] / d;
      ++v.n_relevant_;
      v.detours_.push_back(detour);
      v.worst_detour_ = std::max(v.worst_detour_, detour);
    }
  }
  v.n_stranded_ = n_stranded(connectors, network);
  v.served_ = v.worst_detour_ <= params.max_detour_ && v.n_stranded_ == 0U;
  return v;
}

}  // namespace osr
