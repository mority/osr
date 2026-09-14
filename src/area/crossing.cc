#include "osr/area/crossing.h"

#include <cmath>
#include <utility>

namespace osr {

crossing_drawer::crossing_drawer(area_graph const& graph,
                                 std::vector<area_geometry> areas)
    : graph_{graph} {
  for (auto& a : areas) {
    auto e = std::make_unique<entry>();
    e->geometry_ = std::move(a);
    areas_.push_back(std::move(e));
  }
}

area_geodesics const& crossing_drawer::geodesics(std::size_t const area) const {
  auto& e = *areas_[area];
  std::call_once(e.once_, [&]() {
    // Drawn paths keep kWallClearance off the walls - that is where people
    // walk. What the router paid for the crossing comes from layers 3 and 4,
    // which are fitted to the geodesics without it.
    e.geodesics_ = std::make_unique<area_geodesics>(
        e.geometry_.rings_, e.geometry_.connectors_, e.geometry_.barriers_,
        geodesic_options{.clearance_ = kWallClearance});
  });
  return *e.geodesics_;
}

std::vector<area_crossing> crossing_drawer::crossings(path const& p) const {
  auto const is_hub = [&](node_idx_t const n) { return graph_.is_hub(n); };

  auto out = std::vector<area_crossing>{};
  auto const& s = p.segments_;
  for (auto i = std::size_t{0U}; i < s.size(); ++i) {
    if (is_hub(s[i].from_) || !is_hub(s[i].to_)) {
      continue;  // not entering an area
    }
    auto j = i + 1U;
    while (j < s.size() && is_hub(s[j].from_) && is_hub(s[j].to_)) {
      ++j;  // hub -> neighbouring hub
    }
    if (j >= s.size() || !is_hub(s[j].from_) || is_hub(s[j].to_)) {
      continue;  // a route may not end on a hub; nothing to draw
    }

    auto const area = graph_.area_of(s[i].to_);
    auto c = area_crossing{.first_segment_ = i,
                           .last_segment_ = j,
                           .area_ = area,
                           .level_ = graph_.level_of(s[i].to_)};
    for (auto k = i; k <= j; ++k) {
      c.model_distance_ += s[k].dist_;
    }
    auto const from = graph_.connector_of(area, s[i].from_);
    auto const to = graph_.connector_of(area, s[j].to_);
    if (from.has_value() && to.has_value() && area < areas_.size()) {
      auto const& g = geodesics(area);
      c.geodesic_ = g.path(*from, *to);
      if (auto const d = g.distance(*from, *to);
          d != area_geodesics::kUnreachable) {
        c.geodesic_distance_ = d;
      }
    }
    out.push_back(std::move(c));
    i = j;
  }
  return out;
}

path crossing_drawer::drawn(path const& p) const {
  auto const cs = crossings(p);
  auto out = p;
  out.segments_.clear();

  auto next = begin(cs);
  for (auto i = std::size_t{0U}; i < p.segments_.size(); ++i) {
    if (next != end(cs) && next->first_segment_ == i &&
        !next->geodesic_.empty()) {
      auto seg = path::segment{};
      seg.polyline_ = next->geodesic_;
      seg.from_level_ = seg.to_level_ = next->level_;
      seg.from_ = p.segments_[i].from_;
      seg.to_ = p.segments_[next->last_segment_].to_;
      seg.cost_ = 0U;
      for (auto k = i; k <= next->last_segment_; ++k) {
        seg.cost_ = static_cast<cost_t>(seg.cost_ + p.segments_[k].cost_);
        seg.duration_ += p.segments_[k].duration_;
      }
      seg.dist_ = static_cast<distance_t>(std::lround(
          next->geodesic_distance_ >= 0.0 ? next->geodesic_distance_
                                          : next->model_distance_));
      seg.mode_ = p.segments_[i].mode_;
      out.segments_.push_back(std::move(seg));
      i = next->last_segment_;
      ++next;
      continue;
    }
    if (next != end(cs) && next->first_segment_ == i) {
      ++next;  // no line to draw: keep the segments as routed
    }
    out.segments_.push_back(p.segments_[i]);
  }
  return out;
}

}  // namespace osr
