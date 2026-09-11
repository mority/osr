#include "osr/area/walkable.h"

#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <string>

#include "boost/geometry.hpp"

#include "fmt/core.h"

#include "geo/box.h"
#include "geo/constants.h"

#include "utl/enumerate.h"

#include "osr/util/levels.h"

namespace bg = boost::geometry;

namespace osr {

namespace {

bool one_of(std::optional<std::string_view> const value,
            std::initializer_list<std::string_view> options) {
  return value.has_value() && std::ranges::find(options, *value) != end(options);
}

std::optional<float> number(std::optional<std::string_view> const value) {
  if (!value.has_value()) {
    return std::nullopt;
  }
  auto const s = std::string{*value};
  char* end = nullptr;
  auto const x = std::strtof(s.c_str(), &end);
  return end == s.c_str() ? std::nullopt : std::optional{x};
}

float highest_level(level_bits_t const bits) {
  auto highest = kMinLevel;
  for (auto b = 0U; b != 64U; ++b) {
    if (((bits >> b) & 1U) != 0U) {
      highest = std::max(
          highest, level_t{static_cast<std::uint8_t>(b)}.to_float());
    }
  }
  return highest;
}

// Merging and hole cleanup work in OSM's own fixed-point unit, 1e-7 degrees,
// held as doubles: exact up to 2^53, far beyond the 1.8e9 they reach. Every
// vertex that survives a union is then bit-identical to its node.
using merge_point = bg::model::d2::point_xy<double>;
using merge_polygon = bg::model::polygon<merge_point>;
using merge_multi = bg::model::multi_polygon<merge_polygon>;

// Less overlap than this is two mappers' outlines grazing each other.
constexpr auto kMinMergeOverlap = 1.0;  // m^2

merge_point to_merge(geo::latlng const& p) {
  return {static_cast<double>(std::llround(p.lng() * 1e7)),
          static_cast<double>(std::llround(p.lat() * 1e7))};
}

geo::latlng from_merge(merge_point const& p) {
  return {p.y() / 1e7, p.x() / 1e7};
}

std::uint64_t point_key(merge_point const& p) {
  auto const x = static_cast<std::uint32_t>(static_cast<std::int64_t>(p.x()));
  auto const y = static_cast<std::uint32_t>(static_cast<std::int64_t>(p.y()));
  return (static_cast<std::uint64_t>(x) << 32U) | y;
}

merge_polygon ring_polygon(std::vector<geo::latlng> const& ring) {
  auto poly = merge_polygon{};
  for (auto const& p : ring) {
    poly.outer().push_back(to_merge(p));
  }
  bg::correct(poly);
  return poly;
}

// Nothing for an outline OSM got wrong (self-intersecting, say).
std::optional<merge_polygon> area_polygon(walkable_area const& a) {
  auto poly = merge_polygon{};
  for (auto const [i, r] : utl::enumerate(a.rings_)) {
    auto& out = i == 0U ? poly.outer() : poly.inners().emplace_back();
    for (auto const& p : r.points_) {
      out.push_back(to_merge(p));
    }
  }
  bg::correct(poly);
  return bg::is_valid(poly) ? std::optional{std::move(poly)} : std::nullopt;
}

double area_m2(double const units2, double const lat) {
  auto const m_per_unit = geo::kApproxDistanceLatDegrees / 1e7;
  return std::abs(units2) * m_per_unit * m_per_unit *
         std::cos(lat * geo::kPI / 180.0);
}

}  // namespace

area_levels area_levels::of(tag_lookup const& tags) {
  auto l = area_levels{};
  auto has_level = false;
  for (auto const key : {"level", "indoor:level"}) {
    if (auto const v = tags(key); v.has_value()) {
      has_level = true;
      l.bits_ |= parse_levels(*v);
    }
  }
  l.any_ = !has_level;
  if (auto const layer = tags("layer"); layer.has_value()) {
    l.layer_ = std::atoi(std::string{*layer}.c_str());
  }
  return l;
}

bool area_levels::matches(area_levels const& o, bool const strict) const {
  if (strict && !any_) {
    return !o.any_ && (bits_ & o.bits_) != 0U;
  }
  return any_ || o.any_ || (bits_ & o.bits_) != 0U;
}

bool area_levels::matches_spatially(area_levels const& o,
                                    bool const strict) const {
  return layer_ == o.layer_ && matches(o, strict);
}

bool area_levels::same_as(area_levels const& o) const {
  return any_ == o.any_ && bits_ == o.bits_ && layer_ == o.layer_;
}

void area_levels::merge(area_levels const& o) {
  any_ = any_ || o.any_;
  bits_ |= o.bits_;
}

bool is_walkable_area(tag_lookup const& tags,
                      bool const from_way,
                      walkable_options const& options) {
  // A station covers a whole station on every level; walking across it
  // would teleport between them.
  if (one_of(tags("public_transport"), {"station"}) ||
      one_of(tags("area"), {"no"})) {
    return false;
  }

  auto const is_square = options.place_square_ && one_of(tags("place"), {"square"});
  auto const is_surface =
      one_of(tags("highway"), {"pedestrian", "footway", "path", "living_street",
                               "service", "corridor", "platform", "steps"}) ||
      one_of(tags("public_transport"), {"platform"}) ||
      one_of(tags("railway"), {"platform"}) || is_square;
  if (!is_surface) {
    return false;
  }

  // A closed highway=footway without area=yes is a loop, not a surface. A
  // multipolygon is an area by construction.
  return !from_way || one_of(tags("area"), {"yes"}) || is_square;
}

bool is_linear_way(tag_lookup const& tags) {
  return tags("highway").has_value() && !one_of(tags("area"), {"yes"});
}

blocker blocker_of(tag_lookup const& tags) {
  auto b = blocker{.levels_ = area_levels::of(tags)};
  if (one_of(tags("barrier"), {"yes", "wall", "fence", "hedge",
                               "retaining_wall", "city_wall", "guard_rail"})) {
    b.kind_ = blocker::kind::kBarrier;
    return b;
  }

  auto const building = tags("building");
  // building=roof is a roof on posts - a petrol station, a platform canopy -
  // with open ground underneath.
  if (!building.has_value() || one_of(building, {"no", "roof"})) {
    return b;
  }
  b.kind_ = blocker::kind::kBuilding;
  if (auto const l = number(tags("building:min_level")); l.has_value()) {
    b.min_level_ = *l;
  } else if (auto const h = number(tags("min_height")); h.has_value() && *h > 2.5F) {
    b.min_level_ = 1.F;  // lifted by more than a storey: an upper floor
  }
  return b;
}

bool blocker::stands_in(area_levels const& area, bool const strict) const {
  switch (kind_) {
    case kind::kNone: return false;

    // A fence on a bridge does not block the plaza underneath, so a barrier
    // has to be on the area's layer.
    case kind::kBarrier: return area.matches_spatially(levels_, strict);

    // A building's `layer` orders it against the tunnels and bridges around
    // it - the Park Inn at Alexanderplatz is layer=1, over the station - it
    // does not lift it off the ground. Its level has to match, and an upper
    // floor on columns (building:min_level=1) only blocks the levels it is
    // on.
    case kind::kBuilding:
      if (min_level_ > 0.F &&
          (area.any_ || highest_level(area.bits_) < min_level_)) {
        return false;
      }
      return area.matches(levels_, strict);
  }
  return false;
}

std::vector<std::vector<geo::latlng>> clean_holes(
    std::vector<geo::latlng> const& outer,
    std::vector<std::vector<geo::latlng>> const& holes) {
  auto const outline = ring_polygon(outer);
  if (!bg::is_valid(outline)) {
    return holes;
  }
  try {
    auto merged = merge_multi{};
    for (auto const& h : holes) {
      auto const poly = ring_polygon(h);
      if (!bg::is_valid(poly)) {
        continue;
      }
      auto next = merge_multi{};
      bg::union_(merged, poly, next);
      merged = std::move(next);
    }
    auto clipped = merge_multi{};
    bg::intersection(outline, merged, clipped);

    auto out = std::vector<std::vector<geo::latlng>>{};
    for (auto const& poly : clipped) {
      auto& r = out.emplace_back();
      for (auto const& p : poly.outer()) {
        r.push_back(from_merge(p));
      }
    }
    return out;
  } catch (std::exception const&) {
    return holes;
  }
}

merge_stats merge_overlapping(std::vector<walkable_area>& areas) {
  auto const n = areas.size();
  auto polys = std::vector<std::optional<merge_polygon>>(n);
  auto boxes = std::vector<geo::box>(n);
  for (auto i = std::size_t{0U}; i != n; ++i) {
    polys[i] = area_polygon(areas[i]);
    for (auto const& p : areas[i].rings_.front().points_) {
      boxes[i].extend(p);
    }
  }

  auto root = std::vector<std::size_t>(n);
  for (auto i = std::size_t{0U}; i != n; ++i) {
    root[i] = i;
  }
  auto const find = [&](std::size_t x) {
    while (root[x] != x) {
      x = root[x] = root[root[x]];
    }
    return x;
  };
  for (auto i = std::size_t{0U}; i != n; ++i) {
    for (auto j = i + 1U; j != n; ++j) {
      if (!polys[i].has_value() || !polys[j].has_value() ||
          !areas[i].levels_.same_as(areas[j].levels_) ||
          !boxes[i].overlaps(boxes[j])) {
        continue;
      }
      try {
        auto overlap = merge_multi{};
        bg::intersection(*polys[i], *polys[j], overlap);
        if (area_m2(bg::area(overlap), boxes[i].min_.lat()) >
            kMinMergeOverlap) {
          root[find(i)] = find(j);
        }
      } catch (std::exception const&) {
        // no usable overlap: leave the two apart
      }
    }
  }

  auto groups = std::vector<std::vector<std::size_t>>(n);
  for (auto i = std::size_t{0U}; i != n; ++i) {
    groups[find(i)].push_back(i);
  }

  auto stats = merge_stats{};
  auto next_synthetic_id = std::int64_t{-1};
  auto out = std::vector<walkable_area>{};
  for (auto i = std::size_t{0U}; i != n; ++i) {
    auto const& group = groups[find(i)];
    if (group.size() == 1U) {
      out.push_back(std::move(areas[i]));
      continue;
    }
    if (group.front() != i) {
      continue;  // handled when its group's first member came up
    }

    auto merged = merge_multi{};
    try {
      for (auto const m : group) {
        auto next = merge_multi{};
        bg::union_(merged, *polys[m], next);
        merged = std::move(next);
      }
    } catch (std::exception const&) {
      merged.clear();
    }
    if (merged.empty()) {
      ++stats.n_failed_;
      for (auto const m : group) {
        out.push_back(std::move(areas[m]));
      }
      continue;
    }

    // The merged area is named after its largest member.
    auto largest = group.front();
    auto largest_area = 0.0;
    auto node_ids = hash_map<std::uint64_t, std::int64_t>{};
    for (auto const m : group) {
      if (auto const a = std::abs(bg::area(*polys[m])); a > largest_area) {
        largest_area = a;
        largest = m;
      }
      for (auto const& r : areas[m].rings_) {
        for (auto k = std::size_t{0U}; k != r.ids_.size(); ++k) {
          node_ids.emplace(point_key(to_merge(r.points_[k])), r.ids_[k]);
        }
      }
    }

    auto base = walkable_area{.id_ = areas[largest].id_,
                              .from_way_ = areas[largest].from_way_,
                              .name_ = areas[largest].name_,
                              .levels_ = areas[largest].levels_};
    for (auto const m : group) {
      if (base.name_.empty()) {
        base.name_ = areas[m].name_;
      }
      base.members_.push_back(fmt::format(
          "{}/{}", areas[m].from_way_ ? "way" : "relation", areas[m].id_));
    }

    auto on_outline = hash_set<std::uint64_t>{};
    auto const to_ring = [&](auto const& path) {
      auto r = area_ring{};
      for (auto const& p : path) {  // closed, like the rings osmium assembles
        auto const key = point_key(p);
        on_outline.insert(key);
        auto const it = node_ids.find(key);
        r.ids_.push_back(it != end(node_ids) ? it->second
                                             : next_synthetic_id--);
        r.points_.push_back(from_merge(p));
      }
      if (r.ids_.size() > 1U && r.points_.front() == r.points_.back()) {
        r.ids_.back() = r.ids_.front();
      }
      return r;
    };

    // One area per polygon of the union: normally one, plus any island
    // standing in a hole of it.
    auto pieces = std::vector<std::pair<walkable_area, merge_polygon>>{};
    for (auto const& poly : merged) {
      auto rec = base;
      rec.rings_.push_back(to_ring(poly.outer()));
      for (auto const& hole : poly.inners()) {
        rec.rings_.push_back(to_ring(hole));
      }
      pieces.emplace_back(std::move(rec), poly);
    }

    // The members' outline nodes that ended up inside: a way that met an
    // outline there still enters the merged area there.
    for (auto const m : group) {
      for (auto const& r : areas[m].rings_) {
        for (auto k = std::size_t{0U}; k != r.ids_.size(); ++k) {
          auto const p = to_merge(r.points_[k]);
          if (on_outline.contains(point_key(p))) {
            continue;
          }
          for (auto& [rec, poly] : pieces) {
            if (bg::covered_by(p, poly)) {
              rec.dissolved_.emplace_back(r.ids_[k], r.points_[k]);
              break;
            }
          }
        }
      }
    }

    stats.n_members_ += static_cast<int>(group.size());
    stats.n_merged_ += static_cast<int>(pieces.size());
    for (auto& [rec, poly] : pieces) {
      out.push_back(std::move(rec));
    }
  }
  areas = std::move(out);
  return stats;
}

}  // namespace osr
