#include "osr/area/pipeline.h"

#include <cmath>
#include <algorithm>
#include <functional>
#include <string>
#include <utility>

#include "fmt/core.h"

#include "osmium/area/assembler.hpp"
#include "osmium/area/multipolygon_manager.hpp"
#include "osmium/handler.hpp"
#include "osmium/handler/node_locations_for_ways.hpp"
#include "osmium/index/map/flex_mem.hpp"
#include "osmium/io/pbf_input.hpp"
#include "osmium/osm/area.hpp"
#include "osmium/osm/way.hpp"
#include "osmium/relations/relations_manager.hpp"
#include "osmium/tags/tags_filter.hpp"
#include "osmium/visitor.hpp"

#include "cista/strong.h"

#include "geo/area_db.h"

#include "utl/enumerate.h"

#include "osr/extract/tags.h"

namespace osr {

namespace {

using index_t = osmium::index::map::FlexMem<osmium::unsigned_object_id_type,
                                            osmium::Location>;
using location_handler_t = osmium::handler::NodeLocationsForWays<index_t>;
using area_idx_t = cista::strong<std::uint32_t, struct walkable_area_idx_>;

tag_lookup tags_of(osmium::TagList const& t) {
  return [&t](std::string_view const key) -> std::optional<std::string_view> {
    auto const* const v = t[std::string{key}.c_str()];
    return v == nullptr ? std::nullopt : std::optional<std::string_view>{v};
  };
}

area_levels levels_of(osmium::TagList const& t) {
  return area_levels::of(tags_of(t));
}

bool is_linear(osmium::Way const& w) { return is_linear_way(tags_of(w.tags())); }

struct area_collector : public osmium::handler::Handler {
  explicit area_collector(walkable_options const& o) : options_{o} {}

  void area(osmium::Area const& a) {
    if (!is_walkable_area(tags_of(a.tags()), a.from_way(), options_)) {
      return;
    }

    auto const t = tags{a};
    if (t.has_level_) {
      ++n_with_level_;
    }
    auto rec = walkable_area{.id_ = a.orig_id(),
                             .from_way_ = a.from_way(),
                             .levels_ = levels_of(a.tags())};
    if (auto const* const name = a.tags()["name"]; name != nullptr) {
      rec.name_ = name;
    }

    auto const add_ring = [&](auto const& r) {
      auto& out = rec.rings_.emplace_back();
      for (auto const& n : r) {
        if (!n.location().valid()) {
          continue;
        }
        out.ids_.push_back(n.ref());
        out.points_.push_back({n.lat(), n.lon()});
      }
    };

    // Only the first outer ring is kept: an area with several disjoint
    // outer rings is several areas as far as crossing it goes.
    auto first = true;
    for (auto const& outer : a.outer_rings()) {
      if (!first) {
        break;
      }
      first = false;
      add_ring(outer);
      for (auto const& inner : a.inner_rings(outer)) {
        add_ring(inner);
      }
    }

    if (!rec.rings_.empty() && rec.rings_.front().points_.size() >= 3U) {
      areas_.push_back(std::move(rec));
    }
  }

  walkable_options options_;
  std::vector<walkable_area> areas_;
  int n_with_level_{0};
};

struct linear_collector : public osmium::handler::Handler {
  linear_collector(geo::area_db_lookup<area_idx_t> const& lookup,
                   area_data& data,
                   bool const strict)
      : lookup_{lookup}, data_{data}, strict_{strict} {}

  // Every area this node lies in - strictly inside by the rtree, or on the
  // boundary by being one of its ring nodes.
  std::vector<area_idx_t> const& areas_of(std::int64_t const id,
                                          geo::latlng const& pos) {
    scratch_.clear();
    lookup_.lookup(pos, hits_);
    scratch_.assign(begin(hits_), end(hits_));
    if (auto const it = data_.ring_areas_.find(id);
        it != end(data_.ring_areas_)) {
      for (auto const a : it->second) {
        scratch_.push_back(area_idx_t{a});
      }
    }
    std::ranges::sort(scratch_);
    scratch_.erase(std::ranges::unique(scratch_).begin(), end(scratch_));
    return scratch_;
  }

  area_levels const& levels(area_idx_t const a) const {
    return data_.areas_[to_idx(a)].levels_;
  }

  void node(osmium::Node const& n) {
    if (!n.location().valid()) {
      return;
    }
    auto const t = tags{n};
    // A lift or a marked entrance standing inside an area is an entry to it
    // even though no way of the area passes through it.
    if (t.is_elevator_ || t.is_entrance_) {
      loose_ends_.push_back(loose_end{.id_ = n.id(),
                                      .pos_ = {n.location().lat(),
                                               n.location().lon()},
                                      .levels_ = levels_of(n.tags())});
    }
  }

  void way(osmium::Way const& w) {
    auto const t = tags{w};

    // Barriers and buildings (see blocker_of) are not routable, so they are
    // collected here and then dropped. Containment is the strict rtree test,
    // never ring membership: a fence along the outline is the wall of the
    // area, not an obstacle inside it, and a building mapped as an inner
    // ring has its nodes on the boundary rather than within.
    if (auto const blocker = blocker_of(tags_of(w.tags()));
        blocker.kind_ != blocker::kind::kNone) {
      auto const obstacle = blocker.kind_ == blocker::kind::kBarrier;

      auto pts = std::vector<geo::latlng>{};
      auto inside = std::vector<std::vector<area_idx_t>>{};
      auto touched = std::vector<area_idx_t>{};
      for (auto const& n : w.nodes()) {
        if (!n.location().valid()) {
          continue;
        }
        auto const pos = geo::latlng{n.location().lat(), n.location().lon()};
        lookup_.lookup(pos, hits_);
        auto cur = std::vector<area_idx_t>{};
        for (auto const a : hits_) {
          if (blocker.stands_in(levels(a), strict_)) {
            cur.push_back(a);
            touched.push_back(a);
          }
        }
        pts.push_back(pos);
        inside.push_back(std::move(cur));
      }
      std::ranges::sort(touched);
      touched.erase(std::ranges::unique(touched).begin(), end(touched));

      for (auto const a : touched) {
        auto const in = [&](std::size_t const i) {
          return std::ranges::find(inside[i], a) != end(inside[i]);
        };
        if (obstacle) {
          // Keep a segment if either end is inside: a fence crossing the
          // outline blocks just as much as one wholly within it.
          for (auto i = std::size_t{0U}; i + 1U < pts.size(); ++i) {
            if (in(i) || in(i + 1U)) {
              data_.barriers_[to_idx(a)].push_back({pts[i], pts[i + 1U]});
            }
          }
        } else if (pts.size() >= 4U && pts.front() == pts.back()) {
          data_.buildings_[to_idx(a)].push_back(
              std::vector<geo::latlng>{begin(pts), end(pts) - 1});
        }
      }
    }

    if (!is_linear(w)) {
      return;
    }

    // access=private / access=no: reachable on the map, not in the world.
    auto const restricted =
        t.access_ == override::kBlacklist || t.private_access_;
    if (restricted) {
      ++data_.n_restricted_ways_;
    } else {
      for (auto const& n : w.nodes()) {
        data_.unrestricted_nodes_.insert(n.ref());
      }
    }
    auto const way_levels = levels_of(w.tags());
    if (t.has_level_) {
      ++data_.n_ways_with_level_;
    }
    ++data_.n_ways_;

    for (auto const& n : w.nodes()) {
      data_.node_levels_[n.ref()].merge(way_levels);
    }

    // The stretches of this way that stay inside an area: those are the
    // ways "mapped onto" the area, which the router can already use.
    auto prev_areas = std::vector<area_idx_t>{};
    auto prev_id = std::int64_t{0};
    auto have_prev = false;
    for (auto const& n : w.nodes()) {
      if (!n.location().valid()) {
        have_prev = false;
        continue;
      }
      auto const pos = geo::latlng{n.location().lat(), n.location().lon()};
      auto const& areas = areas_of(n.ref(), pos);
      if (!areas.empty()) {
        data_.node_pos_[n.ref()] = pos;
      }
      if (have_prev) {
        for (auto const a : areas) {
          if (std::ranges::find(prev_areas, a) == end(prev_areas)) {
            continue;
          }
          if (levels(a).matches_spatially(way_levels, strict_)) {
            data_.interior_edges_[to_idx(a)].push_back(
                {.a_ = prev_id, .b_ = n.ref()});
          } else {
            ++data_.n_rejected_edges_;
          }
        }
      }
      prev_areas.assign(begin(areas), end(areas));
      prev_id = n.ref();
      have_prev = true;
    }

    for (auto const* end : {&w.nodes().front(), &w.nodes().back()}) {
      if (end->location().valid()) {
        loose_ends_.push_back(
            loose_end{.id_ = end->ref(),
                      .pos_ = {end->location().lat(), end->location().lon()},
                      .levels_ = way_levels});
      }
    }
  }

  geo::area_db_lookup<area_idx_t> const& lookup_;
  area_data& data_;
  bool strict_;
  std::vector<loose_end> loose_ends_;

private:
  geo::area_db_lookup<area_idx_t>::rtree_results_t hits_;
  std::vector<area_idx_t> scratch_;
};

}  // namespace

area_data collect_areas(std::filesystem::path const& osm,
                        area_options const& options) {
  auto const file = osm.generic_string();
  auto data = area_data{};

  // Areas first: the second pass needs to know which area each node falls
  // in, and that needs the assembled rings.
  auto assembler_config = osmium::area::Assembler::config_type{};
  assembler_config.create_empty_areas = false;
  auto filter = osmium::TagsFilter{false};
  filter.add_rule(true, "highway");
  filter.add_rule(true, "place", "square");
  filter.add_rule(true, "public_transport", "platform");
  filter.add_rule(true, "railway", "platform");
  auto mp_manager = osmium::area::MultipolygonManager<osmium::area::Assembler>{
      assembler_config, filter};
  osmium::relations::read_relations(osmium::io::File{file}, mp_manager);

  auto collector = area_collector{options.walkable_};
  {
    auto index = index_t{};
    auto location_handler = location_handler_t{index};
    location_handler.ignore_errors();
    auto reader = osmium::io::Reader{file, osmium::io::read_meta::no};
    osmium::apply(reader, location_handler,
                  mp_manager.handler([&](osmium::memory::Buffer&& buffer) {
                    osmium::apply(buffer, collector);
                  }));
    reader.close();
  }
  data.areas_ = std::move(collector.areas_);
  data.n_areas_with_level_ = collector.n_with_level_;
  if (options.merge_) {
    data.merged_ = merge_overlapping(data.areas_);
  }

  auto const n = data.areas_.size();
  data.interior_edges_.resize(n);
  data.barriers_.resize(n);
  data.buildings_.resize(n);
  data.loose_ends_.resize(n);
  for (auto const [i, a] : utl::enumerate(data.areas_)) {
    auto seen = hash_set<std::int64_t>{};
    for (auto const& r : a.rings_) {
      for (auto const id : r.ids_) {
        if (seen.insert(id).second) {
          data.ring_areas_[id].push_back(static_cast<std::uint32_t>(i));
        }
      }
    }
  }

  // geo::area_db carries the rings and answers the contains-test.
  auto const db_dir =
      std::filesystem::temp_directory_path() /
      fmt::format("osr-walkable-areas-{:x}", std::hash<std::string>{}(file));
  std::filesystem::remove_all(db_dir);
  std::filesystem::create_directories(db_dir);
  {
    auto storage = geo::area_db_storage<area_idx_t>{
        db_dir, cista::mmap::protection::WRITE};
    for (auto const& a : data.areas_) {
      auto outers = std::vector<std::vector<geo::fixed_latlng>>{};
      auto inners = std::vector<std::vector<std::vector<geo::fixed_latlng>>>{};
      auto& outer = outers.emplace_back();
      auto& inner_group = inners.emplace_back();
      for (auto const& p : a.rings_.front().points_) {
        outer.push_back(geo::fixed_latlng::from_latlng(p));
      }
      for (auto i = std::size_t{1U}; i != a.rings_.size(); ++i) {
        auto& in = inner_group.emplace_back();
        for (auto const& p : a.rings_[i].points_) {
          in.push_back(geo::fixed_latlng::from_latlng(p));
        }
      }
      storage.add_area(outers, inners);
    }
    auto const lookup = geo::area_db_lookup<area_idx_t>{storage};

    auto linear = linear_collector{lookup, data, options.strict_levels_};
    {
      auto index = index_t{};
      auto location_handler = location_handler_t{index};
      location_handler.ignore_errors();
      auto reader = osmium::io::Reader{file, osmium::io::read_meta::no};
      osmium::apply(reader, location_handler, linear);
      reader.close();
    }

    auto hits = geo::area_db_lookup<area_idx_t>::rtree_results_t{};
    for (auto const& le : linear.loose_ends_) {
      lookup.lookup(le.pos_, hits);
      for (auto const a : hits) {
        data.loose_ends_[to_idx(a)].push_back(le);
      }
    }
  }
  std::filesystem::remove_all(db_dir);
  return data;
}

char const* to_str(area_status const s) {
  switch (s) {
    case area_status::kTooFewConnectors: return "too_few_connectors";
    case area_status::kTooBig: return "too_big";
    case area_status::kNoCrossing: return "no_crossing";
    case area_status::kServed: return "served";
    case area_status::kNeedsCells: return "needs_cells";
    case area_status::kMeshed: return "meshed";
    case area_status::kUnreachable: return "unreachable_pair";
  }
  return "?";
}

prepared_area prepare_area(area_data const& data,
                           std::size_t const i,
                           area_options const& o,
                           bool const always_geodesics,
                           bool const build_cells) {
  auto const& a = data.areas_[i];
  auto pa = prepared_area{};
  for (auto const& r : a.rings_) {
    pa.rings_.push_back(r.points_);
  }
  pa.barriers_ = data.barriers_[i];

  // Connectors: where the street network meets the rings, on a level this
  // area is on - and where a neighbouring area shares a ring node.
  auto const& outer = a.rings_.front();
  auto seen = hash_set<std::int64_t>{};
  for (auto const& r : a.rings_) {
    for (auto k = std::size_t{0U}; k != r.ids_.size(); ++k) {
      auto const id = r.ids_[k];

      auto meets_network = false;
      if (auto const it = data.node_levels_.find(id);
          it != end(data.node_levels_)) {
        meets_network = a.levels_.matches(it->second, o.strict_levels_);
      }

      auto meets_area = false;
      if (auto const it = data.ring_areas_.find(id);
          it != end(data.ring_areas_)) {
        auto n_matching = 0;
        for (auto const other : it->second) {
          n_matching +=
              a.levels_.matches(data.areas_[other].levels_, o.strict_levels_)
                  ? 1
                  : 0;
        }
        meets_area = n_matching > 1;
      }

      if ((!meets_network && !meets_area) || !seen.insert(id).second) {
        continue;
      }
      // Reachable on the map but not in the world: every way meeting the
      // area here is access=private or access=no.
      if (meets_network && !meets_area && !data.unrestricted_nodes_.contains(id)) {
        ++pa.counts_.n_restricted_;
        continue;
      }
      ++pa.counts_.n_boundary_;
      pa.connectors_.push_back(
          {.node_ = id,
           .pos_ = r.points_[k],
           .on_ring_ = true,
           .outline_pos_ = &r == &outer ? std::optional{k} : std::nullopt});
    }
  }

  // Outline nodes of areas merged into this one that now lie inside it: a
  // way passing through such a node used to enter there, and still does.
  for (auto const& [id, pos] : a.dissolved_) {
    auto const it = data.node_levels_.find(id);
    if (it == end(data.node_levels_) ||
        !a.levels_.matches(it->second, o.strict_levels_) ||
        !data.unrestricted_nodes_.contains(id) || !seen.insert(id).second) {
      continue;
    }
    ++pa.counts_.n_interior_;
    pa.connectors_.push_back({.node_ = id, .pos_ = pos, .on_ring_ = false});
  }

  // Stairs, lifts and footway stubs that end inside the area instead of on
  // its edge. These carry no ring node, so nothing above can see them.
  for (auto const& le : data.loose_ends_[i]) {
    if (!a.levels_.matches_spatially(le.levels_, o.strict_levels_)) {
      ++pa.counts_.n_rejected_interior_;
      continue;
    }
    auto const duplicate =
        std::ranges::any_of(pa.connectors_, [&](area_connector const& c) {
          return geo::distance(c.pos_, le.pos_) < 0.5;
        });
    if (duplicate) {
      continue;
    }
    ++pa.counts_.n_interior_;
    pa.connectors_.push_back(
        {.node_ = le.id_, .pos_ = le.pos_, .on_ring_ = false});
  }

  if (pa.connectors_.size() < 2U) {
    pa.status_ = area_status::kTooFewConnectors;
    return pa;
  }

  // Holes: the inner rings and the buildings standing in the area, cleaned
  // (see clean_holes).
  {
    auto holes = std::vector<std::vector<geo::latlng>>(
        std::next(begin(pa.rings_)), end(pa.rings_));
    auto const& buildings = data.buildings_[i];
    holes.insert(end(holes), begin(buildings), end(buildings));
    auto outline = std::move(pa.rings_.front());
    pa.rings_.clear();
    pa.rings_.push_back(std::move(outline));
    for (auto& h : clean_holes(pa.rings_.front(), holes)) {
      pa.rings_.push_back(std::move(h));
    }
  }
  for (auto const& r : pa.rings_) {
    pa.n_vertices_ += r.size();
  }
  for (auto const& b : pa.barriers_) {
    pa.n_vertices_ += b.size();
  }
  if (pa.n_vertices_ + pa.connectors_.size() > o.max_vertices_) {
    pa.status_ = area_status::kTooBig;
    return pa;
  }

  // Layer 1 without geometry.
  static auto const kNoInteriorEdges = std::vector<network_edge>{};
  pa.network_ = network_distances(
      a, o.interior_ways_ ? data.interior_edges_[i] : kNoInteriorEdges,
      data.node_pos_, pa.connectors_);
  pa.served_without_geometry_ =
      served_without_geometry(pa.connectors_, pa.network_, o.served_);
  if (pa.served_without_geometry_ && !always_geodesics) {
    pa.status_ = area_status::kServed;
    return pa;
  }

  // Layer 2.
  auto positions = std::vector<geo::latlng>{};
  for (auto const& c : pa.connectors_) {
    positions.push_back(c.pos_);
  }
  auto const& g = pa.geodesics_.emplace(pa.rings_, positions, pa.barriers_);
  auto const k = positions.size();
  pa.geodesic_distances_.resize(k * k);
  auto any_reachable = false;
  for (auto x = std::size_t{0U}; x != k; ++x) {
    for (auto y = std::size_t{0U}; y != k; ++y) {
      auto const d = g.distance(x, y);
      pa.geodesic_distances_[x * k + y] = d;
      any_reachable |= x < y && d != area_geodesics::kUnreachable;
    }
  }
  if (!any_reachable) {
    pa.status_ = area_status::kNoCrossing;
    return pa;
  }

  // Layer 1 with geometry.
  pa.relevance_ =
      relevant_pairs(a, pa.connectors_, pa.geodesic_distances_, o.served_);
  pa.verdict_ = served_by_geodesics(pa.connectors_, pa.network_,
                                    pa.geodesic_distances_,
                                    pa.relevance_.relevant_, o.served_);
  if (pa.verdict_.served_) {
    pa.status_ = area_status::kServed;
    return pa;
  }
  if (!build_cells) {
    pa.status_ = area_status::kNeedsCells;
    return pa;
  }

  // Layer 3.
  pa.cells_ = build_area_cells(
      positions, pa.geodesic_distances_, pa.relevance_.relevant_,
      {.threshold_ = o.cells_threshold_,
       .rings_ = o.shared_borders_ ? &pa.rings_ : nullptr},
      &pa.regions_);
  pa.status_ = pa.cells_.has_value() ? area_status::kMeshed
                                     : area_status::kUnreachable;
  return pa;
}

}  // namespace osr
