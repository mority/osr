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
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <queue>
#include <string>
#include <string_view>
#include <vector>

#include "fmt/core.h"
#include "fmt/format.h"

#include "utl/enumerate.h"
#include "utl/parser/arg_parser.h"
#include "utl/parser/cstr.h"

#include "osmium/area/assembler.hpp"
#include "osmium/area/multipolygon_manager.hpp"
#include "osmium/handler.hpp"
#include "osmium/handler/node_locations_for_ways.hpp"
#include "osmium/index/map/flex_mem.hpp"
#include "osmium/io/pbf_input.hpp"
#include "osmium/io/xml_input.hpp"
#include "osmium/osm/area.hpp"
#include "osmium/osm/way.hpp"
#include "osmium/relations/relations_manager.hpp"
#include "osmium/visitor.hpp"

#include "cista/strong.h"

#include "geo/area_db.h"
#include "geo/box.h"
#include "geo/latlng.h"

#include "osr/area/geodesic.h"
#include "osr/extract/tags.h"
#include "osr/types.h"

namespace {

using index_t = osmium::index::map::FlexMem<osmium::unsigned_object_id_type,
                                            osmium::Location>;
using location_handler_t = osmium::handler::NodeLocationsForWays<index_t>;

using area_idx_t = cista::strong<std::uint32_t, struct area_idx_>;

// Which levels something is on. An object with no `level` tag matches
// everything, following what foot::get_target_level already does with
// kNoLevel - most of OSM is untagged and treating that as "ground floor only"
// would disconnect far more than it fixes.
struct level_set {
  static level_set of(osr::tags const& t, osmium::TagList const& raw) {
    auto const* const l = raw["layer"];
    auto const layer = l == nullptr ? 0 : std::atoi(l);
    return t.has_level_ ? level_set{false, t.level_bits_, layer}
                        : level_set{true, 0U, layer};
  }

  // For things joined only by geometry - a way whose nodes fall inside an
  // area's outline without sharing any of them. A footbridge over a plaza is
  // inside its polygon and connects to nothing, and `layer` is what says so.
  // osr::tags does not parse `layer` at all, so this reads it directly.
  //
  // Deliberately NOT applied where a node is genuinely shared: a mapper who
  // joined a staircase to a plaza at one node meant them to connect, whatever
  // layer the staircase carries.
  bool matches_spatially(level_set const& o, bool const strict) const {
    return layer_ == o.layer_ && matches(o, strict);
  }

  // `*this` is the area. Under the permissive rule an untagged object joins
  // anything, which is right outdoors and dangerous indoors: at Chatelet a
  // footway stub with no level tag would attach to an area three floors below
  // it. The strict rule makes an area that states its level demand the same
  // from whatever connects to it.
  bool matches(level_set const& o, bool const strict) const {
    if (strict && !any_) {
      return !o.any_ && (bits_ & o.bits_) != 0U;
    }
    return any_ || o.any_ || (bits_ & o.bits_) != 0U;
  }

  void merge(level_set const& o) {
    any_ = any_ || o.any_;
    bits_ |= o.bits_;
  }

  bool any_{false};
  osr::level_bits_t bits_{0U};
  int layer_{0};
};

bool g_strict_levels = false;

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

// How much shorter the crossing has to be than the way round for the pair to
// count as one the area model is responsible for.
constexpr auto kRelevantDetour = 1.2;

// How much worse than the geodesic the existing ways may be before the area is
// judged to need meshing of its own.
constexpr auto kServedRatio = 1.2;

// Depth of the placeholder subdivision used to measure error composition.
constexpr auto kMaxSplitDepth = std::size_t{3U};

// The greedy build tests O(n^2) candidate diagonals, each with an O(n)
// visibility check, so it needs a tighter cap than the oracle itself.
constexpr auto kMaxTriangulationVertices = std::size_t{110U};

// Adaptive subdivision needs headroom beyond the fixed-depth experiments.
constexpr auto kAdaptiveMaxDepth = std::size_t{6U};

// How many portals a single cell border may carry. One portal is the
// assumption the whole flat-cost framework was built on; it is only neutral
// for routes already inside the same pair of cells.
std::size_t g_portals_per_border = 1U;

bool one_of_sv(std::string_view const value,
               std::initializer_list<char const*> options) {
  return !value.empty() && std::ranges::any_of(options, [&](char const* o) {
    return value == o;
  });
}

bool one_of(char const* value, std::initializer_list<char const*> options) {
  return value != nullptr && std::ranges::any_of(options, [&](char const* o) {
           return std::strcmp(value, o) == 0;
         });
}

// A way whose interior a pedestrian can cross freely. Deliberately excludes
// `area:highway`-only ways (that scheme draws the road surface alongside the
// linear way that already carries the routing, so meshing it would build a
// second parallel network) and public_transport=station (an area covering the
// whole station would teleport between levels).
bool is_routable_area(osmium::Area const& area) {
  auto const& t = area.tags();
  if (t.has_tag("public_transport", "station") || t.has_tag("area", "no")) {
    return false;
  }

  auto const* const place = t["place"];
  auto const is_surface =
      one_of(t["highway"], {"pedestrian", "footway", "path", "living_street",
                            "service", "corridor", "platform", "steps"}) ||
      one_of(t["public_transport"], {"platform"}) ||
      one_of(t["railway"], {"platform"}) || one_of(place, {"square"});
  if (!is_surface) {
    return false;
  }

  // A closed `highway=footway` without `area=yes` is a loop, not a surface.
  // A multipolygon relation is an area by construction.
  return !area.from_way() || t.has_tag("area", "yes") ||
         one_of(place, {"square"});
}

// The linear network: ways that carry routing in the ordinary way. Their nodes
// are what the areas connect to.
bool is_linear_routable(osmium::Way const& w) {
  return w.tags().has_key("highway") && !w.tags().has_tag("area", "yes");
}

// A point where the linear network stops inside an area rather than on its
// edge - the end of a staircase, a lift shaft, a footway stub. osm-hints.txt
// asks for these explicitly: they are entries the boundary never sees.
struct loose_end {
  osmium::object_id_type id_{0};
  geo::latlng pos_;
  level_set levels_;
};

// An edge of the existing linear network that lies wholly inside one area.
struct interior_edge {
  osmium::object_id_type a_, b_;
};

struct linear_collector : public osmium::handler::Handler {
  linear_collector(geo::area_db_lookup<area_idx_t> const& lookup,
                   osr::hash_map<osmium::object_id_type,
                                 std::vector<area_idx_t>> const& ring_areas,
                   std::vector<level_set> const& area_levels)
      : lookup_{lookup},
        ring_areas_{ring_areas},
        area_levels_{area_levels},
        interior_edges_(area_levels.size()),
        barriers_(area_levels.size()),
        buildings_(area_levels.size()) {}

  // Every area this node lies in - strictly inside by the rtree, or on the
  // boundary by being one of its ring nodes.
  std::vector<area_idx_t> const& areas_of(osmium::object_id_type const id,
                                          geo::latlng const& pos) {
    scratch_.clear();
    lookup_.lookup(pos, hits_);
    scratch_.assign(begin(hits_), end(hits_));
    if (auto const it = ring_areas_.find(id); it != end(ring_areas_)) {
      scratch_.insert(end(scratch_), begin(it->second), end(it->second));
    }
    std::ranges::sort(scratch_);
    scratch_.erase(std::ranges::unique(scratch_).begin(), end(scratch_));
    return scratch_;
  }

  void node(osmium::Node const& n) {
    if (!n.location().valid()) {
      return;
    }
    auto const t = osr::tags{n};
    // A lift or a marked entrance standing inside an area is an entry to it
    // even though no way of the area passes through it.
    if (t.is_elevator_ || t.is_entrance_) {
      loose_ends_.push_back(loose_end{n.id(),
                                      {n.location().lat(), n.location().lon()},
                                      level_set::of(t, n.tags())});
    }
  }

  // Physical obstacles drawn as ways. A fence across a plaza is invisible to
  // the polygon, so both the "already served" test and any visibility graph
  // built later would happily cross it.
  static bool is_obstacle(osr::tags const& t) {
    return one_of_sv(t.barrier_, {"yes", "wall", "fence", "hedge",
                                  "retaining_wall", "city_wall", "guard_rail"});
  }

  void way(osmium::Way const& w) {
    auto const raw_t = osr::tags{w};

    // Barriers and buildings are not routable, so they are counted here and
    // then dropped. Containment is the strict rtree test, never ring
    // membership: a fence along the outline is the wall of the area, not an
    // obstacle inside it, and a building mapped as an inner ring has its nodes
    // on the boundary rather than within.
    if (is_obstacle(raw_t) || w.tags().has_key("building")) {
      auto const obstacle = is_obstacle(raw_t);
      auto const levels = level_set::of(raw_t, w.tags());

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
          if (area_levels_[to_idx(a)].matches_spatially(levels,
                                                        g_strict_levels)) {
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
              barriers_[to_idx(a)].push_back({pts[i], pts[i + 1U]});
            }
          }
        } else if (pts.size() >= 4U && pts.front() == pts.back()) {
          buildings_[to_idx(a)].push_back(
              std::vector<geo::latlng>{begin(pts), end(pts) - 1});
        }
      }

      if (!is_linear_routable(w)) {
        return;
      }
    }

    if (!is_linear_routable(w)) {
      return;
    }
    auto const t = raw_t;

    // access=private / access=no: reachable on the map, not in the world.
    auto const restricted =
        t.access_ == osr::override::kBlacklist || t.private_access_;
    if (restricted) {
      ++n_restricted_ways_;
    } else {
      for (auto const& n : w.nodes()) {
        unrestricted_nodes_.insert(n.ref());
      }
    }
    auto const levels = level_set::of(t, w.tags());
    if (t.has_level_) {
      ++n_with_level_;
    }
    ++n_ways_;

    for (auto const& n : w.nodes()) {
      node_levels_[n.ref()].merge(levels);
    }

    // Record the stretches of this way that stay inside an area: those are
    // the ways "mapped onto" the area, which the router can already use.
    auto prev_areas = std::vector<area_idx_t>{};
    auto prev_id = osmium::object_id_type{0};
    auto have_prev = false;
    for (auto const& n : w.nodes()) {
      if (!n.location().valid()) {
        have_prev = false;
        continue;
      }
      auto const pos = geo::latlng{n.location().lat(), n.location().lon()};
      auto const& areas = areas_of(n.ref(), pos);
      if (!areas.empty()) {
        node_pos_[n.ref()] = pos;
      }
      if (have_prev) {
        for (auto const a : areas) {
          if (std::ranges::find(prev_areas, a) == end(prev_areas)) {
            continue;
          }
          if (area_levels_[to_idx(a)].matches_spatially(levels,
                                                        g_strict_levels)) {
            interior_edges_[to_idx(a)].push_back({prev_id, n.ref()});
          } else {
            ++n_rejected_edges_;
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
            loose_end{end->ref(),
                      {end->location().lat(), end->location().lon()},
                      levels});
      }
    }
  }

  geo::area_db_lookup<area_idx_t> const& lookup_;
  osr::hash_map<osmium::object_id_type, std::vector<area_idx_t>> const&
      ring_areas_;
  std::vector<level_set> const& area_levels_;

  osr::hash_map<osmium::object_id_type, level_set> node_levels_;
  osr::hash_map<osmium::object_id_type, geo::latlng> node_pos_;
  std::vector<loose_end> loose_ends_;
  std::vector<std::vector<interior_edge>> interior_edges_;
  std::vector<std::vector<std::vector<geo::latlng>>> barriers_;
  std::vector<std::vector<std::vector<geo::latlng>>> buildings_;
  osr::hash_set<osmium::object_id_type> unrestricted_nodes_;
  int n_restricted_ways_{0};
  int n_ways_{0};
  int n_with_level_{0};
  int n_rejected_edges_{0};

private:
  geo::area_db_lookup<area_idx_t>::rtree_results_t hits_;
  std::vector<area_idx_t> scratch_;
};

struct ring {
  std::vector<osmium::object_id_type> ids_;
  std::vector<geo::latlng> points_;
};

struct area_record {
  osmium::object_id_type id_{0};
  bool from_way_{false};
  std::string name_;
  level_set levels_;
  std::vector<ring> rings_;
};

struct area_collector : public osmium::handler::Handler {
  void area(osmium::Area const& a) {
    if (!is_routable_area(a)) {
      return;
    }

    auto const t = osr::tags{a};
    if (t.has_level_) {
      ++n_with_level_;
    }
    auto rec = area_record{.id_ = a.orig_id(),
                           .from_way_ = a.from_way(),
                           .levels_ = level_set::of(t, a.tags())};
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

    // Only the first outer ring is kept: an area with several disjoint outer
    // rings is several areas as far as crossing it goes, and splitting them
    // properly is extraction work, not measurement work.
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

  std::vector<area_record> areas_;
  int n_with_level_{0};
};

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

// Shortest distances between the area's connectors over the ways that already
// exist: its ring, plus every linear-network edge lying inside it. A connector
// with no such edge stays unreachable, which is the honest answer - a stub
// ending inside an area cannot be crossed to without crossing the area.
std::vector<double> network_distances(
    area_record const& a,
    std::vector<interior_edge> const& interior,
    osr::hash_map<osmium::object_id_type, geo::latlng> const& node_pos,
    std::vector<osmium::object_id_type> const& connector_ids) {
  auto idx_of = osr::hash_map<osmium::object_id_type, std::size_t>{};
  auto pos = std::vector<geo::latlng>{};
  auto adj = std::vector<std::vector<std::pair<std::size_t, double>>>{};

  auto const add_node = [&](osmium::object_id_type const id,
                            geo::latlng const& p) {
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
  for (auto const& e : interior) {
    auto const pa = node_pos.find(e.a_);
    auto const pb = node_pos.find(e.b_);
    if (pa != end(node_pos) && pb != end(node_pos)) {
      add_edge(add_node(e.a_, pa->second), add_node(e.b_, pb->second));
    }
  }

  auto const k = connector_ids.size();
  auto const inf = std::numeric_limits<double>::infinity();
  auto out = std::vector<double>(k * k, inf);
  auto d = std::vector<double>{};
  using entry = std::pair<double, std::size_t>;
  for (auto i = std::size_t{0U}; i != k; ++i) {
    auto const it = idx_of.find(connector_ids[i]);
    if (it == end(idx_of)) {
      continue;
    }
    d.assign(pos.size(), inf);
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
      if (auto const jt = idx_of.find(connector_ids[j]); jt != end(idx_of)) {
        out[i * k + j] = d[jt->second];
      }
    }
  }
  return out;
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
  ok = ok && g.is_connector_reachable(c.portal_);
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
        auto const d = g.distance(pts[i], pts[j]);
        if (d == osr::area_geodesics::kUnreachable) {
          return {};
        }
        lo = std::min<double>(lo, d);
        hi = std::max<double>(hi, d);
      }
    }
    if (lo == inf) {
      continue;
    }
    auto const c = 0.5 * (lo + hi);
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

    for (auto iter = 0; iter != 4; ++iter) {
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

// Cumulative arc length around the outer ring, computed once per area.
std::vector<double> ring_arc_lengths(ring const& outer) {
  auto cum = std::vector<double>(outer.points_.size() + 1U, 0.0);
  for (auto i = std::size_t{0U}; i != outer.points_.size(); ++i) {
    cum[i + 1U] =
        cum[i] + geo::distance(outer.points_[i],
                               outer.points_[(i + 1U) % outer.points_.size()]);
  }
  return cum;
}

// Shortest way between two connectors that stays on the outer ring - what osr
// does today.
double perimeter_distance(std::vector<double> const& cum,
                          std::size_t const a,
                          std::size_t const b) {
  auto const along = std::abs(cum[a] - cum[b]);
  return std::min(along, cum.back() - along);
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

int main(int argc, char** argv) {
  if (argc < 2) {
    fmt::print(stderr, "usage: {} <osm-file> [more-osm-files...]\n", argv[0]);
    return 1;
  }

  for (auto arg = 1; arg != argc; ++arg) {
    if (std::string_view{argv[arg]} == "--strict-levels") {
      g_strict_levels = true;
      continue;
    }
    auto const path = std::string{argv[arg]};
    fmt::print("\n===== {} =====\n",
               std::filesystem::path{path}.filename().string());

    auto assembler_config = osmium::area::Assembler::config_type{};
    assembler_config.create_empty_areas = false;

    auto filter = osmium::TagsFilter{false};
    filter.add_rule(true, "highway");
    filter.add_rule(true, "place", "square");
    filter.add_rule(true, "public_transport", "platform");
    filter.add_rule(true, "railway", "platform");

    auto mp_manager =
        osmium::area::MultipolygonManager<osmium::area::Assembler>{
            assembler_config, filter};
    osmium::relations::read_relations(osmium::io::File{path}, mp_manager);

    // Areas come first: the linear pass needs to know which area each node
    // falls in, and that needs the assembled rings.
    auto collector = area_collector{};
    {
      auto index = index_t{};
      auto location_handler = location_handler_t{index};
      location_handler.ignore_errors();
      auto reader = osmium::io::Reader{path, osmium::io::read_meta::no};
      osmium::apply(reader, location_handler,
                    mp_manager.handler([&](osmium::memory::Buffer&& buffer) {
                      osmium::apply(buffer, collector);
                    }));
      reader.close();
    }

    auto ring_areas =
        osr::hash_map<osmium::object_id_type, std::vector<area_idx_t>>{};
    auto area_levels = std::vector<level_set>{};
    for (auto const [i, a] : utl::enumerate(collector.areas_)) {
      area_levels.push_back(a.levels_);
      auto seen = osr::hash_set<osmium::object_id_type>{};
      for (auto const& r : a.rings_) {
        for (auto const id : r.ids_) {
          if (seen.insert(id).second) {
            ring_areas[id].push_back(area_idx_t{static_cast<std::uint32_t>(i)});
          }
        }
      }
    }

    // geo::area_db carries the rings and answers the contains-test, which is
    // what osm-hints.txt asks be reused rather than rebuilt.
    auto const db_dir =
        std::filesystem::temp_directory_path() /
        fmt::format("osr-area-stats-{:x}", std::hash<std::string>{}(path));
    std::filesystem::remove_all(db_dir);
    std::filesystem::create_directories(db_dir);
    auto storage = geo::area_db_storage<area_idx_t>{
        db_dir, cista::mmap::protection::WRITE};
    for (auto const& a : collector.areas_) {
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

    auto linear = linear_collector{lookup, ring_areas, area_levels};
    {
      auto index = index_t{};
      auto location_handler = location_handler_t{index};
      location_handler.ignore_errors();
      auto reader = osmium::io::Reader{path, osmium::io::read_meta::no};
      osmium::apply(reader, location_handler, linear);
      reader.close();
    }

    auto loose_by_area =
        std::vector<std::vector<loose_end>>(collector.areas_.size());
    {
      auto hits = geo::area_db_lookup<area_idx_t>::rtree_results_t{};
      for (auto const& le : linear.loose_ends_) {
        lookup.lookup(le.pos_, hits);
        for (auto const a : hits) {
          loose_by_area[to_idx(a)].push_back(le);
        }
      }
    }
    std::filesystem::remove_all(db_dir);

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

    for (auto const [area_i, a] : utl::enumerate(collector.areas_)) {
      auto rings = std::vector<std::vector<geo::latlng>>{};
      auto n_vertices = std::size_t{0U};
      for (auto const& r : a.rings_) {
        rings.push_back(r.points_);
        n_vertices += r.points_.size();
      }

      // Connector positions, plus where each sits on the outer ring so the
      // perimeter route can be measured against the geodesic.
      auto connectors = std::vector<geo::latlng>{};
      auto connector_ids = std::vector<osmium::object_id_type>{};
      auto on_ring = std::vector<char>{};
      auto outer_pos = std::vector<std::size_t>{};
      auto const& outer = a.rings_.front();
      auto seen = osr::hash_set<osmium::object_id_type>{};
      for (auto const& r : a.rings_) {
        for (auto i = std::size_t{0U}; i != r.ids_.size(); ++i) {
          auto const id = r.ids_[i];

          // The linear network meets the boundary here, on a level this area
          // is actually on.
          auto meets_network = false;
          if (auto const it = linear.node_levels_.find(id);
              it != end(linear.node_levels_)) {
            meets_network = a.levels_.matches(it->second, g_strict_levels);
          }

          // Or a neighbouring area shares this node, on a shared level.
          auto meets_area = false;
          if (auto const it = ring_areas.find(id); it != end(ring_areas)) {
            auto n_matching = 0;
            for (auto const other : it->second) {
              n_matching +=
                  a.levels_.matches(area_levels[to_idx(other)], g_strict_levels)
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
          if (meets_network && !meets_area &&
              !linear.unrestricted_nodes_.contains(id)) {
            ++n_restricted_connectors;
            continue;
          }
          ++n_boundary_connectors;
          connector_ids.push_back(id);
          on_ring.push_back(1);
          connectors.push_back(r.points_[i]);
          outer_pos.push_back(
              &r == &outer ? i : std::numeric_limits<std::size_t>::max());
        }
      }

      // Stairs, lifts and footway stubs that end inside the area instead of on
      // its edge. These carry no ring node, so nothing above can see them.
      for (auto const& le : loose_by_area[area_i]) {
        if (!a.levels_.matches_spatially(le.levels_, g_strict_levels)) {
          ++n_rejected_interior;
          continue;
        }
        auto const duplicate =
            std::ranges::any_of(connectors, [&](geo::latlng const& c) {
              return geo::distance(c, le.pos_) < 0.5;
            });
        if (duplicate) {
          continue;
        }
        ++n_interior_connectors;
        connector_ids.push_back(le.id_);
        on_ring.push_back(0);
        connectors.push_back(le.pos_);
        outer_pos.push_back(std::numeric_limits<std::size_t>::max());
      }

      if (connectors.size() < 2U) {
        ++n_no_crossing;
        continue;
      }
      // Epstein & Sack's counting/sampling DP is defined over subpolygons
      // P(i,j) of a SIMPLE polygon, so it does not survive holes - and every
      // barrier and building we just added is a hole.
      n_no_inner_rings += a.rings_.size() == 1U ? 1 : 0;
      n_simple_polygon +=
          (a.rings_.size() == 1U && linear.barriers_[area_i].empty() &&
           linear.buildings_[area_i].empty())
              ? 1
              : 0;

      auto const& barriers = linear.barriers_[area_i];
      auto const& buildings = linear.buildings_[area_i];
      if (!barriers.empty()) {
        ++n_areas_with_barrier;
        n_barrier_edges += static_cast<int>(barriers.size());
      }
      if (!buildings.empty()) {
        ++n_areas_with_building;
        n_buildings_inside += static_cast<int>(buildings.size());
      }

      // A building overlapping the area that was never mapped as an inner ring
      // is still a building: routing has to go round it, so it becomes a hole.
      for (auto const& b : buildings) {
        rings.push_back(b);
        n_vertices += b.size();
      }
      for (auto const& b : barriers) {
        n_vertices += b.size();
      }

      if (n_vertices + connectors.size() > kMaxVertices) {
        ++n_too_big;
        big_vertex_counts.push_back(static_cast<double>(n_vertices));
        continue;
      }

      // Stage 1 (the binary): are any ways mapped into this area at all?
      auto const has_interior_ways = !linear.interior_edges_[area_i].empty();
      n_no_interior_ways += has_interior_ways ? 0 : 1;

      // What the router can already do here: the area's own ring ways, plus
      // any linear way lying inside it.
      auto const network = network_distances(a, linear.interior_edges_[area_i],
                                             linear.node_pos_, connector_ids);

      // Stage 2, cheap half. The geodesic is never shorter than the straight
      // line, so network/straight_line <= T implies network/geodesic <= T.
      // Where that holds for every pair, the area is provably already served
      // and the O(n^3) visibility graph never has to be built for it.
      auto cheap_served = true;
      for (auto i = std::size_t{0U}; i != connectors.size() && cheap_served;
           ++i) {
        for (auto j = i + 1U; j != connectors.size() && cheap_served; ++j) {
          if (on_ring[i] == 0 || on_ring[j] == 0) {
            continue;
          }
          auto const straight = geo::distance(connectors[i], connectors[j]);
          if (straight < 1.0) {
            continue;
          }
          cheap_served =
              network[i * connectors.size() + j] <= kServedRatio * straight;
        }
      }
      n_cheap_served += cheap_served ? 1 : 0;
      // An area with no ways mapped into it can still be adequately served by
      // its own ring - two adjacent entrances, or a thin shape. Sending it
      // straight to subdivision on the binary alone would be wasted work.
      n_no_ways_but_served += (!has_interior_ways && cheap_served) ? 1 : 0;

      auto const g = osr::area_geodesics{rings, connectors, barriers};

      auto all = std::vector<std::size_t>(connectors.size());
      for (auto i = std::size_t{0U}; i != all.size(); ++i) {
        all[i] = i;
      }
      auto const s = compute_pair_stats(g, all);
      if (s.n_pairs_ == 0U) {
        ++n_no_crossing;
        continue;
      }

      auto const arc = ring_arc_lengths(outer);
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

      // A pair the perimeter already serves at near-optimal cost does not
      // need the area edge at all - the ring ways stay in the graph either
      // way. Only pairs the crossing actually shortens have to be modelled,
      // and those are the ones the flat cost should be judged on.
      auto const n_c = connectors.size();
      auto relevant = std::vector<bool>(n_c * n_c, true);
      for (auto i = std::size_t{0U}; i != n_c; ++i) {
        for (auto j = i + 1U; j != n_c; ++j) {
          auto const d = g.distance(i, j);
          if (d == osr::area_geodesics::kUnreachable || d < 1.0 ||
              outer_pos[i] == std::numeric_limits<std::size_t>::max() ||
              outer_pos[j] == std::numeric_limits<std::size_t>::max()) {
            continue;
          }
          auto const ratio =
              perimeter_distance(arc, outer_pos[i], outer_pos[j]) / d;
          detours.push_back(ratio);
          if (ratio < kRelevantDetour) {
            relevant[i * n_c + j] = relevant[j * n_c + i] = false;
          }
        }
      }

      // How well the existing ways already serve the pairs that matter.
      //
      // Only boundary-to-boundary pairs can answer this: those are joined by
      // the ring whatever else is there, so a bad ratio means the interior is
      // genuinely unserved rather than merely outside the subgraph. An
      // interior stub's one edge leaves the area and is not in the subgraph at
      // all, so its distances would all be infinite and say nothing about the
      // area. Stubs are counted separately instead - one that cannot be
      // reached over existing ways is itself a reason to mesh.
      auto worst_network_ratio = 0.0;
      auto n_relevant = 0;
      for (auto i = std::size_t{0U}; i != n_c; ++i) {
        for (auto j = i + 1U; j != n_c; ++j) {
          auto const d = g.distance(i, j);
          if (!relevant[i * n_c + j] ||
              d == osr::area_geodesics::kUnreachable || d < 1.0) {
            continue;
          }
          if (on_ring[i] == 0 || on_ring[j] == 0) {
            continue;
          }
          ++n_relevant;
          auto const r = network[i * n_c + j] / static_cast<double>(d);
          network_ratios.push_back(std::min(r, 1000.0));
          worst_network_ratio = std::max(worst_network_ratio, r);
        }
      }

      auto stranded = 0;
      for (auto i = std::size_t{0U}; i != n_c; ++i) {
        if (on_ring[i] != 0) {
          continue;
        }
        auto reachable = false;
        for (auto j = 0U; j != n_c && !reachable; ++j) {
          reachable = i != j && std::isfinite(network[i * n_c + j]);
        }
        stranded += reachable ? 0 : 1;
      }
      n_stranded_connectors += stranded;
      n_areas_with_stranded += stranded != 0 ? 1 : 0;

      if (n_relevant != 0) {
        worst_network.push_back(std::min(worst_network_ratio, 1000.0));
      }
      // Stays index-aligned with spreads_relevant below.
      needs_mesh.push_back(
          (worst_network_ratio > kServedRatio || stranded != 0) ? 1 : 0);

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
            auto portals_f = std::vector<geo::latlng>{};
            auto root_f = cell_tree{};
            root_f.connectors_ = root.connectors_;
            split_smart(root_f, connectors, rings, portals_f, kAdaptiveMaxDepth,
                        g, relevant, n_c, 0.0);
            auto pts_f = connectors;
            pts_f.insert(end(pts_f), begin(portals_f), end(portals_f));
            if (!portals_f.empty() &&
                n_vertices + pts_f.size() <= kMaxVertices) {
              auto const gf = osr::area_geodesics{rings, pts_f, barriers};
              resolve_portals(root_f, gf);
              auto ok_f = true;
              check_portals(root_f, gf, ok_f);
              if (ok_f) {
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

    fmt::print(
        "areas with a crossing: {}   (skipped, <2 usable connectors: {})\n",
        n_reported, n_no_crossing);
    fmt::print("  level rule: {}\n",
               g_strict_levels ? "strict (a levelled area needs levelled "
                                 "connectors)"
                               : "permissive (untagged joins anything)");
    fmt::print("  level tags: {}/{} areas, {}/{} linear ways\n",
               collector.n_with_level_, collector.areas_.size(),
               linear.n_with_level_, linear.n_ways_);
    if (n_too_big != 0) {
      auto const p = quantiles_of(big_vertex_counts);
      fmt::print(
          "  declined, over {} vertices: {} areas "
          "(p50 {:.0f}, max {:.0f} vertices)\n",
          kMaxVertices, n_too_big, p.p50_, p.max_);
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
        n_restricted_connectors, linear.n_restricted_ways_);
    fmt::print(
        "  rejected by layer (bridge over / tunnel under, no shared "
        "node): {} interior edges, {} interior connectors\n",
        linear.n_rejected_edges_, n_rejected_interior);
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
