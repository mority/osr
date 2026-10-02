#include "osr/routing/reach_cells.h"

#include <cmath>
#include <concepts>
#include <algorithm>
#include <array>
#include <barrier>
#include <future>
#include <mutex>
#include <numbers>
#include <tuple>

#include "utl/enumerate.h"
#include "utl/helpers/algorithm.h"
#include "utl/verify.h"

#include "osr/lookup.h"
#include "osr/routing/dijkstra.h"
#include "osr/routing/with_profile.h"
#include "osr/ways.h"

namespace osr {

namespace {

// Besides the closest way, also start on ways that are at most this much
// further away (e.g. the other carriageway of a dual carriageway).
constexpr auto const kSeedDistanceSlack = 25.0;

constexpr auto const kMetersPerDegreeLat = 111'320.0;

// Profiles whose search state is keyed by `node_idx_t` with one slot per
// (way, direction): the owner of a node is the owner of its cheapest slot.
template <typename P>
concept SlotProfile =
    std::is_same_v<typename P::key, node_idx_t> &&
    requires(
        typename P::node const n, ways::routing const& r, node_idx_t const i) {
      {
        P::entry::get_node(i, std::size_t{0U})
      } -> std::same_as<typename P::node>;
      { P::entry::get_index(n) } -> std::convertible_to<std::size_t>;
      P::entry::storage_t::slot_count(r, i);
    };

template <Profile P>
std::uint64_t slot_key(typename P::node const n) {
  return (static_cast<std::uint64_t>(to_idx(n.get_node())) << 8U) |
         static_cast<std::uint64_t>(P::entry::get_index(n));
}

// Seeds `d` with the start candidates of `src`, calls `on_seed(node, cd)` for
// every start label.
template <Profile P, typename OnSeed>
void seed_source(dijkstra<P>& d,
                 typename P::parameters const& params,
                 ways const& w,
                 lookup const& l,
                 location const& src,
                 cost_t const max,
                 double const max_match_distance,
                 bitvec<node_idx_t> const* blocked,
                 OnSeed&& on_seed) {
  constexpr auto const kDir = direction::kForward;
  auto matches = match_result{};
  l.match<P>(params, src, false, kDir, max_match_distance, blocked, matches);
  auto const m = matches[match_idx_t{0U}];
  for (auto j = std::size_t{0U}; j != m.size(); ++j) {
    if (m.dist_to_way_[j] > m.dist_to_way_[0] + kSeedDistanceSlack) {
      break;
    }
    auto const way = m.way_[j];
    for (auto const& nc : {m.left(j), m.right(j)}) {
      if (!nc.valid() || nc.cost_ >= max) {
        continue;
      }
      auto const start = P::way_cost(
          params, *w.r_, w.timezones_, way, w.r_->way_properties_[way],
          flip(kDir, nc.way_dir_), static_cast<distance_t>(nc.dist_to_node_),
          std::nullopt, duration_t{0}, kDir);
      if (start.cost_ == kInfeasible || start.cost_ >= max) {
        continue;
      }
      P::resolve_start_node(
          *w.r_, way, nc.node_, src.lvl_, kDir, [&](auto const node) {
            auto label = typename P::label{node, start.cost_};
            label.track(label, *w.r_, way, node.get_node(), false);
            d.add_start(label, start.duration_);
            on_seed(node, start);
          });
    }
  }
}

// Cheapest slot of every reached node.
template <Profile P, typename Fn>
void for_each_reached_node(dijkstra<P> const& d, ways const& w, Fn&& fn) {
  for (auto const& [n, e] : d.cost_) {
    if (to_idx(n) >= w.n_nodes()) {
      continue;  // additional node
    }
    auto best = P::node::invalid();
    auto best_cost = kInfeasible;
    auto const n_slots = P::entry::storage_t::slot_count(*w.r_, n);
    for (auto i = std::size_t{0U}; i != n_slots; ++i) {
      auto const x = P::entry::get_node(n, i);
      if (auto const c = e.cost(x); c < best_cost) {
        best = x;
        best_cost = c;
      }
    }
    if (best_cost != kInfeasible) {
      fn(n, e, best, best_cost);
    }
  }
}

template <Profile P>
std::vector<reach_node> compute_owners(typename P::parameters const& params,
                                       ways const& w,
                                       lookup const& l,
                                       std::vector<location> const& sources,
                                       cost_t const max,
                                       double const max_match_distance,
                                       bitvec<node_idx_t> const* blocked) {
  constexpr auto const kDir = direction::kForward;

  struct seed {
    source_idx_t source_;
    cost_and_duration cd_;
  };

  auto d = dijkstra<P>{};
  d.reset({.profile_ = params,
           .w_ = &w,
           .max_ = max,
           .dir_ = kDir,
           .blocked_ = blocked});

  // Start slots that remain without predecessor belong to the source that
  // seeded them with the best (cost, duration), mirroring `entry::update`.
  auto seeds = hash_map<std::uint64_t, seed>{};
  for (auto const [i, src] : utl::enumerate(sources)) {
    auto const source = static_cast<source_idx_t>(i);
    seed_source<P>(
        d, params, w, l, src, max, max_match_distance, blocked,
        [&](typename P::node const node, cost_and_duration const start) {
          auto const [it, inserted] = seeds.emplace(
              slot_key<P>(node), seed{.source_ = source, .cd_ = start});
          if (!inserted && start < it->second.cd_) {
            it->second = seed{.source_ = source, .cd_ = start};
          }
        });
  }

  d.run();
  // The owner of a slot is the owner of the root of its predecessor chain.
  // Memoized per slot (not per node): the cheapest slot of a node can have a
  // different owner than the slot a successor was reached from.
  auto memo = hash_map<std::uint64_t, source_idx_t>{};
  auto stack = std::vector<std::uint64_t>{};
  auto const get_owner = [&](typename P::node x) {
    stack.clear();
    auto owner = kNoSource;
    while (true) {
      auto const key = slot_key<P>(x);
      if (auto const it = memo.find(key); it != end(memo)) {
        owner = it->second;
        break;
      }
      stack.push_back(key);
      auto const pred = d.cost_.at(x.get_key()).pred(x);
      if (!pred.has_value()) {
        auto const it = seeds.find(key);
        owner = it == end(seeds) ? kNoSource : it->second.source_;
        break;
      }
      x = *pred;
    }
    for (auto const key : stack) {
      memo.emplace(key, owner);
    }
    return owner;
  };

  auto result = std::vector<reach_node>{};
  result.reserve(d.cost_.size());
  for_each_reached_node(
      d, w,
      [&](node_idx_t const n, auto const& e, typename P::node const best,
          cost_t const best_cost) {
        auto const owner = get_owner(best);
        if (owner != kNoSource) {
          result.push_back({.node_ = n,
                            .cost_ = best_cost,
                            .duration_ = e.duration(best),
                            .source_ = owner});
        }
      });
  return result;
}

// Best settled cost per search state (node, way, direction) over all sources,
// shared by the per-source searches of `compute_owners_parallel`.
struct shared_best {
  // Returns false if another source has settled this state strictly cheaper:
  // everything reachable from it is then cheaper from that source as well.
  bool claim(std::uint64_t const key, cost_t const cost) {
    auto& s = shards_[ankerl::unordered_dense::hash<std::uint64_t>{}(key) %
                      shards_.size()];
    auto const lock = std::lock_guard{s.mutex_};
    auto const [it, inserted] = s.best_.emplace(key, cost);
    if (inserted || cost <= it->second) {
      it->second = cost;
      return true;
    }
    return false;
  }

  struct shard {
    std::mutex mutex_;
    hash_map<std::uint64_t, cost_t> best_;
  };
  std::vector<shard> shards_ = std::vector<shard>(1024U);
};

// Same result as `compute_owners` (up to ties), one thread per source:
// a search stops expanding a state another source has already settled
// cheaper, so each search mostly explores its own cell. The searches advance
// in lockstep rounds (by cost), so no search runs ahead into the cells of
// others before they had a chance to claim them.
template <Profile P>
std::vector<reach_node> compute_owners_parallel(
    typename P::parameters const& params,
    ways const& w,
    lookup const& l,
    std::vector<location> const& sources,
    cost_t const max,
    double const max_match_distance,
    bitvec<node_idx_t> const* blocked) {
  constexpr auto const kRound = cost_t{60U};

  auto best = shared_best{};
  auto round_end = kRound;  // only changed while all searches wait
  auto sync = std::barrier{static_cast<std::ptrdiff_t>(sources.size()),
                           [&]() noexcept { round_end += kRound; }};

  auto const search = [&](source_idx_t const source) {
    // Leave the lockstep in any case, otherwise the others wait forever.
    struct leave {
      ~leave() { sync_.arrive_and_drop(); }
      decltype(sync)& sync_;
    } const leave_sync{sync};

    auto d = dijkstra<P>{};
    d.reset({.profile_ = params,
             .w_ = &w,
             .max_ = max,
             .dir_ = direction::kForward,
             .blocked_ = blocked});
    seed_source<P>(d, params, w, l, sources[source], max, max_match_distance,
                   blocked, [](auto&&, auto&&) {});
    d.run([&](typename P::label const& x) {
      while (x.cost() >= round_end) {
        sync.arrive_and_wait();
      }
      return best.claim(slot_key<P>(x.get_node()), x.cost());
    });

    auto result = std::vector<reach_node>{};
    result.reserve(d.cost_.size());
    for_each_reached_node(d, w,
                          [&](node_idx_t const n, auto const& e,
                              typename P::node const x, cost_t const cost) {
                            result.push_back({.node_ = n,
                                              .cost_ = cost,
                                              .duration_ = e.duration(x),
                                              .source_ = source});
                          });
    return result;
  };

  auto searches = std::vector<std::future<std::vector<reach_node>>>{};
  for (auto i = std::size_t{0U}; i != sources.size(); ++i) {
    searches.push_back(
        std::async(std::launch::async, search, static_cast<source_idx_t>(i)));
  }

  // Owner of a node = source with the cheapest slot (ties: lower index).
  auto result = std::vector<reach_node>{};
  auto pos = hash_map<node_idx_t, std::size_t>{};
  for (auto& f : searches) {
    for (auto const& n : f.get()) {
      auto const [it, inserted] = pos.emplace(n.node_, result.size());
      if (inserted) {
        result.push_back(n);
      } else if (auto& r = result[it->second];
                 std::tie(n.cost_, n.source_) < std::tie(r.cost_, r.source_)) {
        r = n;
      }
    }
  }
  return result;
}

// Per-source searches for the meeting zones, in lockstep rounds (by cost).
// Tracks every node's latest arrival / total over the sources that reached it
// so far; once all sources reached a node, it bounds the best values. The
// searches stop when the round end passes best + margin + slack for every
// requested zone: all places of the zones are reached by all sources then.
template <Profile P>
std::vector<std::vector<reach_node>> meet_searches(
    typename P::parameters const& params,
    ways const& w,
    lookup const& l,
    std::vector<location> const& sources,
    cost_t const max,
    double const max_match_distance,
    bitvec<node_idx_t> const* blocked,
    meet_bounds const& bounds) {
  constexpr auto const kRound = cost_t{60U};
  constexpr auto const kUnknown = std::numeric_limits<std::uint32_t>::max();
  auto const k = sources.size();
  auto const bounded =
      bounds.latest_margin_.has_value() || bounds.total_margin_.has_value();

  struct acc {
    std::size_t n_;
    std::uint32_t latest_;
    std::uint32_t total_;
  };
  struct shard {
    std::mutex mutex_;
    hash_map<node_idx_t, acc> reached_;
  };
  auto shards = std::vector<shard>(bounded ? 1024U : 0U);
  auto best_mutex = std::mutex{};
  auto best_latest = kUnknown;
  auto best_total = kUnknown;

  // Called the first time source i settles node n (= its cheapest slot).
  auto const reach = [&](node_idx_t const n, std::uint32_t const duration) {
    auto& s = shards[ankerl::unordered_dense::hash<std::uint32_t>{}(to_idx(n)) %
                     shards.size()];
    auto const lock = std::lock_guard{s.mutex_};
    auto& a = s.reached_.emplace(n, acc{0U, 0U, 0U}).first->second;
    ++a.n_;
    a.latest_ = std::max(a.latest_, duration);
    a.total_ += duration;
    if (a.n_ == k) {
      auto const best_lock = std::lock_guard{best_mutex};
      best_latest = std::min(best_latest, a.latest_);
      best_total = std::min(best_total, a.total_);
    }
  };

  // Round bookkeeping, only changed while all searches wait in the barrier.
  auto round_end = kRound;
  auto settled_below = std::numeric_limits<cost_t>::max();  // set on stop
  auto stop = false;
  auto const radius_needed = [&]() -> std::uint64_t {
    auto needed = std::uint64_t{0U};
    auto const add = [&](std::optional<duration_t> const& margin,
                         std::uint32_t const best) {
      if (!margin.has_value()) {
        return;
      }
      needed = std::max(needed, best == kUnknown
                                    ? std::uint64_t{kUnknown}
                                    : std::uint64_t{best} + margin->count() +
                                          kMeetSearchSlack);
    };
    add(bounds.latest_margin_, best_latest);
    add(bounds.total_margin_, best_total);
    return needed;
  };
  auto sync = std::barrier{static_cast<std::ptrdiff_t>(k), [&]() noexcept {
                             // All labels with cost < round_end have been
                             // processed.
                             if (bounded && round_end > radius_needed()) {
                               stop = true;
                               settled_below = round_end;
                             }
                             round_end += kRound;
                           }};

  auto const search = [&](source_idx_t const source) {
    struct leave {
      ~leave() { sync_.arrive_and_drop(); }
      decltype(sync)& sync_;
    } const leave_sync{sync};

    auto d = dijkstra<P>{};
    d.reset({.profile_ = params,
             .w_ = &w,
             .max_ = max,
             .dir_ = direction::kForward,
             .blocked_ = blocked});
    seed_source<P>(d, params, w, l, sources[source], max, max_match_distance,
                   blocked, [](auto&&, auto&&) {});

    auto seen = hash_set<node_idx_t>{};
    d.run([&](typename P::label const& x) {
      while (x.cost() >= round_end) {
        sync.arrive_and_wait();
        if (stop) {
          d.pq_.clear();
          return false;
        }
      }
      auto const node = x.get_node();
      if (bounded && to_idx(node.get_node()) < w.n_nodes() &&
          seen.emplace(node.get_node()).second) {
        reach(node.get_node(),
              d.cost_.at(node.get_key()).duration(node).count());
      }
      return true;
    });

    // Beyond the stopping round, costs are not final (only too high).
    auto result = std::vector<reach_node>{};
    result.reserve(d.cost_.size());
    for_each_reached_node(d, w,
                          [&](node_idx_t const n, auto const& e,
                              typename P::node const x, cost_t const cost) {
                            if (cost < settled_below) {
                              result.push_back({.node_ = n,
                                                .cost_ = cost,
                                                .duration_ = e.duration(x),
                                                .source_ = 0U});
                            }
                          });
    return result;
  };

  auto searches = std::vector<std::future<std::vector<reach_node>>>{};
  for (auto i = std::size_t{0U}; i != k; ++i) {
    searches.push_back(
        std::async(std::launch::async, search, static_cast<source_idx_t>(i)));
  }
  auto per_source = std::vector<std::vector<reach_node>>{};
  for (auto& f : searches) {
    per_source.push_back(f.get());
  }
  return per_source;
}

// Edge directions, counter-clockwise: east, north, west, south.
constexpr auto const kDx = std::array<std::int32_t, 4U>{1, 0, -1, 0};
constexpr auto const kDy = std::array<std::int32_t, 4U>{0, 1, 0, -1};

struct grid_edge {
  grid_vertex from_;
  std::uint8_t dir_;
  bool used_{false};
};

std::uint64_t vertex_key(grid_vertex const v) {
  return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(v.x_)) << 32U) |
         static_cast<std::uint64_t>(static_cast<std::uint32_t>(v.y_));
}

// Twice the signed area: > 0 counter-clockwise.
std::int64_t signed_area2(std::vector<grid_vertex> const& ring) {
  auto a = std::int64_t{0};
  for (auto i = std::size_t{0U}; i != ring.size(); ++i) {
    auto const& p = ring[i];
    auto const& q = ring[(i + 1U) % ring.size()];
    a += static_cast<std::int64_t>(p.x_) * q.y_ -
         static_cast<std::int64_t>(q.x_) * p.y_;
  }
  return a;
}

// Is the center of cell (cx, cy) inside the ring? Works on doubled
// coordinates: the cell center is odd, ring vertices are even, so the point
// is never on the ring.
bool contains_cell(std::vector<grid_vertex> const& ring,
                   std::int32_t const cx,
                   std::int32_t const cy) {
  auto const px = 2 * static_cast<std::int64_t>(cx) + 1;
  auto const py = 2 * static_cast<std::int64_t>(cy) + 1;
  auto inside = false;
  for (auto i = std::size_t{0U}, j = ring.size() - 1U; i != ring.size();
       j = i++) {
    auto const xi = 2 * static_cast<std::int64_t>(ring[i].x_);
    auto const yi = 2 * static_cast<std::int64_t>(ring[i].y_);
    auto const xj = 2 * static_cast<std::int64_t>(ring[j].x_);
    auto const yj = 2 * static_cast<std::int64_t>(ring[j].y_);
    if ((yi > py) != (yj > py)) {
      // x of the crossing, compared without division
      auto const lhs = (px - xi) * (yj - yi);
      auto const rhs = (xj - xi) * (py - yi);
      if ((yj > yi) ? lhs < rhs : lhs > rhs) {
        inside = !inside;
      }
    }
  }
  return inside;
}

template <typename IsPinned>
std::vector<grid_polygon> trace_source(std::vector<grid_edge>& edges,
                                       IsPinned&& is_pinned) {
  auto out = hash_map<std::uint64_t, std::array<std::uint32_t, 2U>>{};
  constexpr auto const kNone = std::numeric_limits<std::uint32_t>::max();
  for (auto i = 0U; i != edges.size(); ++i) {
    auto& slots = out.emplace(vertex_key(edges[i].from_),
                              std::array<std::uint32_t, 2U>{kNone, kNone})
                      .first->second;
    slots[slots[0] == kNone ? 0U : 1U] = i;
  }

  struct ring {
    std::vector<grid_vertex> vertices_;
    grid_edge first_;  // unsimplified, for the hole -> outer ring test
    std::int64_t area2_;
  };

  auto rings = std::vector<ring>{};
  auto const emit = [&](std::span<grid_edge const> loop) {
    // Keep corners and pinned vertices only.
    auto r = ring{.vertices_ = {}, .first_ = loop.front(), .area2_ = 0};
    for (auto i = std::size_t{0U}; i != loop.size(); ++i) {
      auto v = loop[i].from_;
      v.pinned_ = is_pinned(v);
      if (v.pinned_ ||
          loop[(i + loop.size() - 1U) % loop.size()].dir_ != loop[i].dir_) {
        r.vertices_.push_back(v);
      }
    }
    r.area2_ = signed_area2(r.vertices_);
    rings.push_back(std::move(r));
  };

  auto path = std::vector<grid_edge>{};
  auto path_pos = hash_map<std::uint64_t, std::size_t>{};
  for (auto start = 0U; start != edges.size(); ++start) {
    if (edges[start].used_) {
      continue;
    }

    path.clear();
    path_pos.clear();
    auto curr = start;
    while (true) {
      auto& e = edges[curr];
      e.used_ = true;

      // Where four cells around a vertex alternate between this source and
      // others, the boundary touches itself. Cut the loop since the last
      // visit off as a ring of its own, so all rings stay simple (the cut off
      // part becomes a polygon or a hole touching another ring at a point).
      auto const key = vertex_key(e.from_);
      if (auto const it = path_pos.find(key); it != end(path_pos)) {
        auto const k = it->second;
        emit(std::span{path}.subspan(k));
        for (auto i = k; i != path.size(); ++i) {
          path_pos.erase(vertex_key(path[i].from_));
        }
        path.resize(k);
      }
      path_pos.emplace(key, path.size());
      path.push_back(e);

      // At a vertex where two cells of this source only touch diagonally,
      // there are two outgoing edges. Turning left keeps the cells apart.
      auto const to =
          grid_vertex{e.from_.x_ + kDx[e.dir_], e.from_.y_ + kDy[e.dir_]};
      auto const& candidates = out.at(vertex_key(to));
      auto next = kNone;
      for (auto const turn : {1U, 0U, 3U}) {
        auto const want = static_cast<std::uint8_t>((e.dir_ + turn) % 4U);
        for (auto const c : candidates) {
          if (c != kNone && edges[c].dir_ == want) {
            next = c;
            break;
          }
        }
        if (next != kNone) {
          break;
        }
      }
      utl::verify(next != kNone, "trace_reach_grid: open ring");
      if (edges[next].used_) {
        utl::verify(next == start, "trace_reach_grid: ring not closed");
        break;
      }
      curr = next;
    }
    emit(path);
  }

  auto polygons = std::vector<grid_polygon>{};
  auto outer_idx = std::vector<std::size_t>{};
  for (auto const& r : rings) {
    if (r.area2_ > 0) {
      outer_idx.push_back(polygons.size());
      polygons.push_back(grid_polygon{r.vertices_});
    }
  }
  for (auto const& r : rings) {
    if (r.area2_ > 0) {
      continue;
    }
    // The owned cell left of the hole's first edge lies in the polygon the
    // hole belongs to: the smallest outer ring containing it.
    auto const [v, dir, _] = r.first_;
    auto const cx = v.x_ + (dir == 1U || dir == 2U ? -1 : 0);
    auto const cy = v.y_ + (dir == 2U || dir == 3U ? -1 : 0);
    auto best = std::size_t{0U};
    auto best_area = std::numeric_limits<std::int64_t>::max();
    for (auto const i : outer_idx) {
      auto const area = signed_area2(polygons[i].front());
      if (area < best_area && contains_cell(polygons[i].front(), cx, cy)) {
        best = i;
        best_area = area;
      }
    }
    utl::verify(best_area != std::numeric_limits<std::int64_t>::max(),
                "trace_reach_grid: hole without outer ring");
    polygons[best].push_back(r.vertices_);
  }
  return polygons;
}

using label_t = std::uint32_t;
constexpr auto const kNoLabel = std::numeric_limits<label_t>::max();

// A vertex is pinned (kept in place by smoothing) where three or more labels
// meet, or where two labels touch diagonally. Everywhere else, the boundary
// separates exactly two labels, and both rings smooth it the same way.
template <typename GetLabel>
bool is_junction(GetLabel&& get_label, grid_vertex const v) {
  auto const sw = get_label(v.x_ - 1, v.y_ - 1);
  auto const se = get_label(v.x_, v.y_ - 1);
  auto const nw = get_label(v.x_ - 1, v.y_);
  auto const ne = get_label(v.x_, v.y_);
  auto distinct = std::array<label_t, 4U>{sw, se, nw, ne};
  utl::sort(distinct);
  auto const n = static_cast<std::size_t>(std::distance(
      begin(distinct), std::unique(begin(distinct), end(distinct))));
  return n >= 3U || (n == 2U && sw == ne && se == nw);
}

// Traces the cells of every label. `get_label(x, y)` must return kNoLabel
// outside of the grid. Junctions are determined by `get_pin_label`, which
// must be at least as fine as `get_label` (e.g. bands for source outlines, so
// outlines and bands are smoothed the same way).
template <typename GetLabel, typename GetPinLabel>
std::vector<std::vector<grid_polygon>> trace_labels(
    reach_grid const& g,
    GetLabel&& get_label,
    std::size_t const n_labels,
    GetPinLabel&& get_pin_label) {
  // Boundary edges of each label, oriented counter-clockwise around the
  // cells (cell on the left).
  auto edges = std::vector<std::vector<grid_edge>>(n_labels);
  for (auto y = 0; y != static_cast<std::int32_t>(g.ny_); ++y) {
    for (auto x = 0; x != static_cast<std::int32_t>(g.nx_); ++x) {
      auto const o = get_label(x, y);
      if (o == kNoLabel) {
        continue;
      }
      auto& e = edges[o];
      if (get_label(x, y - 1) != o) {
        e.push_back({.from_ = {x, y}, .dir_ = 0U});
      }
      if (get_label(x + 1, y) != o) {
        e.push_back({.from_ = {x + 1, y}, .dir_ = 1U});
      }
      if (get_label(x, y + 1) != o) {
        e.push_back({.from_ = {x + 1, y + 1}, .dir_ = 2U});
      }
      if (get_label(x - 1, y) != o) {
        e.push_back({.from_ = {x, y + 1}, .dir_ = 3U});
      }
    }
  }

  auto result = std::vector<std::vector<grid_polygon>>(n_labels);
  for (auto i = std::size_t{0U}; i != n_labels; ++i) {
    result[i] = trace_source(edges[i], [&](grid_vertex const v) {
      return is_junction(get_pin_label, v);
    });
  }
  return result;
}

}  // namespace

std::vector<reach_node> compute_reach_owners(
    profile_parameters const& params,
    ways const& w,
    lookup const& l,
    search_profile const profile,
    std::vector<location> const& sources,
    cost_t const max,
    double const max_match_distance,
    bitvec<node_idx_t> const* blocked) {
  utl::verify(
      profile == search_profile::kCar || profile == search_profile::kBus ||
          profile == search_profile::kHgv,
      "compute_reach_owners: profile {} not supported", to_str(profile));
  utl::verify(sources.size() < kNoSource,
              "compute_reach_owners: too many sources");
  return with_profile(profile, [&]<Profile P>(P&&) -> std::vector<reach_node> {
    if constexpr (SlotProfile<P>) {
      auto const& pp = std::get<typename P::parameters>(params);
      return sources.size() > 1U
                 ? compute_owners_parallel<P>(pp, w, l, sources, max,
                                              max_match_distance, blocked)
                 : compute_owners<P>(pp, w, l, sources, max, max_match_distance,
                                     blocked);
    } else {
      throw utl::fail("compute_reach_owners: profile {} not supported",
                      to_str(profile));
    }
  });
}

meet_times compute_meet_times(profile_parameters const& params,
                              ways const& w,
                              lookup const& l,
                              search_profile const profile,
                              std::vector<location> const& sources,
                              cost_t const max,
                              double const max_match_distance,
                              bitvec<node_idx_t> const* blocked,
                              meet_bounds const& bounds) {
  auto result = meet_times{};
  if (sources.empty()) {
    return result;
  }
  utl::verify(profile == search_profile::kCar ||
                  profile == search_profile::kBus ||
                  profile == search_profile::kHgv,
              "compute_meet_times: profile {} not supported", to_str(profile));

  // Per node: number of sources that reached it so far, latest arrival,
  // total drive time.
  struct acc {
    std::uint16_t n_;
    std::uint32_t latest_;
    std::uint32_t total_;
  };
  auto reached = hash_map<node_idx_t, acc>{};
  auto const per_source = with_profile(
      profile, [&]<Profile P>(P&&) -> std::vector<std::vector<reach_node>> {
        if constexpr (SlotProfile<P>) {
          return meet_searches<P>(std::get<typename P::parameters>(params), w,
                                  l, sources, max, max_match_distance, blocked,
                                  bounds);
        } else {
          throw utl::fail("compute_meet_times: profile {} not supported",
                          to_str(profile));
        }
      });

  for (auto const [i, nodes] : utl::enumerate(per_source)) {
    for (auto const& n : nodes) {
      auto const d = static_cast<std::uint32_t>(n.duration_.count());
      if (i == 0U) {
        reached.emplace(n.node_, acc{.n_ = 1U, .latest_ = d, .total_ = d});
      } else if (auto const it = reached.find(n.node_);
                 it != end(reached) && it->second.n_ == i) {
        ++it->second.n_;
        it->second.latest_ = std::max(it->second.latest_, d);
        it->second.total_ += d;
      }
    }
  }

  // Ties: the other criterion, then node index.
  auto best_latest = std::optional<meet_node>{};
  auto best_total = std::optional<meet_node>{};
  for (auto const& [n, a] : reached) {
    if (a.n_ != sources.size()) {
      continue;
    }
    auto const m =
        meet_node{.node_ = n, .latest_ = a.latest_, .total_ = a.total_};
    result.nodes_.push_back(m);
    if (!best_latest.has_value() ||
        std::tie(m.latest_, m.total_, m.node_) < std::tie(best_latest->latest_,
                                                          best_latest->total_,
                                                          best_latest->node_)) {
      best_latest = m;
    }
    if (!best_total.has_value() ||
        std::tie(m.total_, m.latest_, m.node_) < std::tie(best_total->total_,
                                                          best_total->latest_,
                                                          best_total->node_)) {
      best_total = m;
    }
  }

  auto const with_per_source = [&](meet_node const& m) {
    auto b = meet_times::best{.node_ = m, .per_source_ = {}};
    for (auto const& nodes : per_source) {
      auto const it = utl::find_if(
          nodes, [&](reach_node const& n) { return n.node_ == m.node_; });
      b.per_source_.push_back(it == end(nodes) ? kMaxDuration : it->duration_);
    }
    return b;
  };
  if (best_latest.has_value()) {
    result.latest_ = with_per_source(*best_latest);
    result.total_ = with_per_source(*best_total);
  }
  return result;
}

reach_grid make_reach_grid(ways const& w,
                           std::span<reach_node const> nodes,
                           std::span<meet_node const> meet,
                           double const cell_size_m,
                           unsigned const fill_steps) {
  auto g = reach_grid{};
  if (nodes.empty()) {
    return g;
  }

  auto min = geo::latlng{90.0, 180.0};
  auto max = geo::latlng{-90.0, -180.0};
  for (auto const& n : nodes) {
    auto const p = w.get_node_pos(n.node_).as_latlng();
    min.lat_ = std::min(min.lat_, p.lat_);
    min.lng_ = std::min(min.lng_, p.lng_);
    max.lat_ = std::max(max.lat_, p.lat_);
    max.lng_ = std::max(max.lng_, p.lng_);
  }

  auto const mid_lat = (min.lat_ + max.lat_) / 2.0;
  auto const cos_lat =
      std::max(0.01, std::cos(mid_lat * std::numbers::pi / 180.0));
  g.dlat_ = cell_size_m / kMetersPerDegreeLat;
  g.dlng_ = cell_size_m / (kMetersPerDegreeLat * cos_lat);

  auto const pad = static_cast<double>(fill_steps + 1U);
  g.origin_ = {min.lat_ - pad * g.dlat_, min.lng_ - pad * g.dlng_};
  g.nx_ = static_cast<std::uint32_t>(
      std::floor((max.lng_ - min.lng_) / g.dlng_) + 1.0 + 2.0 * pad);
  g.ny_ = static_cast<std::uint32_t>(
      std::floor((max.lat_ - min.lat_) / g.dlat_) + 1.0 + 2.0 * pad);
  g.owner_.assign(static_cast<std::size_t>(g.nx_) * g.ny_, kNoSource);
  g.cost_.assign(g.owner_.size(), kInfeasible);
  g.duration_.assign(g.owner_.size(), kMaxDuration);

  for (auto const& n : nodes) {
    auto const p = w.get_node_pos(n.node_).as_latlng();
    auto const x =
        static_cast<std::int32_t>((p.lng_ - g.origin_.lng_) / g.dlng_);
    auto const y =
        static_cast<std::int32_t>((p.lat_ - g.origin_.lat_) / g.dlat_);
    auto const i = g.idx(x, y);
    if (std::tie(n.cost_, n.source_) < std::tie(g.cost_[i], g.owner_[i])) {
      g.cost_[i] = n.cost_;
      g.duration_[i] = n.duration_;
      g.owner_[i] = n.source_;
    }
  }

  if (!meet.empty()) {
    constexpr auto const kNone = std::numeric_limits<std::uint32_t>::max();
    g.meet_latest_.assign(g.owner_.size(), kNone);
    g.meet_total_.assign(g.owner_.size(), kNone);
    for (auto const& m : meet) {
      auto const p = w.get_node_pos(m.node_).as_latlng();
      auto const i =
          g.idx(static_cast<std::int32_t>((p.lng_ - g.origin_.lng_) / g.dlng_),
                static_cast<std::int32_t>((p.lat_ - g.origin_.lat_) / g.dlat_));
      g.meet_latest_[i] = std::min(g.meet_latest_[i], m.latest_);
      g.meet_total_[i] = std::min(g.meet_total_[i], m.total_);
    }
  }

  fill_gaps(g, fill_steps);
  close_holes(g);
  return g;
}

void fill_gaps(reach_grid& g, unsigned const steps) {
  for (auto step = 0U; step != steps; ++step) {
    auto owner = g.owner_;
    auto cost = g.cost_;
    auto duration = g.duration_;
    auto meet_latest = g.meet_latest_;
    auto meet_total = g.meet_total_;
    auto changed = false;
    for (auto y = 0; y != static_cast<std::int32_t>(g.ny_); ++y) {
      for (auto x = 0; x != static_cast<std::int32_t>(g.nx_); ++x) {
        auto const i = g.idx(x, y);
        if (g.owner_[i] != kNoSource) {
          continue;
        }
        for (auto dy = -1; dy <= 1; ++dy) {
          for (auto dx = -1; dx <= 1; ++dx) {
            if (!g.contains(x + dx, y + dy)) {
              continue;
            }
            auto const j = g.idx(x + dx, y + dy);
            if (g.owner_[j] != kNoSource && std::tie(g.cost_[j], g.owner_[j]) <
                                                std::tie(cost[i], owner[i])) {
              cost[i] = g.cost_[j];
              duration[i] = g.duration_[j];
              owner[i] = g.owner_[j];
              if (!meet_latest.empty()) {
                meet_latest[i] = g.meet_latest_[j];
                meet_total[i] = g.meet_total_[j];
              }
              changed = true;
            }
          }
        }
      }
    }
    g.owner_ = std::move(owner);
    g.cost_ = std::move(cost);
    g.duration_ = std::move(duration);
    g.meet_latest_ = std::move(meet_latest);
    g.meet_total_ = std::move(meet_total);
    if (!changed) {
      break;
    }
  }
}

namespace {

auto source_label(reach_grid const& g, std::size_t const n_sources) {
  return [&g, n_sources](std::int32_t const x, std::int32_t const y) {
    auto const o = g.owner(x, y);
    return o == kNoSource || o >= n_sources ? kNoLabel
                                            : static_cast<label_t>(o);
  };
}

auto band_label(reach_grid const& g,
                std::size_t const n_sources,
                duration_t const band_size,
                std::size_t const n_bands) {
  return [&g, n_sources, band_size, n_bands](std::int32_t const x,
                                             std::int32_t const y) {
    auto const o = g.owner(x, y);
    if (o == kNoSource || o >= n_sources) {
      return kNoLabel;
    }
    auto const band =
        std::min(static_cast<std::size_t>(g.duration_[g.idx(x, y)] / band_size),
                 n_bands - 1U);
    return static_cast<label_t>(o * n_bands + band);
  };
}

auto meet_label(reach_grid const& g,
                meet_criterion const c,
                std::uint32_t const limit) {
  auto const& values =
      c == meet_criterion::kLatest ? g.meet_latest_ : g.meet_total_;
  return [&g, &values, limit](std::int32_t const x, std::int32_t const y) {
    return g.contains(x, y) && values[g.idx(x, y)] <= limit ? label_t{0U}
                                                            : kNoLabel;
  };
}

}  // namespace

void close_holes(reach_grid& g) {
  auto const n = g.owner_.size();
  if (n == 0U) {
    return;
  }

  // Empty cells connected to the grid border are outside.
  auto outside = std::vector<bool>(n, false);
  auto queue = std::vector<std::pair<std::int32_t, std::int32_t>>{};
  auto const visit_outside = [&](std::int32_t const x, std::int32_t const y) {
    if (g.contains(x, y) && g.owner_[g.idx(x, y)] == kNoSource &&
        !outside[g.idx(x, y)]) {
      outside[g.idx(x, y)] = true;
      queue.emplace_back(x, y);
    }
  };
  auto const nx = static_cast<std::int32_t>(g.nx_);
  auto const ny = static_cast<std::int32_t>(g.ny_);
  for (auto x = 0; x != nx; ++x) {
    visit_outside(x, 0);
    visit_outside(x, ny - 1);
  }
  for (auto y = 0; y != ny; ++y) {
    visit_outside(0, y);
    visit_outside(nx - 1, y);
  }
  for (auto i = std::size_t{0U}; i != queue.size(); ++i) {
    auto const [x, y] = queue[i];
    visit_outside(x + 1, y);
    visit_outside(x - 1, y);
    visit_outside(x, y + 1);
    visit_outside(x, y - 1);
  }

  // All other empty cells are enclosed: fill them from the nearest owned
  // cell (breadth first, 8-neighbourhood).
  queue.clear();
  for (auto y = 0; y != ny; ++y) {
    for (auto x = 0; x != nx; ++x) {
      if (g.owner_[g.idx(x, y)] != kNoSource) {
        queue.emplace_back(x, y);
      }
    }
  }
  for (auto i = std::size_t{0U}; i != queue.size(); ++i) {
    auto const [x, y] = queue[i];
    auto const from = g.idx(x, y);
    for (auto dy = -1; dy <= 1; ++dy) {
      for (auto dx = -1; dx <= 1; ++dx) {
        if (!g.contains(x + dx, y + dy)) {
          continue;
        }
        auto const to = g.idx(x + dx, y + dy);
        if (g.owner_[to] != kNoSource || outside[to]) {
          continue;
        }
        g.owner_[to] = g.owner_[from];
        g.cost_[to] = g.cost_[from];
        g.duration_[to] = g.duration_[from];
        if (!g.meet_latest_.empty()) {
          g.meet_latest_[to] = g.meet_latest_[from];
          g.meet_total_[to] = g.meet_total_[from];
        }
        queue.emplace_back(x + dx, y + dy);
      }
    }
  }
}

std::vector<std::vector<grid_polygon>> trace_reach_grid(
    reach_grid const& g, std::size_t const n_sources) {
  auto const label = source_label(g, n_sources);
  return trace_labels(g, label, n_sources, label);
}

std::vector<std::vector<grid_polygon>> trace_reach_bands(
    reach_grid const& g,
    std::size_t const n_sources,
    duration_t const band_size,
    std::size_t const n_bands) {
  utl::verify(band_size.count() != 0U && n_bands != 0U,
              "trace_reach_bands: no bands");
  auto const label = band_label(g, n_sources, band_size, n_bands);
  return trace_labels(g, label, n_sources * n_bands, label);
}

std::vector<grid_polygon> trace_meet_zone(reach_grid const& g,
                                          meet_criterion const c,
                                          std::uint32_t const limit) {
  if (g.meet_latest_.empty()) {
    return {};
  }
  auto const label = meet_label(g, c, limit);
  return trace_labels(g, label, 1U, label).front();
}

reach_polygon to_latlng(reach_grid const& g,
                        grid_polygon const& p,
                        unsigned const smooth_iterations) {
  struct pt {
    double x_, y_;
    bool pinned_;
  };

  auto result = reach_polygon{};
  result.reserve(p.size());
  auto curr = std::vector<pt>{};
  auto next = std::vector<pt>{};
  for (auto const& ring : p) {
    curr.clear();
    for (auto const& v : ring) {
      curr.push_back(
          {static_cast<double>(v.x_), static_cast<double>(v.y_), v.pinned_});
    }

    // Chaikin corner cutting: every free vertex is replaced by the points
    // 1/4 towards its neighbours. These only depend on the vertex and its two
    // neighbours, so a border shared by two rings (traversed in opposite
    // directions, with pinned ends) is smoothed identically for both.
    for (auto it = 0U; it != smooth_iterations; ++it) {
      next.clear();
      for (auto i = std::size_t{0U}; i != curr.size(); ++i) {
        auto const& v = curr[i];
        if (v.pinned_) {
          next.push_back(v);
          continue;
        }
        auto const& prev = curr[(i + curr.size() - 1U) % curr.size()];
        auto const& succ = curr[(i + 1U) % curr.size()];
        next.push_back({v.x_ + (prev.x_ - v.x_) / 4.0,
                        v.y_ + (prev.y_ - v.y_) / 4.0, false});
        next.push_back({v.x_ + (succ.x_ - v.x_) / 4.0,
                        v.y_ + (succ.y_ - v.y_) / 4.0, false});
      }
      std::swap(curr, next);
    }

    auto& r = result.emplace_back();
    r.reserve(curr.size() + 1U);
    for (auto const& v : curr) {
      r.push_back(
          {g.origin_.lat_ + v.y_ * g.dlat_, g.origin_.lng_ + v.x_ * g.dlng_});
    }
    r.push_back(r.front());
  }
  return result;
}

reach_cells compute_reach_cells(profile_parameters const& params,
                                ways const& w,
                                lookup const& l,
                                search_profile const profile,
                                std::vector<location> const& sources,
                                reach_cells_params const& p,
                                double const max_match_distance,
                                bitvec<node_idx_t> const* blocked) {
  // The per-source searches for the meeting zone run in parallel to the
  // combined search.
  auto meet_search = p.meet_margin_.has_value() || p.total_margin_.has_value()
                         ? std::async(std::launch::async,
                                      [&]() {
                                        return compute_meet_times(
                                            params, w, l, profile, sources,
                                            p.max_, max_match_distance, blocked,
                                            {.latest_margin_ = p.meet_margin_,
                                             .total_margin_ = p.total_margin_});
                                      })
                         : std::future<meet_times>{};
  auto const nodes = compute_reach_owners(params, w, l, profile, sources,
                                          p.max_, max_match_distance, blocked);
  auto const meet = meet_search.valid() ? meet_search.get() : meet_times{};
  auto const g =
      make_reach_grid(w, nodes, meet.nodes_, p.cell_size_m_, p.fill_steps_);

  auto const n_bands =
      p.band_size_.count() == 0U
          ? std::size_t{0U}
          : std::max(std::size_t{1U}, static_cast<std::size_t>(
                                          (p.max_ + p.band_size_.count() - 1U) /
                                          p.band_size_.count()));

  auto result = reach_cells{};
  result.sources_.resize(sources.size());
  auto const outlines = [&]() {
    auto const label = source_label(g, sources.size());
    // With bands, outlines are pinned at band junctions as well, so they are
    // smoothed exactly like the outer border of the bands.
    return n_bands == 0U ? trace_labels(g, label, sources.size(), label)
                         : trace_labels(g, label, sources.size(),
                                        band_label(g, sources.size(),
                                                   p.band_size_, n_bands));
  }();
  for (auto const [i, polygons] : utl::enumerate(outlines)) {
    for (auto const& poly : polygons) {
      result.sources_[i].push_back(to_latlng(g, poly, p.smooth_iterations_));
    }
  }

  if (n_bands != 0U) {
    for (auto const [label, polygons] : utl::enumerate(
             trace_reach_bands(g, sources.size(), p.band_size_, n_bands))) {
      for (auto const& poly : polygons) {
        result.bands_.push_back(
            {.source_ = static_cast<source_idx_t>(label / n_bands),
             .band_ = static_cast<std::uint16_t>(label % n_bands),
             .polygon_ = to_latlng(g, poly, p.smooth_iterations_)});
      }
    }
  }

  auto const add_meet = [&](std::optional<meet_times::best> const& best,
                            std::optional<duration_t> const& margin,
                            meet_criterion const c) {
    auto out = std::optional<reach_cells::meet>{};
    if (!best.has_value() || !margin.has_value()) {
      return out;
    }
    auto& m = out.emplace(
        reach_cells::meet{.pos_ = w.get_node_pos(best->node_.node_).as_latlng(),
                          .latest_ = best->node_.latest_,
                          .total_ = best->node_.total_,
                          .per_source_ = best->per_source_,
                          .zone_ = {}});
    auto const limit = best->node_.get(c) + margin->count();
    for (auto const& poly : trace_meet_zone(g, c, limit)) {
      m.zone_.push_back(to_latlng(g, poly, p.smooth_iterations_));
    }
    return out;
  };
  result.meet_ =
      add_meet(meet.latest_, p.meet_margin_, meet_criterion::kLatest);
  result.total_meet_ =
      add_meet(meet.total_, p.total_margin_, meet_criterion::kTotal);
  return result;
}

}  // namespace osr
