#include "osr/routing/route.h"

#include <cstdint>

#include <algorithm>
#include <optional>

#include <map>
#include <memory>
#include <mutex>
#include <shared_mutex>

#include "boost/thread/tss.hpp"

#include "utl/concat.h"
#include "utl/enumerate.h"
#include "utl/helpers/algorithm.h"
#include "utl/to_vec.h"
#include "utl/verify.h"

#include "osr/elevation_storage.h"
#include "osr/lookup.h"
#include "osr/routing/astar.h"
#include "osr/routing/bidirectional.h"
#include "osr/routing/cch/cch.h"
#include "osr/routing/cch/customize.h"
#include "osr/routing/cch/query.h"
#include "osr/routing/cch/rphast.h"
#include "osr/routing/cch/unpack.h"
#include "osr/routing/dijkstra.h"
#include "osr/routing/path_reconstruction.h"
#include "osr/routing/profiles/bike.h"
#include "osr/routing/profiles/bike_sharing.h"
#include "osr/routing/profiles/car.h"
#include "osr/routing/profiles/car_parking.h"
#include "osr/routing/profiles/car_sharing.h"
#include "osr/routing/profiles/foot.h"
#include "osr/routing/sharing_data.h"
#include "osr/routing/with_profile.h"
#include "osr/util/infinite.h"
#include "osr/util/reverse.h"

namespace osr {

constexpr auto const kMaxMatchingDistanceSquaredRatio = 9.0;
constexpr auto const kBottomKDefinitelyConsidered = 5;
constexpr auto const kMinCostSettled = cost_t{900};

template <Profile P>
bidirectional<P>& get_bidirectional() {
  static auto s = boost::thread_specific_ptr<bidirectional<P>>{};
  if (s.get() == nullptr) {
    s.reset(new bidirectional<P>{});
  }
  return *s.get();
}

template <Profile P>
dijkstra<P>& get_dijkstra() {
  static auto s = boost::thread_specific_ptr<dijkstra<P>>{};
  if (s.get() == nullptr) {
    s.reset(new dijkstra<P>{});
  }
  return *s.get();
}

template <Profile P>
astar<P>& get_astar() {
  static auto s = boost::thread_specific_ptr<astar<P>>{};
  if (s.get() == nullptr) {
    s.reset(new astar<P>{});
  }
  return *s.get();
}

duration_t sum_segment_durations(std::vector<path::segment> const& segments,
                                 duration_t total = duration_t{0}) {
  for (auto const& segment : segments) {
    total = clamp_add_duration(total, segment.duration_);
  }
  return total;
}

routing_algorithm to_algorithm(std::string_view s) {
  switch (cista::hash(s)) {
    case cista::hash("dijkstra"): return routing_algorithm::kDijkstra;
    case cista::hash("bidirectional"): return routing_algorithm::kAStarBi;
    case cista::hash("cch"): return routing_algorithm::kCCH;
  }
  throw utl::fail("unknown routing algorithm: {}", s);
}

template <Profile P>
path reconstruct_bi(typename P::parameters const& params,
                    ways const& w,
                    lookup const& l,
                    bitvec<node_idx_t> const* blocked,
                    sharing_data const* sharing,
                    elevation_storage const* elevations,
                    bidirectional<P> const& b,
                    location const& from,
                    location const& to,
                    way_idx_t const start_way,
                    candidate_node const& start_left,
                    candidate_node const& start_right,
                    way_idx_t const dest_way,
                    candidate_node const& dest_left,
                    candidate_node const& dest_right,
                    cost_t const cost,
                    direction const dir) {
  auto forward_n = b.meet_point_1_;

  // TODO subtract meetpoint node cost?

  auto forward_segments = std::vector<path::segment>{};
  auto forward_dist = 0.0;

  while (true) {
    auto const& e = b.cost1_.at(forward_n.get_key());
    auto const pred = e.pred(forward_n);
    if (pred.has_value()) {
      auto const pred_duration = b.cost1_.at(pred->get_key()).duration(*pred);
      auto const expected_cost = static_cast<cost_t>(
          e.cost(forward_n) - b.template get_cost<direction::kForward>(*pred));
      forward_dist +=
          add_path<P>(params, w, *w.r_, blocked, sharing, elevations, *pred,
                      forward_n, pred_duration, {}, expected_cost,
                      clamp_sub_duration(e.duration(forward_n), pred_duration),
                      forward_segments, dir);
    } else {
      break;
    }
    forward_n = *pred;
  }

  auto const& start_node_candidate =
      forward_n.get_node() == start_left.node_ ? start_left : start_right;

  forward_segments.push_back(
      {.polyline_ = l.get_node_candidate_path(
           start_way, start_node_candidate.node_, start_node_candidate.way_dir_,
           false, from),
       .from_level_ = start_node_candidate.lvl_,
       .to_level_ = start_node_candidate.lvl_,
       .from_ = dir == direction::kBackward ? forward_n.get_node()
                                            : node_idx_t::invalid(),
       .to_ = dir == direction::kForward ? forward_n.get_node()
                                         : node_idx_t::invalid(),

       .way_ = way_idx_t::invalid(),
       .cost_ = start_node_candidate.cost_,
       .duration_ = duration_from_cost(start_node_candidate.cost_),
       .dist_ = static_cast<distance_t>(start_node_candidate.dist_to_node_),
       .mode_ = forward_n.get_mode()});

  auto backward_segments = std::vector<path::segment>{};
  auto backward_n = b.meet_point_2_;
  auto backward_dist = 0.0;

  while (true) {
    auto const& e = b.cost2_.at(backward_n.get_key());
    auto const pred = e.pred(backward_n);
    if (pred.has_value()) {
      auto const expected_cost =
          static_cast<cost_t>(e.cost(backward_n) -
                              b.template get_cost<direction::kBackward>(*pred));
      auto const curr_duration = e.duration(backward_n);
      auto const pred_duration = b.cost2_.at(pred->get_key()).duration(*pred);
      auto const expected_duration =
          clamp_sub_duration(curr_duration, pred_duration);
      backward_dist +=
          add_path<P>(params, w, *w.r_, blocked, sharing, elevations, *pred,
                      backward_n, pred_duration, {}, expected_cost,
                      expected_duration, backward_segments, opposite(dir));
    } else {
      break;
    }
    backward_n = *pred;
  }

  auto const& dest_node_candidate =
      backward_n.get_node() == dest_left.node_ ? dest_left : dest_right;

  backward_segments.push_back(
      {.polyline_ =
           l.get_node_candidate_path(dest_way, dest_node_candidate.node_,
                                     dest_node_candidate.way_dir_, true, to),
       .from_level_ = dest_node_candidate.lvl_,
       .to_level_ = dest_node_candidate.lvl_,
       .from_ = dir == direction::kForward ? backward_n.get_node()
                                           : node_idx_t::invalid(),
       .to_ = dir == direction::kBackward ? backward_n.get_node()
                                          : node_idx_t::invalid(),
       .way_ = way_idx_t::invalid(),
       .cost_ = dest_node_candidate.cost_,
       .duration_ = duration_from_cost(dest_node_candidate.cost_),
       .dist_ = static_cast<distance_t>(dest_node_candidate.dist_to_node_),
       .mode_ = backward_n.get_mode()});

  if (dir == direction::kForward) {
    std::reverse(forward_segments.begin(), forward_segments.end());
  } else {
    std::reverse(backward_segments.begin(), backward_segments.end());
  }
  forward_segments.insert(forward_segments.end(), backward_segments.begin(),
                          backward_segments.end());

  auto total_dist = start_node_candidate.dist_to_node_ + forward_dist +
                    backward_dist + dest_node_candidate.dist_to_node_;

  auto path_elevation = elevation_storage::elevation{};
  for (auto const& segment : forward_segments) {
    path_elevation += segment.elevation_;
  }
  auto p =
      path{.cost_ = cost,
           .duration_ = sum_segment_durations(forward_segments, duration_t{0}),
           .dist_ = total_dist,
           .elevation_ = path_elevation,
           .segments_ = forward_segments};

  b.cost2_.at(backward_n.get_key()).write(backward_n, p);
  return p;
}

template <Profile P, typename Search>
path reconstruct(typename P::parameters const& params,
                 ways const& w,
                 lookup const& l,
                 bitvec<node_idx_t> const* blocked,
                 sharing_data const* sharing,
                 elevation_storage const* elevations,
                 Search const& search,
                 location const& from,
                 location const& to,
                 way_idx_t const start_way,
                 candidate_node const& start_left,
                 candidate_node const& start_right,
                 way_idx_t const dest_way,
                 candidate_node const& dest_nc,
                 typename P::node const dest_node,
                 cost_t const cost,
                 direction const dir,
                 std::optional<routing_time_t> const start_time) {

  auto n = dest_node;
  auto segments = std::vector<path::segment>{
      {.polyline_ =
           l.get_node_candidate_path(dest_way, dest_nc.node_, dest_nc.way_dir_,
                                     dir == direction::kForward, to),
       .from_level_ = dest_nc.lvl_,
       .to_level_ = dest_nc.lvl_,
       .from_ =
           dir == direction::kForward ? n.get_node() : node_idx_t::invalid(),
       .to_ =
           dir == direction::kBackward ? n.get_node() : node_idx_t::invalid(),
       .way_ = way_idx_t::invalid(),
       .cost_ = dest_nc.cost_,
       .duration_ = duration_from_cost(dest_nc.cost_),
       .dist_ = static_cast<distance_t>(dest_nc.dist_to_node_),
       .mode_ = dest_node.get_mode()}};
  auto dist = 0.0;
  while (true) {
    auto const& e = search.cost_.at(n.get_key());
    auto const pred = e.pred(n);
    if (pred.has_value()) {
      auto const pred_duration =
          search.cost_.at(pred->get_key()).duration(*pred);
      auto const expected_cost =
          static_cast<cost_t>(e.cost(n) - search.get_cost(*pred));
      dist += add_path<P>(params, w, *w.r_, blocked, sharing, elevations, *pred,
                          n, pred_duration, start_time, expected_cost,
                          clamp_sub_duration(e.duration(n), pred_duration),
                          segments, dir);
    } else {
      break;
    }
    n = *pred;
  }

  auto const& start_nc =
      n.get_node() == start_left.node_ ? start_left : start_right;
  segments.push_back(
      {.polyline_ = l.get_node_candidate_path(
           start_way, start_nc.node_, start_nc.way_dir_,
           dir == direction::kBackward, from),
       .from_level_ = start_nc.lvl_,
       .to_level_ = start_nc.lvl_,
       .from_ =
           dir == direction::kBackward ? n.get_node() : node_idx_t::invalid(),
       .to_ = dir == direction::kForward ? n.get_node() : node_idx_t::invalid(),
       .way_ = way_idx_t::invalid(),
       .cost_ = start_nc.cost_,
       .duration_ = duration_from_cost(start_nc.cost_),
       .dist_ = static_cast<distance_t>(start_nc.dist_to_node_),
       .mode_ = n.get_mode()});
  if (dir == direction::kForward) {
    std::reverse(begin(segments), end(segments));
  }
  auto path_elevation = elevation_storage::elevation{};
  for (auto const& segment : segments) {
    path_elevation += segment.elevation_;
  }
  auto p = path{.cost_ = cost,
                .duration_ = sum_segment_durations(segments),
                .dist_ = start_nc.dist_to_node_ + dist + dest_nc.dist_to_node_,
                .elevation_ = path_elevation,
                .segments_ = segments};
  search.cost_.at(dest_node.get_key()).write(dest_node, p);
  return p;
}

bool component_seen(ways const& w,
                    match_view_t const& matches,
                    size_t match_idx,
                    unsigned times = 1) {
  auto this_component = w.r_->way_component_[matches.way_[match_idx]];
  for (auto j = 0U; j < match_idx; ++j) {
    if (w.r_->way_component_[matches.way_[j]] == this_component) {
      if (--times == 0) {
        return true;
      }
    }
  }
  return false;
}

template <Profile P, typename Search>
std::optional<std::tuple<candidate_node, way_idx_t, typename P::node, path>>
best_candidate(typename P::parameters const& params,
               ways const& w,
               Search& search,
               level_t const lvl,
               match_view_t const& m,
               cost_t const max,
               direction const dir,
               std::optional<routing_time_t> const start_time,
               bool should_continue,
               way_idx_t const start_way,
               double const limit_squared_max_matching_distance) {
  auto best_cost = path{.cost_ = std::numeric_limits<cost_t>::max(),
                        .duration_ = kMaxDuration};
  auto best_node = P::node::invalid();
  auto best = candidate_node{};
  auto have_best = false;

  auto const get_best = [&](way_idx_t const dest_way, candidate_node const& x) {
    P::resolve_all(*w.r_, x.node_, lvl, [&](auto&& node) {
      auto const target_cost = search.get_cost(node);
      if (target_cost == kInfeasible || target_cost > best_cost.cost_) {
        return;
      }

      auto const target_duration =
          search.cost_.at(node.get_key()).duration(node);
      if (!P::is_dest_reachable(params, *w.r_, w.timezones_, node, dest_way,
                                flip(opposite(dir), x.way_dir_), dir,
                                start_time, target_duration)) {
        return;
      }

      auto const dest_way_cost = P::way_cost(
          params, *w.r_, w.timezones_, dest_way,
          w.r_->way_properties_[dest_way], flip(opposite(dir), x.way_dir_),
          static_cast<distance_t>(x.dist_to_node_), start_time, target_duration,
          dir);
      if (dest_way_cost.cost_ == kInfeasible) {
        return;
      }

      auto const total_cost = target_cost + dest_way_cost.cost_;
      auto const total_duration =
          clamp_add_duration(target_duration, dest_way_cost.duration_);
      if (total_cost < best_cost.cost_ ||
          (total_cost == best_cost.cost_ &&
           total_duration < best_cost.duration_)) {
        best_node = node;
        best = x;
        have_best = true;
        best_cost.cost_ = static_cast<cost_t>(total_cost);
        best_cost.duration_ = total_duration;
      }
    });
  };

  auto const start_component = w.r_->way_component_[start_way];
  auto component_seen_ctr = 0;
  auto const n = m.size();
  for (auto j = std::size_t{0U}; j != n; ++j) {
    auto const dest_way = m.way_[j];
    if (start_component != w.r_->way_component_[dest_way]) {
      continue;
    }
    if (!should_continue && ++component_seen_ctr > 1) {
      break;
    }
    if (std::pow(m.dist_to_way_[j], 2) > limit_squared_max_matching_distance &&
        j > kBottomKDefinitelyConsidered) {
      break;
    }

    for (auto const& x : {m.left(j), m.right(j)}) {
      if (x.valid()) {
        get_best(dest_way, x);
      }
    }

    if (have_best) {
      return best_cost.cost_ < max ? std::optional{std::tuple{
                                         best, dest_way, best_node, best_cost}}
                                   : std::nullopt;
    }
  }
  return std::nullopt;
}

std::optional<path> try_direct(osr::location const& from,
                               osr::location const& to) {
  auto const dist = geo::distance(from.pos_, to.pos_);
  if (dist < 8.0) {
    return std::optional{path{
        .cost_ = 60U,
        .duration_ = duration_from_cost(60U),
        .dist_ = dist,
        .segments_ = {path::segment{.polyline_ = {from.pos_, to.pos_},
                                    .from_level_ = from.lvl_,
                                    .to_level_ = to.lvl_,
                                    .from_ = node_idx_t::invalid(),
                                    .to_ = node_idx_t::invalid(),
                                    .way_ = way_idx_t::invalid(),
                                    .cost_ = 60U,
                                    .duration_ = duration_from_cost(60U),
                                    .dist_ = static_cast<distance_t>(dist)}},
        .uses_elevator_ = false}};
  } else {
    return std::nullopt;
  }
}

template <Profile P>
std::optional<path> route_bidirectional(typename P::parameters const& params,
                                        ways const& w,
                                        lookup const& l,
                                        bidirectional<P>& b,
                                        location const& from,
                                        location const& to,
                                        match_view_t const& from_match,
                                        match_view_t const& to_match,
                                        cost_t const max,
                                        direction const dir,
                                        bitvec<node_idx_t> const* blocked,
                                        sharing_data const* sharing,
                                        elevation_storage const* elevations) {
  if (auto const direct = try_direct(from, to); direct.has_value()) {
    return *direct;
  }

  b.reset(params, std::max(kMinCostSettled, max), from, to);
  if (!b.search_bounds_valid_) {
    return std::nullopt;
  }

  auto const limit_squared_max_matching_distance =
      geo::approx_squared_distance(from.pos_, to.pos_,
                                   b.distance_lon_degrees_) /
      kMaxMatchingDistanceSquaredRatio;

  for (auto i = std::size_t{0U}; i != from_match.size(); ++i) {
    if (b.max_reached_1_ && component_seen(w, from_match, i)) {
      continue;
    }
    auto const start_way = from_match.way_[i];
    auto const start_left = from_match.left(i);
    auto const start_right = from_match.right(i);
    for (auto const* nc : {&start_left, &start_right}) {
      if (nc->valid() && nc->cost_ < max) {
        auto const start_cost = P::way_cost(
            params, *w.r_, w.timezones_, start_way,
            w.r_->way_properties_[start_way], flip(dir, nc->way_dir_),
            static_cast<distance_t>(nc->dist_to_node_), {}, duration_t{0}, dir);
        if (start_cost.cost_ == kInfeasible || start_cost.cost_ >= max) {
          continue;
        }
        P::resolve_start_node(
            *w.r_, start_way, nc->node_, from.lvl_, dir, [&](auto const node) {
              auto label = typename P::label{node, start_cost.cost_};
              label.track(label, *w.r_, start_way, node.get_node(), false);
              b.add_start(params, w, label, sharing, start_cost.duration_);
            });
      }
    }
    if (b.pq1_.empty()) {
      continue;
    }
    for (auto j = std::size_t{0U}; j != to_match.size(); ++j) {
      auto const end_way = to_match.way_[j];
      if (w.r_->way_component_[start_way] != w.r_->way_component_[end_way]) {
        continue;
      }
      if (b.max_reached_2_ && component_seen(w, to_match, j)) {
        continue;
      }
      if (std::pow(to_match.dist_to_way_[j], 2) >
              limit_squared_max_matching_distance &&
          j > kBottomKDefinitelyConsidered) {
        break;
      }
      auto const end_left = to_match.left(j);
      auto const end_right = to_match.right(j);
      for (auto const* nc : {&end_left, &end_right}) {
        if (nc->valid() && nc->cost_ < max) {
          P::resolve_start_node(
              *w.r_, end_way, nc->node_, to.lvl_, opposite(dir),
              [&](auto const node) {
                auto label = typename P::label{node, nc->cost_};
                label.track(label, *w.r_, end_way, node.get_node(), false);
                b.add_end(params, w, label, sharing);
              });
        }
      }
      if (b.pq2_.empty()) {
        continue;
      }
      auto const should_continue =
          b.run(params, w, *w.r_, std::max(kMinCostSettled, max), blocked,
                sharing, elevations, dir);

      if (b.meet_point_1_.get_node() == node_idx_t::invalid()) {
        if (should_continue) {
          continue;
        }
        return std::nullopt;
      }

      auto const cost = b.get_cost_to_mp(b.meet_point_1_, b.meet_point_2_);

      if (cost >= max) {
        return std::nullopt;
      }

      return reconstruct_bi(params, w, l, blocked, sharing, elevations, b, from,
                            to, start_way, start_left, start_right, end_way,
                            end_left, end_right, cost, dir);
    }
    b.pq1_.clear();
    b.pq2_.clear();
    b.cost2_.clear();
    b.max_reached_2_ = false;
  }
  return std::nullopt;
}

template <Profile P>
std::optional<path> route_dijkstra(
    typename P::parameters const& params,
    ways const& w,
    lookup const& l,
    dijkstra<P>& d,
    location const& from,
    location const& to,
    match_view_t const& from_match,
    match_view_t const& to_match,
    cost_t const max,
    direction const dir,
    std::optional<routing_time_t> const start_time,
    bitvec<node_idx_t> const* blocked,
    sharing_data const* sharing,
    elevation_storage const* elevations) {
  if (auto const direct = try_direct(from, to); direct.has_value()) {
    return *direct;
  }

  auto const limit_squared_max_matching_distance =
      std::pow(geo::distance(from.pos_, to.pos_), 2) /
      kMaxMatchingDistanceSquaredRatio;

  d.reset(std::max(kMinCostSettled, max));
  auto should_continue = true;
  for (auto i = std::size_t{0U}; i != from_match.size(); ++i) {
    if (!should_continue && component_seen(w, from_match, i)) {
      continue;
    }
    auto const start_way = from_match.way_[i];
    auto const start_left = from_match.left(i);
    auto const start_right = from_match.right(i);
    auto const same_component = [&] {
      for (auto k = std::size_t{0U}; k != to_match.size(); ++k) {
        if (w.r_->way_component_[start_way] ==
            w.r_->way_component_[to_match.way_[k]]) {
          return true;
        }
      }
      return false;
    }();
    if (!same_component) {
      continue;
    }

    for (auto const* nc : {&start_left, &start_right}) {
      if (nc->valid() && nc->cost_ < max) {
        auto const start_cost = P::way_cost(
            params, *w.r_, w.timezones_, start_way,
            w.r_->way_properties_[start_way], flip(dir, nc->way_dir_),
            static_cast<distance_t>(nc->dist_to_node_), start_time,
            duration_t{0}, dir);
        if (start_cost.cost_ == kInfeasible || start_cost.cost_ >= max) {
          continue;
        }
        P::resolve_start_node(
            *w.r_, start_way, nc->node_, from.lvl_, dir, [&](auto const node) {
              d.add_start(w, {node, start_cost.cost_}, start_cost.duration_);
            });
      }
    }

    if (d.pq_.empty()) {
      continue;
    }

    should_continue = d.run(params, w, *w.r_, std::max(kMinCostSettled, max),
                            start_time, blocked, sharing, elevations, dir) &&
                      should_continue;

    auto const c = best_candidate<P>(params, w, d, to.lvl_, to_match, max, dir,
                                     start_time, should_continue, start_way,
                                     limit_squared_max_matching_distance);
    if (c.has_value()) {
      auto const [nc, wc, node, p] = *c;
      return reconstruct<P>(params, w, l, blocked, sharing, elevations, d, from,
                            to, start_way, start_left, start_right, wc, nc,
                            node, p.cost_, dir, start_time);
    }
  }

  return std::nullopt;
}

template <Profile P>
std::optional<path> route_astar(typename P::parameters const& params,
                                ways const& w,
                                lookup const& l,
                                astar<P>& a,
                                location const& from,
                                location const& to,
                                match_view_t const& from_match,
                                match_view_t const& to_match,
                                cost_t const max,
                                direction const dir,
                                std::optional<routing_time_t> const start_time,
                                bitvec<node_idx_t> const* blocked,
                                sharing_data const* sharing,
                                elevation_storage const* elevations) {
  if (auto const direct = try_direct(from, to); direct.has_value()) {
    return *direct;
  }

  auto const limit_squared_max_matching_distance =
      std::pow(geo::distance(from.pos_, to.pos_), 2) /
      kMaxMatchingDistanceSquaredRatio;

  a.reset(std::max(kMinCostSettled, max), from, to);
  auto should_continue = true;
  for (auto i = std::size_t{0U}; i != from_match.size(); ++i) {
    if (!should_continue && component_seen(w, from_match, i)) {
      continue;
    }
    auto const start_way = from_match.way_[i];
    auto const start_left = from_match.left(i);
    auto const start_right = from_match.right(i);
    auto const same_component = [&] {
      for (auto k = std::size_t{0U}; k != to_match.size(); ++k) {
        if (w.r_->way_component_[start_way] ==
            w.r_->way_component_[to_match.way_[k]]) {
          return true;
        }
      }
      return false;
    }();
    if (!same_component) {
      continue;
    }

    a.reset(std::max(kMinCostSettled, max), from, to);
    auto component_seen_ctr = 0;
    for (auto j = std::size_t{0U}; j != to_match.size(); ++j) {
      auto const end_way = to_match.way_[j];
      if (w.r_->way_component_[start_way] != w.r_->way_component_[end_way]) {
        continue;
      }
      if (!should_continue && ++component_seen_ctr > 1) {
        continue;
      }
      if (std::pow(to_match.dist_to_way_[j], 2) >
              limit_squared_max_matching_distance &&
          j > kBottomKDefinitelyConsidered) {
        break;
      }

      auto const end_left = to_match.left(j);
      auto const end_right = to_match.right(j);
      for (auto const* nc : {&end_left, &end_right}) {
        if (nc->valid() && nc->cost_ < max) {
          P::resolve_all(*w.r_, nc->node_, to.lvl_, [&](auto const node) {
            if (!P::is_dest_reachable(params, *w.r_, w.timezones_, node,
                                      end_way,
                                      flip(opposite(dir), nc->way_dir_), dir,
                                      start_time, duration_t{0})) {
              return;
            }
            a.add_destination(params, w, sharing, node);
          });
        }
      }
    }

    if (a.destinations_.empty()) {
      continue;
    }

    for (auto const* nc : {&start_left, &start_right}) {
      if (nc->valid() && nc->cost_ < max) {
        auto const start_cost = P::way_cost(
            params, *w.r_, w.timezones_, start_way,
            w.r_->way_properties_[start_way], flip(dir, nc->way_dir_),
            static_cast<distance_t>(nc->dist_to_node_), start_time,
            duration_t{0}, dir);
        if (start_cost.cost_ == kInfeasible || start_cost.cost_ >= max) {
          continue;
        }
        P::resolve_start_node(
            *w.r_, start_way, nc->node_, from.lvl_, dir, [&](auto const node) {
              a.add_start(params, w, sharing,
                          typename P::label{node, start_cost.cost_},
                          start_cost.duration_);
            });
      }
    }

    if (a.pq_.empty()) {
      continue;
    }

    should_continue = a.run(params, w, *w.r_, std::max(kMinCostSettled, max),
                            start_time, blocked, sharing, elevations, dir) &&
                      should_continue;

    auto const c = best_candidate<P>(params, w, a, to.lvl_, to_match, max, dir,
                                     start_time, should_continue, start_way,
                                     limit_squared_max_matching_distance);
    if (c.has_value()) {
      auto const [nc, wc, node, p] = *c;
      return reconstruct<P>(params, w, l, blocked, sharing, elevations, a, from,
                            to, start_way, start_left, start_right, wc, nc,
                            node, p.cost_, dir, start_time);
    }
  }

  return std::nullopt;
}


// ---------------------------------------------------------------------------
// customizable contraction hierarchies
// ---------------------------------------------------------------------------

namespace {

struct cch_registry {
  cch const* get(std::filesystem::path const& p) {
    auto const key = p.generic_string();

    {
      auto const lock = std::shared_lock{cch_m_};
      if (auto const it = cch_.find(key); it != end(cch_)) {
        return it->second == nullptr ? nullptr : it->second->get();
      }
    }

    auto const lock = std::unique_lock{cch_m_};
    if (auto const it = cch_.find(key); it != end(cch_)) {
      return it->second == nullptr ? nullptr : it->second->get();
    }
    auto& slot = cch_[key];
    if (cch::exists(p)) {
      slot = std::make_unique<cista::wrapped<cch>>(cch::read(p));
      return slot->get();
    }
    return nullptr;
  }

  // The customization is expensive (seconds on a country sized graph), so it
  // must not run under the registry wide lock: that would block every other
  // path / profile for its whole duration. Only the lookup of the entry is
  // globally synchronized, the build itself is serialized per entry.
  //
  // The metric is handed out as a shared pointer, not as a reference: a query
  // that runs with different parameters may replace the entry's metric while
  // this one is still using it.
  template <typename Fn>
  std::shared_ptr<cch_metric const> metric(std::filesystem::path const& p,
                                           search_profile const profile,
                                           std::string_view const params,
                                           Fn&& build) {
    auto* e = static_cast<metric_entry*>(nullptr);
    {
      auto const lock = std::lock_guard{metrics_m_};
      auto& slot = metrics_[std::pair{p.generic_string(), profile}];
      if (slot == nullptr) {
        slot = std::make_unique<metric_entry>();
      }
      e = slot.get();
    }

    {
      auto const lock = std::shared_lock{e->m_};
      if (e->metric_ != nullptr && e->params_ == params) {
        return e->metric_;
      }
    }

    auto const lock = std::unique_lock{e->m_};
    if (e->metric_ != nullptr && e->params_ == params) {
      return e->metric_;
    }

    // A metric stored next to the graph is mapped instead of rebuilt, but
    // only if it was customized for exactly these parameters. Everything
    // else falls back to customizing here.
    if (profile == search_profile::kCar && cch_metric::exists(p)) {
      auto mapped =
          std::make_shared<cista::wrapped<cch_metric>>(cch_metric::read(p));
      if (mapped->get()->matches(params)) {
        e->metric_ =
            std::shared_ptr<cch_metric const>{mapped, mapped->get()};
        e->params_ = std::string{params};
        return e->metric_;
      }
    }

    auto m = std::make_shared<cch_metric>();
    build(*m);
    e->metric_ = std::move(m);
    e->params_ = std::string{params};
    return e->metric_;
  }

  struct metric_entry {
    std::shared_mutex m_;
    std::string params_;
    std::shared_ptr<cch_metric const> metric_;
  };

  std::shared_mutex cch_m_;
  std::mutex metrics_m_;
  std::map<std::string, std::unique_ptr<cista::wrapped<cch>>> cch_;
  std::map<std::pair<std::string, search_profile>,
           std::unique_ptr<metric_entry>>
      metrics_;
};

cch_registry& get_cch_registry() {
  static auto r = cch_registry{};
  return r;
}

}  // namespace

template <WayAwareProfile P>
cch_search<P>& get_cch_search() {
  static auto s = boost::thread_specific_ptr<cch_search<P>>{};
  if (s.get() == nullptr) {
    s.reset(new cch_search<P>{});
  }
  return *s.get();
}

template <Profile P>
rphast<P>& get_rphast() {
  static auto s = boost::thread_specific_ptr<rphast<P>>{};
  if (s.get() == nullptr) {
    s.reset(new rphast<P>{});
  }
  return *s.get();
}

template <WayAwareProfile P>
path reconstruct_cch(typename P::parameters const& params,
                     ways const& w,
                     lookup const& l,
                     cch const& c,
                     cch_metric const& m,
                     cch_search<P> const& s,
                     location const& from,
                     location const& to,
                     way_idx_t const start_way,
                     candidate_node const& start_left,
                     candidate_node const& start_right,
                     way_idx_t const dest_way,
                     candidate_node const& dest_left,
                     candidate_node const& dest_right,
                     direction const dir) {
  auto const& r = *w.r_;

  // 1. collect the arcs of the up-down path
  struct arc {
    cch_slot_idx_t slot_;
    cch_entry_idx_t entry_;
    bool up_;
  };

  auto arcs = std::vector<arc>{};

  auto rank = s.meet_rank_;
  auto port = s.meet_f_port_;
  while (true) {
    auto const& e = s.f_.at(to_idx(rank));
    auto const slot = e.slot_[port];
    if (slot == cch::kNoSlot) {
      break;
    }
    auto const up = e.up_[port] != 0U;
    arcs.emplace_back(arc{slot, e.arc_[port], up});
    port = e.pred_port_[port];
    if (slot != cch::kLoopSlot) {
      rank = up ? c.tail(slot) : c.adj_head_[slot];
    }
  }
  std::reverse(begin(arcs), end(arcs));

  auto const start_rank = rank;
  auto const start_port = port;

  rank = s.meet_rank_;
  port = s.meet_f_port_;
  if (s.meet_b_port_ != kMaxPorts) {
    port = s.meet_b_port_;
    while (true) {
      auto const& e = s.b_.at(to_idx(rank));
      auto const slot = e.slot_[port];
      utl::verify(slot != cch::kNoSlot, "cch: broken backward chain");
      arcs.emplace_back(arc{slot, e.arc_[port], false});
      auto const done = e.seed_pred_[port] != 0U;
      auto const pred_port = e.pred_port_[port];
      if (slot != cch::kLoopSlot) {
        rank = c.tail(slot);
      }
      port = pred_port;
      if (done) {
        break;
      }
    }
  }

  auto const dest_rank = rank;

  // the last mile may only be reachable after a turn around at the target
  if (auto const it = s.target_at_.find(to_idx(dest_rank));
      it != end(s.target_at_)) {
    for (auto p = port; it->second.loop_[p] != cch_search<P>::kNoLoop;
         p = it->second.next_[p]) {
      arcs.emplace_back(arc{cch::kLoopSlot, it->second.loop_[p], false});
    }
  }

  // 2. unpack into original edges
  auto edges = std::vector<cch_path_edge>{};
  auto unpacker = cch_unpacker<P>{params, w, c, m};
  for (auto const& a : arcs) {
    if (a.slot_ == cch::kLoopSlot) {
      unpacker.unpack_loop(a.entry_, edges);
    } else {
      unpacker.unpack(a.slot_, a.entry_, a.up_, edges);
    }
  }

  // 3. turn the edges into profile states + per step costs
  auto const start_node = c.order_[start_rank];
  auto const dest_node = c.order_[dest_rank];

  auto const& start_nc =
      start_node == start_left.node_ ? start_left : start_right;
  auto const& dest_nc = dest_node == dest_left.node_ ? dest_left : dest_right;

  auto const start_cost =
      P::way_cost(params, r, w.timezones_, start_way,
                  r.way_properties_[start_way], flip(dir, start_nc.way_dir_),
                  static_cast<distance_t>(start_nc.dist_to_node_), std::nullopt,
                  duration_t{0}, dir)
          .cost_;

  auto states = std::vector<typename P::node>{typename P::node{
      start_node, port_way_pos(start_port), port_dir(start_port)}};
  auto cum = std::vector<cost_t>{start_cost};
  auto steps = std::vector<cost_t>{0U};

  auto prev_port = start_port;
  for (auto const& e : edges) {
    auto const turn = cch_turn_cost<P>(params, r, w.timezones_, e.from_,
                                       prev_port, e.from_port_);
    utl::verify(turn != kInfeasible, "cch: infeasible turn while unpacking");
    auto const step = clamp_cost(static_cast<std::uint64_t>(turn) + e.cost_);
    states.emplace_back(typename P::node{e.to_, port_way_pos(e.to_port_),
                                         port_dir(e.to_port_)});
    steps.emplace_back(step);
    cum.emplace_back(clamp_cost(static_cast<std::uint64_t>(cum.back()) + step));
    prev_port = e.to_port_;
  }

  utl::verify(states.back().n_ == dest_node,
              "cch: unpacked path does not end at the destination");

  // 4. build the path
  auto segments = std::vector<path::segment>{
      {.polyline_ = l.get_node_candidate_path(dest_way, dest_nc.node_,
                                              dest_nc.way_dir_,
                                              dir == direction::kForward, to),
       .from_level_ = dest_nc.lvl_,
       .to_level_ = dest_nc.lvl_,
       .from_ = dir == direction::kForward ? states.back().get_node()
                                           : node_idx_t::invalid(),
       .to_ = dir == direction::kBackward ? states.back().get_node()
                                          : node_idx_t::invalid(),
       .way_ = way_idx_t::invalid(),
       .cost_ = dest_nc.cost_,
       .duration_ = duration_from_cost(dest_nc.cost_),
       .dist_ = static_cast<distance_t>(dest_nc.dist_to_node_),
       .mode_ = states.back().get_mode()}};

  auto dist = 0.0;
  for (auto k = states.size(); k-- > 1U;) {
    dist += add_path<P>(
        params, w, r, nullptr, nullptr, nullptr, states[k - 1U], states[k],
        duration_from_cost(cum[k - 1U]), std::nullopt, steps[k],
        clamp_sub_duration(duration_from_cost(cum[k]),
                           duration_from_cost(cum[k - 1U])),
        segments, dir);
  }

  segments.push_back(
      {.polyline_ = l.get_node_candidate_path(start_way, start_nc.node_,
                                              start_nc.way_dir_,
                                              dir == direction::kBackward,
                                              from),
       .from_level_ = start_nc.lvl_,
       .to_level_ = start_nc.lvl_,
       .from_ = dir == direction::kBackward ? states.front().get_node()
                                            : node_idx_t::invalid(),
       .to_ = dir == direction::kForward ? states.front().get_node()
                                         : node_idx_t::invalid(),
       .way_ = way_idx_t::invalid(),
       .cost_ = start_nc.cost_,
       .duration_ = duration_from_cost(start_nc.cost_),
       .dist_ = static_cast<distance_t>(start_nc.dist_to_node_),
       .mode_ = states.front().get_mode()});

  if (dir == direction::kForward) {
    std::reverse(begin(segments), end(segments));
  }

  return path{.cost_ = s.best(),
              .duration_ = sum_segment_durations(segments),
              .dist_ = start_nc.dist_to_node_ + dist + dest_nc.dist_to_node_,
              .segments_ = segments};
}

template <WayAwareProfile P>
std::optional<path> route_cch(typename P::parameters const& params,
                              ways const& w,
                              lookup const& l,
                              cch const& c,
                              cch_metric const& m,
                              cch_search<P>& s,
                              location const& from,
                              location const& to,
                              match_view_t const& from_match,
                              match_view_t const& to_match,
                              cost_t const max,
                              direction const dir) {
  if (auto const direct = try_direct(from, to); direct.has_value()) {
    return *direct;
  }

  auto const& r = *w.r_;
  auto const limit_squared_max_matching_distance =
      std::pow(geo::distance(from.pos_, to.pos_), 2) /
      kMaxMatchingDistanceSquaredRatio;

  s.clear();
  for (auto i = std::size_t{0U}; i != from_match.size(); ++i) {
    auto const start_way = from_match.way_[i];
    auto const start_left = from_match.left(i);
    auto const start_right = from_match.right(i);

    auto const same_component = [&] {
      for (auto k = std::size_t{0U}; k != to_match.size(); ++k) {
        if (r.way_component_[start_way] ==
            r.way_component_[to_match.way_[k]]) {
          return true;
        }
      }
      return false;
    }();
    if (!same_component) {
      continue;
    }

    for (auto const* nc : {&start_left, &start_right}) {
      if (nc->valid() && nc->cost_ < max) {
        auto const start_cost = P::way_cost(
            params, r, w.timezones_, start_way,
            r.way_properties_[start_way], flip(dir, nc->way_dir_),
            static_cast<distance_t>(nc->dist_to_node_), std::nullopt,
            duration_t{0}, dir);
        if (start_cost.cost_ == kInfeasible || start_cost.cost_ >= max) {
          continue;
        }
        P::resolve_start_node(
            r, start_way, nc->node_, from.lvl_, dir, [&](auto const node) {
              s.add_start(c, node.n_, make_port(node.way_, node.dir_),
                          start_cost.cost_);
            });
      }
    }

    if (s.starts_.empty()) {
      continue;
    }

    for (auto j = std::size_t{0U}; j != to_match.size(); ++j) {
      auto const dest_way = to_match.way_[j];
      if (r.way_component_[start_way] != r.way_component_[dest_way]) {
        continue;
      }
      if (std::pow(to_match.dist_to_way_[j], 2) >
              limit_squared_max_matching_distance &&
          j > kBottomKDefinitelyConsidered) {
        break;
      }

      auto const dest_left = to_match.left(j);
      auto const dest_right = to_match.right(j);

      s.targets_.clear();
      for (auto const* x : {&dest_left, &dest_right}) {
        if (!x->valid()) {
          continue;
        }
        auto const way_dir = flip(opposite(dir), x->way_dir_);
        auto const dest_way_cost =
            P::way_cost(params, r, w.timezones_, dest_way,
                        r.way_properties_[dest_way], way_dir,
                        static_cast<distance_t>(x->dist_to_node_),
                        std::nullopt, duration_t{0}, dir);
        if (dest_way_cost.cost_ == kInfeasible) {
          continue;
        }
        P::resolve_all(r, x->node_, to.lvl_, [&](auto const node) {
          if (!P::is_dest_reachable(params, r, w.timezones_, node, dest_way,
                                    way_dir, dir, std::nullopt,
                                    duration_t{0})) {
            return;
          }
          s.add_target(c, node.n_, make_port(node.way_, node.dir_),
                       dest_way_cost.cost_);
        });
      }

      if (s.targets_.empty()) {
        continue;
      }

      s.run(params, w, c, m, std::max(kMinCostSettled, max));

      if (!s.found()) {
        continue;
      }
      if (s.best() >= max) {
        return std::nullopt;
      }

      return reconstruct_cch<P>(params, w, l, c, m, s, from, to, start_way,
                                start_left, start_right, dest_way, dest_left,
                                dest_right, dir);
    }
  }

  return std::nullopt;
}

// One-to-many over the CCH via RPHAST. The targets are selected once and the
// upward search from the source runs once; a single sweep over the selected set
// then yields the distance to every target at the same time.
//
// Distances only: reconstruction is deliberately not done here. The caller
// learns which targets matter (a station on an optimal journey) only afterwards
// and re-routes point to point for those few, rather than paying to unpack
// hundreds of paths that will be discarded.
template <Profile P>
std::vector<std::optional<path>> route_rphast(
    typename P::parameters const& params,
    ways const& w,
    cch const& c,
    cch_metric const& m,
    rphast<P>& rp,
    location const& from,
    std::vector<location> const& to,
    match_view_t const& from_match,
    match_result const& to_match,
    cost_t const max,
    direction const dir) {
  auto result = std::vector<std::optional<path>>{};
  result.resize(to_match.size());
  if (from_match.empty()) {
    return result;
  }

  auto const& r = *w.r_;
  rp.clear();

  // Selection only needs the set of target nodes, so it is done once over all
  // candidates of all targets and stays valid for every source afterwards.
  for (auto k = std::size_t{0U}; k != to_match.size(); ++k) {
    auto const tm = to_match[match_idx_t{static_cast<match_idx_t::value_t>(k)}];
    for (auto j = std::size_t{0U}; j != tm.size(); ++j) {
      auto const dest_way = tm.way_[j];
      auto const dest_left = tm.left(j);
      auto const dest_right = tm.right(j);
      for (auto const* x : {&dest_left, &dest_right}) {
        if (!x->valid()) {
          continue;
        }
        auto const way_dir = flip(opposite(dir), x->way_dir_);
        auto const dc = P::way_cost(
            params, r, w.timezones_, dest_way, r.way_properties_[dest_way],
            way_dir, static_cast<distance_t>(x->dist_to_node_), std::nullopt,
            duration_t{0}, dir);
        if (dc.cost_ == kInfeasible || dc.cost_ >= max) {
          continue;
        }
        P::resolve_all(r, x->node_, to[k].lvl_, [&](auto const node) {
          if (!P::is_dest_reachable(params, r, w.timezones_, node, dest_way,
                                    way_dir, dir, std::nullopt,
                                    duration_t{0})) {
            return;
          }
          rp.add_target(c, node.n_);
        });
      }
    }
  }
  rp.select(r, c, /* pack */ false, dir == direction::kBackward);

  auto const distance_lng_degrees = geo::approx_distance_lng_degrees(from.pos_);
  auto found = std::size_t{0U};
  auto should_continue = true;

  // Start candidates are added one at a time and the sweep re-run, mirroring
  // how the reference Dijkstra widens its start set: the nearest match is tried
  // first and a further one only enters if some target is still unresolved.
  // Pooling every candidate up front instead would let a target be reached from
  // a start the reference never considers, and answer a different question --
  // in practice a strictly better one, which is exactly why it has to be
  // suppressed here rather than left in.
  for (auto i = std::size_t{0U}; i != from_match.size(); ++i) {
    if (!should_continue && component_seen(w, from_match, i)) {
      continue;
    }
    auto const start_way = from_match.way_[i];
    auto const start_left = from_match.left(i);
    auto const start_right = from_match.right(i);
    for (auto const* nc : {&start_left, &start_right}) {
      if (!nc->valid() || nc->cost_ >= max) {
        continue;
      }
      auto const sc = P::way_cost(
          params, r, w.timezones_, start_way, r.way_properties_[start_way],
          flip(dir, nc->way_dir_),
          static_cast<distance_t>(nc->dist_to_node_), std::nullopt,
          duration_t{0}, dir);
      if (sc.cost_ == kInfeasible || sc.cost_ >= max) {
        continue;
      }
      P::resolve_start_node(
          r, start_way, nc->node_, from.lvl_, dir, [&](auto const node) {
            rp.add_start(c, node.n_, make_port(node.way_, node.dir_),
                         sc.cost_);
          });
    }
    if (rp.no_starts()) {
      continue;
    }

    should_continue =
        rp.run(params, w, c, m, std::max(kMinCostSettled, max)) &&
        should_continue;

    auto const start_component = r.way_component_[start_way];

    for (auto k = std::size_t{0U}; k != to_match.size(); ++k) {
      if (result[k].has_value()) {
        continue;
      }
      auto const tm =
          to_match[match_idx_t{static_cast<match_idx_t::value_t>(k)}];
      auto const& t = to[k];

      if (auto const direct = try_direct(from, t); direct.has_value()) {
        result[k] = direct;
        ++found;
        continue;
      }

      auto const limit_squared_max_matching_distance =
          geo::approx_squared_distance(from.pos_, t.pos_,
                                       distance_lng_degrees) /
          kMaxMatchingDistanceSquaredRatio;
      if (std::pow(from_match.dist_to_way_[i], 2) >
              limit_squared_max_matching_distance &&
          i > kBottomKDefinitelyConsidered) {
        continue;
      }

      // Which of a location's matched ways may be used is a matching question,
      // not a routing one, and is resolved exactly as `best_candidate` does it
      // for the reference: candidates in order, bounded by how far off the
      // straight line the match sits, and the first way that yields any
      // candidate decides.
      auto component_seen_ctr = 0;
      for (auto j = std::size_t{0U}; j != tm.size(); ++j) {
        auto const dest_way = tm.way_[j];
        if (start_component != r.way_component_[dest_way]) {
          continue;
        }
        // Once the search has run into the cost bound, a second way in the same
        // component is not going to be reachable either. The reference gives up
        // here and so must this, or it would report routes the reference does
        // not -- which it otherwise does, being the more complete of the two.
        if (!should_continue && ++component_seen_ctr > 1) {
          break;
        }
        if (std::pow(tm.dist_to_way_[j], 2) >
                limit_squared_max_matching_distance &&
            j > kBottomKDefinitelyConsidered) {
          break;
        }

        auto best = kInfeasible;
        for (auto const& x : {tm.left(j), tm.right(j)}) {
          if (!x.valid()) {
            continue;
          }
          auto const way_dir = flip(opposite(dir), x.way_dir_);
          auto const dc = P::way_cost(
              params, r, w.timezones_, dest_way, r.way_properties_[dest_way],
              way_dir, static_cast<distance_t>(x.dist_to_node_), std::nullopt,
              duration_t{0}, dir);
          if (dc.cost_ == kInfeasible) {
            continue;
          }
          P::resolve_all(r, x.node_, t.lvl_, [&](auto const node) {
            if (!P::is_dest_reachable(params, r, w.timezones_, node, dest_way,
                                      way_dir, dir, std::nullopt,
                                      duration_t{0})) {
              return;
            }
            auto const d = rp.get(c, node.n_, make_port(node.way_, node.dir_));
            if (d == kInfeasible) {
              return;
            }
            auto const total =
                clamp_cost(static_cast<std::uint64_t>(d) + dc.cost_);
            if (total < best) {
              best = total;
            }
          });
        }

        if (best != kInfeasible) {
          if (best < max) {
            // The hierarchy carries cost only -- the metric has no duration
            // array -- so the duration reported here is the cost, which for
            // this profile is a conservative bound: cost includes penalties
            // (u-turns and the like) that are not time, so it is never below
            // the true duration. Measured over 2,896 reference paths it was
            // exact for 98% and never underestimated, with a worst case of
            // +1,176s. An exact duration needs a second customized metric;
            // callers that need one for a specific target re-route point to
            // point, which reconstructs and reports it exactly.
            result[k] = path{.cost_ = best,
                             .duration_ = duration_from_cost(best)};
            ++found;
          }
          break;
        }
      }
    }

    if (found == result.size()) {
      break;
    }
  }

  return result;
}

template <Profile P>
std::vector<std::optional<path>> route(
    typename P::parameters const& params,
    ways const& w,
    lookup const& l,
    dijkstra<P>& d,
    location const& from,
    std::vector<location> const& to,
    match_view_t const& from_match,
    match_result const& to_match,
    cost_t const max,
    direction const dir,
    std::optional<routing_time_t> const start_time,
    bitvec<node_idx_t> const* blocked,
    sharing_data const* sharing,
    elevation_storage const* elevations,
    std::function<bool(path const&)> const& do_reconstruct) {
  auto result = std::vector<std::optional<path>>{};
  result.resize(to_match.size());

  if (from_match.empty()) {
    return result;
  }

  auto const distance_lng_degrees = geo::approx_distance_lng_degrees(from.pos_);

  d.reset(std::max(kMinCostSettled, max));
  auto should_continue = true;
  for (auto i = std::size_t{0U}; i != from_match.size(); ++i) {
    if (!should_continue && component_seen(w, from_match, i)) {
      continue;
    }
    auto const start_way = from_match.way_[i];
    auto const start_left = from_match.left(i);
    auto const start_right = from_match.right(i);
    for (auto const* nc : {&start_left, &start_right}) {
      if (nc->valid() && nc->cost_ < max) {
        auto const start_cost = P::way_cost(
            params, *w.r_, w.timezones_, start_way,
            w.r_->way_properties_[start_way], flip(dir, nc->way_dir_),
            static_cast<distance_t>(nc->dist_to_node_), start_time,
            duration_t{0}, dir);
        if (start_cost.cost_ == kInfeasible || start_cost.cost_ >= max) {
          continue;
        }
        P::resolve_start_node(
            *w.r_, start_way, nc->node_, from.lvl_, dir, [&](auto const node) {
              auto label = typename P::label{node, start_cost.cost_};
              label.track(label, *w.r_, start_way, node.get_node(), false);
              d.add_start(w, label, start_cost.duration_);
            });
      }
    }

    should_continue = d.run(params, w, *w.r_, std::max(kMinCostSettled, max),
                            start_time, blocked, sharing, elevations, dir) &&
                      should_continue;

    auto found = 0U;
    for (auto k = std::size_t{0U}; k != result.size(); ++k) {
      auto const m =
          to_match[match_idx_t{static_cast<match_idx_t::value_t>(k)}];
      auto const& t = to[k];
      auto& r = result[k];
      if (r.has_value()) {
        ++found;
      } else if (auto const direct = try_direct(from, t); direct.has_value()) {
        r = direct;
      } else {
        auto const limit_squared_max_matching_distance =
            geo::approx_squared_distance(from.pos_, t.pos_,
                                         distance_lng_degrees) /
            kMaxMatchingDistanceSquaredRatio;
        if (std::pow(from_match.dist_to_way_[i], 2) >
                limit_squared_max_matching_distance &&
            i > kBottomKDefinitelyConsidered) {
          continue;
        }

        auto const c = best_candidate<P>(params, w, d, t.lvl_, m, max, dir,
                                         start_time, should_continue, start_way,
                                         limit_squared_max_matching_distance);
        if (c.has_value()) {
          auto [nc, wc, n, p] = *c;
          d.cost_.at(n.get_key()).write(n, p);
          if (do_reconstruct(p)) {
            p = reconstruct<P>(params, w, l, blocked, sharing, elevations, d,
                               from, t, start_way, start_left, start_right, wc,
                               nc, n, p.cost_, dir, start_time);
            p.uses_elevator_ = true;
          }
          r = std::make_optional(p);
          ++found;
        }
      }
    }

    if (found == result.size()) {
      return result;
    }
  }

  return result;
}

std::optional<path> route_bidirectional(profile_parameters const& params,
                                        ways const& w,
                                        lookup const& l,
                                        search_profile const profile,
                                        location const& from,
                                        location const& to,
                                        cost_t const max,
                                        direction const dir,
                                        double const max_match_distance,
                                        bitvec<node_idx_t> const* blocked,
                                        sharing_data const* sharing,
                                        elevation_storage const* elevations) {
  return with_profile(profile, [&]<Profile P>(P&&) -> std::optional<path> {
    auto const& pp = std::get<typename P::parameters>(params);
    auto from_m = match_result{};
    l.complete_match<P>(pp, from, false, dir, max_match_distance, blocked,
                        std::nullopt, {}, from_m);
    auto to_m = match_result{};
    l.complete_match<P>(pp, to, true, dir, max_match_distance, blocked,
                        std::nullopt, {}, to_m);
    auto const from_match = from_m[match_idx_t{0U}];
    auto const to_match = to_m[match_idx_t{0U}];

    if (from_match.empty() || to_match.empty()) {
      return std::nullopt;
    }

    return route_bidirectional(pp, w, l, get_bidirectional<P>(), from, to,
                               from_match, to_match, max, dir, blocked, sharing,
                               elevations);
  });
}

std::vector<std::optional<path>> route(
    profile_parameters const& params,
    ways const& w,
    lookup const& l,
    search_profile const profile,
    location const& from,
    std::vector<location> const& to,
    cost_t const max,
    direction const dir,
    double const max_match_distance,
    bitvec<node_idx_t> const* blocked,
    sharing_data const* sharing,
    elevation_storage const* elevations,
    std::function<bool(path const&)> const& do_reconstruct,
    std::optional<routing_time_t> const start_time) {
  return with_profile(
      profile, [&]<Profile P>(P&&) -> std::vector<std::optional<path>> {
        auto const& pp = std::get<typename P::parameters>(params);
        auto from_m = match_result{};
        l.match<P>(pp, from, false, dir, max_match_distance, blocked, from_m,
                   start_time);
        auto const from_match = from_m[match_idx_t{0U}];
        if (from_match.empty()) {
          return std::vector<std::optional<path>>(to.size());
        }
        auto to_match = match_result{};
        for (auto const& x : to) {
          l.match<P>(pp, x, true, dir, max_match_distance, blocked, to_match,
                     start_time);
        }
        return route(pp, w, l, get_dijkstra<P>(), from, to, from_match,
                     to_match, max, dir, start_time, blocked, sharing,
                     elevations, do_reconstruct);
      });
}

std::optional<path> route_dijkstra(
    profile_parameters const& params,
    ways const& w,
    lookup const& l,
    search_profile const profile,
    location const& from,
    location const& to,
    cost_t const max,
    direction const dir,
    double const max_match_distance,
    bitvec<node_idx_t> const* blocked,
    sharing_data const* sharing,
    elevation_storage const* elevations,
    std::optional<routing_time_t> const start_time) {
  return with_profile(profile, [&]<Profile P>(P&&) -> std::optional<path> {
    auto const& pp = std::get<typename P::parameters>(params);
    auto from_m = match_result{};
    l.complete_match<P>(pp, from, false, dir, max_match_distance, blocked,
                        start_time, {}, from_m);
    auto to_m = match_result{};
    l.complete_match<P>(pp, to, true, dir, max_match_distance, blocked,
                        start_time, {}, to_m);
    auto const from_match = from_m[match_idx_t{0U}];
    auto const to_match = to_m[match_idx_t{0U}];

    if (from_match.empty() || to_match.empty()) {
      return std::nullopt;
    }

    return route_dijkstra(pp, w, l, get_dijkstra<P>(), from, to, from_match,
                          to_match, max, dir, start_time, blocked, sharing,
                          elevations);
  });
}

std::optional<path> route_astar(
    profile_parameters const& params,
    ways const& w,
    lookup const& l,
    search_profile const profile,
    location const& from,
    location const& to,
    cost_t const max,
    direction const dir,
    double const max_match_distance,
    bitvec<node_idx_t> const* blocked,
    sharing_data const* sharing,
    elevation_storage const* elevations,
    std::optional<routing_time_t> const start_time) {
  return with_profile(profile, [&]<Profile P>(P&&) -> std::optional<path> {
    auto const& pp = std::get<typename P::parameters>(params);
    auto from_m = match_result{};
    l.complete_match<P>(pp, from, false, dir, max_match_distance, blocked,
                        start_time, {}, from_m);
    auto to_m = match_result{};
    l.complete_match<P>(pp, to, true, dir, max_match_distance, blocked,
                        start_time, {}, to_m);
    auto const from_match = from_m[match_idx_t{0U}];
    auto const to_match = to_m[match_idx_t{0U}];

    if (from_match.empty() || to_match.empty()) {
      return std::nullopt;
    }

    return route_astar(pp, w, l, get_astar<P>(), from, to, from_match, to_match,
                       max, dir, start_time, blocked, sharing, elevations);
  });
}

std::vector<std::optional<path>> route(
    profile_parameters const& params,
    ways const& w,
    lookup const& l,
    search_profile const profile,
    location const& from,
    std::vector<location> const& to,
    match_view_t const& from_match,
    match_result const& to_match,
    cost_t const max,
    direction const dir,
    bitvec<node_idx_t> const* blocked,
    sharing_data const* sharing,
    elevation_storage const* elevations,
    std::function<bool(path const&)> const& do_reconstruct,
    std::optional<routing_time_t> const start_time,
    routing_algorithm const algo) {
  if (from_match.empty()) {
    return std::vector<std::optional<path>>(to.size());
  }
  return with_profile(
      profile, [&]<Profile P>(P&&) -> std::vector<std::optional<path>> {
        if (algo == routing_algorithm::kCCH) {
          if constexpr (WayAwareProfile<P> &&
                        requires(typename P::parameters const& pp) {
                          pp.uturn_penalty_;
                        }) {
            auto const* c = get_cch_registry().get(w.p_);
            utl::verify(c != nullptr, "no cch found in {}", w.p_);
            auto const& pp = std::get<typename P::parameters>(params);
            auto const blob = params_blob(pp);
            auto const m = get_cch_registry().metric(
                w.p_, profile, blob,
                [&](cch_metric& out) { customize<P>(pp, w, *c, out); });
            return route_rphast<P>(pp, w, *c, *m, get_rphast<P>(), from, to,
                                   from_match, to_match, max, dir);
          } else {
            throw utl::fail("cch not supported for profile {}",
                            to_str(profile));
          }
        }
        return route(std::get<typename P::parameters>(params), w, l,
                     get_dijkstra<P>(), from, to, from_match, to_match, max,
                     dir, start_time, blocked, sharing, elevations,
                     do_reconstruct);
      });
}


namespace {

bool cch_supported(search_profile const p,
                   direction const dir,
                   std::optional<routing_time_t> const& start_time,
                   bitvec<node_idx_t> const* blocked,
                   sharing_data const* sharing) {
  return dir == direction::kForward && !start_time.has_value() &&
         blocked == nullptr && sharing == nullptr &&
         (p == search_profile::kCar || p == search_profile::kBus ||
          p == search_profile::kHgv);
}

}  // namespace

std::optional<path> route_cch(profile_parameters const& params,
                              ways const& w,
                              lookup const& l,
                              search_profile const profile,
                              location const& from,
                              location const& to,
                              match_view_t const& from_match,
                              match_view_t const& to_match,
                              cost_t const max,
                              direction const dir) {
  auto const* c = get_cch_registry().get(w.p_);
  utl::verify(c != nullptr, "no cch found in {}", w.p_);

  return with_profile(profile, [&]<Profile P>(P&&) -> std::optional<path> {
    if constexpr (WayAwareProfile<P> &&
                  requires(typename P::parameters const& pp) {
                    pp.uturn_penalty_;
                  }) {
      auto const& pp = std::get<typename P::parameters>(params);
      auto const blob = params_blob(pp);
      auto const m = get_cch_registry().metric(
          w.p_, profile, blob,
          [&](cch_metric& out) { customize<P>(pp, w, *c, out); });
      return route_cch<P>(pp, w, l, *c, *m, get_cch_search<P>(), from, to,
                          from_match, to_match, max, dir);
    } else {
      throw utl::fail("cch not supported for profile {}", to_str(profile));
    }
  });
}

std::optional<path> route_cch(profile_parameters const& params,
                              ways const& w,
                              lookup const& l,
                              search_profile const profile,
                              location const& from,
                              location const& to,
                              cost_t const max,
                              direction const dir,
                              double const max_match_distance) {
  return with_profile(profile, [&]<Profile P>(P&&) -> std::optional<path> {
    auto const& pp = std::get<typename P::parameters>(params);
    auto from_m = match_result{};
    l.complete_match<P>(pp, from, false, dir, max_match_distance, nullptr,
                        std::nullopt, {}, from_m);
    auto to_m = match_result{};
    l.complete_match<P>(pp, to, true, dir, max_match_distance, nullptr,
                        std::nullopt, {}, to_m);
    return route_cch(params, w, l, profile, from, to, from_m[match_idx_t{0U}],
                     to_m[match_idx_t{0U}], max, dir);
  });
}

std::optional<path> route(profile_parameters const& params,
                          ways const& w,
                          lookup const& l,
                          search_profile const profile,
                          location const& from,
                          location const& to,
                          match_view_t const& from_match,
                          match_view_t const& to_match,
                          cost_t const max,
                          direction const dir,
                          bitvec<node_idx_t> const* blocked,
                          sharing_data const* sharing,
                          elevation_storage const* elevations,
                          routing_algorithm algo,
                          std::optional<routing_time_t> const start_time) {
  if (from_match.empty() || to_match.empty()) {
    return std::nullopt;
  }

  if (profile == search_profile::kBikeSharing ||
      profile == search_profile::kCarSharing) {
    algo = routing_algorithm::kDijkstra;  // TODO
  }
  if (algo == routing_algorithm::kCCH &&
      (!cch_supported(profile, dir, start_time, blocked, sharing) ||
       get_cch_registry().get(w.p_) == nullptr)) {
    algo = routing_algorithm::kDijkstra;
  } else if (profile == search_profile::kHgv &&
             algo != routing_algorithm::kCCH) {
    algo = routing_algorithm::kDijkstra;  // TODO
  }

  switch (algo) {
    case routing_algorithm::kDijkstra:
      return with_profile(profile, [&]<Profile P>(P&&) {
        return route_dijkstra(std::get<typename P::parameters>(params), w, l,
                              get_dijkstra<P>(), from, to, from_match, to_match,
                              max, dir, start_time, blocked, sharing,
                              elevations);
      });
    case routing_algorithm::kCCH:
      return route_cch(params, w, l, profile, from, to, from_match, to_match,
                       max, dir);
    case routing_algorithm::kAStarBi:
      return with_profile(profile, [&]<Profile P>(P&&) {
        auto const& pp = std::get<typename P::parameters>(params);
        if constexpr (requires { P::kExactBidirectional; }) {
          if constexpr (!P::kExactBidirectional) {
            return route_astar(pp, w, l, get_astar<P>(), from, to, from_match,
                               to_match, max, dir, start_time, blocked, sharing,
                               elevations);
          }
        }
        auto result = route_bidirectional(pp, w, l, get_bidirectional<P>(),
                                          from, to, from_match, to_match, max,
                                          dir, blocked, sharing, elevations);
        if constexpr (requires(typename P::node const n) {
                        P::bidirectional_meet_cost(pp, *w.r_, n, n);
                      }) {
          if (!result.has_value()) {
            return route_dijkstra(pp, w, l, get_dijkstra<P>(), from, to,
                                  from_match, to_match, max, dir, start_time,
                                  blocked, sharing, elevations);
          }
        }
        return result;
      });
  }
  throw utl::fail("not implemented");
}

std::optional<path> route(profile_parameters const& params,
                          ways const& w,
                          lookup const& l,
                          search_profile const profile,
                          location const& from,
                          location const& to,
                          cost_t const max,
                          direction const dir,
                          double const max_match_distance,
                          bitvec<node_idx_t> const* blocked,
                          sharing_data const* sharing,
                          elevation_storage const* elevations,
                          routing_algorithm algo,
                          std::optional<routing_time_t> const start_time) {
  if (profile == search_profile::kBikeSharing ||
      profile == search_profile::kCarSharing ||
      profile == search_profile::kCarParkingWheelchair ||
      profile == search_profile::kCarParking) {
    algo = routing_algorithm::kDijkstra;  // TODO
  }
  if (algo == routing_algorithm::kCCH &&
      (!cch_supported(profile, dir, start_time, blocked, sharing) ||
       get_cch_registry().get(w.p_) == nullptr)) {
    algo = routing_algorithm::kDijkstra;
  } else if (profile == search_profile::kHgv &&
             algo != routing_algorithm::kCCH) {
    algo = routing_algorithm::kDijkstra;  // TODO
  }
  switch (algo) {
    case routing_algorithm::kDijkstra:
      return route_dijkstra(params, w, l, profile, from, to, max, dir,
                            max_match_distance, blocked, sharing, elevations,
                            start_time);
    case routing_algorithm::kCCH:
      return route_cch(params, w, l, profile, from, to, max, dir,
                       max_match_distance);
    case routing_algorithm::kAStarBi:
      return route_bidirectional(params, w, l, profile, from, to, max, dir,
                                 max_match_distance, blocked, sharing,
                                 elevations);
  }
  throw utl::fail("not implemented");
}

}  // namespace osr
