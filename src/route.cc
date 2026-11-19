#include "osr/routing/route.h"

#include <cstdint>

#include <cmath>
#include <algorithm>
#include <memory>
#include <optional>
#include <ranges>
#include <span>
#include <utility>

#include "utl/erase_if.h"
#include "utl/helpers/algorithm.h"
#include "utl/to_vec.h"
#include "utl/verify.h"

#include "osr/elevation_storage.h"
#include "osr/lookup.h"
#include "osr/routing/astar.h"
#include "osr/routing/bidirectional.h"
#include "osr/routing/cost_search_limit.h"
#include "osr/routing/dijkstra.h"
#include "osr/routing/path_reconstruction.h"
#include "osr/routing/profiles/bike.h"
#include "osr/routing/profiles/bike_sharing.h"
#include "osr/routing/profiles/car.h"
#include "osr/routing/profiles/car_parking.h"
#include "osr/routing/profiles/car_sharing.h"
#include "osr/routing/profiles/foot.h"
#include "osr/routing/profiles/hgv.h"
#include "osr/routing/reachability.h"
#include "osr/routing/sharing_data.h"
#include "osr/routing/with_profile.h"
#include "osr/util/infinite.h"
#include "osr/util/reverse.h"

namespace osr {

constexpr auto const kMatchingPenaltyFactor = 20.0;

// Budget contribution of matching penalties per endpoint. Candidates with a
// higher penalty remain eligible; routes through them may exceed the budget.
constexpr auto const kMatchingPenaltyCap = cost_t{300U};

duration_t to_duration_limit(std::chrono::seconds const duration) {
  utl::verify(duration.count() >= 0 && duration.count() < kMaxDuration.count(),
              "duration limit must be between 0 and {} seconds",
              kMaxDuration.count() - 1U);
  return duration_t{static_cast<duration_t::rep>(duration.count())};
}

routing_algorithm to_algorithm(std::string_view s) {
  switch (cista::hash(s)) {
    case cista::hash("dijkstra"): return routing_algorithm::kDijkstra;
    case cista::hash("bidirectional"): return routing_algorithm::kAStarBi;
  }
  throw utl::fail("unknown routing algorithm: {}", s);
}

bool requires_dijkstra(search_profile const p) {
  switch (p) {
    case search_profile::kBikeSharing:
    case search_profile::kCarSharing:
    case search_profile::kCarParking:
    case search_profile::kCarParkingWheelchair:
    case search_profile::kCarDropOff:
    case search_profile::kCarDropOffWheelchair:
    case search_profile::kHgv: return true;
    default: return false;
  }
}

struct endpoint_candidate {
  way_idx_t way_;
  candidate_node node_;
  endpoint_state_t state_;
  double graph_distance_;  // projected point -> node
  cost_t matching_penalty_;
};

template <Profile P>
struct endpoint_root {
  typename P::node node_;
  endpoint_candidate endpoint_;
  cost_and_duration connection_;
};

cost_t matching_penalty(double const distance, double const closest_distance) {
  auto const additional_distance = std::max(0.0, distance - closest_distance);
  auto const penalty =
      std::round(distance + additional_distance * kMatchingPenaltyFactor);
  return !(penalty < static_cast<double>(kInfeasible))
             ? kInfeasible
             : static_cast<cost_t>(penalty);
}

template <typename Fn>
void for_each_endpoint_candidate(match_view_t const& matches, Fn&& fn) {
  for (auto i = std::size_t{0U}; i != matches.size(); ++i) {
    for (auto const& node : {matches.left(i), matches.right(i)}) {
      if (node.valid()) {
        auto const distance_to_way =
            static_cast<double>(matches.dist_to_way_[i]);
        fn(endpoint_candidate{.way_ = matches.way_[i],
                              .node_ = node,
                              .state_ = matches.state_[i],
                              .graph_distance_ = std::max(
                                  0.0, node.dist_to_node_ - distance_to_way),
                              .matching_penalty_ = matching_penalty(
                                  distance_to_way, matches.penalty_ref_)});
      }
    }
  }
}

template <Profile P>
cost_and_duration get_endpoint_connection(
    typename P::parameters const& params,
    ways const& w,
    endpoint_candidate const& endpoint,
    typename P::node const node,
    route_end const end,
    endpoint_role const role,
    std::optional<routing_time_t> const start_time,
    duration_t const current_duration,
    direction const search_dir) {
  auto const way_dir = flip(travel_dir_of(end), endpoint.node_.way_dir_);
  if (role == endpoint_role::kRoot &&
      !P::endpoint_root_allowed(params, node, way_dir)) {
    return infeasible_cost_and_duration();
  }
  auto const get_way_cost = [&](double const distance) {
    return P::endpoint_way_cost(params, *w.r_, w.timezones_, node,
                                endpoint.way_,
                                w.r_->way_properties_[endpoint.way_], way_dir,
                                static_cast<distance_t>(distance), start_time,
                                current_duration, search_dir);
  };
  auto total = clamp_add(
      cost_and_duration{
          .cost_ = get_way_cost(endpoint.graph_distance_).cost_,
          .duration_ = get_way_cost(endpoint.node_.dist_to_node_).duration_},
      endpoint.matching_penalty_, duration_t{0U});
  // Node costs are charged when arriving at a node. For the origin, that is
  // the connector from the query position.
  if (end == route_end::kOrigin) {
    total = clamp_add(
        total, P::endpoint_node_cost(params, node,
                                     w.r_->node_properties_[node.get_node()]));
  }
  if (role == endpoint_role::kGoal) {
    total = clamp_add(
        total, P::endpoint_transition_cost(params, *w.r_, w.timezones_, node,
                                           endpoint.way_, way_dir, search_dir,
                                           start_time, current_duration));
  }
  return total;
}

path::segment make_endpoint_segment(lookup const& l,
                                    location const& location,
                                    endpoint_candidate const& endpoint,
                                    cost_and_duration const connection,
                                    node_idx_t const graph_node,
                                    mode const mode,
                                    endpoint_role const role,
                                    bool const reverse,
                                    direction const dir) {
  auto const graph_node_at_from = role == endpoint_role::kRoot
                                      ? dir == direction::kBackward
                                      : dir == direction::kForward;
  return {.polyline_ = l.get_node_candidate_path(
              endpoint.way_, endpoint.node_.node_, endpoint.node_.way_dir_,
              reverse, location),
          .from_level_ = endpoint.node_.lvl_,
          .to_level_ = endpoint.node_.lvl_,
          .from_ = graph_node_at_from ? graph_node : node_idx_t::invalid(),
          .to_ = graph_node_at_from ? node_idx_t::invalid() : graph_node,
          .way_ = way_idx_t::invalid(),
          .cost_ = connection.cost_,
          .duration_ = connection.duration_,
          .dist_ = static_cast<distance_t>(endpoint.graph_distance_),
          .mode_ = mode};
}

template <Profile P, typename CostMap>
endpoint_root<P> const& get_endpoint_root(
    CostMap const& costs,
    std::span<endpoint_root<P> const> const roots,
    typename P::node const node) {
  auto const& entry = costs.at(node.get_key());
  auto const cd = cost_and_duration{.cost_ = entry.cost(node),
                                    .duration_ = entry.duration(node)};
  auto const root = std::find_if(begin(roots), end(roots), [&](auto const& r) {
    return r.node_ == node && r.connection_ == cd;
  });
  utl::verify(root != end(roots), "no exact endpoint root for settled state");
  return *root;
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
                    std::span<endpoint_root<P> const> const starts,
                    std::span<endpoint_root<P> const> const destinations,
                    cost_t const cost,
                    duration_t const duration,
                    direction const dir) {
  auto forward_n = b.meet_point_1_;

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

  auto const& start = get_endpoint_root<P>(b.cost1_, starts, forward_n);
  forward_segments.push_back(make_endpoint_segment(
      l, from, start.endpoint_, start.connection_, forward_n.get_node(),
      forward_n.get_mode(), endpoint_role::kRoot, dir == direction::kBackward,
      dir));

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

  auto const& destination =
      get_endpoint_root<P>(b.cost2_, destinations, backward_n);
  backward_segments.push_back(make_endpoint_segment(
      l, to, destination.endpoint_, destination.connection_,
      backward_n.get_node(), backward_n.get_mode(), endpoint_role::kGoal,
      dir == direction::kForward, dir));

  // Neither search half includes the turn joining the two meeting states,
  // so we add it to the segment leaving the meeting node in travel direction.
  auto& outgoing = dir == direction::kForward ? backward_segments.front()
                                              : forward_segments.front();
  outgoing.cost_ = clamp_cost(static_cast<std::uint64_t>(outgoing.cost_) +
                              b.best_transition_.cost_);
  outgoing.duration_ =
      clamp_add_duration(outgoing.duration_, b.best_transition_.duration_);

  if (dir == direction::kForward) {
    std::reverse(forward_segments.begin(), forward_segments.end());
  } else {
    std::reverse(backward_segments.begin(), backward_segments.end());
    forward_segments.swap(backward_segments);
  }
  forward_segments.insert(forward_segments.end(), backward_segments.begin(),
                          backward_segments.end());

  auto total_dist = start.endpoint_.graph_distance_ + forward_dist +
                    backward_dist + destination.endpoint_.graph_distance_;

  auto path_elevation = elevation_storage::elevation{};
  for (auto const& segment : forward_segments) {
    path_elevation += segment.elevation_;
  }
  auto p = path{.cost_ = cost,
                .duration_ = duration,
                .dist_ = total_dist,
                .elevation_ = path_elevation,
                .segments_ = forward_segments};

  b.cost1_.at(b.meet_point_1_.get_key()).write(b.meet_point_1_, p);
  auto const uses_elevator = p.uses_elevator_;
  b.cost2_.at(b.meet_point_2_.get_key()).write(b.meet_point_2_, p);
  p.uses_elevator_ = p.uses_elevator_ || uses_elevator;
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
                 std::span<endpoint_root<P> const> const starts,
                 endpoint_candidate const& destination,
                 cost_and_duration const destination_connection,
                 typename P::node const dest_node,
                 cost_and_duration const total,
                 direction const dir,
                 std::optional<routing_time_t> const start_time) {
  auto n = dest_node;
  auto segments = std::vector<path::segment>{make_endpoint_segment(
      l, to, destination, destination_connection, n.get_node(),
      dest_node.get_mode(), endpoint_role::kGoal, dir == direction::kForward,
      dir)};
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

  auto const& start = get_endpoint_root<P>(search.cost_, starts, n);
  segments.push_back(make_endpoint_segment(
      l, from, start.endpoint_, start.connection_, n.get_node(), n.get_mode(),
      endpoint_role::kRoot, dir == direction::kBackward, dir));
  if (dir == direction::kForward) {
    std::reverse(begin(segments), end(segments));
  }
  auto path_elevation = elevation_storage::elevation{};
  for (auto const& segment : segments) {
    path_elevation += segment.elevation_;
  }
  auto p = path{.cost_ = total.cost_,
                .duration_ = total.duration_,
                .dist_ = start.endpoint_.graph_distance_ + dist +
                         destination.graph_distance_,
                .elevation_ = path_elevation,
                .segments_ = segments};
  search.cost_.at(dest_node.get_key()).write(dest_node, p);
  return p;
}

// Turns every match candidate into the search roots it can be entered from,
// with an initial cost covering the virtual edge from the query position to
// the graph node.
template <Profile P, typename AddFn>
std::vector<endpoint_root<P>> add_endpoint_roots(
    typename P::parameters const& params,
    ways const& w,
    match_view_t const& matches,
    route_end const end,
    direction const dir,
    std::optional<routing_time_t> const start_time,
    cost_t const max,
    AddFn&& add) {
  auto roots = std::vector<endpoint_root<P>>{};
  for_each_endpoint_candidate(matches, [&](endpoint_candidate const& endpoint) {
    P::resolve_endpoint(
        *w.r_, endpoint.way_, endpoint.node_.node_, matches.lvl_, end,
        endpoint_role::kRoot, endpoint.state_, [&](auto const node) {
          auto const connection = get_endpoint_connection<P>(
              params, w, endpoint, node, end, endpoint_role::kRoot, start_time,
              duration_t{0U}, dir);
          if (!connection.feasible() || connection.cost_ >= max) {
            return;
          }
          auto label = typename P::label{node, connection.cost_};
          label.track(label, *w.r_, endpoint.way_, node.get_node(), false);
          roots.emplace_back(endpoint_root<P>{node, endpoint, connection});
          add(std::move(label), connection.duration_);
        });
  });
  return roots;
}

template <Profile P>
struct destination_candidate {
  endpoint_candidate endpoint_;
  typename P::node node_;
  cost_and_duration connection_;
  cost_and_duration total_;
};

template <Profile P, typename Search>
std::optional<destination_candidate<P>> best_candidate(
    typename P::parameters const& params,
    ways const& w,
    Search const& search,
    match_view_t const& matches,
    cost_t const max,
    direction const dir,
    std::optional<routing_time_t> const start_time) {
  auto best = std::optional<destination_candidate<P>>{};
  auto const end = route_end_of(opposite(dir));
  for_each_endpoint_candidate(matches, [&](endpoint_candidate const& endpoint) {
    auto const& candidate_node = endpoint.node_;
    auto const consider = [&](auto const node) {
      auto const target_cost = search.get_cost(node);
      if (target_cost == kInfeasible) {
        return;
      }
      auto const target_duration =
          search.cost_.at(node.get_key()).duration(node);
      auto const connection = get_endpoint_connection<P>(
          params, w, endpoint, node, end, endpoint_role::kGoal, start_time,
          target_duration, dir);
      if (!connection.feasible()) {
        return;
      }
      auto const total = clamp_add(
          cost_and_duration{.cost_ = target_cost, .duration_ = target_duration},
          connection);
      if (!total.feasible()) {
        return;
      }
      if (total.cost_ >= max) {
        return;
      }
      if (!best.has_value() || total < best->total_) {
        best = destination_candidate<P>{endpoint, node, connection, total};
      }
    };
    P::resolve_endpoint(*w.r_, endpoint.way_, candidate_node.node_,
                        matches.lvl_, end, endpoint_role::kGoal,
                        endpoint.state_, consider);
  });
  return best;
}

template <typename P>
concept single_mode_profile = requires { P::node::get_mode(); };

// Along way from the projection a to the projection b.
geo::polyline stretch_polyline(ways const& w,
                               way_idx_t const way,
                               lookup::way_stretch const& a,
                               lookup::way_stretch const& b) {
  auto polyline = geo::polyline{a.projection_};
  auto const way_polyline = w.way_polylines_[way];
  if (b.offset_ >= a.offset_) {
    for (auto k = a.segment_idx_ + 1U; k <= b.segment_idx_; ++k) {
      polyline.push_back(way_polyline[k]);
    }
  } else {
    for (auto k = a.segment_idx_; k > b.segment_idx_; --k) {
      polyline.push_back(way_polyline[k]);
    }
  }
  polyline.push_back(b.projection_);
  return polyline;
}

cost_and_duration total_of(path const& p) {
  return {.cost_ = p.cost_, .duration_ = p.duration_};
}

// With both ends between the same two routing nodes of a way, the search can
// only connect them via one of these nodes. This returns the direct piece
// (if below max).
template <Profile P>
std::optional<path> same_stretch_path(
    typename P::parameters const& params,
    ways const& w,
    lookup const& l,
    location const& from,
    location const& to,
    match_view_t const& from_match,
    match_view_t const& to_match,
    direction const dir,
    std::optional<routing_time_t> const start_time,
    cost_t const max) {
  if constexpr (!single_mode_profile<P>) {
    return std::nullopt;
  } else {
    auto const fwd = dir == direction::kForward;
    auto const& origin = fwd ? from : to;
    auto const& destination = fwd ? to : from;
    auto const& origin_match = fwd ? from_match : to_match;
    auto const& destination_match = fwd ? to_match : from_match;
    auto const on_level = [](way_properties const& props, location const& x) {
      return !x.lvl_.has_level() || is_on_level(props, x.lvl_);
    };

    auto best = std::optional<path>{};
    for (auto i = std::size_t{0U}; i != origin_match.size(); ++i) {
      auto const way = origin_match.way_[i];
      auto const it = utl::find(destination_match.way_, way);
      auto const props = w.r_->way_properties_[way];
      if (it == end(destination_match.way_) || !on_level(props, origin) ||
          !on_level(props, destination)) {
        continue;
      }
      auto const a = l.get_way_stretch(way, origin.pos_);
      auto const b = l.get_way_stretch(way, destination.pos_);
      if (a.candidate_.left_.node_ != b.candidate_.left_.node_ ||
          a.candidate_.right_.node_ != b.candidate_.right_.node_) {
        continue;
      }

      auto const along = std::abs(b.offset_ - a.offset_);
      auto const origin_dist =
          static_cast<double>(origin_match.dist_to_way_[i]);
      auto const destination_dist = static_cast<double>(
          destination_match.dist_to_way_[static_cast<std::size_t>(
              std::distance(begin(destination_match.way_), it))]);
      auto const get_way_cost = [&](double const distance) {
        return P::way_cost(
            params, *w.r_, w.timezones_, way, props,
            b.offset_ >= a.offset_ ? direction::kForward : direction::kBackward,
            static_cast<distance_t>(std::round(distance)), start_time,
            duration_t{0U}, dir);
      };
      auto const total = clamp_add(
          cost_and_duration{
              .cost_ = get_way_cost(along).cost_,
              .duration_ = get_way_cost(along + origin_dist + destination_dist)
                               .duration_},
          clamp_cost(static_cast<std::uint64_t>(matching_penalty(
                         origin_dist, origin_match.penalty_ref_)) +
                     matching_penalty(destination_dist,
                                      destination_match.penalty_ref_)),
          duration_t{0U});
      if (!total.feasible() || total.cost_ >= max ||
          (best.has_value() && total_of(*best) <= total)) {
        continue;
      }

      best = path{.cost_ = total.cost_,
                  .duration_ = total.duration_,
                  .dist_ = along,
                  .segments_ = {path::segment{
                      .polyline_ = stretch_polyline(w, way, a, b),
                      .from_level_ = props.from_level(),
                      .to_level_ = props.to_level(),
                      .way_ = way,
                      .cost_ = total.cost_,
                      .duration_ = total.duration_,
                      .dist_ = static_cast<distance_t>(std::round(along)),
                      .mode_ = P::node::get_mode()}}};
    }
    return best;
  }
}

std::optional<path> cheaper(std::optional<path> a, std::optional<path> b) {
  if (!a.has_value()) {
    return b;
  }
  if (!b.has_value()) {
    return a;
  }
  return total_of(*b) < total_of(*a) ? std::move(b) : std::move(a);
}

std::optional<path> try_direct(osr::location from,
                               osr::location to,
                               direction const dir) {
  auto const dist = geo::distance(from.pos_, to.pos_);
  if (dist >= 8.0) {
    return std::nullopt;
  }
  if (dir == direction::kBackward) {
    std::swap(from, to);
  }
  return std::optional{
      path{.cost_ = 60U,
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
}

cost_t max_matching_penalty(match_view_t const& matches) {
  auto max = cost_t{0U};
  for_each_endpoint_candidate(matches, [&](endpoint_candidate const& e) {
    if (e.matching_penalty_ != kInfeasible) {
      max = std::max(max, e.matching_penalty_);
    }
  });
  return max;
}

cost_t endpoint_cost_limit(cost_t const profile_limit,
                           cost_t const from_penalty,
                           cost_t const to_penalty) {
  auto const limit = static_cast<std::uint64_t>(profile_limit) +
                     std::min(from_penalty, kMatchingPenaltyCap) +
                     std::min(to_penalty, kMatchingPenaltyCap) + 1U;
  return static_cast<cost_t>(
      std::min<std::uint64_t>(limit, kMaxDurationSearchCost));
}

// Applied only to the final choice: a route rejected here must not make the
// caller fall back to another search that could return a different route.
std::optional<path> within_duration_limit(std::optional<path> p,
                                          duration_t const max_duration) {
  return p.has_value() && p->duration_ <= max_duration ? std::move(p)
                                                       : std::nullopt;
}

// How much the bidirectional heuristic may overestimate the connector to a
// match d meters away, which costs at least its matching penalty (>= d -
// 0.5).
template <Profile P>
cost_t heuristic_connector_slack(typename P::parameters const& params,
                                 match_view_t const& matches) {
  auto slack = 0.0;
  for (auto i = std::size_t{0U}; i != matches.size(); ++i) {
    auto const d = static_cast<double>(matches.dist_to_way_[i]);
    slack = std::max(slack, P::lower_bound_heuristic(params, d) - d + 0.5);
  }
  return clamp_cost(static_cast<std::uint64_t>(std::ceil(slack)));
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
  b.reset({.profile_ = params,
           .w_ = &w,
           .max_ = max,
           .dir_ = dir,
           .blocked_ = blocked,
           .sharing_ = sharing,
           .elevations_ = elevations,
           .start_loc_ = from,
           .end_loc_ = to},
          heuristic_connector_slack<P>(params, from_match),
          heuristic_connector_slack<P>(params, to_match));
  if (!b.search_bounds_valid_) {
    return std::nullopt;
  }

  auto const starts = add_endpoint_roots<P>(
      params, w, from_match, route_end_of(dir), dir, std::nullopt, max,
      [&](auto&& label, duration_t const duration) {
        b.add_start(std::forward<decltype(label)>(label), duration);
      });
  auto const destinations = add_endpoint_roots<P>(
      params, w, to_match, route_end_of(opposite(dir)), dir, std::nullopt, max,
      [&](auto&& label, duration_t const duration) {
        b.add_end(std::forward<decltype(label)>(label), duration);
      });
  if (starts.empty() || destinations.empty() || b.pq1_.empty() ||
      b.pq2_.empty()) {
    return std::nullopt;
  }

  b.run();
  if (b.meet_point_1_.get_node() == node_idx_t::invalid()) {
    return std::nullopt;
  }
  auto const cost = b.best_cost_;
  if (cost >= max) {
    return std::nullopt;
  }
  return reconstruct_bi<P>(params, w, l, blocked, sharing, elevations, b, from,
                           to, starts, destinations, cost, b.best_duration(),
                           dir);
}

template <Profile P>
std::optional<path> route_dijkstra(
    typename P::parameters const& params,
    ways const& w,
    lookup const& l,
    dijkstra<P, false>& d,
    location const& from,
    location const& to,
    match_view_t const& from_match,
    match_view_t const& to_match,
    cost_t const max,
    duration_t const max_duration,
    direction const dir,
    std::optional<routing_time_t> const start_time,
    bitvec<node_idx_t> const* blocked,
    sharing_data const* sharing,
    elevation_storage const* elevations,
    std::optional<path> stretch) {
  d.reset({.profile_ = params,
           .w_ = &w,
           .max_ = max,
           .max_duration_ = max_duration,
           .dir_ = dir,
           .start_time_ = start_time,
           .blocked_ = blocked,
           .sharing_ = sharing,
           .elevations_ = elevations,
           .start_loc_ = from,
           .end_loc_ = to});
  auto const starts = add_endpoint_roots<P>(
      params, w, from_match, route_end_of(dir), dir, start_time, max,
      [&](auto&& label, duration_t const duration) {
        d.add_start(std::forward<decltype(label)>(label), duration);
      });
  if (d.pq_.empty()) {
    return within_duration_limit(std::move(stretch), max_duration);
  }
  d.run();
  auto const find_candidate = [&]() {
    return best_candidate<P>(params, w, d, to_match, max, dir, start_time);
  };
  auto candidate = find_candidate();
  // Settle through the complete cost of the candidate and the direct piece:
  // a cheaper route through an over-duration state may still beat them.
  auto settle_cost = cost_t{0U};
  if (candidate.has_value() && candidate->total_.duration_ <= max_duration) {
    settle_cost = candidate->total_.cost_;
  }
  if (stretch.has_value() && stretch->duration_ <= max_duration) {
    settle_cost = std::max(settle_cost, stretch->cost_);
  }
  if (d.settle_up_to(settle_cost)) {
    candidate = find_candidate();
  }
  if (stretch.has_value() &&
      (!candidate.has_value() || total_of(*stretch) < candidate->total_)) {
    return within_duration_limit(std::move(stretch), max_duration);
  }
  if (!candidate.has_value() || candidate->total_.duration_ > max_duration) {
    return std::nullopt;
  }
  return reconstruct<P>(params, w, l, blocked, sharing, elevations, d, from, to,
                        starts, candidate->endpoint_, candidate->connection_,
                        candidate->node_, candidate->total_, dir, start_time);
}

template <Profile P>
std::optional<path> route_astar(typename P::parameters const& params,
                                ways const& w,
                                lookup const& l,
                                astar<P, false>& a,
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
  a.reset({.profile_ = params,
           .w_ = &w,
           .max_ = max,
           .dir_ = dir,
           .start_time_ = start_time,
           .blocked_ = blocked,
           .sharing_ = sharing,
           .elevations_ = elevations,
           .start_loc_ = from,
           .end_loc_ = to});
  for_each_endpoint_candidate(
      to_match, [&](endpoint_candidate const& endpoint) {
        auto const& candidate_node = endpoint.node_;
        auto const add = [&](auto const node) { a.add_destination(node); };
        P::resolve_endpoint(*w.r_, endpoint.way_, candidate_node.node_,
                            to_match.lvl_, route_end_of(opposite(dir)),
                            endpoint_role::kGoal, endpoint.state_, add);
      });
  if (a.destinations_.empty()) {
    return std::nullopt;
  }
  auto const starts = add_endpoint_roots<P>(
      params, w, from_match, route_end_of(dir), dir, start_time, max,
      [&](auto&& label, duration_t const duration) {
        a.add_start(std::forward<decltype(label)>(label), duration);
      });
  if (a.pq_.empty()) {
    return std::nullopt;
  }
  a.run();
  auto const candidate =
      best_candidate<P>(params, w, a, to_match, max, dir, start_time);
  if (!candidate.has_value()) {
    return std::nullopt;
  }
  return reconstruct<P>(params, w, l, blocked, sharing, elevations, a, from, to,
                        starts, candidate->endpoint_, candidate->connection_,
                        candidate->node_, candidate->total_, dir, start_time);
}

// Keeps a finished one-to-many search alive so that paths to individual
// destinations can be reconstructed later on demand.
template <Profile P>
struct one_to_many_state_impl final : public one_to_many_state {
  explicit one_to_many_state_impl(std::vector<location> to)
      : to_{std::move(to)}, candidates_(to_.size()) {}

  std::vector<std::optional<path>> const& results() const override {
    return results_;
  }

  std::optional<path> reconstruct(ways const& w,
                                  lookup const& l,
                                  std::size_t const k,
                                  sharing_data const* sharing) override {
    if (k >= results_.size() || !results_[k].has_value()) {
      return std::nullopt;
    }
    if (k >= candidates_.size() || !candidates_[k].has_value()) {
      return results_[k];  // direct path (from ~ to), nothing to reconstruct
    }
    auto const& c = *candidates_[k];
    auto const& sp = d_.params_;
    return osr::reconstruct<P>(sp.profile_, w, l, sp.blocked_, sharing,
                               sp.elevations_, d_, sp.start_loc_, to_[k],
                               starts_, c.endpoint_, c.connection_, c.node_,
                               c.total_, sp.dir_, sp.start_time_);
  }

  // The search owns everything it ran with (profile parameters, blocked,
  // direction, start time, the start location, ...) - see `osr::search_params`.
  // `sharing_data` is the exception: it is only valid while the search runs,
  // so `reconstruct()` takes it from the caller.
  dijkstra<P> d_;
  std::vector<location> to_;
  std::vector<endpoint_root<P>> starts_;
  std::vector<std::optional<destination_candidate<P>>> candidates_;
  std::vector<std::optional<path>> results_;
};

// Target indices per color, for colors with at least two members.
std::vector<std::vector<std::size_t>> color_groups(
    std::span<color_idx_t const> colors) {
  auto order = std::vector<std::size_t>{};
  for (auto i = std::size_t{0U}; i != colors.size(); ++i) {
    if (colors[i] != color_idx_t::invalid()) {
      order.push_back(i);
    }
  }
  std::ranges::stable_sort(order, {},
                           [&](std::size_t const i) { return colors[i]; });
  auto groups = std::vector<std::vector<std::size_t>>{};
  for (auto a = begin(order); a != end(order);) {
    auto const b = std::find_if(a, end(order), [&](std::size_t const i) {
      return colors[i] != colors[*a];
    });
    if (std::distance(a, b) > 1) {
      groups.emplace_back(a, b);
    }
    a = b;
  }
  return groups;
}

// Cost a member may take on top of a found member of its color: the search
// limit for the straight-line gap between them at the slowest speed.
template <Profile P>
cost_t color_budget(typename P::parameters const& params,
                    profile_parameters const& pp,
                    location const& found,
                    location const& member) {
  auto const gap =
      P::upper_bound_heuristic(params, geo::distance(found.pos_, member.pos_));
  return cost_search_limit(
      pp, clamp_duration(static_cast<std::uint64_t>(std::ceil(gap))));
}

// Continues a one-to-many search for the members of colors with a found
// member, each up to the cheapest found member's cost + `color_budget`.
template <Profile P>
void settle_colors(
    typename P::parameters const& params,
    ways const& w,
    lookup const& l,
    dijkstra<P, false>& d,
    location const& from,
    std::vector<location> const& to,
    match_view_t const& from_match,
    match_result const& to_match,
    direction const dir,
    std::optional<routing_time_t> const start_time,
    bitvec<node_idx_t> const* blocked,
    sharing_data const* sharing,
    elevation_storage const* elevations,
    std::function<bool(path const&)> const& do_reconstruct,
    std::vector<std::vector<std::size_t>> const& groups,
    std::vector<cost_t> const& destination_limits,
    std::span<endpoint_root<P> const> const starts,
    std::vector<std::optional<destination_candidate<P>>>& candidates,
    std::vector<std::optional<path>>& result) {
  constexpr auto const kStep = cost_t{60U};

  struct member {
    std::size_t i_;
    cost_t budget_;
    std::optional<path> stretch_;
  };

  auto const match_of = [&](std::size_t const i) {
    return to_match[match_idx_t{static_cast<match_idx_t::value_t>(i)}];
  };

  auto const pp = profile_parameters{params};
  auto pending = std::vector<member>{};
  for (auto const& g : groups) {
    for (auto const i : g) {
      if (result[i].has_value() || destination_limits[i] == 0U ||
          match_of(i).empty()) {
        continue;
      }
      auto budget = kInfeasible;
      for (auto const found : g) {
        if (result[found].has_value()) {
          budget = std::min(
              budget,
              clamp_cost(static_cast<std::uint64_t>(result[found]->cost_) +
                         color_budget<P>(params, pp, to[found], to[i])));
        }
      }
      if (budget == kInfeasible) {
        continue;  // no member of this color found
      }
      budget = std::min(budget, static_cast<cost_t>(d.params_.max_ - 1U));
      pending.push_back(
          member{i, budget,
                 same_stretch_path<P>(params, w, l, from, to[i], from_match,
                                      match_of(i), dir, start_time,
                                      clamp_cost(budget + std::uint64_t{1U}))});
    }
  }

  auto const finish = [&](member& m,
                          std::optional<destination_candidate<P>> const& c) {
    if (m.stretch_.has_value() &&
        (!c.has_value() || total_of(*m.stretch_) < c->total_)) {
      result[m.i_] = std::move(m.stretch_);
      candidates[m.i_].reset();
      return;
    }
    auto p = path{.cost_ = c->total_.cost_, .duration_ = c->total_.duration_};
    d.cost_.at(c->node_.get_key()).write(c->node_, p);
    if (do_reconstruct(p)) {
      p = reconstruct<P>(params, w, l, blocked, sharing, elevations, d, from,
                         to[m.i_], starts, c->endpoint_, c->connection_,
                         c->node_, c->total_, dir, start_time);
    }
    result[m.i_] = std::move(p);
    candidates[m.i_] = c;
  };

  // All states with a cost <= settled are final. A member is final once its
  // best total (candidate or same-stretch path) is <= settled.
  auto settled = d.settle_horizon_;
  while (!pending.empty()) {
    auto next = kInfeasible;
    utl::erase_if(pending, [&](member& m) {
      auto const c = best_candidate<P>(
          params, w, d, match_of(m.i_),
          clamp_cost(m.budget_ + std::uint64_t{1U}), dir, start_time);
      auto const total =
          std::min(c.has_value() ? c->total_.cost_ : kInfeasible,
                   m.stretch_.has_value() ? m.stretch_->cost_ : kInfeasible);
      if (total != kInfeasible && total <= settled) {
        finish(m, c);
        return true;
      }
      if (m.budget_ <= settled) {
        return true;  // not within its budget
      }
      next = std::min(
          next, total != kInfeasible
                    ? total
                    : std::min(m.budget_,
                               clamp_cost(static_cast<std::uint64_t>(settled) +
                                          kStep)));
      return false;
    });
    if (!pending.empty()) {
      settled = !d.settle_up_to(next) && d.pq_.empty() ? kInfeasible : next;
    }
  }
}

template <Profile P>
std::vector<std::optional<path>> route(
    typename P::parameters const& params,
    ways const& w,
    lookup const& l,
    dijkstra<P, false>& d,
    location const& from,
    std::vector<location> const& to,
    match_view_t const& from_match,
    match_result const& to_match,
    duration_t const max_duration,
    direction const dir,
    std::optional<routing_time_t> const start_time,
    bitvec<node_idx_t> const* blocked,
    sharing_data const* sharing,
    elevation_storage const* elevations,
    std::function<bool(path const&)> const& do_reconstruct,
    one_to_many_state_impl<P>* const state = nullptr,
    std::span<color_idx_t const> const colors = {}) {
  utl::verify(colors.empty() || colors.size() == to.size(),
              "route: {} colors for {} destinations", colors.size(), to.size());
  auto result = std::vector<std::optional<path>>(to_match.size());

  if (from_match.empty() || to_match.empty()) {
    return result;
  }

  auto const* components = reachability_components<P>(w, sharing);
  auto const profile_limit = cost_search_limit(params, max_duration);
  auto const from_penalty = max_matching_penalty(from_match);
  auto const destination_limits = utl::to_vec(
      std::views::iota(match_idx_t::value_t{0U},
                       static_cast<match_idx_t::value_t>(to_match.size())),
      [&](match_idx_t::value_t const i) {
        auto const to_i = to_match[match_idx_t{i}];
        if (components != nullptr &&
            !is_connected(*components, from_match, to_i)) {
          return cost_t{0U};
        }
        return endpoint_cost_limit(profile_limit, from_penalty,
                                   max_matching_penalty(to_i));
      });
  auto const from_connected_m =
      connected_candidates<P>(w, sharing, from_match, to_match.way_);
  auto const from_connected = from_connected_m[match_idx_t{0U}];
  auto const max = std::ranges::max(destination_limits);

  // Labels past `max` stay in the queue for the second phase (colors), up to
  // the largest budget a member can get.
  auto const groups = color_groups(colors);
  auto color_slack = cost_t{0U};
  if (!groups.empty()) {
    auto const pp = profile_parameters{params};
    for (auto const& g : groups) {
      for (auto const a : g) {
        for (auto const b : g) {
          color_slack =
              std::max(color_slack, color_budget<P>(params, pp, to[a], to[b]));
        }
      }
    }
  }

  d.reset({.profile_ = params,
           .w_ = &w,
           .max_ = static_cast<cost_t>(std::min(
               static_cast<std::uint64_t>(max) + color_slack,
               std::max(static_cast<std::uint64_t>(max),
                        static_cast<std::uint64_t>(kMaxDurationSearchCost)))),
           .settle_max_ = max,
           .max_duration_ = max_duration,
           .dir_ = dir,
           .start_time_ = start_time,
           .blocked_ = blocked,
           .sharing_ = sharing,
           .elevations_ = elevations,
           .start_loc_ = from});
  auto local_starts = std::vector<endpoint_root<P>>{};
  auto& starts = state == nullptr ? local_starts : state->starts_;
  starts = add_endpoint_roots<P>(
      params, w, from_connected, route_end_of(dir), dir, start_time, max,
      [&](auto&& label, duration_t const duration) {
        d.add_start(std::forward<decltype(label)>(label), duration);
      });
  d.run();
  auto const find_candidate = [&](std::size_t const i) {
    return best_candidate<P>(
        params, w, d,
        to_match[match_idx_t{static_cast<match_idx_t::value_t>(i)}],
        destination_limits[i], dir, start_time);
  };
  auto local_candidates =
      std::vector<std::optional<destination_candidate<P>>>{};
  auto& candidates = state == nullptr ? local_candidates : state->candidates_;
  candidates.resize(result.size());
  auto stretches = std::vector<std::optional<path>>(result.size());
  auto settle_cost = cost_t{0U};
  for (auto i = std::size_t{0U}; i != result.size(); ++i) {
    candidates[i] = find_candidate(i);
    if (candidates[i].has_value() &&
        candidates[i]->total_.duration_ <= max_duration) {
      settle_cost = std::max(settle_cost, candidates[i]->total_.cost_);
    }
    if (try_direct(from, to[i], dir).has_value()) {
      continue;
    }
    stretches[i] = same_stretch_path<P>(
        params, w, l, from, to[i], from_match,
        to_match[match_idx_t{static_cast<match_idx_t::value_t>(i)}], dir,
        start_time, destination_limits[i]);
    if (stretches[i].has_value() && stretches[i]->duration_ <= max_duration) {
      settle_cost = std::max(settle_cost, stretches[i]->cost_);
    }
  }
  // Finalize all candidates (and settle up to the direct pieces) before
  // reconstructing paths from the shared search.
  if (d.settle_up_to(settle_cost)) {
    for (auto i = std::size_t{0U}; i != result.size(); ++i) {
      candidates[i] = find_candidate(i);
    }
  }
  for (auto i = std::size_t{0U}; i != result.size(); ++i) {
    if (auto direct = try_direct(from, to[i], dir); direct.has_value()) {
      result[i] = within_duration_limit(std::move(direct), max_duration);
      candidates[i].reset();
      continue;
    }
    auto const& candidate = candidates[i];
    if (auto& stretch = stretches[i];
        stretch.has_value() &&
        (!candidate.has_value() || total_of(*stretch) < candidate->total_)) {
      result[i] = within_duration_limit(std::move(stretch), max_duration);
      candidates[i].reset();
      continue;
    }
    if (!candidate.has_value() || candidate->total_.duration_ > max_duration) {
      continue;
    }
    auto p = path{.cost_ = candidate->total_.cost_,
                  .duration_ = candidate->total_.duration_};
    d.cost_.at(candidate->node_.get_key()).write(candidate->node_, p);
    if (do_reconstruct(p)) {
      p = reconstruct<P>(params, w, l, blocked, sharing, elevations, d, from,
                         to[i], starts, candidate->endpoint_,
                         candidate->connection_, candidate->node_,
                         candidate->total_, dir, start_time);
    }
    result[i] = std::move(p);
  }

  if (!groups.empty()) {
    settle_colors<P>(params, w, l, d, from, to, from_match, to_match, dir,
                     start_time, blocked, sharing, elevations, do_reconstruct,
                     groups, destination_limits, starts, candidates, result);
  }
  return result;
}

std::optional<path> route_bidirectional(profile_parameters const& params,
                                        ways const& w,
                                        lookup const& l,
                                        search_profile const profile,
                                        location const& from,
                                        location const& to,
                                        std::chrono::seconds const max_duration,
                                        direction const dir,
                                        double const max_match_distance,
                                        bitvec<node_idx_t> const* blocked,
                                        sharing_data const* sharing,
                                        elevation_storage const* elevations) {
  return route(params, w, l, profile, from, to, max_duration, dir,
               max_match_distance, blocked, sharing, elevations,
               routing_algorithm::kAStarBi, std::nullopt);
}

std::vector<std::optional<path>> route(
    profile_parameters const& params,
    ways const& w,
    lookup const& l,
    search_profile const profile,
    location const& from,
    std::vector<location> const& to,
    std::chrono::seconds const max_duration,
    direction const dir,
    double const max_match_distance,
    bitvec<node_idx_t> const* blocked,
    sharing_data const* sharing,
    elevation_storage const* elevations,
    std::function<bool(path const&)> const& do_reconstruct,
    std::optional<routing_time_t> const start_time) {
  auto const duration_limit = to_duration_limit(max_duration);
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
        auto d = dijkstra<P>{};
        return route(pp, w, l, d, from, to, from_match, to_match,
                     duration_limit, dir, start_time, blocked, sharing,
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
    std::chrono::seconds const max_duration,
    direction const dir,
    double const max_match_distance,
    bitvec<node_idx_t> const* blocked,
    sharing_data const* sharing,
    elevation_storage const* elevations,
    std::optional<routing_time_t> const start_time) {
  return route(params, w, l, profile, from, to, max_duration, dir,
               max_match_distance, blocked, sharing, elevations,
               routing_algorithm::kDijkstra, start_time);
}

std::optional<path> route_astar(
    profile_parameters const& params,
    ways const& w,
    lookup const& l,
    search_profile const profile,
    location const& from,
    location const& to,
    std::chrono::seconds const max_duration,
    direction const dir,
    double const max_match_distance,
    bitvec<node_idx_t> const* blocked,
    sharing_data const* sharing,
    elevation_storage const* elevations,
    std::optional<routing_time_t> const start_time) {
  auto const duration_limit = to_duration_limit(max_duration);
  return with_profile(profile, [&]<Profile P>(P&&) -> std::optional<path> {
    auto const& pp = std::get<typename P::parameters>(params);
    auto from_m = match_result{};
    l.match<P>(pp, from, false, dir, max_match_distance, blocked, from_m,
               start_time);
    auto to_m = match_result{};
    l.match<P>(pp, to, true, dir, max_match_distance, blocked, to_m,
               start_time);
    auto const from_match = from_m[match_idx_t{0U}];
    auto const to_match = to_m[match_idx_t{0U}];

    if (from_match.empty() || to_match.empty()) {
      return std::nullopt;
    }

    if (auto direct = try_direct(from, to, dir); direct.has_value()) {
      return within_duration_limit(std::move(direct), duration_limit);
    }
    auto const max = endpoint_cost_limit(
        cost_search_limit(params, duration_limit),
        max_matching_penalty(from_match), max_matching_penalty(to_match));
    auto const from_connected_m =
        connected_candidates<P>(w, sharing, from_match, to_match.way_);
    auto const to_connected_m =
        connected_candidates<P>(w, sharing, to_match, from_match.way_);
    auto const from_connected = from_connected_m[match_idx_t{0U}];
    auto const to_connected = to_connected_m[match_idx_t{0U}];
    if (from_connected.empty() || to_connected.empty()) {
      return std::nullopt;
    }
    auto a = astar<P>{};
    return within_duration_limit(
        cheaper(route_astar(pp, w, l, a, from, to, from_connected, to_connected,
                            max, dir, start_time, blocked, sharing, elevations),
                same_stretch_path<P>(pp, w, l, from, to, from_connected,
                                     to_connected, dir, start_time, max)),
        duration_limit);
  });
}

std::unique_ptr<one_to_many_state> route_one_to_many(
    profile_parameters const& params,
    ways const& w,
    lookup const& l,
    search_profile const profile,
    location const& from,
    std::vector<location> const& to,
    match_view_t const& from_match,
    match_result const& to_match,
    std::chrono::seconds const max_duration,
    direction const dir,
    bitvec<node_idx_t> const* blocked,
    sharing_data const* sharing,
    elevation_storage const* elevations,
    std::function<bool(path const&)> const& do_reconstruct,
    std::optional<routing_time_t> const start_time,
    std::span<color_idx_t const> const colors) {
  auto const duration_limit = to_duration_limit(max_duration);
  return with_profile(
      profile, [&]<Profile P>(P&&) -> std::unique_ptr<one_to_many_state> {
        auto s = std::make_unique<one_to_many_state_impl<P>>(to);
        if (from_match.empty()) {
          s->results_.resize(to.size());
          return s;
        }
        s->results_ = route(std::get<typename P::parameters>(params), w, l,
                            s->d_, from, s->to_, from_match, to_match,
                            duration_limit, dir, start_time, blocked, sharing,
                            elevations, do_reconstruct, s.get(), colors);
        return s;
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
                          std::chrono::seconds const max_duration,
                          direction const dir,
                          bitvec<node_idx_t> const* blocked,
                          sharing_data const* sharing,
                          elevation_storage const* elevations,
                          routing_algorithm algo,
                          std::optional<routing_time_t> const start_time) {
  auto const duration_limit = to_duration_limit(max_duration);
  if (from_match.empty() || to_match.empty()) {
    return std::nullopt;
  }

  if (requires_dijkstra(profile)) {
    algo = routing_algorithm::kDijkstra;
  }

  auto p = with_profile(profile, [&]<Profile P>(P&&) -> std::optional<path> {
    auto const& pp = std::get<typename P::parameters>(params);
    if (auto direct = try_direct(from, to, dir); direct.has_value()) {
      return direct;
    }
    auto const max = endpoint_cost_limit(
        cost_search_limit(params, duration_limit),
        max_matching_penalty(from_match), max_matching_penalty(to_match));
    auto const from_connected_m =
        connected_candidates<P>(w, sharing, from_match, to_match.way_);
    auto const to_connected_m =
        connected_candidates<P>(w, sharing, to_match, from_match.way_);
    auto const from_connected = from_connected_m[match_idx_t{0U}];
    auto const to_connected = to_connected_m[match_idx_t{0U}];
    if (from_connected.empty() || to_connected.empty()) {
      return std::nullopt;
    }
    auto stretch = same_stretch_path<P>(pp, w, l, from, to, from_connected,
                                        to_connected, dir, start_time, max);
    switch (algo) {
      case routing_algorithm::kAStarBi: {
        auto b = bidirectional<P>{};
        return cheaper(route_bidirectional(pp, w, l, b, from, to,
                                           from_connected, to_connected, max,
                                           dir, blocked, sharing, elevations),
                       std::move(stretch));
      }
      case routing_algorithm::kDijkstra: {
        auto d = dijkstra<P>{};
        return route_dijkstra(pp, w, l, d, from, to, from_connected,
                              to_connected, max, duration_limit, dir,
                              start_time, blocked, sharing, elevations,
                              std::move(stretch));
      }
    }
    throw utl::fail("not implemented");
  });
  return within_duration_limit(std::move(p), duration_limit);
}

std::optional<path> route(profile_parameters const& params,
                          ways const& w,
                          lookup const& l,
                          search_profile const profile,
                          location const& from,
                          location const& to,
                          std::chrono::seconds const max_duration,
                          direction const dir,
                          double const max_match_distance,
                          bitvec<node_idx_t> const* blocked,
                          sharing_data const* sharing,
                          elevation_storage const* elevations,
                          routing_algorithm algo,
                          std::optional<routing_time_t> const start_time) {
  return with_profile(profile, [&]<Profile P>(P&&) {
    auto const& pp = std::get<typename P::parameters>(params);
    auto from_matches = match_result{};
    auto to_matches = match_result{};
    l.match<P>(pp, from, false, dir, max_match_distance, blocked, from_matches,
               start_time);
    l.match<P>(pp, to, true, dir, max_match_distance, blocked, to_matches,
               start_time);
    return route(params, w, l, profile, from, to, from_matches[match_idx_t{0U}],
                 to_matches[match_idx_t{0U}], max_duration, dir, blocked,
                 sharing, elevations, algo, start_time);
  });
}

}  // namespace osr
