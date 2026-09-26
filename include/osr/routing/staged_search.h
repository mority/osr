#pragma once

// Staged search (prototype).
//
// A staged search runs a sequence of single-profile Dijkstra searches
// ("stages") that are connected by transitions. Each stage is one of the plain
// profiles (foot, bike, car, ...) with its own cost table; a transition seeds a
// later stage from the settled nodes of an earlier stage, either through the
// additional nodes of a `sharing_data` (walk to a station, ride away from it)
// or directly at a node (free-floating pick-up / drop-off).
//
// The stages form a DAG in index order, so seeding stage j only from stages
// < j and running the stages in order yields the same costs as one joint
// search over the product state space (the state machine of bike_sharing /
// car_sharing is exactly such a DAG: initial foot -> vehicle -> trailing foot).
//
// The point of the decomposition: the foot stages can be shared between any
// number of rental products, only the vehicle stages are per product. The
// caller (e.g. motis) builds the stage list in a loop over its products.
//
// Lifetime: the search keeps pointers to the `sharing_data` of its
// transitions; they have to outlive the search (also for reconstruction).
//
// Known differences to bike_sharing / car_sharing:
//  - Whether a direct switch checks the target profile's node cost at the
//    switch node is a per-transition flag (`check_target_node_`).
//    bike_sharing checks it for the backward search only, car_sharing never
//    does. `add_rental_transitions()` follows bike_sharing.
//  - Reconstruction renders a direct switch as a zero-length segment that
//    carries the switch penalty, followed by the first edge after the switch;
//    the joint profiles fold the penalty into that edge.

#include <cassert>
#include <cmath>
#include <cstdint>

#include <algorithm>
#include <functional>
#include <map>
#include <memory>
#include <optional>
#include <vector>

#include "utl/helpers/algorithm.h"
#include "utl/verify.h"

#include "osr/elevation_storage.h"
#include "osr/location.h"
#include "osr/lookup.h"
#include "osr/routing/additional_edge.h"
#include "osr/routing/dijkstra.h"
#include "osr/routing/path.h"
#include "osr/routing/path_reconstruction.h"
#include "osr/routing/profile.h"
#include "osr/routing/sharing_data.h"
#include "osr/types.h"
#include "osr/ways.h"

namespace osr {

// Way properties used to cost the virtual edges between a real node and an
// additional node (station / free-floating vehicle). Accessible for every
// mode so that any profile can be on either side of a transition.
constexpr auto const kStageAdditionalWayProperties =
    way_properties{.is_foot_accessible_ = true,
                   .is_bike_accessible_ = true,
                   .is_car_accessible_ = true,
                   .is_destination_ = false,
                   .is_oneway_car_ = false,
                   .is_oneway_bike_ = false,
                   .is_elevator_ = false,
                   .is_steps_ = false,
                   .speed_limit_ = speed_limit::kmh_10,
                   .is_platform_ = 0,
                   .is_parking_ = false,
                   .is_ramp_ = false,
                   .is_sidewalk_separate_ = false,
                   .motor_vehicle_no_ = false,
                   .from_level_ = 0,
                   .has_toll_ = false,
                   .is_big_street_ = false,
                   .to_level_ = 0,
                   .is_bus_accessible_ = false,
                   .in_route_ = false,
                   .is_railway_accessible_ = false,
                   .is_oneway_bus_psv_ = false,
                   .is_incline_down_ = false,
                   .is_bus_accessible_with_penalty_ = false,
                   .is_ferry_accessible_ = false,
                   .is_railway_accessible_with_penalty_ = false,
                   .has_hgv_info_ = false,
                   .has_conditionals_ = false,
                   .is_in_low_emission_zone_ = false,
                   .is_detour_ = false,
                   .is_oneway_reverse_ = false};

// Profile-independent description of a profile node. Enough to rebuild the
// node of any profile with `P::create_node()`; fields a profile does not have
// keep their defaults.
struct generic_node {
  friend bool operator==(generic_node, generic_node) = default;

  node_idx_t n_{node_idx_t::invalid()};
  level_t lvl_{kNoLevel};
  way_pos_t way_{0U};
  direction dir_{direction::kForward};
};

template <Profile P>
generic_node to_generic(typename P::node const& n) {
  auto g = generic_node{.n_ = n.get_node()};
  if constexpr (requires { n.lvl_; }) {
    g.lvl_ = n.lvl_;
  }
  if constexpr (requires { n.way_; }) {
    g.way_ = n.way_;
  }
  if (auto const d = n.get_direction(); d.has_value()) {
    g.dir_ = *d;
  }
  return g;
}

template <Profile P>
typename P::node from_generic(generic_node const g) {
  return P::create_node(g.n_, g.lvl_, g.way_, g.dir_);
}

// Wraps a plain profile for use as a stage: identical search state, but the
// expansion is restricted to nodes allowed by `through_allowed_` (the
// geofence of a rental product). Additional nodes never appear inside a
// stage, they are handled by the transitions.
template <Profile P>
struct stage_profile {
  using inner_t = P;

  struct parameters {
    using profile_t = stage_profile<P>;
    typename P::parameters inner_{};
    bitvec<node_idx_t> const* through_allowed_{nullptr};
  };

  using key = typename P::key;
  using node = typename P::node;
  using label = typename P::label;
  using entry = typename P::entry;
  using hash = typename P::hash;

  static node create_node(node_idx_t const n,
                          level_t const lvl,
                          way_pos_t const way,
                          direction const dir) {
    return P::create_node(n, lvl, way, dir);
  }

  template <typename Fn>
  static void resolve_start_node(ways::routing const& w,
                                 way_idx_t const way,
                                 node_idx_t const n,
                                 level_t const lvl,
                                 direction const search_dir,
                                 Fn&& f) {
    P::resolve_start_node(w, way, n, lvl, search_dir, std::forward<Fn>(f));
  }

  template <typename Fn>
  static void resolve_all(ways::routing const& w,
                          node_idx_t const n,
                          level_t const lvl,
                          Fn&& f) {
    P::resolve_all(w, n, lvl, std::forward<Fn>(f));
  }

  template <direction SearchDir, bool WithBlocked, typename Fn>
  static void adjacent(parameters const& params,
                       ways::routing const& w,
                       timezone_cache_t const& timezones,
                       node const n,
                       duration_t const current_duration,
                       std::optional<routing_time_t> const start_time,
                       bitvec<node_idx_t> const* blocked,
                       sharing_data const*,
                       elevation_storage const* elevations,
                       Fn&& fn) {
    P::template adjacent<SearchDir, WithBlocked>(
        params.inner_, w, timezones, n, current_duration, start_time, blocked,
        nullptr, elevations,
        [&](node const neighbor, std::uint32_t const cost,
            duration_t const duration, distance_t const dist,
            way_idx_t const way, std::uint16_t const from,
            std::uint16_t const to,
            elevation_storage::elevation const elevation, bool const track) {
          if (is_allowed(params.through_allowed_, neighbor.get_node())) {
            fn(neighbor, cost, duration, dist, way, from, to, elevation, track);
          }
        });
  }

  static bool is_dest_reachable(parameters const& params,
                                ways::routing const& w,
                                timezone_cache_t const& timezones,
                                node const n,
                                way_idx_t const way,
                                direction const way_dir,
                                direction const search_dir,
                                std::optional<routing_time_t> const start_time,
                                duration_t const current_duration) {
    return P::is_dest_reachable(params.inner_, w, timezones, n, way, way_dir,
                                search_dir, start_time, current_duration);
  }

  static constexpr cost_and_duration way_cost(
      parameters const& params,
      ways::routing const& w,
      timezone_cache_t const& timezones,
      way_idx_t const way,
      way_properties const& e,
      direction const dir,
      distance_t const dist,
      std::optional<routing_time_t> const start_time,
      duration_t const current_duration,
      direction const search_dir) {
    return P::way_cost(params.inner_, w, timezones, way, e, dir, dist,
                       start_time, current_duration, search_dir);
  }

  static constexpr cost_and_duration node_cost(parameters const& params,
                                               node_properties const n) {
    return P::node_cost(params.inner_, n);
  }

  static constexpr double lower_bound_heuristic(parameters const& params,
                                                double const dist) {
    return P::lower_bound_heuristic(params.inner_, dist);
  }

  static constexpr double upper_bound_heuristic(parameters const& params,
                                                double const dist) {
    return P::upper_bound_heuristic(params.inner_, dist);
  }

  static constexpr node get_reverse(node const n) { return P::get_reverse(n); }
};

// Where a seed of a stage came from. Needed to continue path reconstruction
// in the previous stage once the predecessor chain of a stage ends.
//
// Station transition: from_node_ -> additional_ (mode of from_stage_)
//                     -> seed node (mode of the seeded stage, incl. penalty).
// Direct switch:      from_node_ is the switch node; the seed is a neighbor
//                     of the switch node in the seeded stage, reached via
//                     penalty + one edge (like the joint profiles, which
//                     never create the label at the switch node itself).
struct seed_origin {
  bool via_additional_node() const noexcept {
    return additional_ != node_idx_t::invalid();
  }

  std::size_t from_stage_{};
  generic_node from_node_{};  // node in `from_stage_` (always a real node)

  node_idx_t additional_{node_idx_t::invalid()};
  sharing_data const* sharing_{nullptr};  // coordinates of additional_
  cost_t enter_cost_{0U};
  distance_t enter_dist_{0U};
  cost_t leave_cost_{0U};
  distance_t leave_dist_{0U};

  generic_node switch_node_{};  // node of the seeded stage at from_node_.n_
  cost_t penalty_{0U};
  cost_t edge_cost_{0U};  // switch_node_ -> seed
  duration_t edge_duration_{0};
  duration_t switch_duration_{0};  // duration at the switch node
};

struct stage_transition {
  std::size_t from_{};
  std::size_t to_{};

  // Additional nodes (stations, free-floating vehicles) and their edges. The
  // switch happens at the additional node: walk/ride from a real node to the
  // additional node in the mode of `from_`, continue from the additional node
  // in the mode of `to_`. nullptr = no station transitions.
  sharing_data const* sharing_{nullptr};

  // Where the switch is allowed: checked at the additional node for station
  // transitions and at the real node for direct switches. nullptr = anywhere.
  bitvec<node_idx_t> const* allowed_{nullptr};

  // Allow switching at any settled real node of `from_` that is allowed
  // (free-floating zones).
  bool direct_switch_{false};

  // Direct switch only at nodes where the node cost of the target profile is
  // feasible.
  bool check_target_node_{false};

  cost_t penalty_{0U};
};

struct stage_search_params {
  ways const* w_{nullptr};
  cost_t max_{0U};
  direction dir_{direction::kForward};
  std::optional<routing_time_t> start_time_{};
  bitvec<node_idx_t> const* blocked_{nullptr};
  elevation_storage const* elevations_{nullptr};
  location start_loc_{};
};

// Destination hit in one stage (see `best_candidate()` in route.cc).
struct stage_dest_hit {
  cost_t cost_{kInfeasible};
  duration_t duration_{kMaxDuration};
  candidate_node nc_{};
  way_idx_t way_{way_idx_t::invalid()};
  generic_node node_{};
};

struct stage_base {
  virtual ~stage_base() = default;

  virtual mode get_mode() const = 0;
  virtual void reset(stage_search_params const&) = 0;
  virtual void add_start_candidate(way_idx_t start_way,
                                   candidate_node const&,
                                   level_t lvl) = 0;
  virtual bool has_pending() const = 0;
  virtual bool run() = 0;

  // Cheapest label of any node variant (level / direction / way) at `n`.
  virtual std::optional<std::pair<generic_node, cost_t>> min_cost_at(
      node_idx_t n) const = 0;
  virtual cost_t cost_of(generic_node) const = 0;
  virtual bool node_feasible(node_idx_t n) const = 0;
  virtual cost_and_duration additional_edge_cost(distance_t) const = 0;

  // Adds a seed at `n` with cost `c` if that improves the node. Records the
  // origin for reconstruction.
  virtual bool seed(node_idx_t n,
                    level_t lvl,
                    cost_t c,
                    seed_origin const&) = 0;

  // Direct switch at `n` with cost `c` (before the penalty): seeds every
  // neighbor of `n` in this stage with c + penalty + edge cost. `origin` has
  // the from-stage fields set. Returns true if the cost limit was hit.
  virtual bool seed_adjacent(node_idx_t n,
                             level_t lvl,
                             cost_t c,
                             cost_t penalty,
                             seed_origin origin) = 0;

  virtual std::optional<seed_origin> origin_of(generic_node) const = 0;

  // Every settled real node of this stage.
  virtual void for_each_settled(
      std::function<void(generic_node, cost_t)> const&) const = 0;

  virtual std::optional<stage_dest_hit> find_dest(
      location const& to,
      match_view_t const& m,
      way_idx_t start_way,
      double limit_squared_max_matching_distance,
      bool should_continue) const = 0;

  // Walks the predecessor chain from `from` to the root of this stage and
  // appends the segments (in search order, i.e. from destination towards the
  // start). Returns the root.
  virtual generic_node walk_preds(generic_node from,
                                  std::vector<path::segment>&,
                                  double& dist) const = 0;

  // Appends the segment for the edge `from` -> `to` of this stage's profile.
  virtual void add_edge_segment(generic_node from,
                                generic_node to,
                                cost_t expected_cost,
                                duration_t expected_duration,
                                duration_t from_duration,
                                std::vector<path::segment>&,
                                double& dist) const = 0;
};

template <Profile P>
struct stage_impl final : public stage_base {
  using sp = stage_profile<P>;
  using node = typename P::node;
  using label = typename P::label;
  using key = typename P::key;

  struct node_less {
    bool operator()(node const& a, node const& b) const { return a < b; }
  };

  stage_impl(typename P::parameters const& p,
             bitvec<node_idx_t> const* through_allowed)
      : params_{.inner_ = p, .through_allowed_ = through_allowed} {}

  mode get_mode() const override {
    return node::invalid().get_mode();  // static for all plain profiles
  }

  void reset(stage_search_params const& p) override {
    sp_ = p;
    d_.reset({.profile_ = params_,
              .w_ = p.w_,
              .max_ = p.max_,
              .dir_ = p.dir_,
              .start_time_ = p.start_time_,
              .blocked_ = p.blocked_,
              .sharing_ = nullptr,
              .elevations_ = p.elevations_,
              .start_loc_ = p.start_loc_});
    origins_.clear();
  }

  void add_start_candidate(way_idx_t const start_way,
                           candidate_node const& nc,
                           level_t const lvl) override {
    auto const& w = *sp_.w_;
    if (!nc.valid() || nc.cost_ >= sp_.max_) {
      return;
    }
    auto const start_cost = P::way_cost(
        params_.inner_, *w.r_, w.timezones_, start_way,
        w.r_->way_properties_[start_way], flip(sp_.dir_, nc.way_dir_),
        static_cast<distance_t>(nc.dist_to_node_), sp_.start_time_,
        duration_t{0}, sp_.dir_);
    if (start_cost.cost_ == kInfeasible || start_cost.cost_ >= sp_.max_) {
      return;
    }
    P::resolve_start_node(*w.r_, start_way, nc.node_, lvl, sp_.dir_,
                          [&](node const n) {
                            auto l = label{n, start_cost.cost_};
                            l.track(l, *w.r_, start_way, n.get_node(), false);
                            d_.add_start(l, start_cost.duration_);
                          });
  }

  bool has_pending() const override { return !d_.pq_.empty(); }

  bool run() override { return d_.run(); }

  std::optional<std::pair<generic_node, cost_t>> min_cost_at(
      node_idx_t const n) const override {
    auto best = std::optional<std::pair<generic_node, cost_t>>{};
    P::resolve_all(*sp_.w_->r_, n, kNoLevel, [&](node const x) {
      auto const c = d_.get_cost(x);
      if (c != kInfeasible && (!best.has_value() || c < best->second)) {
        best = std::pair{to_generic<P>(x), c};
      }
    });
    return best;
  }

  cost_t cost_of(generic_node const g) const override {
    return d_.get_cost(from_generic<P>(g));
  }

  bool node_feasible(node_idx_t const n) const override {
    return P::node_cost(params_.inner_, sp_.w_->r_->node_properties_[n])
               .cost_ != kInfeasible;
  }

  cost_and_duration additional_edge_cost(distance_t const dist) const override {
    auto const& w = *sp_.w_;
    return P::way_cost(params_.inner_, *w.r_, w.timezones_,
                       way_idx_t::invalid(), kStageAdditionalWayProperties,
                       direction::kForward, dist, sp_.start_time_,
                       duration_t{0}, sp_.dir_);
  }

  bool seed(node_idx_t const n,
            level_t const lvl,
            cost_t const c,
            seed_origin const& o) override {
    return seed_node(P::create_node(n, lvl, way_pos_t{0U}, direction::kForward),
                     c, duration_from_cost(c), o);
  }

  bool seed_adjacent(node_idx_t const n,
                     level_t const lvl,
                     cost_t const c,
                     cost_t const penalty,
                     seed_origin origin) override {
    auto const& w = *sp_.w_;
    auto const x = P::create_node(n, lvl, way_pos_t{0U}, direction::kForward);
    auto const x_duration = duration_from_cost(c);
    auto max_reached = false;

    origin.switch_node_ = to_generic<P>(x);
    origin.penalty_ = penalty;
    origin.switch_duration_ = x_duration;

    auto const expand = [&]<direction SearchDir, bool WithBlocked>() {
      P::template adjacent<SearchDir, WithBlocked>(
          params_.inner_, *w.r_, w.timezones_, x, x_duration, sp_.start_time_,
          sp_.blocked_, nullptr, sp_.elevations_,
          [&](node const neighbor, std::uint32_t const cost,
              duration_t const duration, distance_t, way_idx_t, std::uint16_t,
              std::uint16_t, elevation_storage::elevation, bool) {
            if (!is_allowed(params_.through_allowed_, neighbor.get_node())) {
              return;
            }
            auto const total = static_cast<std::uint64_t>(c) + penalty + cost;
            if (total >= sp_.max_) {
              max_reached = true;
              return;
            }
            origin.edge_cost_ = cost;
            origin.edge_duration_ = duration;
            seed_node(
                neighbor, static_cast<cost_t>(total),
                clamp_add_duration(clamp_add_duration(x_duration, duration),
                                   duration_from_cost(penalty)),
                origin);
          });
    };
    if (sp_.dir_ == direction::kForward) {
      sp_.blocked_ == nullptr
          ? expand.template operator()<direction::kForward, false>()
          : expand.template operator()<direction::kForward, true>();
    } else {
      sp_.blocked_ == nullptr
          ? expand.template operator()<direction::kBackward, false>()
          : expand.template operator()<direction::kBackward, true>();
    }
    return max_reached;
  }

  std::optional<seed_origin> origin_of(generic_node const g) const override {
    auto const it = origins_.find(from_generic<P>(g));
    return it == end(origins_) ? std::nullopt : std::optional{it->second};
  }

  void for_each_settled(
      std::function<void(generic_node, cost_t)> const& f) const override {
    auto const& r = *sp_.w_->r_;
    for (auto const& [k, e] : d_.cost_) {
      if constexpr (std::is_same_v<key, node>) {
        auto const c = e.cost(k);
        if (c != kInfeasible) {
          f(to_generic<P>(k), c);
        }
      } else {
        static_assert(std::is_same_v<key, node_idx_t>);
        P::resolve_all(r, k, kNoLevel, [&](node const x) {
          auto const c = e.cost(x);
          if (c != kInfeasible) {
            f(to_generic<P>(x), c);
          }
        });
      }
    }
  }

  std::optional<stage_dest_hit> find_dest(
      location const& to,
      match_view_t const& m,
      way_idx_t const start_way,
      double const limit_squared_max_matching_distance,
      bool const should_continue) const override {
    auto const& w = *sp_.w_;
    auto const dir = sp_.dir_;
    auto best = std::optional<stage_dest_hit>{};

    auto const get_best = [&](way_idx_t const dest_way,
                              candidate_node const& x) {
      P::resolve_all(*w.r_, x.node_, to.lvl_, [&](node const n) {
        auto const target_cost = d_.get_cost(n);
        if (target_cost == kInfeasible ||
            (best.has_value() && target_cost > best->cost_)) {
          return;
        }
        auto const target_duration = d_.cost_.at(n.get_key()).duration(n);
        if (!P::is_dest_reachable(params_.inner_, *w.r_, w.timezones_, n,
                                  dest_way, flip(opposite(dir), x.way_dir_),
                                  dir, sp_.start_time_, target_duration)) {
          return;
        }
        auto const dest_way_cost = P::way_cost(
            params_.inner_, *w.r_, w.timezones_, dest_way,
            w.r_->way_properties_[dest_way], flip(opposite(dir), x.way_dir_),
            static_cast<distance_t>(x.dist_to_node_), sp_.start_time_,
            target_duration, dir);
        if (dest_way_cost.cost_ == kInfeasible) {
          return;
        }
        auto const total_cost =
            static_cast<std::uint64_t>(target_cost) + dest_way_cost.cost_;
        auto const total_duration =
            clamp_add_duration(target_duration, dest_way_cost.duration_);
        if (!best.has_value() || total_cost < best->cost_ ||
            (total_cost == best->cost_ && total_duration < best->duration_)) {
          best = stage_dest_hit{.cost_ = static_cast<cost_t>(total_cost),
                                .duration_ = total_duration,
                                .nc_ = x,
                                .way_ = dest_way,
                                .node_ = to_generic<P>(n)};
        }
      });
    };

    auto const start_component = w.r_->way_component_[start_way];
    auto component_seen_ctr = 0;
    for (auto j = std::size_t{0U}; j != m.size(); ++j) {
      auto const dest_way = m.way_[j];
      if (start_component != w.r_->way_component_[dest_way]) {
        continue;
      }
      if (!should_continue && ++component_seen_ctr > 1) {
        break;
      }
      if (std::pow(m.dist_to_way_[j], 2) >
              limit_squared_max_matching_distance &&
          j > kBottomKDefinitelyConsidered) {
        break;
      }
      for (auto const& x : {m.left(j), m.right(j)}) {
        if (x.valid()) {
          get_best(dest_way, x);
        }
      }
      if (best.has_value()) {
        return best->cost_ < sp_.max_ ? best : std::nullopt;
      }
    }
    return std::nullopt;
  }

  generic_node walk_preds(generic_node const from,
                          std::vector<path::segment>& segments,
                          double& dist) const override {
    auto const& w = *sp_.w_;
    auto n = from_generic<P>(from);
    while (true) {
      auto const& e = d_.cost_.at(n.get_key());
      auto const pred = e.pred(n);
      if (!pred.has_value()) {
        break;
      }
      auto const pred_duration = d_.cost_.at(pred->get_key()).duration(*pred);
      auto const expected_cost =
          static_cast<cost_t>(e.cost(n) - d_.get_cost(*pred));
      dist += add_path<P>(
          params_.inner_, w, *w.r_, sp_.blocked_, nullptr, sp_.elevations_,
          *pred, n, pred_duration, sp_.start_time_, expected_cost,
          clamp_sub_duration(e.duration(n), pred_duration), segments, sp_.dir_);
      n = *pred;
    }
    return to_generic<P>(n);
  }

  void add_edge_segment(generic_node const from,
                        generic_node const to,
                        cost_t const expected_cost,
                        duration_t const expected_duration,
                        duration_t const from_duration,
                        std::vector<path::segment>& segments,
                        double& dist) const override {
    auto const& w = *sp_.w_;
    dist += add_path<P>(params_.inner_, w, *w.r_, sp_.blocked_, nullptr,
                        sp_.elevations_, from_generic<P>(from),
                        from_generic<P>(to), from_duration, sp_.start_time_,
                        expected_cost, expected_duration, segments, sp_.dir_);
  }

  static constexpr auto const kBottomKDefinitelyConsidered = 5;

private:
  bool seed_node(node const x,
                 cost_t const c,
                 duration_t const duration,
                 seed_origin const& o) {
    if (c >= sp_.max_ || c >= d_.get_cost(x)) {
      return false;
    }
    d_.add_start(label{x, c}, duration);
    origins_[x] = o;
    return true;
  }

  typename sp::parameters params_;
  stage_search_params sp_;
  dijkstra<sp> d_;
  std::map<node, seed_origin, node_less> origins_;
};

struct staged_search {
  static constexpr auto const kMinCostSettled = cost_t{900};
  static constexpr auto const kMaxMatchingDistanceSquaredRatio = 9.0;
  static constexpr auto const kBottomKDefinitelyConsidered = 5U;

  struct dest_candidate {
    std::size_t stage_;
    stage_dest_hit hit_;
  };

  template <Profile P>
  std::size_t add_stage(typename P::parameters const& p,
                        bitvec<node_idx_t> const* through_allowed = nullptr,
                        bool const terminal = false) {
    stages_.emplace_back(std::make_unique<stage_impl<P>>(p, through_allowed));
    terminal_.push_back(terminal);
    return stages_.size() - 1U;
  }

  void add_transition(stage_transition const t) {
    utl::verify(t.from_ < t.to_ && t.to_ < stages_.size(),
                "staged_search: transition {} -> {} not in stage order ({} "
                "stages)",
                t.from_, t.to_, stages_.size());
    transitions_.push_back(t);
  }

  std::vector<std::optional<path>> const& results() const { return results_; }

  void run(ways const& w,
           location const& from,
           std::vector<location> const& to,
           match_view_t const& from_match,
           match_result const& to_match,
           cost_t const max,
           direction const dir,
           bitvec<node_idx_t> const* blocked = nullptr,
           elevation_storage const* elevations = nullptr,
           std::optional<routing_time_t> const start_time = std::nullopt) {
    utl::verify(!stages_.empty(), "staged_search: no stages");

    w_ = &w;
    from_ = from;
    to_ = to;
    dir_ = dir;
    from_match_ = match_result{};
    from_match_.start(from_match.lvl_);
    for (auto j = std::size_t{0U}; j != from_match.size(); ++j) {
      from_match_.add(from_match.dist_to_way_[j], from_match.way_[j],
                      from_match.nodes_[j]);
    }
    from_match_.finish();

    results_.assign(to_match.size(), std::nullopt);
    dests_.assign(to_match.size(), std::nullopt);
    if (from_match.empty()) {
      return;
    }

    auto const params =
        stage_search_params{.w_ = &w,
                            .max_ = std::max(kMinCostSettled, max),
                            .dir_ = dir,
                            .start_time_ = start_time,
                            .blocked_ = blocked,
                            .elevations_ = elevations,
                            .start_loc_ = from};
    max_ = params.max_;
    for (auto& s : stages_) {
      s->reset(params);
    }

    auto const distance_lng_degrees =
        geo::approx_distance_lng_degrees(from.pos_);
    auto found = std::size_t{0U};
    max_reached_ = false;
    for (auto i = std::size_t{0U}; i != from_match.size(); ++i) {
      if (max_reached_ && component_seen(w, from_match, i)) {
        continue;
      }
      auto const start_way = from_match.way_[i];
      stages_[0]->add_start_candidate(start_way, from_match.left(i),
                                      from_match.lvl_);
      stages_[0]->add_start_candidate(start_way, from_match.right(i),
                                      from_match.lvl_);
      if (!stages_[0]->has_pending()) {
        continue;
      }

      // Stages are in topological order: everything that can seed stage j
      // has finished before j runs. Re-running after new start candidates is
      // plain Dijkstra with additional sources (labels only improve).
      max_reached_ |= !stages_[0]->run();
      for (auto j = std::size_t{1U}; j != stages_.size(); ++j) {
        for (auto const& t : transitions_) {
          if (t.to_ == j) {
            apply(t);
          }
        }
        max_reached_ |= !stages_[j]->run();
      }

      for (auto k = std::size_t{0U}; k != results_.size(); ++k) {
        if (results_[k].has_value()) {
          continue;
        }
        auto const& t = to[k];
        if (auto const direct = try_direct(from, t); direct.has_value()) {
          results_[k] = direct;
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
        auto const m =
            to_match[match_idx_t{static_cast<match_idx_t::value_t>(k)}];
        auto best = std::optional<dest_candidate>{};
        for (auto s = std::size_t{0U}; s != stages_.size(); ++s) {
          if (!terminal_[s]) {
            continue;
          }
          auto const hit = stages_[s]->find_dest(
              t, m, start_way, limit_squared_max_matching_distance,
              !max_reached_);
          if (hit.has_value() &&
              (!best.has_value() || hit->cost_ < best->hit_.cost_)) {
            best = dest_candidate{s, *hit};
          }
        }
        if (best.has_value()) {
          dests_[k] = best;
          results_[k] = path{.cost_ = best->hit_.cost_,
                             .duration_ = best->hit_.duration_};
          ++found;
        }
      }

      if (found == results_.size()) {
        return;
      }
    }
  }

  std::optional<path> reconstruct(lookup const& l, std::size_t const k) const {
    if (k >= results_.size() || !results_[k].has_value()) {
      return std::nullopt;
    }
    if (!dests_[k].has_value()) {
      return results_[k];  // direct path
    }

    auto const& w = *w_;
    auto const dir = dir_;
    auto const& dc = *dests_[k];
    auto const& hit = dc.hit_;
    auto const* stage = stages_[dc.stage_].get();

    auto segments = std::vector<path::segment>{
        {.polyline_ = l.get_node_candidate_path(
             hit.way_, hit.nc_.node_, hit.nc_.way_dir_,
             dir == direction::kForward, to_[k]),
         .from_level_ = hit.nc_.lvl_,
         .to_level_ = hit.nc_.lvl_,
         .from_ =
             dir == direction::kForward ? hit.node_.n_ : node_idx_t::invalid(),
         .to_ =
             dir == direction::kBackward ? hit.node_.n_ : node_idx_t::invalid(),
         .way_ = way_idx_t::invalid(),
         .cost_ = hit.nc_.cost_,
         .duration_ = duration_from_cost(hit.nc_.cost_),
         .dist_ = static_cast<distance_t>(hit.nc_.dist_to_node_),
         .mode_ = stage->get_mode()}};
    auto dist = 0.0;

    auto const node_pos = [&](node_idx_t const n,
                              sharing_data const* sharing) -> geo::latlng {
      if (sharing != nullptr && sharing->is_additional_node(n)) {
        return sharing->get_additional_node_coordinates(n);
      }
      return w.get_node_pos(n).as_latlng();
    };

    // Virtual edge pred -> curr (in search order) as one segment, like the
    // additional-edge branch of `add_path()`.
    auto const add_virtual = [&](node_idx_t const pred, node_idx_t const curr,
                                 cost_t const cost, distance_t const d,
                                 mode const m, sharing_data const* sharing) {
      auto& s = segments.emplace_back();
      s.from_ = dir == direction::kBackward ? curr : pred;
      s.to_ = dir == direction::kBackward ? pred : curr;
      s.polyline_ = {node_pos(s.from_, sharing), node_pos(s.to_, sharing)};
      s.from_level_ = level_t{0.0F};
      s.to_level_ = level_t{0.0F};
      s.way_ = way_idx_t::invalid();
      s.cost_ = cost;
      s.duration_ = duration_from_cost(cost);
      s.dist_ = d;
      s.mode_ = m;
      dist += d;
    };

    auto g = hit.node_;
    while (true) {
      auto const root = stage->walk_preds(g, segments, dist);
      auto const origin = stage->origin_of(root);
      if (!origin.has_value()) {
        g = root;
        break;
      }
      auto const* from_stage = stages_[origin->from_stage_].get();
      if (origin->via_additional_node()) {
        add_virtual(origin->additional_, root.n_, origin->leave_cost_,
                    origin->leave_dist_, stage->get_mode(), origin->sharing_);
        add_virtual(origin->from_node_.n_, origin->additional_,
                    origin->enter_cost_, origin->enter_dist_,
                    from_stage->get_mode(), origin->sharing_);
      } else {
        stage->add_edge_segment(origin->switch_node_, root, origin->edge_cost_,
                                origin->edge_duration_,
                                origin->switch_duration_, segments, dist);
        add_virtual(origin->switch_node_.n_, origin->switch_node_.n_,
                    origin->penalty_, distance_t{0U}, stage->get_mode(),
                    nullptr);
      }
      stage = from_stage;
      g = origin->from_node_;
    }

    // Start candidate: the from_match entry the root was seeded from.
    auto const root_cost = stage->cost_of(g);
    auto const fm = from_match_[match_idx_t{0U}];
    auto const is_start = [&](auto const& n) {
      return n.node_ == g.n_ && n.cost_ == root_cost;
    };
    auto const it = utl::find_if(fm.nodes_, [&](auto const& n) {
      return is_start(n.left_) || is_start(n.right_);
    });
    if (it == end(fm.nodes_)) {
      assert(false);
      return std::nullopt;
    }
    auto const start_idx =
        static_cast<std::size_t>(std::distance(begin(fm.nodes_), it));
    auto const start_way = fm.way_[start_idx];
    auto const start_nc =
        is_start(it->left_) ? fm.left(start_idx) : fm.right(start_idx);
    segments.push_back(
        {.polyline_ = l.get_node_candidate_path(
             start_way, start_nc.node_, start_nc.way_dir_,
             dir == direction::kBackward, from_),
         .from_level_ = start_nc.lvl_,
         .to_level_ = start_nc.lvl_,
         .from_ = dir == direction::kBackward ? g.n_ : node_idx_t::invalid(),
         .to_ = dir == direction::kForward ? g.n_ : node_idx_t::invalid(),
         .way_ = way_idx_t::invalid(),
         .cost_ = start_nc.cost_,
         .duration_ = duration_from_cost(start_nc.cost_),
         .dist_ = static_cast<distance_t>(start_nc.dist_to_node_),
         .mode_ = stage->get_mode()});

    if (dir == direction::kForward) {
      std::reverse(begin(segments), end(segments));
    }

    auto duration = duration_t{0};
    auto elevation = elevation_storage::elevation{};
    for (auto const& s : segments) {
      duration = clamp_add_duration(duration, s.duration_);
      elevation += s.elevation_;
    }
    return path{.cost_ = results_[k]->cost_,
                .duration_ = duration,
                .dist_ = start_nc.dist_to_node_ + dist + hit.nc_.dist_to_node_,
                .elevation_ = elevation,
                .segments_ = std::move(segments)};
  }

private:
  void apply(stage_transition const& t) {
    auto& a = *stages_[t.from_];
    auto& b = *stages_[t.to_];

    if (t.sharing_ != nullptr) {
      auto const& sharing = *t.sharing_;
      for (auto const& [s, edges] : sharing.additional_edges_) {
        if (sharing.is_additional_node(s)) {
          continue;
        }
        auto const best = a.min_cost_at(s);
        if (!best.has_value()) {
          continue;
        }
        auto const& [s_node, s_cost] = *best;
        for (auto const& ae : edges) {
          auto const add = ae.to_;
          if (!sharing.is_additional_node(add) ||
              !is_allowed(t.allowed_, add)) {
            continue;
          }
          auto const enter = a.additional_edge_cost(ae.distance_);
          if (!enter.feasible()) {
            continue;
          }
          auto const at_add = static_cast<std::uint64_t>(s_cost) + enter.cost_;
          if (at_add >= max_) {
            max_reached_ = true;
            continue;
          }
          auto const it = sharing.additional_edges_.find(add);
          if (it == end(sharing.additional_edges_)) {
            continue;
          }
          for (auto const& leave_edge : it->second) {
            auto const leave = clamp_add(
                b.additional_edge_cost(leave_edge.distance_), t.penalty_);
            if (!leave.feasible()) {
              continue;
            }
            auto const total = at_add + leave.cost_;
            if (total >= max_) {
              max_reached_ = true;
              continue;
            }
            b.seed(leave_edge.to_, s_node.lvl_, static_cast<cost_t>(total),
                   seed_origin{.from_stage_ = t.from_,
                               .from_node_ = s_node,
                               .additional_ = add,
                               .sharing_ = &sharing,
                               .enter_cost_ = enter.cost_,
                               .enter_dist_ = ae.distance_,
                               .leave_cost_ = leave.cost_,
                               .leave_dist_ = leave_edge.distance_});
          }
        }
      }
    }

    if (t.direct_switch_) {
      a.for_each_settled([&](generic_node const g, cost_t const c) {
        if (!is_allowed(t.allowed_, g.n_) ||
            (t.check_target_node_ && !b.node_feasible(g.n_))) {
          return;
        }
        max_reached_ |= b.seed_adjacent(
            g.n_, g.lvl_, c, t.penalty_,
            seed_origin{.from_stage_ = t.from_, .from_node_ = g});
      });
    }
  }

  static bool component_seen(ways const& w,
                             match_view_t const& matches,
                             std::size_t const match_idx) {
    auto const this_component = w.r_->way_component_[matches.way_[match_idx]];
    for (auto j = std::size_t{0U}; j != match_idx; ++j) {
      if (w.r_->way_component_[matches.way_[j]] == this_component) {
        return true;
      }
    }
    return false;
  }

  static std::optional<path> try_direct(location const& from,
                                        location const& to) {
    auto const dist = geo::distance(from.pos_, to.pos_);
    if (dist >= 8.0) {
      return std::nullopt;
    }
    return path{
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
                                    .dist_ = static_cast<distance_t>(dist)}}};
  }

  std::vector<std::unique_ptr<stage_base>> stages_;
  std::vector<bool> terminal_;
  std::vector<stage_transition> transitions_;

  ways const* w_{nullptr};
  location from_{};
  std::vector<location> to_;
  match_result from_match_;
  direction dir_{direction::kForward};
  cost_t max_{0U};
  bool max_reached_{false};
  std::vector<std::optional<path>> results_;
  std::vector<std::optional<dest_candidate>> dests_;
};

// Stage layout equivalent to bike_sharing / car_sharing for any number of
// rental products sharing the two foot stages:
//
//   forward:   foot(initial) -> vehicle_p -> foot(trailing)
//   backward:  foot(trailing) -> vehicle_p -> foot(initial)
//
// Usage: `first = add_stage<Foot>(..., terminal=true)`, then one
// `add_stage<Vehicle>(params, sharing.through_allowed_)` per product, then
// `last = add_stage<Foot>(..., terminal=true)`, then `add_rental_transitions`
// per product. `direct_pickup` = the vehicle may be picked up anywhere in
// the start zone, not only at additional nodes (car_sharing: yes,
// bike_sharing: no).
inline void add_rental_transitions(staged_search& s,
                                   std::size_t const first_foot,
                                   std::size_t const vehicle,
                                   std::size_t const last_foot,
                                   sharing_data const& sharing,
                                   direction const dir,
                                   bool const direct_pickup,
                                   cost_t const start_penalty = cost_t{30U},
                                   cost_t const end_penalty = cost_t{30U}) {
  if (dir == direction::kForward) {
    s.add_transition({.from_ = first_foot,
                      .to_ = vehicle,
                      .sharing_ = &sharing,
                      .allowed_ = sharing.start_allowed_,
                      .direct_switch_ = direct_pickup,
                      .penalty_ = start_penalty});
    s.add_transition({.from_ = vehicle,
                      .to_ = last_foot,
                      .sharing_ = &sharing,
                      .allowed_ = sharing.end_allowed_,
                      .direct_switch_ = true,
                      .penalty_ = end_penalty});
  } else {
    s.add_transition({.from_ = first_foot,
                      .to_ = vehicle,
                      .sharing_ = &sharing,
                      .allowed_ = sharing.end_allowed_,
                      .direct_switch_ = true,
                      .check_target_node_ = true,  // as bike_sharing
                      .penalty_ = end_penalty});
    // Note: car_sharing uses the end penalty for the direct switch here.
    s.add_transition({.from_ = vehicle,
                      .to_ = last_foot,
                      .sharing_ = &sharing,
                      .allowed_ = sharing.start_allowed_,
                      .direct_switch_ = direct_pickup,
                      .penalty_ = start_penalty});
  }
}

}  // namespace osr
