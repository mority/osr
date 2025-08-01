#pragma once

#include <bitset>

#include "boost/json.hpp"

#include "utl/helpers/algorithm.h"

#include "osr/elevation_storage.h"
#include "osr/routing/mode.h"
#include "osr/routing/profiles/car.h"
#include "osr/routing/profiles/foot.h"
#include "osr/routing/route.h"
#include "osr/ways.h"

namespace osr {

struct sharing_data;

template <bool IsWheelchair, bool UseParking = true>
struct car_parking {
  using footp = foot<IsWheelchair>;

  static constexpr auto const kSwitchPenalty = cost_t{200U};
  static constexpr auto const kMaxMatchDistance = car::kMaxMatchDistance;

  using key = node_idx_t;

  enum class node_type : std::uint8_t { kCar, kFoot, kInvalid };

  static constexpr std::string_view node_type_to_str(node_type const type) {
    switch (type) {
      case node_type::kCar: return "car";
      case node_type::kFoot: return "foot";
      case node_type::kInvalid: return "invalid";
    }
    std::unreachable();
  }

  struct node {
    friend bool operator==(node const& a, node const& b) {
      auto const is_zero = [](level_t const l) {
        return l == kNoLevel || l == level_t{0.F};
      };
      return a.n_ == b.n_ && a.type_ == b.type_ && a.dir_ == b.dir_ &&
             a.way_ == b.way_ &&
             (a.lvl_ == b.lvl_ || (is_zero(a.lvl_) && is_zero(b.lvl_)));
    }

    boost::json::object geojson_properties(ways const& w) const {
      auto properties =
          boost::json::object{{"osm_node_id", to_idx(w.node_to_osm_[n_])},
                              {"level", lvl_.to_float()},
                              {"type", node_type_to_str(type_)}};
      if (is_car_node()) {
        properties.emplace("direction", to_str(dir_));
      }
      return properties;
    }

    std::ostream& print(std::ostream& out, ways const& w) const {
      return out << "(node=" << w.node_to_osm_[n_] << ", level=" << lvl_
                 << ", dir=" << to_str(dir_)
                 << ", way=" << w.way_osm_idx_[w.r_->node_ways_[n_][way_]]
                 << ", type=" << node_type_to_str(type_) << ")";
    }

    static constexpr node invalid() noexcept { return node{}; }
    constexpr node_idx_t get_node() const noexcept { return n_; }
    constexpr node_idx_t get_key() const noexcept { return n_; }

    constexpr mode get_mode() const noexcept {
      return is_car_node() ? mode::kCar : mode::kFoot;
    }

    constexpr bool is_car_node() const noexcept {
      return type_ == node_type::kCar;
    }

    constexpr bool is_foot_node() const noexcept {
      return type_ == node_type::kFoot;
    }

    constexpr bool is_invalid_node() const noexcept {
      return type_ == node_type::kInvalid;
    }

    node_idx_t n_{node_idx_t::invalid()};
    node_type type_{node_type::kInvalid};
    level_t lvl_;
    direction dir_;
    way_pos_t way_;
  };

  struct label {
    label(node const n, cost_t const c)
        : n_{n.n_},
          cost_{c},
          type_{n.type_},
          lvl_{n.lvl_},
          dir_{n.dir_},
          way_(n.way_) {}

    constexpr node get_node() const noexcept {
      return {
          .n_ = n_, .type_ = type_, .lvl_ = lvl_, .dir_ = dir_, .way_ = way_};
    }

    constexpr cost_t cost() const noexcept { return cost_; }

    void track(
        label const&, ways::routing const&, way_idx_t, node_idx_t, bool) {}

    node_idx_t n_;
    cost_t cost_;
    node_type type_;
    level_t lvl_;
    direction dir_;
    way_pos_t way_;
  };

  struct entry {
    static constexpr auto const kMaxWays = way_pos_t{16U};
    static constexpr auto const kN = kMaxWays * 2U + 1 /* FWD+BWD + foot */;

    entry() { utl::fill(cost_, kInfeasible); }

    constexpr std::optional<node> pred(node const n) const noexcept {
      auto const idx = get_index(n);
      return pred_[idx] == node_idx_t::invalid()
                 ? std::nullopt
                 : std::optional{node{.n_ = pred_[idx],
                                      .type_ = to_node_type(pred_type_[idx]),
                                      .lvl_ = pred_lvl_[idx],
                                      .dir_ = to_dir(pred_dir_[idx]),
                                      .way_ = pred_way_[idx]}};
    }

    constexpr cost_t cost(node const n) const noexcept {
      return cost_[get_index(n)];
    }

    constexpr bool update(label const,
                          node const n,
                          cost_t const c,
                          node const pred) noexcept {
      auto const idx = get_index(n);
      if (c < cost_[idx]) {
        cost_[idx] = c;
        pred_[idx] = pred.n_;
        pred_lvl_[idx] = pred.lvl_;
        pred_type_[idx] = to_bool(pred.type_);
        pred_way_[idx] = pred.way_;
        pred_dir_[idx] = to_bool(pred.dir_);
        return true;
      }
      return false;
    }

    static constexpr std::size_t get_index(node const n) {
      return n.is_foot_node()
                 ? 0U
                 : 1U + (n.dir_ == direction::kForward ? 0U : 1U) * kMaxWays +
                       n.way_;
    }

    static constexpr direction to_dir(bool const b) {
      return b ? direction::kBackward : direction::kForward;
    }

    static constexpr bool to_bool(direction const d) {
      return d == direction::kBackward;
    }

    static constexpr node_type to_node_type(bool const b) {
      return b ? node_type::kFoot : node_type::kCar;
    }

    static constexpr bool to_bool(node_type const t) {
      return t == node_type::kFoot;
    }

    void write(node, path&) const {}

    std::array<node_idx_t, kN> pred_;
    std::array<cost_t, kN> cost_;
    std::array<way_pos_t, kN> pred_way_;
    std::array<level_t, kN> pred_lvl_;
    std::bitset<kN> pred_dir_;
    std::bitset<kN> pred_type_;
    std::bitset<kN> pred_parking_;
  };

  struct hash {
    using is_avalanching = void;
    auto operator()(key const n) const noexcept -> std::uint64_t {
      using namespace ankerl::unordered_dense::detail;
      return wyhash::hash(static_cast<std::uint64_t>(to_idx(n)));
    }
  };

  static car::node to_car(node const n) {
    return {.n_ = n.n_, .way_ = n.way_, .dir_ = n.dir_};
  }

  static footp::node to_foot(node const n) {
    return {.n_ = n.n_, .lvl_ = n.lvl_};
  }

  static node to_node(car::node const n, level_t const lvl) {
    return {.n_ = n.n_,
            .type_ = node_type::kCar,
            .lvl_ = lvl,
            .dir_ = n.dir_,
            .way_ = n.way_};
  }

  static node to_node(footp::node const n) {
    return {.n_ = n.n_,
            .type_ = node_type::kFoot,
            .lvl_ = n.lvl_,
            .dir_ = direction::kForward,
            .way_ = 0};
  }

  template <typename Fn>
  static void resolve_all(ways::routing const& w,
                          node_idx_t const n,
                          level_t const lvl,
                          Fn&& f) {
    footp::resolve_all(
        w, n, lvl, [&](footp::node const neighbor) { f(to_node(neighbor)); });
    car::resolve_all(w, n, lvl, [&](car::node const neighbor) {
      auto const p = w.way_properties_[w.node_ways_[n][neighbor.way_]];
      auto const node_level = lvl == kNoLevel ? p.from_level() : lvl;
      f(to_node(neighbor, node_level));
    });
  }

  template <direction SearchDir, bool WithBlocked, typename Fn>
  static void adjacent(ways::routing const& w,
                       node const n,
                       bitvec<node_idx_t> const* blocked,
                       sharing_data const*,
                       elevation_storage const* elevations,
                       Fn&& fn) {
    static constexpr auto const kFwd = SearchDir == direction::kForward;
    static constexpr auto const kBwd = SearchDir == direction::kBackward;

    auto const is_parking =
        !UseParking || w.node_properties_[n.n_].is_parking() ||
        utl::any_of(w.node_ways_[n.n_], [&](way_idx_t const way) {
          return w.way_properties_[way].is_parking();
        });

    if (n.is_foot_node() || (kFwd && n.is_car_node() && is_parking)) {
      footp::template adjacent<SearchDir, WithBlocked>(
          w, to_foot(n), blocked, nullptr, elevations,
          [&](footp::node const neighbor, std::uint32_t const cost,
              distance_t const dist, way_idx_t const way,
              std::uint16_t const from, std::uint16_t const to,
              elevation_storage::elevation const elevation, bool) {
            fn(to_node(neighbor),
               cost + (n.is_foot_node() ? 0 : kSwitchPenalty), dist, way, from,
               to, elevation, false);
          });
    }

    if (n.is_car_node() || (kBwd && n.is_foot_node() && is_parking)) {
      car::template adjacent<SearchDir, WithBlocked>(
          w, to_car(n), blocked, nullptr, elevations,
          [&](car::node const neighbor, std::uint32_t const cost,
              distance_t const dist, way_idx_t const way,
              std::uint16_t const from, std::uint16_t const to,
              elevation_storage::elevation const elevation, bool) {
            auto const way_prop = w.way_properties_[way];
            fn(to_node(neighbor, way_prop.from_level()),
               cost + (n.is_car_node() ? 0 : kSwitchPenalty), dist, way, from,
               to, elevation, false);
          });
    }
  }

  template <typename Fn>
  static void resolve_start_node(ways::routing const& w,
                                 way_idx_t const way,
                                 node_idx_t const n,
                                 level_t lvl,
                                 direction search_dir,
                                 Fn&& f) {
    auto const way_properties = w.way_properties_[way];
    search_dir == direction::kForward
        ? car::resolve_start_node(
              w, way, n, lvl, search_dir,
              [&](car::node const cn) {
                auto const node_level =
                    lvl == kNoLevel ? way_properties.from_level() : lvl;
                f(to_node(cn, node_level));
              })
        : footp::resolve_start_node(
              w, way, n, lvl, search_dir,
              [&](footp::node const fn) { f(to_node(fn)); });
  }

  static bool is_dest_reachable(ways::routing const& w,
                                node const n,
                                way_idx_t const way,
                                direction const way_dir,
                                direction const search_dir) {
    return !UseParking || w.way_properties_[way].is_parking() ||
           (search_dir == direction::kForward
                ? n.is_foot_node() &&
                      footp::is_dest_reachable(w, to_foot(n), way, way_dir,
                                               search_dir)
                : n.is_car_node() &&
                      car::is_dest_reachable(w, to_car(n), way, way_dir,
                                             search_dir));
  }

  static constexpr cost_t way_cost(way_properties const& e,
                                   direction const dir,
                                   std::uint16_t const dist) {
    return footp::way_cost(e, dir, dist);
  }

  static constexpr cost_t node_cost(node_properties const n) {
    return footp::node_cost(n);
  }

  static constexpr double heuristic(double dist) {
    return car::heuristic(dist);
  }

  static constexpr node get_reverse(node n) {
    return {n.n_, n.type_, n.lvl_, opposite(n.dir_), n.way_};
  }
};

}  // namespace osr
