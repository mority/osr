#pragma once

#include <optional>

#include "osr/routing/cch/cch.h"
#include "osr/routing/profile.h"
#include "osr/routing/profiles/car.h"
#include "osr/routing/profiles/common.h"
#include "osr/types.h"

namespace osr {

// `for_each_adjacent_node` selects the turn restriction set with a template
// parameter that is not part of the profile itself, so it has to be recovered
// here.
template <typename P>
struct cch_profile_traits {
  static constexpr auto const kIsBus = false;
};

template <>
struct cch_profile_traits<generic_car<true>> {
  static constexpr auto const kIsBus = true;
};

// Cost of the turn from the incoming port `in` to the outgoing port `out` at
// node `n`, `kInfeasible` if the turn is forbidden. This is exactly what
// `for_each_adjacent_node` adds on top of way cost and node cost, which is why
// turns never have to be materialized as nodes: shortcuts remember the ports of
// their first and last edge and the turn is evaluated on the fly.
template <WayAwareProfile P>
cost_t cch_turn_cost(typename P::parameters const& params,
                     ways::routing const& r,
                     timezone_cache_t const& timezones,
                     node_idx_t const n,
                     port_t const in,
                     port_t const out) {
  auto const from = port_way_pos(in);
  auto const to = port_way_pos(out);

  if (is_profile_turn_restricted<P, direction::kForward,
                                 cch_profile_traits<P>::kIsBus>(
          params, r, timezones, n, from, to, std::nullopt, duration_t{0},
          direction::kForward)) {
    return kInfeasible;
  }

  if (is_uturn(in, out)) {
    return params.uturn_penalty_;
  }

  return P::turn_cost(
      params, r.get_turn_angle(n, from, port_dir(in), to, port_dir(out)));
}

}  // namespace osr
