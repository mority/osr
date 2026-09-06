#include "osr/routing/cch/order.h"

#include <numeric>
#include <thread>

#include "utl/verify.h"

#include "inertialflowcutter/run.h"

#include <algorithm>
#include <array>
#include <limits>
#include <numeric>

#include "utl/progress_tracker.h"

namespace osr {


std::vector<std::uint32_t> compute_inertial_flow_cutter_order(
    nd_graph const& g, unsigned const n_threads) {
  auto pt = utl::get_active_progress_tracker_or_activate("osr");
  pt->status("CCH Order (IFC)").reset_bounds().out_bounds(2.0F, 6.0F).in_high(1);

  // Every undirected edge once: the cutter adds the reverse arcs itself.
  auto tail = std::vector<unsigned>{};
  auto head = std::vector<unsigned>{};
  tail.reserve(g.adj_.size() / 2U);
  head.reserve(g.adj_.size() / 2U);
  for (auto u = std::uint32_t{0U}; u != g.n(); ++u) {
    for (auto const v : g.adj(u)) {
      if (u < v) {
        tail.push_back(u);
        head.push_back(v);
      }
    }
  }

  // `0` means "all cores" everywhere else in the build, but TBB's global
  // control rejects a parallelism of zero, so resolve it here.
  auto const threads = n_threads == 0U
                           ? std::max(1U, std::thread::hardware_concurrency())
                           : n_threads;

  auto order = std::vector<std::uint32_t>(g.n());
  ifc::run_inertial_flow_cutter(
      static_cast<int>(threads), static_cast<int>(g.n()), tail, head,
      [&](int const i) {
        return std::pair<double, double>{
            static_cast<double>(g.x_[static_cast<std::size_t>(i)]),
            static_cast<double>(g.y_[static_cast<std::size_t>(i)])};
      },
      [&](int const node, int const position) {
        order[static_cast<std::size_t>(position)] =
            static_cast<std::uint32_t>(node);
      },
      // The graph above is simple by construction: `build_cch` sorts and
      // uniques the adjacency of every node and drops `v == u`, and the loop
      // here emits each undirected edge once. The cutter would otherwise strip
      // multi arcs and loops itself, which is two more passes over every arc.
      //
      // Four cutters instead of the default eight: measured on a planet build
      // that is 27:31 of ordering down to 19:30, against 1.3% on the median
      // query. Of every parameter on this path it buys the most preprocessing
      // time per unit of query cost, and the rate falls off steeply beyond it.
      ifc::options{.already_simple = true, .geo_pos_cutters = 4});
  pt->update_monotonic(1U);
  return order;
}

}  // namespace osr
