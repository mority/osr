#pragma once

#include <cstdint>

#include <span>
#include <vector>

namespace osr {

// Undirected graph the contraction order is computed on. Nodes are numbered
// 0..n-1 and carry a position, which the flow cutter uses to seed its cuts.
struct nd_graph {
  std::uint32_t n() const {
    return static_cast<std::uint32_t>(ofs_.size() - 1U);
  }

  std::span<std::uint32_t const> adj(std::uint32_t const v) const {
    return {&adj_[ofs_[v]], ofs_[v + 1U] - ofs_[v]};
  }

  std::vector<std::uint64_t> ofs_;
  std::vector<std::uint32_t> adj_;
  std::vector<std::int32_t> x_, y_;
};

// The contraction order: the returned vector holds the nodes in the order they
// are contracted in, i.e. index == rank, with separators on top of the
// hierarchy. InertialFlowCutter produces this directly -- it reports the
// position of every node in the order it computes.
std::vector<std::uint32_t> compute_inertial_flow_cutter_order(nd_graph const&,
                                                              unsigned n_threads);

}  // namespace osr
