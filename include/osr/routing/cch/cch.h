#pragma once

#include <cinttypes>

#include <algorithm>
#include <filesystem>
#include <limits>
#include <span>

#include "cista/memory_holder.h"

#include "utl/zip.h"

#include "osr/types.h"
#include "osr/ways.h"

namespace osr {

// Rank of a node in the contraction order. Nodes are contracted in increasing
// rank order, i.e. rank 0 is contracted first and therefore is the "lowest"
// node of the hierarchy.
using cch_rank_t = cista::strong<std::uint32_t, struct cch_rank_>;

// Index of an entry in the up-/downward arc entry arrays. 64 bit: a planet
// sized hierarchy passes two billion entries during contraction.
using cch_entry_idx_t = std::uint64_t;

// Index of an adjacency slot, i.e. of an undirected pair {lower, higher}.
using cch_slot_idx_t = std::uint32_t;

// A "port" identifies one of the directed original edges incident to a node:
// the way slot (position in `node_ways_[n]`) plus the direction in which that
// way is traversed. Turn costs and turn restrictions are a function of the
// incoming and the outgoing port at a node, which is exactly what the routing
// state `(node, way_pos, dir)` of the way aware profiles encodes.
using port_t = std::uint8_t;

constexpr auto const kMaxPorts = port_t{32U};

constexpr port_t make_port(way_pos_t const way_pos, direction const dir) {
  return static_cast<port_t>((way_pos << 1U) |
                             (dir == direction::kBackward ? 1U : 0U));
}

constexpr way_pos_t port_way_pos(port_t const p) {
  return static_cast<way_pos_t>(p >> 1U);
}

constexpr direction port_dir(port_t const p) {
  return (p & 1U) != 0U ? direction::kBackward : direction::kForward;
}

// A u-turn keeps the way slot and reverses the direction of travel, so it flips
// exactly the direction bit of the port.
constexpr bool is_uturn(port_t const in, port_t const out) {
  return (in ^ out) == 1U;
}

constexpr port_t n_ports(ways::routing const& r, node_idx_t const n) {
  return static_cast<port_t>(r.node_ways_[n].size() << 1U);
}

// One (directed) arc of the contraction hierarchy. The pair of nodes it
// connects is given by the adjacency slot it belongs to, the two ports say
// which original edge the represented path starts / ends with. This is how turn
// restrictions and turn costs stay implicit: instead of blowing up the graph to
// one node per turn we remember with every shortcut which turn we leave the
// tail node with and which turn we arrive at the head node with.
struct cch_entry {
  friend bool operator==(cch_entry, cch_entry) = default;
  friend auto operator<=>(cch_entry, cch_entry) = default;

  port_t entry_;  // port of the first edge at the tail node
  port_t exit_;  // port of the last edge at the head node
};

static_assert(sizeof(cch_entry) == 2U);

// Metric independent part of the customizable contraction hierarchy.
//
// Layout: every undirected pair {l, h} with rank(l) < rank(h) is stored exactly
// once, in the adjacency list of the lower node `l` ("slot"). Both search
// directions of a CH query need exactly these slots:
//   * the forward search relaxes upward arcs l -> h,
//   * the backward search relaxes downward arcs h -> l (traversed backwards).
// Per slot there are two lists of entries: `up` for l -> h and `dn` for h -> l.
struct cch {
  static constexpr auto const kMode =
      cista::mode::WITH_INTEGRITY | cista::mode::WITH_STATIC_VERSION;

  static cista::wrapped<cch> read(std::filesystem::path const&);
  void write(std::filesystem::path const&) const;
  static bool exists(std::filesystem::path const&);

  cch_rank_t::value_t n_ranks() const {
    return static_cast<cch_rank_t::value_t>(order_.size());
  }

  cch_slot_idx_t n_slots() const {
    return static_cast<cch_slot_idx_t>(adj_head_.size());
  }

  bool contains(node_idx_t const n) const {
    return n < rank_.size() && rank_[n] != cch_rank_t::invalid();
  }

  std::span<cch_rank_t const> upper(cch_rank_t const r) const {
    return {&adj_head_[adj_ofs_[r]], adj_ofs_[r + 1U] - adj_ofs_[r]};
  }

  cch_slot_idx_t upper_begin(cch_rank_t const r) const { return adj_ofs_[r]; }
  cch_slot_idx_t upper_end(cch_rank_t const r) const {
    return adj_ofs_[r + 1U];
  }

  std::span<std::uint32_t const> lower_slots(cch_rank_t const r) const {
    return {&lower_slot_[lower_ofs_[r]], lower_ofs_[r + 1U] - lower_ofs_[r]};
  }

  std::span<cch_entry const> up_entries(cch_slot_idx_t const s) const {
    return {&up_[up_ofs_[s]], up_ofs_[s + 1U] - up_ofs_[s]};
  }

  std::span<cch_entry const> dn_entries(cch_slot_idx_t const s) const {
    return {&dn_[dn_ofs_[s]], dn_ofs_[s + 1U] - dn_ofs_[s]};
  }

  // Self loops: a path that leaves a node and comes back to it. They are the
  // shortcuts of "turn around at a lower ranked node", which is the only way
  // through an intersection whose direct turn is forbidden. Without them the
  // hierarchy would not preserve all shortest paths.
  std::span<cch_entry const> loop_entries(cch_rank_t const r) const {
    return {&loop_[loop_ofs_[r]], loop_ofs_[r + 1U] - loop_ofs_[r]};
  }

  cch_entry_idx_t loop_begin(cch_rank_t const r) const { return loop_ofs_[r]; }

  // Finds the slot of the pair {l, h} with rank(l) < rank(h). Returns
  // `kNoSlot` if the two nodes are not adjacent in the hierarchy.
  static constexpr auto const kNoSlot =
      std::numeric_limits<cch_slot_idx_t>::max();
  // marks a predecessor arc that is a self loop instead of an adjacency slot
  static constexpr auto const kLoopSlot = kNoSlot - 1U;

  // Lower node of a slot. `adj_ofs_` is the CSR over the slots indexed by the
  // lower node, so the tail is the last rank whose range starts at or below
  // the slot. Deriving it costs a binary search but saves 4 bytes per slot,
  // which is 566 MB on Germany.
  cch_rank_t tail(cch_slot_idx_t const s) const {
    auto const it = std::upper_bound(begin(adj_ofs_), end(adj_ofs_), s);
    return cch_rank_t{static_cast<cch_rank_t::value_t>(
        std::distance(begin(adj_ofs_), it) - 1)};
  }

  cch_slot_idx_t find_slot(cch_rank_t const l, cch_rank_t const h) const {
    auto const from = begin(adj_head_) + adj_ofs_[l];
    auto const to = begin(adj_head_) + adj_ofs_[l + 1U];
    auto const it = std::lower_bound(from, to, h);
    return (it == to || *it != h)
               ? kNoSlot
               : static_cast<cch_slot_idx_t>(std::distance(
                     begin(adj_head_), it));
  }

  // Finds the index of the entry (entry_port, exit_port) within the up-/down
  // entry list of a slot. Entries are sorted, so this is a binary search.
  static cch_entry_idx_t find_entry(std::span<cch_entry const> entries,
                                    cch_entry const e) {
    auto const it = std::lower_bound(begin(entries), end(entries), e);
    return (it == end(entries) || *it != e)
               ? std::numeric_limits<cch_entry_idx_t>::max()
               : static_cast<cch_entry_idx_t>(
                     std::distance(begin(entries), it));
  }

  // node -> rank (invalid for nodes that are not part of the hierarchy)
  vec_map<node_idx_t, cch_rank_t> rank_;
  // rank -> node
  vec_map<cch_rank_t, node_idx_t> order_;

  // CSR over the slots, indexed by the *lower* node of each pair.
  vec_map<cch_rank_t, cch_slot_idx_t> adj_ofs_;  // size: n_ranks + 1
  vec<cch_rank_t> adj_head_;  // higher ranked node of the slot
  // the lower ranked node is not stored, see `tail()`

  // CSR over the entries of every slot.
  vec<cch_entry_idx_t> up_ofs_;  // size: n_slots + 1
  vec<cch_entry_idx_t> dn_ofs_;  // size: n_slots + 1
  vec64<cch_entry> up_;
  vec64<cch_entry> dn_;

  vec_map<cch_rank_t, cch_entry_idx_t> loop_ofs_;  // size: n_ranks + 1
  vec64<cch_entry> loop_;

  // Entries that an original edge maps to. Shortcut entries never have to be
  // checked against the original graph while unpacking.
  bitvec64 up_is_edge_, dn_is_edge_, loop_is_edge_;

  // Transposed adjacency: for every node the slots in which it is the *higher*
  // node. Needed to enumerate lower triangles during customization / unpacking.
  vec_map<cch_rank_t, cch_slot_idx_t> lower_ofs_;  // size: n_ranks + 1
  vec<cch_slot_idx_t> lower_slot_;
};

// Enumerates all original directed edges leaving `u`.
//
// fn(target_node, way, tail_port, head_port, distance, from_idx, to_idx)
template <typename Fn>
void for_each_edge(ways::routing const& r, node_idx_t const u, Fn&& fn) {
  auto way_pos = way_pos_t{0U};
  for (auto const [way, i] :
       utl::zip_unchecked(r.node_ways_[u], r.node_in_way_idx_[u])) {
    auto const expand = [&](direction const way_dir, std::uint16_t const from,
                            std::uint16_t const to) {
      auto const v = r.way_nodes_[way][to];
      fn(v, way, make_port(way_pos, way_dir),
         make_port(r.get_way_pos(v, way, to), way_dir),
         r.get_way_node_distance(way, std::min(from, to)), from, to);
    };

    if (i != 0U) {
      expand(direction::kBackward, i, static_cast<std::uint16_t>(i - 1U));
    }
    if (i != r.way_nodes_[way].size() - 1U) {
      expand(direction::kForward, i, static_cast<std::uint16_t>(i + 1U));
    }

    ++way_pos;
  }
}

}  // namespace osr
