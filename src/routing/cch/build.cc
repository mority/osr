#include "osr/routing/cch/build.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <condition_variable>
#include <memory>
#include <mutex>
#include <thread>
#include <span>
#include <vector>

#include "utl/helpers/algorithm.h"
#include "utl/parallel_for.h"
#include "utl/progress_tracker.h"
#include "utl/timer.h"
#include "utl/verify.h"

#include "osr/routing/cch/order.h"

namespace osr {

bool is_cch_way(way_properties const& p) {
  return p.is_car_accessible() || p.is_bus_accessible() ||
         p.is_bus_accessible_with_penalty();
}

bool is_restricted_for_all(ways::routing const& r,
                           node_idx_t const n,
                           way_pos_t const from,
                           way_pos_t const to) {
  if (!r.node_is_restricted_[n]) {
    return false;
  }
  return utl::any_of(r.node_restrictions_[n], [&](restriction const& x) {
    return x.from_ == from && x.to_ == to && x.applies_to_default_ &&
           x.applies_to_bus_ && x.applies_to_hgv_ &&
           x.condition_set_ == conditional_condition_set_idx_t::invalid();
  });
}

namespace {

// One entry of the growing adjacency during contraction. Every arc is stored
// exactly once, namely at its lower ranked end point: when the higher ranked
// end point is contracted the arc is already gone.
struct arc_rec {
  // Sorting these records is the hottest part of the contraction, so the
  // members are packed into one integer instead of being compared one by one.
  // The fields do not overlap and sit in declaration order, so comparing the
  // keys is exactly the lexicographic order of the members.
  constexpr std::uint64_t key() const {
    return (static_cast<std::uint64_t>(other_) << 24U) |
           (static_cast<std::uint64_t>(out_) << 16U) |
           (static_cast<std::uint64_t>(p_) << 8U) |
           static_cast<std::uint64_t>(q_);
  }

  friend bool operator==(arc_rec const a, arc_rec const b) {
    return a.key() == b.key();
  }
  friend auto operator<=>(arc_rec const a, arc_rec const b) {
    return a.key() <=> b.key();
  }

  cch_rank_t::value_t other_;
  std::uint8_t out_;  // 1: this -> other, 0: other -> this
  port_t p_;  // entry port at the tail node
  port_t q_;  // exit port at the head node
};

// One neighbor's entries, grouped by the port at the *far* node.
//
// A record the contraction pushes only keeps the two far ports -- the ports at
// the node being contracted are consumed by the turn check and then thrown
// away. Enumerating the entry lists pairwise therefore produces the same
// record over and over: on Germany 27.6 billion records collapse to 3.0
// billion distinct ones, and the nine duplicates out of ten are pushed,
// radix sorted, merged and dropped again. Grouping first turns the inner
// product into one record per pair of far ports.
//
// `mask_` holds the ports at the contracted node that go with `port_`: for an
// outgoing group the departure ports as they are, for an incoming group the
// *closure* of the arrival ports under `allowed_`, so that testing a pair is a
// single and. `wide_` marks a group that contains a port beyond `kMaxPorts`,
// which has no bit in either mask -- such a group is treated as reachable,
// which can only add arcs, never lose one.
struct port_group {
  port_t port_;
  bool wide_;
  std::uint32_t mask_;
};

// Everything the contraction of one node needs besides its own arc list. One
// per thread on a wide level, one per node on a narrow one -- the buffers are
// reused, so a node pays no allocation for them.
struct node_scratch {
  std::vector<cch_rank_t::value_t> neighbors_;
  std::vector<std::vector<cch_entry>> in_of_, out_of_;
  std::vector<port_group> in_grp_, out_grp_;
  std::vector<std::uint32_t> in_ofs_, out_ofs_;

  // Scratch for grouping an outgoing entry list by its far port: those lists
  // are sorted by the *near* port, so the groups are not contiguous. Indexed
  // by the full `port_t` range so that a node with more than `kMaxPorts / 2`
  // ways cannot run off the end. Bit 32 marks a group as wide, bit 33 marks it
  // as used, so that a group whose mask is empty is still found.
  static constexpr auto const kSeenBit = std::uint64_t{1U} << 33U;
  static constexpr auto const kWideBit = std::uint64_t{1U} << 32U;
  std::array<std::uint64_t, 256U> acc_{};
  std::vector<port_t> touched_;

  std::array<std::uint32_t, kMaxPorts> reach_{};
  std::array<std::uint32_t, kMaxPorts> allowed_{};
  bool allowed_all_{false};

  std::span<port_group const> in_group(std::size_t const k) const {
    return {in_grp_.data() + in_ofs_[k], in_ofs_[k + 1U] - in_ofs_[k]};
  }
  std::span<port_group const> out_group(std::size_t const k) const {
    return {out_grp_.data() + out_ofs_[k], out_ofs_[k + 1U] - out_ofs_[k]};
  }
};

// The finished row of one node, in the layout `write_csr` wants it: a header
// per slot, then that node's up entries and dn entries, both grouped by slot.
//
// Contracting by level means a node's arcs have to outlive the node, because
// the CSR comes out in rank order and a level is not a range of ranks. Keeping
// the raw `arc_rec`s for that costs eight bytes a record plus a vector header a
// node -- 54 GB on the planet, which is more memory than is left once the order
// has run. The same content re-encoded is two bytes an entry in one allocation.
struct slot_head {
  cch_rank_t::value_t head_;
  std::uint16_t n_up_, n_dn_;
};

static_assert(sizeof(slot_head) == 8U);

// Stripe lock over the target lists. A target is written by one node at a
// time, but two nodes of the same level can share one, so `process_target`
// holds the stripe of its target for the whole call -- one acquisition per
// (node, neighbor) pair rather than one per record.
constexpr auto const kLockStripes = std::size_t{4096U};

// Workers for the contraction, created once instead of per node.
//
// The parallel region is entered once per dense node, and a node's work is
// far too small to pay for creating and joining threads every time: at 32
// threads that overhead was large enough to make the contraction slower than
// at 8. The workers therefore stay alive for the whole contraction and are
// handed one node's targets at a time.
struct worker_pool {
  explicit worker_pool(unsigned const n) {
    threads_.reserve(n);
    for (auto i = 0U; i != n; ++i) {
      threads_.emplace_back([this]() { loop(); });
    }
  }

  worker_pool(worker_pool const&) = delete;
  worker_pool& operator=(worker_pool const&) = delete;
  worker_pool(worker_pool&&) = delete;
  worker_pool& operator=(worker_pool&&) = delete;

  ~worker_pool() {
    {
      auto const lock = std::lock_guard{m_};
      stop_ = true;
      ++generation_;
    }
    cv_.notify_all();
    for (auto& t : threads_) {
      t.join();
    }
  }

  // Calls `job(i)` for every i < count and returns when all of them are done.
  // The calling thread takes part, so the pool holds one worker less than the
  // requested thread count -- the caller would only block otherwise.
  template <typename Fn>
  void run(std::size_t const count, Fn const& job) {
    {
      auto const lock = std::lock_guard{m_};
      job_ = [](void const* ctx, std::size_t const i) {
        (*static_cast<Fn const*>(ctx))(i);
      };
      job_ctx_ = &job;
      count_ = count;
      next_.store(0U, std::memory_order_relaxed);
      pending_ = static_cast<unsigned>(threads_.size());
      ++generation_;
    }
    cv_.notify_all();

    work();

    auto lock = std::unique_lock{m_};
    done_cv_.wait(lock, [this]() { return pending_ == 0U; });
  }

private:
  void work() {
    for (auto i = next_.fetch_add(1U, std::memory_order_relaxed); i < count_;
         i = next_.fetch_add(1U, std::memory_order_relaxed)) {
      job_(job_ctx_, i);
    }
  }

  void loop() {
    auto seen = std::uint64_t{0U};
    for (;;) {
      {
        auto lock = std::unique_lock{m_};
        cv_.wait(lock, [&]() { return stop_ || generation_ != seen; });
        seen = generation_;
        if (stop_) {
          return;
        }
      }

      work();

      {
        auto const lock = std::lock_guard{m_};
        if (--pending_ == 0U) {
          done_cv_.notify_one();
        }
      }
    }
  }

  std::mutex m_;
  std::condition_variable cv_, done_cv_;
  std::vector<std::thread> threads_;

  // written under `m_` before the generation is bumped, read by the workers
  // after they wake up, i.e. after they acquired `m_`
  void (*job_)(void const*, std::size_t){nullptr};
  void const* job_ctx_{nullptr};
  std::size_t count_{0U};

  std::atomic<std::size_t> next_{0U};
  std::uint64_t generation_{0U};
  unsigned pending_{0U};
  bool stop_{false};
};

// Groups the ranks by their level in the elimination tree.
//
// The parent of a node is the lowest ranked of its upper neighbours in the
// *filled* graph, which Liu's algorithm derives from the original graph alone:
// walking the ranks in order and merging every lower neighbour's component
// into the current node makes the root of a component the node that will
// absorb it. That needs no fill-in, only a union find with path compression.
//
// The level of a node is its height in that tree. A node's lower neighbours in
// the filled graph are exactly the tree descendants it is adjacent to, so
// their height is smaller -- the height of a node is one more than the largest
// height among its lower neighbours, which is the same level `customize()`
// computes on the finished hierarchy.
std::vector<std::vector<cch_rank_t::value_t>> compute_levels(
    nd_graph const& g, std::vector<std::uint32_t> const& order) {
  auto const n = static_cast<cch_rank_t::value_t>(order.size());
  constexpr auto const kNoParent = std::numeric_limits<std::uint32_t>::max();

  auto rank_of = std::vector<std::uint32_t>(n);
  for (auto i = cch_rank_t::value_t{0U}; i != n; ++i) {
    rank_of[order[i]] = i;
  }

  auto ancestor = std::vector<std::uint32_t>(n);
  auto parent = std::vector<std::uint32_t>(n, kNoParent);
  for (auto i = cch_rank_t::value_t{0U}; i != n; ++i) {
    ancestor[i] = i;
    for (auto const l : g.adj(order[i])) {
      auto const j = rank_of[l];
      if (j >= i) {
        continue;
      }
      auto root = j;
      while (ancestor[root] != root) {
        root = ancestor[root];
      }
      for (auto x = j; ancestor[x] != root;) {
        auto const next = ancestor[x];
        ancestor[x] = root;
        x = next;
      }
      if (root != i) {
        parent[root] = i;
        ancestor[root] = i;
      }
    }
  }

  // Pushed forward rather than pulled: the parent of a node has a larger rank,
  // so walking the ranks in order reaches a node only after all of its
  // children are final.
  auto height = std::vector<std::uint32_t>(n, 0U);
  auto max_height = std::uint32_t{0U};
  for (auto i = cch_rank_t::value_t{0U}; i != n; ++i) {
    max_height = std::max(max_height, height[i]);
    if (parent[i] != kNoParent) {
      height[parent[i]] = std::max(height[parent[i]], height[i] + 1U);
    }
  }

  auto by_level =
      std::vector<std::vector<cch_rank_t::value_t>>(max_height + 1U);
  auto cnt = std::vector<std::uint32_t>(max_height + 1U, 0U);
  for (auto const h : height) {
    ++cnt[h];
  }
  for (auto h = std::uint32_t{0U}; h != by_level.size(); ++h) {
    by_level[h].reserve(cnt[h]);
  }
  for (auto i = cch_rank_t::value_t{0U}; i != n; ++i) {
    by_level[height[i]].emplace_back(i);
  }
  return by_level;
}

struct contractor {
  contractor(ways const& w, cch& c, unsigned const n_threads)
      : r_{*w.r_},
        c_{c},
        n_threads_{n_threads == 0U
                       ? std::max(1U, std::thread::hardware_concurrency())
                       : n_threads},
        pool_{n_threads_ - 1U} {}

  bool turn_ok(node_idx_t const n, port_t const in, port_t const out) const {
    return !is_restricted_for_all(r_, n, port_way_pos(in), port_way_pos(out));
  }

  template <typename T>
  void dedup(std::vector<T>& v) const {
    std::sort(begin(v), end(v));
    v.erase(std::unique(begin(v), end(v)), end(v));
  }

  // Below this an LSD pass costs more than it saves.
  static constexpr auto const kRadixMin = std::size_t{2048U};

  // Number of byte digits of `arc_rec::key()` (56 significant bits).
  static constexpr auto const kRadixPasses = 7U;

  static std::uint8_t digit(arc_rec const a, unsigned const pass) {
    return static_cast<std::uint8_t>((a.key() >> (pass * 8U)) & 0xFFU);
  }

  // LSD radix sort of `v[from, end)` by the packed key. The key is an integer,
  // so the records do not have to be compared at all: seven counting passes
  // put them in order. Passes whose digit is the same for every record are
  // skipped, which is common for the high bytes of the neighbor id.
  void radix_sort(std::vector<arc_rec>& v, std::size_t const from) {
    auto const n = v.size() - from;

    auto hist = std::array<std::array<std::uint32_t, 256U>, kRadixPasses>{};
    for (auto i = from; i != v.size(); ++i) {
      for (auto p = 0U; p != kRadixPasses; ++p) {
        ++hist[p][digit(v[i], p)];
      }
    }

    static thread_local auto tmp = std::vector<arc_rec>{};
    tmp.resize(n);
    auto* src = v.data() + from;
    auto* dst = tmp.data();

    for (auto p = 0U; p != kRadixPasses; ++p) {
      if (hist[p][digit(src[0], p)] == n) {
        continue;  // constant digit, nothing to do
      }
      auto ofs = std::array<std::uint32_t, 256U>{};
      auto sum = std::uint32_t{0U};
      for (auto b = std::size_t{0U}; b != 256U; ++b) {
        ofs[b] = sum;
        sum += hist[p][b];
      }
      for (auto i = std::size_t{0U}; i != n; ++i) {
        dst[ofs[digit(src[i], p)]++] = src[i];
      }
      std::swap(src, dst);
    }

    if (src != v.data() + from) {
      std::copy(src, src + n, v.data() + from);
    }
  }

  // Same as `dedup`, but for a vector whose first `clean` elements are already
  // sorted and unique: only the appended tail has to be ordered, the rest is a
  // linear merge. Re-sorting the whole vector every time it doubles is what
  // used to dominate the build.
  void dedup_tail(std::vector<arc_rec>& v, std::size_t const clean) {
    auto const n = v.size() - clean;
    if (n == 0U) {
      return;
    }
    if (n < kRadixMin) {
      std::sort(begin(v) + static_cast<std::ptrdiff_t>(clean), end(v));
    } else {
      radix_sort(v, clean);
    }
    if (clean != 0U) {
      std::inplace_merge(
          begin(v), begin(v) + static_cast<std::ptrdiff_t>(clean), end(v));
    }
    v.erase(std::unique(begin(v), end(v)), end(v));
  }

  void push_loop(cch_rank_t::value_t const at, cch_entry const e) {
    auto& v = loops_[at];
    v.emplace_back(e);
    if (v.size() > 64U) {
      dedup(v);
    }
  }

  // Reachability between the *arriving* ports of `n` using its self loops:
  // bit q of reach[i] is set if we can arrive at `n` with port q after having
  // arrived with port i and driven one or more loops.
  void build_reach(node_scratch& s,
                   node_idx_t const n,
                   std::span<cch_entry const> loops,
                   port_t const ports) {
    auto& reach_ = s.reach_;
    reach_.fill(0U);
    for (auto i = port_t{0U}; i != ports; ++i) {
      reach_[i] = std::uint32_t{1U} << i;
    }
    if (loops.empty()) {
      return;
    }

    auto step = std::array<std::uint32_t, kMaxPorts>{};
    for (auto const& e : loops) {
      for (auto i = port_t{0U}; i != ports; ++i) {
        if (turn_ok(n, i, e.entry_) && e.exit_ < ports) {
          step[i] |= std::uint32_t{1U} << e.exit_;
        }
      }
    }

    for (auto changed = true; changed;) {
      changed = false;
      for (auto i = port_t{0U}; i != ports; ++i) {
        auto next = reach_[i];
        auto rest = reach_[i];
        while (rest != 0U) {
          auto const k = static_cast<port_t>(std::countr_zero(rest));
          rest &= rest - 1U;
          next |= step[k];
        }
        if (next != reach_[i]) {
          reach_[i] = next;
          changed = true;
        }
      }
    }
  }

  // `allowed_[i]`: the departure ports that can be reached from arrival port
  // `i`, driving zero or more of the node's self loops on the way. This is the
  // per pair `turn_reachable` test of the old contraction hoisted out of the
  // record loop: it is evaluated once per node instead of 27.6 billion times.
  //
  // A node without a turn restriction permits every turn, so its table is all
  // ones and neither the loop closure nor `is_restricted_for_all` is needed --
  // that is 99.2% of the nodes.
  void build_allowed(node_scratch& s,
                     node_idx_t const n,
                     std::span<cch_entry const> loops,
                     port_t const ports) {
    auto& allowed_ = s.allowed_;
    if (!r_.node_is_restricted_[n]) {
      if (!s.allowed_all_) {
        allowed_.fill(std::numeric_limits<std::uint32_t>::max());
        s.allowed_all_ = true;
      }
      return;
    }
    s.allowed_all_ = false;

    build_reach(s, n, loops, ports);

    auto ok = std::array<std::uint32_t, kMaxPorts>{};
    for (auto k = port_t{0U}; k != ports; ++k) {
      for (auto o = port_t{0U}; o != ports; ++o) {
        if (turn_ok(n, k, o)) {
          ok[k] |= std::uint32_t{1U} << o;
        }
      }
    }

    allowed_.fill(0U);
    for (auto i = port_t{0U}; i != ports; ++i) {
      auto rest = s.reach_[i];
      while (rest != 0U) {
        auto const k = static_cast<port_t>(std::countr_zero(rest));
        rest &= rest - 1U;
        allowed_[i] |= ok[k];
      }
    }
  }

  static bool reachable(port_group const& in, port_group const& out) {
    return in.wide_ || out.wide_ || (in.mask_ & out.mask_) != 0U;
  }

  void push(cch_rank_t::value_t const at, arc_rec const a) {
    auto& v = adj_[at];
    v.emplace_back(a);
    auto& clean = clean_size_[at];
    if (v.size() > std::max(std::size_t{32U}, 2U * clean)) {
      dedup_tail(v, clean);
      clean = v.size();
    }
  }

  // Prepares one node: makes its arc list final and groups it by neighbor.
  // No global state is touched, so every node of a level can run this at the
  // same time -- a node is never a neighbor of another node of its own level,
  // so nothing can still be pushed into its list.
  void prepare(cch_rank_t::value_t const i, node_scratch& s) {
    auto const node = c_.order_[cch_rank_t{i}];

    auto& v = adj_[i];
    dedup_tail(v, clean_size_[i]);
    clean_size_[i] = v.size();

    auto& lp = loops_[i];
    dedup(lp);
    auto const ports = n_capped_ports(r_, node);
    build_allowed(s, node, std::span<cch_entry const>{lp}, ports);

    // group the arc records by neighbor. The inner lists are reused instead
    // of rebuilt: at 15 million nodes their allocation is measurable.
    //
    // The row this node contributes to the CSR is counted on the way: the
    // records are already being walked here, and `write_csr` would otherwise
    // have to walk all 1.8 GB of them once more just to size its arrays.
    s.neighbors_.clear();
    auto n_up = cch_entry_idx_t{0U};
    auto n_dn = cch_entry_idx_t{0U};
    auto k = std::size_t{0U};
    for (auto it = begin(v); it != end(v); ++k) {
      auto const other = it->other_;
      s.neighbors_.emplace_back(other);
      if (s.in_of_.size() == k) {
        s.in_of_.emplace_back();
        s.out_of_.emplace_back();
      }
      auto& in = s.in_of_[k];
      auto& out = s.out_of_[k];
      in.clear();
      out.clear();
      for (; it != end(v) && it->other_ == other; ++it) {
        (it->out_ != 0U ? out : in).emplace_back(cch_entry{it->p_, it->q_});
      }
      n_up += out.size();
      n_dn += in.size();
    }
    auto const rank = cch_rank_t{i};
    c_.adj_ofs_[rank + 1U] = static_cast<cch_slot_idx_t>(k);
    c_.loop_ofs_[rank + 1U] = lp.size();
    up_base_[i + 1U] = n_up;
    dn_base_[i + 1U] = n_dn;

    // The arc records are dead from here on -- everything the contraction still
    // needs sits in the scratch, and everything the layout needs is the row.
    // Writing it now is what keeps the retained set to the row rather than to
    // the records it was built from.
    auto const bytes = k * sizeof(slot_head) +
                       static_cast<std::size_t>(n_up + n_dn) * sizeof(cch_entry);
    auto buf = std::make_unique_for_overwrite<std::byte[]>(bytes);
    auto* const heads = reinterpret_cast<slot_head*>(buf.get());
    auto* up = reinterpret_cast<cch_entry*>(buf.get() + k * sizeof(slot_head));
    auto* dn = up + n_up;
    for (auto j = std::size_t{0U}; j != k; ++j) {
      auto const& out = s.out_of_[j];
      auto const& in = s.in_of_[j];
      heads[j] = slot_head{s.neighbors_[j],
                           static_cast<std::uint16_t>(out.size()),
                           static_cast<std::uint16_t>(in.size())};
      up = std::copy(begin(out), end(out), up);
      dn = std::copy(begin(in), end(in), dn);
    }
    row_[i] = std::move(buf);

    v.clear();
    v.shrink_to_fit();

    // group both entry lists of every neighbor by their far port
    s.in_grp_.clear();
    s.out_grp_.clear();
    s.in_ofs_.clear();
    s.out_ofs_.clear();
    for (auto j = std::size_t{0U}; j != s.neighbors_.size(); ++j) {
      s.in_ofs_.emplace_back(static_cast<std::uint32_t>(s.in_grp_.size()));
      s.out_ofs_.emplace_back(static_cast<std::uint32_t>(s.out_grp_.size()));

      // incoming: sorted by the far port, so the groups are runs. The mask is
      // the closure of the arrival ports, which makes the pair test an and
      // against the departure ports of an outgoing group.
      auto const& in = s.in_of_[j];
      for (auto it = begin(in); it != end(in);) {
        auto const far = it->entry_;
        auto mask = std::uint32_t{0U};
        auto wide = false;
        for (; it != end(in) && it->entry_ == far; ++it) {
          if (it->exit_ < kMaxPorts) {
            mask |= s.allowed_[it->exit_];
          } else {
            wide = true;
          }
        }
        s.in_grp_.emplace_back(port_group{far, wide, mask});
      }

      // outgoing: sorted by the near port, so the groups are collected in
      // `acc_` and read back in the order they were first seen
      s.touched_.clear();
      for (auto const& e : s.out_of_[j]) {
        auto& a = s.acc_[e.exit_];
        if (a == 0U) {
          s.touched_.emplace_back(e.exit_);
        }
        a |= node_scratch::kSeenBit |
             (e.entry_ < kMaxPorts ? std::uint64_t{1U} << e.entry_
                                   : node_scratch::kWideBit);
      }
      for (auto const far : s.touched_) {
        auto const a = s.acc_[far];
        s.acc_[far] = 0U;
        s.out_grp_.emplace_back(
            port_group{far, (a & node_scratch::kWideBit) != 0U,
                       static_cast<std::uint32_t>(a & 0xFFFFFFFFU)});
      }
    }
    s.in_ofs_.emplace_back(static_cast<std::uint32_t>(s.in_grp_.size()));
    s.out_ofs_.emplace_back(static_cast<std::uint32_t>(s.out_grp_.size()));
  }

  // Connects the neighbors of a contracted node that meet at `k`. `k == b`
  // creates a self loop, i.e. the shortcut of turning around at this node.
  //
  // Every arc is pushed to the lower ranked of its two ends, so the pairs are
  // enumerated per *target* instead of as a full matrix: the neighbor `k`
  // collects the pairs whose lower end it is. All of a call's writes therefore
  // go to one target list, which one stripe lock covers.
  void process_target(node_scratch const& s, std::size_t const k) {
    auto const& neighbors = s.neighbors_;
    auto const t = neighbors[k];
    auto const in_k = s.in_group(k);
    auto const out_k = s.out_group(k);
    if (in_k.empty() && out_k.empty()) {
      return;
    }

    auto const lock = std::lock_guard{target_locks_[t % kLockStripes]};

    // pairs (k, b): this node is the tail, so it owns them if t < y
    if (!in_k.empty()) {
      for (auto b = std::size_t{0U}; b != neighbors.size(); ++b) {
        auto const out_b = s.out_group(b);
        if (out_b.empty()) {
          continue;
        }
        auto const y = neighbors[b];
        if (b != k && t >= y) {
          continue;  // owned by `b`
        }
        for (auto const& g1 : in_k) {
          for (auto const& g2 : out_b) {
            if (!reachable(g1, g2)) {
              continue;
            }
            if (b == k) {
              push_loop(t, cch_entry{g1.port_, g2.port_});
            } else {
              push(t, arc_rec{y, 1U, g1.port_, g2.port_});
            }
          }
        }
      }
    }

    // pairs (a, k): this node is the head, so it owns them if t < x
    if (!out_k.empty()) {
      for (auto a = std::size_t{0U}; a != neighbors.size(); ++a) {
        if (a == k) {
          continue;
        }
        auto const in_a = s.in_group(a);
        if (in_a.empty()) {
          continue;
        }
        auto const x = neighbors[a];
        if (x <= t) {
          continue;  // owned by `a`
        }
        for (auto const& g1 : in_a) {
          for (auto const& g2 : out_k) {
            if (!reachable(g1, g2)) {
              continue;
            }
            push(t, arc_rec{x, 0U, g1.port_, g2.port_});
          }
        }
      }
    }
  }

  // Contracts the hierarchy level by level.
  //
  // The contraction order is a valid order, but so is any order that finishes
  // a level before it starts the next one: a node only ever pushes into the
  // list of a node on a larger level, so within a level the nodes are
  // independent of each other. That is what makes the *whole* per node cost
  // parallel -- deduplicating the arc list, grouping it and enumerating the
  // pairs -- and not just the pairs of the few dense nodes.
  //
  // Most of the work sits in the few, very dense levels at the top, and those
  // cannot be spread over their nodes. They are spread over the (node,
  // neighbor) pairs instead, which needs the whole level's grouping to be held
  // at once -- hence `kLevelParallelCutoff`, which `customize()` shares.
  void contract(std::vector<std::vector<cch_rank_t::value_t>> const& by_level) {
    auto const n = c_.n_ranks();

    auto pt = utl::get_active_progress_tracker_or_activate("osr");
    pt->status("CCH Contraction")
        .reset_bounds()
        .out_bounds(8.0F, 90.0F)
        .in_high(n);

    auto done = std::size_t{0U};
    for (auto const& lvl : by_level) {
      if (n_threads_ == 1U || lvl.size() == 1U) {
        static thread_local auto s = node_scratch{};
        for (auto const rank : lvl) {
          prepare(rank, s);
          for (auto k = std::size_t{0U}; k != s.neighbors_.size(); ++k) {
            process_target(s, k);
          }
        }
      } else if (lvl.size() > kLevelParallelCutoff) {
        // wide level: one job per node, every worker reusing its scratch
        pool_.run(lvl.size(), [&](std::size_t const i) {
          static thread_local auto s = node_scratch{};
          prepare(lvl[i], s);
          for (auto k = std::size_t{0U}; k != s.neighbors_.size(); ++k) {
            process_target(s, k);
          }
        });
      } else {
        // narrow level: group every node first, then one job per neighbor
        if (level_scratch_.size() < lvl.size()) {
          level_scratch_.resize(lvl.size());
        }
        pool_.run(lvl.size(), [&](std::size_t const i) {
          prepare(lvl[i], level_scratch_[i]);
        });

        slot_jobs_.clear();
        for (auto i = std::uint32_t{0U}; i != lvl.size(); ++i) {
          for (auto k = std::uint32_t{0U};
               k != level_scratch_[i].neighbors_.size(); ++k) {
            slot_jobs_.emplace_back(i, k);
          }
        }
        pool_.run(slot_jobs_.size(), [&](std::size_t const j) {
          auto const [i, k] = slot_jobs_[j];
          process_target(level_scratch_[i], k);
        });
      }

      done += lvl.size();
      pt->update_monotonic(done);
    }

    write_csr();
  }

  // Writes the finished arc lists out as the hierarchy's CSR.
  //
  // This cannot happen while the levels run: the rows have to come out in rank
  // order and a level is not a range of ranks. Counting the rows first makes
  // the write itself parallel again, so the only sequential part left is the
  // prefix over the counts.
  void write_csr() {
    auto const n = c_.n_ranks();

    auto pt = utl::get_active_progress_tracker_or_activate("osr");
    pt->status("CCH Layout").reset_bounds().out_bounds(90.0F, 97.0F).in_high(2);

    // The counts are already in place, `prepare` filled them in. The four
    // prefixes do not depend on each other, so they run as four jobs.
    pool_.run(4U, [&](std::size_t const which) {
      switch (which) {
        case 0U:
          for (auto i = cch_rank_t::value_t{0U}; i != n; ++i) {
            c_.adj_ofs_[cch_rank_t{i} + 1U] += c_.adj_ofs_[cch_rank_t{i}];
          }
          break;
        case 1U:
          for (auto i = cch_rank_t::value_t{0U}; i != n; ++i) {
            c_.loop_ofs_[cch_rank_t{i} + 1U] += c_.loop_ofs_[cch_rank_t{i}];
          }
          break;
        case 2U:
          for (auto i = std::size_t{0U}; i != n; ++i) {
            up_base_[i + 1U] += up_base_[i];
          }
          break;
        default:
          for (auto i = std::size_t{0U}; i != n; ++i) {
            dn_base_[i + 1U] += dn_base_[i];
          }
          break;
      }
    });
    pt->update_monotonic(1U);

    auto const n_slots = c_.adj_ofs_[cch_rank_t{n}];
    c_.adj_head_.resize(n_slots);
    c_.up_ofs_.resize(n_slots + 1U);
    c_.dn_ofs_.resize(n_slots + 1U);
    c_.up_.resize(up_base_[n]);
    c_.dn_.resize(dn_base_[n]);
    c_.loop_.resize(c_.loop_ofs_[cch_rank_t{n}]);
    c_.up_ofs_[n_slots] = up_base_[n];
    c_.dn_ofs_[n_slots] = dn_base_[n];

    pool_.run(n, [&](std::size_t const i) {
      auto const rank = cch_rank_t{static_cast<cch_rank_t::value_t>(i)};
      auto const n_nb = c_.adj_ofs_[rank + 1U] - c_.adj_ofs_[rank];
      auto const* const heads =
          reinterpret_cast<slot_head const*>(row_[i].get());
      auto const* src_up = reinterpret_cast<cch_entry const*>(
          row_[i].get() + n_nb * sizeof(slot_head));
      auto const* src_dn = src_up + (up_base_[i + 1U] - up_base_[i]);

      auto slot = c_.adj_ofs_[rank];
      auto up = up_base_[i];
      auto dn = dn_base_[i];
      for (auto j = cch_slot_idx_t{0U}; j != n_nb; ++j) {
        c_.adj_head_[slot] = cch_rank_t{heads[j].head_};
        c_.up_ofs_[slot] = up;
        c_.dn_ofs_[slot] = dn;
        std::copy(src_up, src_up + heads[j].n_up_, c_.up_.data() + up);
        std::copy(src_dn, src_dn + heads[j].n_dn_, c_.dn_.data() + dn);
        src_up += heads[j].n_up_;
        src_dn += heads[j].n_dn_;
        up += heads[j].n_up_;
        dn += heads[j].n_dn_;
        ++slot;
      }
      row_[i].reset();

      auto& lp = loops_[i];
      std::copy(begin(lp), end(lp), begin(c_.loop_) + static_cast<std::ptrdiff_t>(
                                        c_.loop_ofs_[rank]));
      lp.clear();
      lp.shrink_to_fit();
    });
    pt->update_monotonic(2U);
  }

  void build_transpose() {
    auto const n = c_.n_ranks();
    c_.lower_ofs_.resize(n + 1U, cch_slot_idx_t{0U});
    for (auto const h : c_.adj_head_) {
      ++c_.lower_ofs_[h];
    }
    auto sum = cch_slot_idx_t{0U};
    for (auto i = cch_rank_t::value_t{0U}; i != n; ++i) {
      auto const cnt = c_.lower_ofs_[cch_rank_t{i}];
      c_.lower_ofs_[cch_rank_t{i}] = sum;
      sum += cnt;
    }
    c_.lower_ofs_[cch_rank_t{n}] = sum;

    c_.lower_slot_.resize(sum);
    auto next = std::vector<cch_slot_idx_t>{};
    next.resize(n);
    for (auto i = cch_rank_t::value_t{0U}; i != n; ++i) {
      next[i] = c_.lower_ofs_[cch_rank_t{i}];
    }
    auto const n_slots = static_cast<cch_slot_idx_t>(c_.adj_head_.size());
    for (auto s = cch_slot_idx_t{0U}; s != n_slots; ++s) {
      c_.lower_slot_[next[to_idx(c_.adj_head_[s])]++] = s;
    }
  }

  ways::routing const& r_;
  cch& c_;
  unsigned n_threads_;
  worker_pool pool_;
  std::vector<std::vector<arc_rec>> adj_;
  std::vector<std::vector<cch_entry>> loops_;
  std::vector<std::size_t> clean_size_;
  std::array<std::mutex, kLockStripes> target_locks_;

  // row sizes of the CSR, filled by `prepare` and prefixed by `write_csr`
  std::vector<cch_entry_idx_t> up_base_, dn_base_;
  // the finished rows, one packed buffer per rank
  std::vector<std::unique_ptr<std::byte[]>> row_;

  // narrow levels: one scratch per node of the level, one job per neighbor
  std::vector<node_scratch> level_scratch_;
  std::vector<std::pair<std::uint32_t, std::uint32_t>> slot_jobs_;
};

}  // namespace

cista::wrapped<cch> build_cch(ways const& w, unsigned const n_threads) {
  auto const t = utl::scoped_timer{"build cch"};

  auto c = cista::wrapped<cch>{cista::raw::make_unique<cch>()};
  auto const& r = *w.r_;
  auto const n_nodes = w.n_nodes();

  auto pt = utl::get_active_progress_tracker_or_activate("osr");

  // ---------------------------------------------------------------------
  // 1. determine the sub graph
  // ---------------------------------------------------------------------
  pt->status("CCH Sub Graph").reset_bounds().out_bounds(0.0F, 1.0F).in_high(
      n_nodes);

  // Chunk size of the node parallel passes. Large enough that the per chunk
  // bookkeeping disappears, small enough that the last chunk cannot hold up
  // the whole pass.
  static constexpr auto const kChunk = std::uint32_t{1U << 16U};
  auto const n_chunks_of = [](std::uint32_t const count) {
    return static_cast<std::size_t>((count + kChunk - 1U) / kChunk);
  };

  auto constexpr kNoLocal = std::numeric_limits<std::uint32_t>::max();
  auto local = vec_map<node_idx_t, std::uint32_t>{};
  local.resize(n_nodes, kNoLocal);
  auto nodes = std::vector<node_idx_t>{};
  {
    // Whether a node belongs to the sub graph is a read only test, so it runs
    // over chunks in parallel. Only the numbering is sequential, and a second
    // pass over the flags produces it without any synchronization: the first
    // pass counts per chunk, the prefix over those counts gives every chunk
    // the index its first node gets.
    auto const n_chunks = n_chunks_of(n_nodes);
    auto ofs = std::vector<std::uint32_t>(n_chunks + 1U, 0U);
    auto const chunk_end = [&](std::size_t const ch) {
      return std::min(n_nodes,
                      static_cast<node_idx_t::value_t>((ch + 1U) * kChunk));
    };

    utl::parallel_for_run(
        n_chunks,
        [&](std::size_t const ch) {
          auto cnt = std::uint32_t{0U};
          for (auto i = static_cast<node_idx_t::value_t>(ch * kChunk);
               i != chunk_end(ch); ++i) {
            auto const n = node_idx_t{i};
            if (utl::any_of(r.node_ways_[n], [&](way_idx_t const way) {
                  return is_cch_way(r.way_properties_[way]);
                })) {
              local[n] = 0U;  // marked, numbered by the second pass
              ++cnt;
            }
          }
          ofs[ch + 1U] = cnt;
        },
        utl::noop_progress_update{}, utl::parallel_error_strategy::QUIT_EXEC,
        n_threads);

    for (auto ch = std::size_t{0U}; ch != n_chunks; ++ch) {
      ofs[ch + 1U] += ofs[ch];
    }
    nodes.resize(ofs[n_chunks]);

    utl::parallel_for_run(
        n_chunks,
        [&](std::size_t const ch) {
          auto idx = ofs[ch];
          for (auto i = static_cast<node_idx_t::value_t>(ch * kChunk);
               i != chunk_end(ch); ++i) {
            auto const n = node_idx_t{i};
            if (local[n] != kNoLocal) {
              local[n] = idx;
              nodes[idx] = n;
              ++idx;
            }
          }
        },
        utl::noop_progress_update{}, utl::parallel_error_strategy::QUIT_EXEC,
        n_threads);
    pt->update_monotonic(n_nodes);
  }

  utl::verify(!nodes.empty(), "cch: empty sub graph");

  // ---------------------------------------------------------------------
  // 2. contraction order (nested dissection)
  // ---------------------------------------------------------------------
  auto g = nd_graph{};
  {
    pt->status("CCH Graph").reset_bounds().out_bounds(1.0F, 2.0F).in_high(
        nodes.size());
    auto const n = static_cast<std::uint32_t>(nodes.size());
    g.ofs_.resize(n + 1U, 0U);
    g.x_.resize(n);
    g.y_.resize(n);

    // Same shape as the sub graph pass, except that the rows are built once
    // into a buffer per chunk and only copied into place afterwards -- the
    // adjacency of a node is not known before it has been deduplicated, so
    // counting first would mean building every row twice.
    auto const n_chunks = n_chunks_of(n);
    auto const chunk_end = [&](std::size_t const ch) {
      return std::min(n, static_cast<std::uint32_t>((ch + 1U) * kChunk));
    };
    auto bufs = std::vector<std::vector<std::uint32_t>>(n_chunks);

    utl::parallel_for_run(
        n_chunks,
        [&](std::size_t const ch) {
          auto& adj = bufs[ch];
          for (auto u = static_cast<std::uint32_t>(ch * kChunk);
               u != chunk_end(ch); ++u) {
            auto const node = nodes[u];
            g.x_[u] = w.get_node_pos(node).lat_;
            g.y_[u] = w.get_node_pos(node).lng_;
            auto const before = adj.size();
            for_each_edge(r, node,
                          [&](node_idx_t const v, way_idx_t const way, port_t,
                              port_t, distance_t, std::uint16_t,
                              std::uint16_t) {
                            if (!is_cch_way(r.way_properties_[way])) {
                              return;
                            }
                            auto const l = local[v];
                            if (l != kNoLocal && l != u) {
                              adj.emplace_back(l);
                            }
                          });
            std::sort(begin(adj) + static_cast<std::ptrdiff_t>(before),
                      end(adj));
            adj.erase(
                std::unique(begin(adj) + static_cast<std::ptrdiff_t>(before),
                            end(adj)),
                end(adj));
            g.ofs_[u + 1U] = adj.size() - before;  // degree, summed below
          }
        },
        utl::noop_progress_update{}, utl::parallel_error_strategy::QUIT_EXEC,
        n_threads);

    for (auto u = std::uint32_t{0U}; u != n; ++u) {
      g.ofs_[u + 1U] += g.ofs_[u];
    }
    g.adj_.resize(g.ofs_[n]);

    utl::parallel_for_run(
        n_chunks,
        [&](std::size_t const ch) {
          auto& adj = bufs[ch];
          std::copy(begin(adj), end(adj),
                    begin(g.adj_) + static_cast<std::ptrdiff_t>(
                                        g.ofs_[ch * kChunk]));
          adj.clear();
          adj.shrink_to_fit();
        },
        utl::noop_progress_update{}, utl::parallel_error_strategy::QUIT_EXEC,
        n_threads);
    pt->update_monotonic(n);
  }

  auto const order = compute_inertial_flow_cutter_order(g, n_threads);
  utl::verify(order.size() == nodes.size(), "cch: bad order");

  // ---------------------------------------------------------------------
  // 3. rank <-> node mapping
  // ---------------------------------------------------------------------
  auto const n_ranks = static_cast<cch_rank_t::value_t>(order.size());
  c->rank_.resize(n_nodes, cch_rank_t::invalid());
  c->order_.resize(n_ranks);
  for (auto i = cch_rank_t::value_t{0U}; i != n_ranks; ++i) {
    auto const node = nodes[order[i]];
    c->order_[cch_rank_t{i}] = node;
    c->rank_[node] = cch_rank_t{i};
  }

  // ---------------------------------------------------------------------
  // 4. elimination tree levels
  // ---------------------------------------------------------------------
  auto const levels = compute_levels(g, order);
  g = nd_graph{};  // not needed any more, and it is not small

  // ---------------------------------------------------------------------
  // 5. contraction
  // ---------------------------------------------------------------------
  auto ctr = contractor{w, *c, n_threads};
  ctr.adj_.resize(n_ranks);
  ctr.loops_.resize(n_ranks);
  ctr.clean_size_.resize(n_ranks, 0U);
  ctr.up_base_.resize(n_ranks + 1U, 0U);
  ctr.dn_base_.resize(n_ranks + 1U, 0U);
  ctr.row_.resize(n_ranks);
  c->adj_ofs_.resize(n_ranks + 1U, cch_slot_idx_t{0U});
  c->loop_ofs_.resize(n_ranks + 1U, cch_entry_idx_t{0U});

  pt->status("CCH Seed").reset_bounds().out_bounds(6.0F, 8.0F).in_high(n_ranks);

  // Runs over the ranks in parallel. `push` and `push_loop` are only safe
  // without a lock as long as every target list is written by one thread, so
  // a rank seeds exactly the records that belong to *itself*: the edge to a
  // higher ranked neighbour, and its reverse, which the neighbour would
  // otherwise have contributed.
  //
  // `for_each_edge` reports the two directions of an original edge from their
  // respective tail, so the reverse of `node --way--> v` is the edge that `v`
  // reports for the same pair of way positions, travelled the other way round.
  // Its ports are therefore the ports of this edge with the direction bit
  // flipped -- for `v` that is the same way slot, for `node` the way slot has
  // to be looked up, because a node can sit on the same way more than once.
  utl::parallel_for_run(
      n_ranks,
      [&](std::size_t const rank) {
        auto const i = static_cast<cch_rank_t::value_t>(rank);
        auto const node = c->order_[cch_rank_t{i}];
        for_each_edge(
            r, node,
            [&](node_idx_t const v, way_idx_t const way,
                port_t const tail_port, port_t const head_port, distance_t,
                std::uint16_t const from, std::uint16_t) {
              if (!is_cch_way(r.way_properties_[way]) || !c->contains(v)) {
                return;
              }
              // a closed way can lead from a node back to itself: that is an
              // original self loop, and both of its directions are reported
              // here anyway
              if (v == node) {
                ctr.push_loop(i, cch_entry{tail_port, head_port});
                return;
              }
              auto const other = to_idx(c->rank_[v]);
              if (i >= other) {
                return;  // both records belong to `other`
              }
              ctr.push(i, arc_rec{other, 1U, tail_port, head_port});

              auto const rev_dir = port_dir(tail_port) == direction::kForward
                                       ? direction::kBackward
                                       : direction::kForward;
              ctr.push(i, arc_rec{other, 0U,
                                  static_cast<port_t>(head_port ^ 1U),
                                  make_port(r.get_way_pos(node, way, from),
                                            rev_dir)});
            });
      },
      utl::noop_progress_update{}, utl::parallel_error_strategy::QUIT_EXEC,
      n_threads);
  pt->update_monotonic(n_ranks);

  ctr.contract(levels);
  ctr.build_transpose();

  // remember which entries an original edge maps to
  c->up_is_edge_.resize(c->up_.size());
  c->dn_is_edge_.resize(c->dn_.size());
  c->loop_is_edge_.resize(c->loop_.size());
  pt->status("CCH Edge Flags")
      .reset_bounds()
      .out_bounds(97.0F, 100.0F)
      .in_high(n_ranks);

  // Runs over the ranks in parallel. Two ranks never set the same bit, but
  // the bits of neighbouring ranks share a block, so the writes are atomic.
  utl::parallel_for_run(
      n_ranks,
      [&](std::size_t const i) {
        auto const rank = cch_rank_t{static_cast<cch_rank_t::value_t>(i)};
        auto const node = c->order_[rank];
        for_each_edge(
            r, node,
            [&](node_idx_t const v, way_idx_t const way, port_t const tail_port,
                port_t const head_port, distance_t, std::uint16_t,
                std::uint16_t) {
              if (!is_cch_way(r.way_properties_[way]) || !c->contains(v)) {
                return;
              }
              auto const ref =
                  c->find_edge(rank, v, cch_entry{tail_port, head_port});
              switch (ref.kind_) {
                case cch_arc::kUp: c->up_is_edge_.set<true>(ref.idx_); break;
                case cch_arc::kDn: c->dn_is_edge_.set<true>(ref.idx_); break;
                case cch_arc::kLoop:
                  c->loop_is_edge_.set<true>(ref.idx_);
                  break;
                case cch_arc::kNone: break;
              }
            });
      },
      utl::noop_progress_update{}, utl::parallel_error_strategy::QUIT_EXEC,
      n_threads);
  pt->update_monotonic(n_ranks);

  return c;
}

}  // namespace osr
