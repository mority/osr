#include "osr/routing/cch/build.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <condition_variable>
#include <mutex>
#include <thread>
#include <span>
#include <vector>

#include "utl/helpers/algorithm.h"
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

struct contractor {
  contractor(ways const& w, cch& c, unsigned const n_threads)
      : r_{*w.r_},
        c_{c},
        n_threads_{n_threads == 0U
                       ? std::max(1U, std::thread::hardware_concurrency())
                       : n_threads},
        pool_{n_threads_ - 1U} {}

  // Below this degree the threads cost more than the pairs they split.
  static constexpr auto const kMinParallelDegree = std::size_t{32U};

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
  void build_reach(node_idx_t const n,
                   std::span<cch_entry const> loops,
                   port_t const ports) {
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

  bool turn_reachable(node_idx_t const n,
                      port_t const in,
                      port_t const out) const {
    auto rest = reach_[in];
    while (rest != 0U) {
      auto const k = static_cast<port_t>(std::countr_zero(rest));
      rest &= rest - 1U;
      if (turn_ok(n, k, out)) {
        return true;
      }
    }
    return false;
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

  void contract() {
    auto const n = c_.n_ranks();

    auto pt = utl::get_active_progress_tracker_or_activate("osr");
    // The bar advances by rank, which is the only counter available here.
    // Be aware that it is optimistic: the work per node grows with its degree
    // and the highest ranked nodes are by far the densest, so the last few
    // percent of the ranks are most of the time.
    pt->status("CCH Contraction")
        .reset_bounds()
        .out_bounds(8.0F, 97.0F)
        .in_high(n);

    c_.adj_ofs_.resize(n + 1U);
    c_.loop_ofs_.resize(n + 1U);
    c_.up_ofs_.emplace_back(0U);
    c_.dn_ofs_.emplace_back(0U);

    auto in_of = std::vector<std::vector<cch_entry>>{};
    auto out_of = std::vector<std::vector<cch_entry>>{};
    auto neighbors = std::vector<cch_rank_t::value_t>{};

    for (auto i = cch_rank_t::value_t{0U}; i != n; ++i) {
      auto const r = cch_rank_t{i};
      auto const node = c_.order_[r];
      c_.adj_ofs_[r] = static_cast<cch_slot_idx_t>(c_.adj_head_.size());
      c_.loop_ofs_[r] = static_cast<cch_entry_idx_t>(c_.loop_.size());

      auto& v = adj_[i];
      dedup_tail(v, clean_size_[i]);

      auto& lp = loops_[i];
      dedup(lp);
      for (auto const& e : lp) {
        c_.loop_.emplace_back(e);
      }
      auto const ports = std::min(n_ports(r_, node), kMaxPorts);
      build_reach(node, std::span<cch_entry const>{lp}, ports);
      lp.clear();
      lp.shrink_to_fit();

      // group the arc records by neighbor
      neighbors.clear();
      in_of.clear();
      out_of.clear();
      for (auto it = begin(v); it != end(v);) {
        auto const other = it->other_;
        neighbors.emplace_back(other);
        in_of.emplace_back();
        out_of.emplace_back();
        for (; it != end(v) && it->other_ == other; ++it) {
          (it->out_ != 0U ? out_of : in_of).back().push_back(
              cch_entry{it->p_, it->q_});
        }
      }
      v.clear();
      v.shrink_to_fit();
      clean_size_[i] = 0U;

      // write the final CSR row of this node
      for (auto k = std::size_t{0U}; k != neighbors.size(); ++k) {
        c_.adj_head_.emplace_back(cch_rank_t{neighbors[k]});
        for (auto const& e : out_of[k]) {
          c_.up_.emplace_back(e);
        }
        for (auto const& e : in_of[k]) {
          c_.dn_.emplace_back(e);
        }
        c_.up_ofs_.emplace_back(static_cast<cch_entry_idx_t>(c_.up_.size()));
        c_.dn_ofs_.emplace_back(static_cast<cch_entry_idx_t>(c_.dn_.size()));
      }

      // contract: connect all pairs of remaining neighbors. `a == b` creates
      // a self loop, i.e. the shortcut of turning around at this node.
      //
      // Every arc is pushed to the lower ranked of its two ends, so the pairs
      // are enumerated per *target* instead of as a full matrix: the neighbor
      // `k` collects the pairs whose lower end it is. Each target is then
      // touched by exactly one task, which is what makes this parallel
      // without any locking -- `adj_`, `loops_` and `clean_size_` are all
      // indexed by the target.
      auto const process_target = [&](std::size_t const k) {
        auto const t = neighbors[k];

        // pairs (k, b): this node is the tail, so it owns them if t < y
        if (!in_of[k].empty()) {
          for (auto b = std::size_t{0U}; b != neighbors.size(); ++b) {
            if (out_of[b].empty()) {
              continue;
            }
            auto const y = neighbors[b];
            if (b != k && t >= y) {
              continue;  // owned by `b`
            }
            for (auto const& e1 : in_of[k]) {
              for (auto const& e2 : out_of[b]) {
                if (!turn_reachable(node, e1.exit_, e2.entry_)) {
                  continue;
                }
                if (b == k) {
                  push_loop(t, cch_entry{e1.entry_, e2.exit_});
                } else {
                  push(t, arc_rec{y, 1U, e1.entry_, e2.exit_});
                }
              }
            }
          }
        }

        // pairs (a, k): this node is the head, so it owns them if t < x
        if (!out_of[k].empty()) {
          for (auto a = std::size_t{0U}; a != neighbors.size(); ++a) {
            if (a == k || in_of[a].empty()) {
              continue;
            }
            auto const x = neighbors[a];
            if (x <= t) {
              continue;  // owned by `a`
            }
            for (auto const& e1 : in_of[a]) {
              for (auto const& e2 : out_of[k]) {
                if (!turn_reachable(node, e1.exit_, e2.entry_)) {
                  continue;
                }
                push(t, arc_rec{x, 0U, e1.entry_, e2.exit_});
              }
            }
          }
        }
      };

      // The work of a node grows with the square of its degree, so only the
      // dense nodes are worth spreading over threads -- and those are exactly
      // the ones that dominate the contraction.
      if (n_threads_ > 1U && neighbors.size() >= kMinParallelDegree) {
        pool_.run(neighbors.size(), process_target);
      } else {
        for (auto k = std::size_t{0U}; k != neighbors.size(); ++k) {
          process_target(k);
        }
      }

      pt->update_monotonic(i);
    }

    c_.adj_ofs_[cch_rank_t{n}] =
        static_cast<cch_slot_idx_t>(c_.adj_head_.size());
    c_.loop_ofs_[cch_rank_t{n}] =
        static_cast<cch_entry_idx_t>(c_.loop_.size());
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
  std::array<std::uint32_t, kMaxPorts> reach_{};
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
  auto local = vec_map<node_idx_t, std::uint32_t>{};
  local.resize(n_nodes, std::numeric_limits<std::uint32_t>::max());
  auto nodes = std::vector<node_idx_t>{};
  for (auto i = node_idx_t{0U}; i != node_idx_t{n_nodes}; ++i) {
    auto const ways_of_node = r.node_ways_[i];
    auto const relevant = utl::any_of(ways_of_node, [&](way_idx_t const way) {
      return is_cch_way(r.way_properties_[way]);
    });
    if (relevant) {
      local[i] = static_cast<std::uint32_t>(nodes.size());
      nodes.emplace_back(i);
    }
    pt->update_monotonic(to_idx(i));
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
    auto adj = std::vector<std::uint32_t>{};
    for (auto u = std::uint32_t{0U}; u != n; ++u) {
      auto const node = nodes[u];
      g.x_[u] = w.get_node_pos(node).lat_;
      g.y_[u] = w.get_node_pos(node).lng_;
      auto const before = adj.size();
      for_each_edge(r, node,
                    [&](node_idx_t const v, way_idx_t const way, port_t,
                        port_t, distance_t, std::uint16_t, std::uint16_t) {
                      if (!is_cch_way(r.way_properties_[way])) {
                        return;
                      }
                      auto const l = local[v];
                      if (l != std::numeric_limits<std::uint32_t>::max() &&
                          l != u) {
                        adj.emplace_back(l);
                      }
                    });
      std::sort(begin(adj) + static_cast<std::ptrdiff_t>(before), end(adj));
      adj.erase(std::unique(begin(adj) + static_cast<std::ptrdiff_t>(before),
                            end(adj)),
                end(adj));
      g.ofs_[u + 1U] = adj.size();
      pt->update_monotonic(u);
    }
    g.adj_ = std::move(adj);
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
  // 4. contraction
  // ---------------------------------------------------------------------
  auto ctr = contractor{w, *c, n_threads};
  ctr.adj_.resize(n_ranks);
  ctr.loops_.resize(n_ranks);
  ctr.clean_size_.resize(n_ranks, 0U);

  pt->status("CCH Seed").reset_bounds().out_bounds(6.0F, 8.0F).in_high(n_ranks);
  for (auto i = cch_rank_t::value_t{0U}; i != n_ranks; ++i) {
    auto const node = c->order_[cch_rank_t{i}];
    for_each_edge(r, node,
                  [&](node_idx_t const v, way_idx_t const way,
                      port_t const tail_port, port_t const head_port,
                      distance_t, std::uint16_t, std::uint16_t) {
                    if (!is_cch_way(r.way_properties_[way]) ||
                        !c->contains(v)) {
                      return;
                    }
                    // a closed way can lead from a node back to itself: that
                    // is an original self loop
                    if (v == node) {
                      ctr.push_loop(i, cch_entry{tail_port, head_port});
                      return;
                    }
                    auto const other = to_idx(c->rank_[v]);
                    if (i < other) {
                      ctr.push(i, arc_rec{other, 1U, tail_port, head_port});
                    } else {
                      ctr.push(other, arc_rec{i, 0U, tail_port, head_port});
                    }
                  });
    pt->update_monotonic(i);
  }

  ctr.contract();
  ctr.build_transpose();

  // remember which entries an original edge maps to
  c->up_is_edge_.resize(c->up_.size());
  c->dn_is_edge_.resize(c->dn_.size());
  c->loop_is_edge_.resize(c->loop_.size());
  pt->status("CCH Edge Flags")
      .reset_bounds()
      .out_bounds(97.0F, 100.0F)
      .in_high(n_ranks);
  for (auto i = cch_rank_t::value_t{0U}; i != n_ranks; ++i) {
    auto const rank = cch_rank_t{i};
    auto const node = c->order_[rank];
    for_each_edge(
        r, node,
        [&](node_idx_t const v, way_idx_t const way, port_t const tail_port,
            port_t const head_port, distance_t, std::uint16_t,
            std::uint16_t) {
          if (!is_cch_way(r.way_properties_[way]) || !c->contains(v)) {
            return;
          }
          auto const e = cch_entry{tail_port, head_port};
          if (v == node) {
            auto const idx = cch::find_entry(c->loop_entries(rank), e);
            if (idx != std::numeric_limits<cch_entry_idx_t>::max()) {
              c->loop_is_edge_.set(c->loop_begin(rank) + idx);
            }
            return;
          }
          auto const other = c->rank_[v];
          auto const up = rank < other;
          auto const slot =
              up ? c->find_slot(rank, other) : c->find_slot(other, rank);
          if (slot == cch::kNoSlot) {
            return;
          }
          auto const idx = cch::find_entry(
              up ? c->up_entries(slot) : c->dn_entries(slot), e);
          if (idx == std::numeric_limits<cch_entry_idx_t>::max()) {
            return;
          }
          if (up) {
            c->up_is_edge_.set(c->up_ofs_[slot] + idx);
          } else {
            c->dn_is_edge_.set(c->dn_ofs_[slot] + idx);
          }
        });
    pt->update_monotonic(i);
  }

  return c;
}

}  // namespace osr
