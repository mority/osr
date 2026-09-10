#pragma once

#include <algorithm>
#include <cstdint>
#include <thread>
#include <utility>
#include <vector>

#include "oneapi/tbb/blocked_range.h"
#include "oneapi/tbb/parallel_for.h"
#include "oneapi/tbb/task_arena.h"

#include "osr/routing/cch/cch.h"

namespace osr {

// A fixed set of workers, held for the whole contraction resp. customization.
//
// The parallel region is entered once or twice per level and there are
// thousands of levels, so the workers have to outlive a level: creating and
// joining a set of threads per region costs more than the level's own work --
// at 32 threads that overhead alone made the contraction slower than at 8.
struct cch_arena {
  explicit cch_arena(unsigned const n_threads)
      : arena_{static_cast<int>(
            n_threads == 0U ? std::max(1U, std::thread::hardware_concurrency())
                            : n_threads)} {}

  // Calls `job(i)` for every i < count and returns once all of them are done.
  //
  // Grain size one with the simple partitioner, i.e. ranges that split down to
  // the single index: the cost of a node varies by orders of magnitude, so a
  // range that cannot be split further would let one dense node hold up the
  // whole level.
  template <typename Fn>
  void run(std::size_t const count, Fn const& job) {
    if (count == 1U) {
      job(std::size_t{0U});  // most levels hold a single node
      return;
    }
    arena_.execute([&]() {
      tbb::parallel_for(
          tbb::blocked_range<std::size_t>{std::size_t{0U}, count, 1U},
          [&](tbb::blocked_range<std::size_t> const& r) {
            for (auto i = r.begin(); i != r.end(); ++i) {
              job(i);
            }
          },
          tbb::simple_partitioner{});
    });
  }

private:
  tbb::task_arena arena_;
};

// Walks the elimination tree level by level.
//
// The contraction order is a valid order, but so is any order that finishes a
// level before it starts the next one: a node only ever touches state that
// belongs to a node on a larger level, so within a level the nodes are
// independent of each other. The contraction and the customization both walk
// the levels this way, which is why they share `kLevelParallelCutoff`.
//
// Most of the work sits in the few, very dense levels at the top: on Hamburg,
// levels with fewer than 64 nodes hold 83% of it. Those cannot be spread over
// their nodes, so they are spread over the (node, slot) pairs of the whole
// level instead -- distinct slots of a node touch distinct arcs, so that is as
// race free as the split by node. It needs the scratch of every node of the
// level at once, which is why wide levels -- cheap, and there can be tens of
// thousands of them -- keep the split by node.
//
//   prepare(scratch, rank)          makes `rank` ready, filling its scratch
//   slots(scratch, rank) -> [b, e)  half open slot range of `rank`
//   job(scratch, rank, slot)        the work of one (node, slot) pair
//   progress(n_nodes_of_level)      called once per finished level
template <typename Rank,
          typename Scratch,
          typename PrepareFn,
          typename SlotsFn,
          typename JobFn,
          typename ProgressFn>
void for_each_level(cch_arena& arena,
                    std::vector<std::vector<Rank>> const& by_level,
                    std::vector<Scratch>& level_scratch,
                    PrepareFn const& prepare,
                    SlotsFn const& slots,
                    JobFn const& job,
                    ProgressFn const& progress) {
  using slot_t = decltype(slots(std::declval<Scratch const&>(),
                                std::declval<Rank>())
                              .first);
  auto slot_jobs = std::vector<std::pair<std::uint32_t, slot_t>>{};

  for (auto const& lvl : by_level) {
    if (lvl.size() > kLevelParallelCutoff) {
      // wide level: one job per node, every worker reusing its own scratch
      arena.run(lvl.size(), [&](std::size_t const i) {
        static thread_local auto s = Scratch{};
        prepare(s, lvl[i]);
        auto const [b, e] = slots(s, lvl[i]);
        for (auto k = b; k != e; ++k) {
          job(s, lvl[i], k);
        }
      });
    } else {
      // narrow level: prepare every node first, then one job per slot
      if (level_scratch.size() < lvl.size()) {
        level_scratch.resize(lvl.size());
      }
      arena.run(lvl.size(), [&](std::size_t const i) {
        prepare(level_scratch[i], lvl[i]);
      });

      slot_jobs.clear();
      for (auto i = std::uint32_t{0U}; i != lvl.size(); ++i) {
        auto const [b, e] = slots(level_scratch[i], lvl[i]);
        for (auto k = b; k != e; ++k) {
          slot_jobs.emplace_back(i, k);
        }
      }
      arena.run(slot_jobs.size(), [&](std::size_t const j) {
        auto const [i, k] = slot_jobs[j];
        job(level_scratch[i], lvl[i], k);
      });
    }

    progress(lvl.size());
  }
}

}  // namespace osr
