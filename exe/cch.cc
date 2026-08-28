#ifdef _WIN32
#include "windows.h"
#endif

#include <filesystem>
#include <iostream>

#include "fmt/core.h"
#include "fmt/std.h"

#include "conf/options_parser.h"

#include "utl/progress_tracker.h"
#include "utl/timer.h"

#include "osr/lookup.h"
#include "osr/routing/cch/build.h"
#include "osr/routing/cch/cch.h"
#include "osr/routing/cch/customize.h"
#include "osr/routing/profiles/car.h"
#include "osr/ways.h"

namespace fs = std::filesystem;
using namespace osr;

struct settings : public conf::configuration {
  settings() : configuration{"Options"} {
    param(data_, "data,d", "osr data directory");
    param(customize_, "customize,c", "run a customization for the car profile");
    param(threads_, "threads,t", "build/customization threads (0 = all cores)");
  }

  fs::path data_{"osr"};
  bool customize_{true};
  unsigned threads_{0U};
};

int main(int ac, char const** av) {
  auto opt = settings{};
  auto parser = conf::options_parser{{&opt}};
  parser.read_command_line_args(ac, av);
  if (parser.help()) {
    parser.print_help(std::cout);
    return 0;
  }
  parser.print_used(std::cout);

  if (!fs::is_directory(opt.data_)) {
    fmt::println("data directory {} not found", opt.data_);
    return 1;
  }

  utl::activate_progress_tracker("osr");
  auto const silencer = utl::global_progress_bars{false};

  auto const w = ways{opt.data_, cista::mmap::protection::READ};
  auto const c = build_cch(w, opt.threads_);
  c->write(opt.data_);

  auto n_loops = std::size_t{0U};
  for (auto i = cch_rank_t::value_t{0U}; i != c->n_ranks(); ++i) {
    n_loops += c->loop_entries(cch_rank_t{i}).size();
  }
  fmt::println(
      "cch: {} nodes, {} arcs, {} up entries, {} dn entries, {} self loops",
      c->n_ranks(), c->n_slots(), c->up_.size(), c->dn_.size(), n_loops);

  if (opt.customize_) {
    auto m = cch_metric{};
    {
      auto const t = utl::scoped_timer{"customize car"};
      customize<car>(car::parameters{}, w, *c, m, opt.threads_);
    }
    auto n_inf = std::size_t{0U};
    for (auto const x : m.up_) {
      n_inf += (x == cch_metric::kWeightInfeasible ? 1U : 0U);
    }
    for (auto const x : m.dn_) {
      n_inf += (x == cch_metric::kWeightInfeasible ? 1U : 0U);
    }
    fmt::println("customization: {} infeasible entries, {} triangles applied",
                 n_inf, m.applied_);
    m.write(opt.data_);
    fmt::println("wrote {}", cch_metric::file(opt.data_));
  }
}
