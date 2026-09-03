// Reports the shape of a hierarchy that is already on disk, without rebuilding
// it: node and arc counts, entry counts, and the share of the graph that
// entered the hierarchy.
#include <filesystem>
#include <iostream>

#include "fmt/core.h"
#include "fmt/std.h"

#include "conf/options_parser.h"

#include "osr/routing/cch/cch.h"
#include "osr/routing/cch/customize.h"
#include "osr/ways.h"

namespace fs = std::filesystem;
using namespace osr;

struct settings : public conf::configuration {
  settings() : configuration{"Options"} {
    param(data_, "data,d", "osr data directory");
    param(name_, "name,n", "label for the row");
  }
  fs::path data_{"osr"};
  std::string name_{"network"};
};

int main(int ac, char const** av) {
  auto opt = settings{};
  auto parser = conf::options_parser{{&opt}};
  parser.read_command_line_args(ac, av);
  if (parser.help()) {
    parser.print_help(std::cout);
    return 0;
  }
  if (!cch::exists(opt.data_)) {
    fmt::println("no cch in {}", opt.data_);
    return 1;
  }
  auto const w = ways{opt.data_, cista::mmap::protection::READ};
  auto const c = cch::read(opt.data_);

  auto n_loops = std::size_t{0U};
  for (auto i = cch_rank_t::value_t{0U}; i != c->n_ranks(); ++i) {
    n_loops += c->loop_entries(cch_rank_t{i}).size();
  }
  auto const entries = c->up_.size() + c->dn_.size() + n_loops;
  auto const graph_nodes = static_cast<std::size_t>(w.n_nodes());
  auto const cch_nodes = static_cast<std::size_t>(c->n_ranks());

  fmt::println(
      "INFO {} | graph_nodes {} | cch_nodes {} | share {:.4f} | slots {} | "
      "up {} | dn {} | loops {} | entries {} | entries_per_node {:.2f}",
      opt.name_, graph_nodes, cch_nodes,
      graph_nodes == 0 ? 0.0
                       : static_cast<double>(cch_nodes) /
                             static_cast<double>(graph_nodes),
      c->n_slots(), c->up_.size(), c->dn_.size(), n_loops, entries,
      cch_nodes == 0 ? 0.0
                     : static_cast<double>(entries) /
                           static_cast<double>(cch_nodes));
  return 0;
}
