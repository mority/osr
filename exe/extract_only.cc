// Extraction without the contraction hierarchy.
//
// `osr-extract` on this branch builds and customizes the CCH as part of
// extraction, so it cannot answer the question a reader actually has: what does
// adding the hierarchy cost on top of the preprocessing I already run? This
// produces the other half of that comparison. Wall time and peak RSS come from
// /usr/bin/time around the process, so both halves are measured the same way.
#include <filesystem>
#include <iostream>

#include "fmt/core.h"
#include "fmt/std.h"

#include "conf/options_parser.h"

#include "osr/extract/extract.h"

namespace fs = std::filesystem;

struct settings : public conf::configuration {
  settings() : configuration{"Options"} {
    param(in_, "in,i", "input pbf");
    param(out_, "out,o", "output directory");
  }
  fs::path in_, out_{"osr"};
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
  osr::extract_ways(true, opt.in_, opt.out_, fs::path{});
  return 0;
}
