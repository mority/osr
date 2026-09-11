#pragma once

#include <algorithm>
#include <string_view>

#include "utl/parser/arg_parser.h"
#include "utl/parser/cstr.h"

#include "osr/types.h"

namespace osr {

// A `level` value - "0", "-1", "0;1", "1.5" - as the set of levels it names,
// clamped to the range level_t can hold.
inline level_bits_t parse_levels(std::string_view const value) {
  auto bits = level_bits_t{0U};
  auto s = utl::cstr{value};
  while (s) {
    auto l = 0.0F;
    utl::parse_arg(s, l);
    auto const lvl = level_t{std::clamp(l, kMinLevel, kMaxLevel)};
    bits |= (static_cast<level_bits_t>(1) << to_idx(lvl));
    if (s) {
      ++s;
    }
  }
  return bits;
}

}  // namespace osr
