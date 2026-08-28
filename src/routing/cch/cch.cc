#include "osr/routing/cch/cch.h"

#include "cista/serialization.h"

namespace osr {

cista::wrapped<cch> cch::read(std::filesystem::path const& p) {
  return cista::read<cch, kMode>(p / "cch.bin");
}

void cch::write(std::filesystem::path const& p) const {
  cista::write<kMode>(p / "cch.bin", *this);
}

bool cch::exists(std::filesystem::path const& p) {
  return std::filesystem::is_regular_file(p / "cch.bin");
}

}  // namespace osr
