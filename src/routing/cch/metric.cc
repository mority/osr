#include "osr/routing/cch/customize.h"

#include "cista/serialization.h"

namespace osr {

std::filesystem::path cch_metric::file(std::filesystem::path const& p) {
  return p / "cch_metric_car.bin";
}

bool cch_metric::exists(std::filesystem::path const& p) {
  return std::filesystem::is_regular_file(file(p));
}

cista::wrapped<cch_metric> cch_metric::read(std::filesystem::path const& p) {
  return cista::read<cch_metric, kMode>(file(p));
}

void cch_metric::write(std::filesystem::path const& p) const {
  cista::write<kMode>(file(p), *this);
}

}  // namespace osr
