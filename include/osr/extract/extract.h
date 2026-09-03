#include <filesystem>

namespace osr {

// Extraction only, without the contraction hierarchy that `extract` goes on to
// build. Exposed so preprocessing can be measured phase by phase.
void extract_ways(bool with_platforms,
                  std::filesystem::path const& in,
                  std::filesystem::path const& out,
                  std::filesystem::path const& elevation_dir);

void extract(bool with_platforms,
             std::filesystem::path const& in,
             std::filesystem::path const& out,
             std::filesystem::path const& elevation_dir);

}  // namespace osr