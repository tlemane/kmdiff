#include <string>
#include <filesystem>
#include <kmdiff/accumulator.hpp>

namespace fs = std::filesystem;

namespace kmdiff {

  bool partitions_exist(const std::string& prefix, size_t nb_partitions, const std::string& dir)
  {
    for (std::size_t i = 0; i < nb_partitions; ++i)
    {
      if (!fs::exists(fmt::format(prefix, dir, i)))
        return false;
    }
    return true;
  }

} // end of namespace kmdiff
