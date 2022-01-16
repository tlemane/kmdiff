#pragma once

#include <string>

namespace kmdiff {

  enum class CorrectionType
  {
    NOTHING,
    BONFERRONI,
    BENJAMINI
  };

  std::string correction_type_str(const CorrectionType type);
}

