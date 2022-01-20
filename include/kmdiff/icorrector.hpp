#pragma once

#include <memory>
#include <string>
#include <kmdiff/correction.hpp>

namespace kmdiff {

  class ICorrector
  {
    public:
      ICorrector() = default;
      virtual CorrectionType type() = 0;
      virtual std::string str_type() = 0;
      virtual bool apply(double pvalue) = 0;
  };

  using corrector_t = std::shared_ptr<ICorrector>;
} // end of namespace kmdiff

