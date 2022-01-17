#pragma once

#include <string>

namespace kmdiff {

  class ICorrector
  {
    public:
      ICorrector = default;

      virtual void configure(const std::string& config) = 0;
      virtual bool apply(double pvalue) = 0;
  };

} // end of namespace kmdiff

