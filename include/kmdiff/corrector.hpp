#pragma once

#include <kmdiff/icorrector.hpp>

namespace kmdiff {

  class Bonferroni : public ICorrector
  {
    public:
      Bonferroni(double theshold, size_t total);
      bool apply(double pvalue) override;
    private:
      size_t m_total {0};
      double m_threshold {0.0};
  };

  class Benjamini : public ICorrector
  {
    public:
      Benjamini(double fdr, size_t total);
      bool apply(double pvalue) override;
    private:
      size_t m_total {0};
      size_t m_rank {1};
      double m_fdr {0.0};
  };

} // end of namespace kmdiff

