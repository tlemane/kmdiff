#pragma once

#include <kmdiff/icorrector.hpp>

namespace kmdiff {

  class bonferroni : public ICorrector
  {
    public:
      bonferroni(double threshold, std::size_t total);
      bool apply(double pvalue) override;
      CorrectionType type() override;
      std::string str_type() override;
    private:
      std::size_t m_total {0};
      double m_threshold {0.0};
  };

  class benjamini : public ICorrector
  {
    public:
      benjamini(double fdr, std::size_t total);
      bool apply(double pvalue) override;
      CorrectionType type() override;
      std::string str_type() override;

    private:
      std::size_t m_total {0};
      std::size_t m_rank {1};
      double m_fdr {0.0};
  };

  class holm : public ICorrector
  {
    public:
      holm(double threshold, std::size_t total);
      bool apply(double pvalue) override;
      CorrectionType type() override;
      std::string str_type() override;

    private:
      std::size_t m_total {0};
      double m_threshold {0};
  };

  class sidak : public ICorrector
  {
    public:
      sidak(double threshold, std::size_t total);
      bool apply(double pvalue) override;
      CorrectionType type() override;
      std::string str_type() override;

    private:
      std::size_t m_total {0};
      double m_threshold {0.0};
  };

  class basic_threshold  : public ICorrector
  {
    public:
      basic_threshold(double threshold);
      bool apply(double pvalue) override;
      CorrectionType type() override;
      std::string str_type() override;

    private:
      double m_threshold {0};
  };

  std::shared_ptr<ICorrector> make_corrector(CorrectionType type, double threshold, std::size_t kmers);

} // end of namespace kmdiff

