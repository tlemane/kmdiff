#include <kmdiff/corrector.hpp>
#include <spdlog/spdlog.h>

namespace kmdiff {

  bonferroni::bonferroni(double threshold, std::size_t total)
    : m_total(total), m_threshold(threshold) {}

  bool bonferroni::apply(double pvalue)
  {
    return pvalue < (m_threshold / m_total);
  }

  CorrectionType bonferroni::type()
  {
    return CorrectionType::BONFERRONI;
  }

  std::string bonferroni::str_type()
  {
    return "bonferroni";
  }

  benjamini::benjamini(double fdr, std::size_t total)
    : m_total(total), m_fdr(fdr) {}

  bool benjamini::apply(double pvalue)
  {
    if (pvalue < ((m_rank / static_cast<double>(m_total)) * m_fdr))
    {
      m_rank++;
      return true;
    }
    return false;
  }

  CorrectionType benjamini::type()
  {
    return CorrectionType::BENJAMINI;
  }

  std::string benjamini::str_type()
  {
    return "benjamini";
  }

  basic_threshold::basic_threshold(double threshold)
    : m_threshold(threshold) {}

  bool basic_threshold::apply(double pvalue)
  {
    return pvalue < m_threshold;
  }

  CorrectionType basic_threshold::type()
  {
    return CorrectionType::NOTHING;
  }

  std::string basic_threshold::str_type()
  {
    return "threshold";
  }

} // end of namespace kmdiff

