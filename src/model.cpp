#include <kmdiff/model.hpp>

namespace kmdiff
{
std::string correction_type_str(const CorrectionType type)
{
  switch (type)
  {
    case CorrectionType::NOTHING:
      return "NOTHING";
      break;
    case CorrectionType::BONFERRONI:
      return "BONFERRONI";
      break;
    default:
      return "";
      break;
  }
}

bool ICorrector::apply(double p_value) { return true; }

Bonferroni::Bonferroni(double threshold, size_t total)
    : m_threshold(threshold), m_total(static_cast<double>(total))
{
}

bool Bonferroni::apply(double p_value) { return p_value < m_threshold / m_total; }

};  // namespace kmdiff