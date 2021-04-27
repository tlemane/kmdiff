/*****************************************************************************
 *   kmdiff
 *   Authors: T. Lemane
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

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
    case CorrectionType::BENJAMINI:
      return "BENJAMINI";
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
  this->m_type = CorrectionType::BONFERRONI;
}

bool Bonferroni::apply(double p_value) { return p_value < m_threshold / m_total; }

BenjaminiHochberg::BenjaminiHochberg(double fdr, size_t total)
  : m_fdr(fdr), m_total(static_cast<double>(total))
{
  this->m_type = CorrectionType::BENJAMINI;
}

bool BenjaminiHochberg::apply(double p_value)
{
  if (p_value < (m_rank/m_total)*m_fdr)
  {
    m_rank++;
    return true;
  }
  return false;
}

};  // namespace kmdiff