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