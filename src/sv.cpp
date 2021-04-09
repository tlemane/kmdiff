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

#include <kmdiff/sv.hpp>

namespace kmdiff
{
SV::SV(
    const std::string& chr,
    const std::string& type,
    const std::string& seq,
    size_t start,
    size_t end)
    : m_chr(chr), m_type(type), m_seq(seq), m_start(start), m_end(end)
{
}

std::string SV::to_bed_entry()
{
  std::stringstream ss;
  ss << m_chr << '\t';
  ss << std::to_string(m_start) << '\t';
  ss << std::to_string(m_end) << '\t';
  ss << m_type << '\t';
  ss << m_seq << '\t';
  ss << '0';
  return ss.str();
}

bool SV::operator==(const SV& rhs) const
{
  return m_chr == rhs.m_chr && m_type == rhs.m_type && m_seq == rhs.m_seq &&
         m_start == rhs.m_start && m_end == rhs.m_end;
}

bool SV::operator<(const SV& rhs) const { return m_start < rhs.m_start; }

};  // namespace kmdiff