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
#pragma once

// std
#include <sstream>
#include <string>

// ext
#include <xxhash.h>

namespace kmdiff
{
class SV
{
  friend struct std::hash<SV>;

 public:
  std::string m_chr;
  std::string m_type;
  std::string m_seq;
  size_t m_start;
  size_t m_end;

 public:
  SV(const std::string& chr,
     const std::string& type,
     const std::string& seq,
     size_t start,
     size_t end);

  std::string to_bed_entry();

  bool operator==(const SV& rhs) const;
  bool operator<(const SV& rhs) const;
};

};  // end of namespace kmdiff

template <>
struct std::hash<kmdiff::SV>
{
  uint64_t operator()(const kmdiff::SV& sv) const noexcept
  {
    std::string s =
        sv.m_chr + sv.m_type + sv.m_seq + std::to_string(sv.m_start) + std::to_string(sv.m_end);
    return static_cast<uint64_t>(XXH64(s.c_str(), s.size(), 0));
  }
};