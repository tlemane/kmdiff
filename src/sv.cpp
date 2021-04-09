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