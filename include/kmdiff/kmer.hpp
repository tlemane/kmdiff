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

#include <cstddef>
#include <cstdint>
#include <vector>
#include <memory>
#include <fstream>

#include <xxhash.h>

#include <kmtricks/kmer.hpp>

namespace kmdiff {

  enum class Significance
  {
    CONTROL,
    CASE,
    NO
  };

  inline char significance_to_char(Significance sign)
  {
    switch (sign)
    {
      case Significance::CONTROL:
        return '-';
      case Significance::CASE:
        return '+';
      case Significance::NO:
        return '$';
      default:
        return '?';
    }
  }

  template <size_t MAX_K>
  class KmerSign
  {
    friend struct std::hash<KmerSign<MAX_K>>;
   public:
  #ifdef WITH_POPSTRAT
    KmerSign(km::Kmer<MAX_K>&& kmer, double pvalue, Significance sign,
             std::vector<double>& counts_ratio, double mean_control = 0, double mean_case = 0)
        : m_kmer(std::move(kmer)), m_pvalue(pvalue), m_sign(sign), m_counts_ratio(counts_ratio),
          m_mean_control(mean_control), m_mean_case(mean_case)
    {
    }
  #else
    KmerSign(km::Kmer<MAX_K>&& kmer, double pvalue, Significance sign,
             double mean_control = 0, double mean_case = 0)
        : m_kmer(std::move(kmer)), m_pvalue(pvalue), m_sign(sign),
          m_mean_control(mean_control), m_mean_case(mean_case)
    {
    }
  #endif

    KmerSign() {}

    void set_pval(double pvalue)
    {
      m_pvalue = pvalue;
    }

    std::string to_string() const
    {
      return m_kmer.to_string();
    }

    void load(std::shared_ptr<std::istream> stream, size_t size)
    {
      m_kmer.load(*stream);
      stream->read(reinterpret_cast<char*>(&m_pvalue), sizeof(m_pvalue));
      stream->read(reinterpret_cast<char*>(&m_sign), sizeof(m_sign));
      stream->read(reinterpret_cast<char*>(&m_mean_control), sizeof(m_mean_control));
      stream->read(reinterpret_cast<char*>(&m_mean_case), sizeof(m_mean_case));
  #ifdef WITH_POPSTRAT
      std::uint16_t size_ = 0;
      stream->read(reinterpret_cast<char*>(&size_), sizeof(size_));
      m_counts_ratio.resize(size_, 0);
      stream->read(reinterpret_cast<char*>(m_counts_ratio.data()),
                   size_*sizeof(double));
  #endif
    }

    void dump(std::shared_ptr<std::ostream> stream)
    {
      m_kmer.dump(*stream);
      stream->write(reinterpret_cast<char*>(&m_pvalue), sizeof(m_pvalue));
      stream->write(reinterpret_cast<char*>(&m_sign), sizeof(m_sign));
      stream->write(reinterpret_cast<char*>(&m_mean_control), sizeof(m_mean_control));
      stream->write(reinterpret_cast<char*>(&m_mean_case), sizeof(m_mean_case));
  #ifdef WITH_POPSTRAT
      std::uint16_t size = m_counts_ratio.size();
      stream->write(reinterpret_cast<char*>(&size), sizeof(size));
      stream->write(reinterpret_cast<char*>(m_counts_ratio.data()),
                   size*sizeof(double));
  #endif
    }

    bool operator==(const KmerSign<MAX_K>& rhs) const { return m_kmer == rhs.m_kmer; }

   public:
    km::Kmer<MAX_K> m_kmer;
    double m_pvalue{0};
    Significance m_sign{Significance::NO};
  #ifdef WITH_POPSTRAT
    std::vector<double> m_counts_ratio;
  #endif
    double m_mean_control;
    double m_mean_case;
  };

  template<size_t MAX_K>
  inline bool operator<(const KmerSign<MAX_K>& lhs, const KmerSign<MAX_K>& rhs)
  {
    return lhs.m_pvalue > rhs.m_pvalue;
  }

} // end of namespace kmdiff

template <size_t MAX_K>
struct std::hash<km::Kmer<MAX_K>>
{
  uint64_t operator()(const km::Kmer<MAX_K>& kmer) const noexcept
  {
    return static_cast<uint64_t>(XXH64(kmer.m_data8, sizeof(kmer.m_data8), 0));
  }
};

template <size_t MAX_K>
struct std::hash<kmdiff::KmerSign<MAX_K>>
{
  uint64_t operator()(const kmdiff::KmerSign<MAX_K>& kmer) const noexcept
  {
    return static_cast<uint64_t>(XXH64(kmer.m_kmer.m_data8, sizeof(kmer.m_kmer.m_data8), 0));
  }
};
