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
#include <algorithm>
#include <array>
#include <bitset>
#include <cstddef>
#include <cstdint>
#include <kmdiff/utils.hpp>
#include <sstream>
#include <vector>

// ext
#include <xxhash.h>
#define _KM_LIB_INCLUDE_
#include <kmtricks/utilities.hpp>

namespace kmdiff
{
const uint8_t bToN[] = {'A', 'C', 'T', 'G'};
const uint8_t revN[] = {'T', 'G', 'A', 'C'};
const uint8_t NToB[256] = {
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 0, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 0, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};


inline std::string str_rev_comp(const std::string& s)
{
  std::string rev;
  for (auto it=s.rbegin(); it!=s.rend(); it++)
    rev.push_back(revN[NToB[*it]]);
  return rev;
}

template <size_t MAX_K>
class Kmer
{
  friend struct std::hash<Kmer<MAX_K>>;
  using data_ptr64 = const uint64_t*;
  using data_ptr8 = const uint8_t*;

 private:
  const static uint16_t m_max_data{(MAX_K + 31) / 32};
  union
  {
    uint64_t m_data[m_max_data];
    uint8_t m_data8[m_max_data * 8];
  };
  inline static size_t m_kmer_size;
  inline static uint16_t m_n_data;

 public:
  Kmer() {}

  Kmer(const std::string& str_kmer)
  {
    set_k(str_kmer.size());
    set_from_str(str_kmer);
  }

  template <typename T>
  Kmer(const T* data, size_t size, size_t kmer_size)
  {
    set_from_data_2(data, size, kmer_size);
  }

  void set_from_str(const std::string& str_kmer)
  {
    std::fill(std::begin(m_data), std::end(m_data), 0);

    int i = 0;
    for (auto& c : str_kmer)
    {
      m_data[i / 32] <<= 2;
      m_data[i / 32] |= NToB[c];
      i++;
    }
  }

  void dump(const std::shared_ptr<std::ostream>& stream)
  {
    stream->write(reinterpret_cast<char*>(m_data8), sizeof(m_data8));
    return;
  }

  void load(const std::shared_ptr<std::istream>& stream, size_t kmer_size)
  {
    set_k(kmer_size);
    stream->read(reinterpret_cast<char*>(m_data8), sizeof(m_data8));
  }

  template <typename T>
  void from_km(T* data, size_t kmer_size)
  {
    set_k(kmer_size);
    if constexpr (sizeof(T) > 8)
    {
      m_data[1] = static_cast<uint64_t>(*data);
      m_data[0] = static_cast<uint64_t>((*data) >> 64);
    }
    else
    {
      m_data[0] = static_cast<uint64_t>(*data);
    }
  }

  void set_k(size_t kmer_size)
  {
    m_kmer_size = kmer_size;
    m_n_data = (m_kmer_size + 31) / 32;
  }

  data_ptr64 get_data64() const { return m_data; }

  data_ptr8 get_data8() const { return m_data8; }

  std::string to_bit_string() const
  {
    std::stringstream ss;
    for (uint16_t i = 0; i < m_max_data; i++)
      ss << std::bitset<sizeof(uint64_t) * 8>(m_data[i]).to_string();
    return ss.str();
  }

  std::string to_bit_string_real() const { return ""; }

  bool operator>(const Kmer<MAX_K>& k) const
  {
    for (int i = 0; i < m_max_data; i++)
    {
      if (m_data[i] > k.m_data[i]) return true;
      if (m_data[i] < k.m_data[i]) return false;
    }
    return false;
  }

  bool operator<(const Kmer<MAX_K>& k) const
  {
    for (int i = 0; i < m_max_data; i++)
    {
      if (m_data[i] < k.m_data[i]) return true;
      if (m_data[i] > k.m_data[i]) return false;
    }
    return false;
  }

  bool operator==(const Kmer<MAX_K>& k) const
  {
    for (int i = 0; i < m_max_data; i++)
      if (!(m_data[i] == k.m_data[i])) return false;
    return true;
  }

  std::string to_string() const
  {
    char tmp[m_kmer_size + 1];
    int shift = 0;

    for (int k = m_kmer_size - 1; k >= 0; k--)
    {
      tmp[k] = bToN[(m_data[k / 32] >> shift) & 3];

      if ((k % 32) == 0)
        shift = 0;
      else
        shift += 2;
    }
    tmp[m_kmer_size] = '\0';
    return tmp;
  }
};

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
  KmerSign(Kmer<MAX_K>&& kmer, double pvalue, Significance sign,
           std::vector<double>& counts_ratio, double mean_control = 0, double mean_case = 0)
      : m_kmer(std::move(kmer)), m_pvalue(pvalue), m_sign(sign), m_counts_ratio(counts_ratio),
        m_mean_control(mean_control), m_mean_case(mean_case)
  {
  }
#else
  KmerSign(Kmer<MAX_K>&& kmer, double pvalue, Significance sign,
           double mean_control = 0, double mean_case = 0)
      : m_kmer(std::move(kmer)), m_pvalue(pvalue), m_sign(sign),
        m_mean_control(mean_control), m_mean_case(mean_case)
  {
  }
#endif

  KmerSign() {}

  void load(std::shared_ptr<std::istream> stream, size_t size)
  {
    m_kmer.load(stream, size);
    stream->read(reinterpret_cast<char*>(&m_pvalue), sizeof(m_pvalue));
    stream->read(reinterpret_cast<char*>(&m_sign), sizeof(m_sign));
    stream->read(reinterpret_cast<char*>(&m_mean_control), sizeof(m_mean_control));
    stream->read(reinterpret_cast<char*>(&m_mean_case), sizeof(m_mean_case));
#ifdef WITH_POPSTRAT
    size_t size_ = 0;
    stream->read(reinterpret_cast<char*>(&size_), sizeof(size_));
    m_counts_ratio.resize(size_, 0);
    stream->read(reinterpret_cast<char*>(m_counts_ratio.data()),
                 size_*sizeof(double));
#endif
  }

  void dump(std::shared_ptr<std::ostream> stream)
  {
    m_kmer.dump(stream);
    stream->write(reinterpret_cast<char*>(&m_pvalue), sizeof(m_pvalue));
    stream->write(reinterpret_cast<char*>(&m_sign), sizeof(m_sign));
    stream->write(reinterpret_cast<char*>(&m_mean_control), sizeof(m_mean_control));
    stream->write(reinterpret_cast<char*>(&m_mean_case), sizeof(m_mean_case));
#ifdef WITH_POPSTRAT
    size_t size = m_counts_ratio.size();
    stream->write(reinterpret_cast<char*>(&size), sizeof(size));
    stream->write(reinterpret_cast<char*>(m_counts_ratio.data()),
                 size*sizeof(double));
#endif
  }

  bool operator==(const KmerSign<MAX_K>& rhs) const { return m_kmer == rhs.m_kmer; }

 public:
  Kmer<MAX_K> m_kmer;
  double m_pvalue{0};
  double m_mean_control;
  double m_mean_case;
  Significance m_sign{Significance::NO};
#ifdef WITH_POPSTRAT
  std::vector<double> m_counts_ratio;
  double m_corrected {0};
#endif
};

template<size_t MAX_K>
inline bool operator<(const KmerSign<MAX_K>& lhs, const KmerSign<MAX_K>& rhs)
{
  return lhs.m_pvalue > rhs.m_pvalue;
}

};  // namespace kmdiff

template <size_t MAX_K>
struct std::hash<kmdiff::Kmer<MAX_K>>
{
  uint64_t operator()(const kmdiff::Kmer<MAX_K>& kmer) const noexcept
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
