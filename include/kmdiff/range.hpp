#pragma once

#include <vector>

namespace kmdiff {

  template <typename T>
  class Range
  {
    std::vector<T>& m_data;
    size_t m_start;
    size_t m_size;

   public:
    Range(std::vector<T>& data, size_t start, size_t size)
        : m_data(data), m_start(start), m_size(size) {}

    auto begin() const
    {
      auto start = m_data.begin();
      start += m_start;
      return start;
    }

    auto end() const
    {
      auto end = m_data.begin();
      end += m_start + m_size;
      return end;
    }

    size_t size() const { return m_size; }

    const T& operator[](size_t index) const { return m_data[m_start + index];}
  };

} // end of namespace kmdiff

