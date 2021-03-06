#pragma once

#include <cmath>
#include <cstdint>
#include <vector>

namespace kmdiff {

  class LogFactorialTable
  {
  public:
    LogFactorialTable(size_t size);

    double operator[](size_t i)
    {
      if (i < m_size) return m_table[i];
      return log_factorial(i);
    }

  private:
    double log_factorial(size_t k);

  private:
    std::vector<double> m_table;
    size_t m_size;
  };

} // end of namespace kmdiff

