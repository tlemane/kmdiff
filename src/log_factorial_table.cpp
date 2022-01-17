#include <kmdiff/log_factorial_table.hpp>

namespace kmdiff {

  LogFactorialTable::LogFactorialTable(size_t size)
    : m_size(size)
  {
    m_table.reserve(m_size);
    for (int i=0; i<size; i++)
      m_table.push_back(log_factorial(i));
  }

  double LogFactorialTable::log_factorial(int k)
  {
    double res = 0;
    while (k > 1)
    {
      res += log(k);
      k--;
    }
    return res;
  }

} // end of namespace kmdiff
