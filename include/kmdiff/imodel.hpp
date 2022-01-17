#pragma once

#include <limits>
#include <tuple>
#include <string>

#include <kmdiff/range.hpp>
#include <kmdiff/kmer.hpp>

#define KMTRICKS_PUBLIC
#include <kmtricks/utils.hpp>

namespace kmdiff {

  using pvalue_t = double;

  using model_ret_t = std::tuple<pvalue_t, Significance, double, double>;

  constexpr size_t maxc8 = std::numeric_limits<uint8_t>::max();
  constexpr size_t maxc16 = std::numeric_limits<uint16_t>::max();
  constexpr size_t maxc32 = std::numeric_limits<uint32_t>::max();

  template<std::size_t MAX_C>
  class IModel
  {
    public:
      using count_type = typename km::selectC<MAX_C>::type;
      using range_type = Range<count_type>;

      IModel() = default;

      virtual ~IModel() {}

      virtual void configure(const std::string& config) = 0;

      virtual model_ret_t process(const range_type& controls, const range_type& cases) = 0;

      template<typename Iterable>
      static double mean(const Iterable& iter)
      {
        return std::accumulate(iter.begin(), iter.end(), 0.0) / iter.size();
      }

      template<typename Iterable>
      static std::tuple<double, size_t> mean_count(const Iterable& iter)
      {
        auto [s, p] = sum_count(iter);
        return std::make_tuple(s / static_cast<double>(iter.size()), p);
      }

      template<typename Iterable>
      static std::tuple<double, size_t> sum_count(const Iterable& iter)
      {
        double s = 0;
        size_t pos = 0;
        for (auto& e : iter)
        {
          s += e;
          if (e > 0)
            pos++;
        }
        return std::make_tuple(s, pos);
      }

      template<typename Iterable>
      static double sd(const Iterable& iter)
      {
        double m = mean(iter);
        double s = std::inner_product(iter.begin(), iter.end(), iter.begin(), 0.0);
        return std::sqrt(s / iter.size() - m * m);
      }
  };

} // end of namespace kmdiff

