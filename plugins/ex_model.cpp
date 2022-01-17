#include <kmdiff/imodel.hpp>

template<std::size_t MAX_C>
class ExModel : public kmdiff::IModel<MAX_C>
{
  using base = kmdiff::IModel<MAX_C>;
  using count_type = typename base::count_type;

  public:
    ExModel() = default;

    virtual ~ExModel() {}

    void configure(const std::string& config) override {}

    kmdiff::model_ret_t process(const kmdiff::Range<count_type>& controls,
                                const kmdiff::Range<count_type>& cases) override
    {
      double mean_controls = base::mean(controls);
      double mean_cases = base::mean(cases);

      double pvalue = 0.000000000001;

      return std::make_tuple(pvalue,  kmdiff::Significance::CONTROL, mean_controls, mean_cases);
    }

};

extern "C" std::string plugin_name() { return "ExModel"; }
extern "C" ExModel<kmdiff::maxc8>* create8() { return new ExModel<kmdiff::maxc8>(); }
extern "C" ExModel<kmdiff::maxc16>* create16() { return new ExModel<kmdiff::maxc16>(); }
extern "C" ExModel<kmdiff::maxc32>* create32() { return new ExModel<kmdiff::maxc32>(); }
