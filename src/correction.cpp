#include <string>
#include <kmdiff/correction.hpp>

namespace kmdiff {

  std::string correction_type_str(const CorrectionType type)
  {
    switch (type)
    {
      case CorrectionType::NOTHING:
        return "NOTHING";
        break;
      case CorrectionType::BONFERRONI:
        return "BONFERRONI";
        break;
      case CorrectionType::BENJAMINI:
        return "BENJAMINI";
        break;
      default:
        return "";
        break;
    }
  }

}
