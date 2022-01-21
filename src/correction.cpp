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
      case CorrectionType::SIDAK:
        return "SIDAK";
        break;
      case CorrectionType::HOLM:
        return "HOLM";
        break;
      default:
        return "";
        break;
    }
  }

} // end of namespace kmdiff

