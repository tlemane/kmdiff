#pragma once

#include <indicators/progress_bar.hpp>

namespace kmdiff {
  inline indicators::ProgressBar* get_progress_bar(const std::string& name,
                                       size_t size,
                                       size_t width,
                                       indicators::Color color,
                                       bool time)
  {
    using namespace indicators;
    return new indicators::ProgressBar {
      option::Start{" ["},
      option::Fill{"■"},
      option::Lead{">"},
      option::Remainder{"·"},
      option::End{"]"},
      option::ForegroundColor{color},
      option::ShowPercentage{true},
      option::ShowElapsedTime{true},
      option::ShowRemainingTime{false},
      option::PrefixText{name},
      option::FontStyles{std::vector<FontStyle>{FontStyle::bold}},
      option::MaxProgress{size},
      option::BarWidth{width},
      option::Stream{std::cerr}
    };
  }
} // end of namespace kmdiff

