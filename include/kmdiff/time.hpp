#pragma once

#include <string>
#include <chrono>
#include <functional>

namespace kmdiff {

  class Timer
  {
    using time_point_t = std::chrono::time_point<std::chrono::steady_clock>;
    using days = std::chrono::duration<int, std::ratio<86400>>;

   public:
    Timer();

    template <typename Unit>
    auto elapsed()
    {
      if (m_running) end();
      return std::chrono::duration_cast<Unit>(m_end_time - m_start_time);
    }

    template <typename Unit>
    static auto time_it(std::function<void()> func)
    {
      auto start = std::chrono::steady_clock::now();
      func();
      auto end = std::chrono::steady_clock::now();
      return std::chrono::duration_cast<Unit>(end - start).count();
    }

    void reset();

    std::string formatted();

   private:
    void start();
    void end();

   private:
    time_point_t m_start_time{};
    time_point_t m_end_time{};
    bool m_running{false};
  };

} // end of namespace kmdiff

