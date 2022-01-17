#include <kmdiff/time.hpp>

#include <sstream>
#include <iomanip>

namespace kmdiff {

  Timer::Timer() { start(); }

  void Timer::reset()
  {
    m_running = false;
    m_start_time = time_point_t{};
    m_end_time = time_point_t{};
    start();
  }

  void Timer::start()
  {
    m_running = true;
    m_start_time = std::chrono::steady_clock::now();
  }

  void Timer::end()
  {
    m_running = false;
    m_end_time = std::chrono::steady_clock::now();
  }

  std::string Timer::formatted()
  {
    if (m_running) end();

    std::chrono::seconds seconds(std::chrono::duration_cast<std::chrono::seconds>(m_end_time - m_start_time));
    auto d = std::chrono::duration_cast<days>(seconds); seconds -= d;
    auto h = std::chrono::duration_cast<std::chrono::hours>(seconds); seconds -= h;
    auto m = std::chrono::duration_cast<std::chrono::minutes>(seconds); seconds -= m;
    auto s = std::chrono::duration_cast<std::chrono::seconds>(seconds);

    std::stringstream ss; ss.fill('0');
    if (d.count())
      ss << std::setw(2) << d.count() << "d";
    if (h.count())
      ss << std::setw(2) << h.count() << "h";
    if (m.count())
      ss << std::setw(2) << m.count() << "m";
    ss << std::setw(2) << s.count() << "s";
    return ss.str();
  }

} // end of namespace kmdiff
