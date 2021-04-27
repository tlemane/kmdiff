/*****************************************************************************
 *   kmdiff
 *   Authors: T. Lemane
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

#pragma once

// std
#include <sys/types.h>
#include <sys/wait.h>

#include <array>
#include <chrono>
#include <filesystem>
#include <iterator>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>
#include <unistd.h>
#include <sys/resource.h>

#if __APPLE__
#include <mach-o/dyld.h>
#include <mach/mach.h>
#endif

// ext
#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include <bcli/bcli.hpp>

// int
#include <kmdiff/exceptions.hpp>

namespace fs = std::filesystem;

namespace kmdiff
{
enum class VerbosityLevel
{
  DEBUG,
  INFO,
  WARNING,
  ERROR
};

std::string command_exists(const std::string& dir, const std::string& cmd);

std::string get_binary_dir();

std::string get_uname_sr();

VerbosityLevel str_to_verbosity_level(const std::string& str_level);

void set_verbosity_level(const std::string& level);

int exec_external_cmd(const std::string& cmd, const std::string& args);

std::string random_dna_seq(size_t size);

size_t get_current_rss();

size_t get_peak_rss();

std::string& str_to_upper(std::string& s);

template <template<typename> typename Container, typename T>
void destroy_container(Container<T>& container)
{
  Container<T>().swap(container);
}

template <typename T>
class Range
{
  std::vector<T>& m_data;
  size_t m_start;
  size_t m_size;

 public:
  Range(std::vector<T>& data, size_t start, size_t size)
      : m_data(data), m_start(start), m_size(size)
  {
  }

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
};

template <typename T>
std::vector<T> slice(const std::vector<T>& v, int i, int j)
{
  if (j > v.size() - 1) j = v.size() - 1;
  return std::vector<T>(v.cbegin() + i, v.cbegin() + j + 1);
}

class Timer
{
  using time_point_t = std::chrono::time_point<std::chrono::steady_clock>;

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

 private:
  void start();
  void end();

 private:
  time_point_t m_start_time{};
  time_point_t m_end_time{};
  bool m_running{false};
};

template <typename T>
struct has_insert
{
  template <typename U>
  static uint8_t test(decltype(&U::insert));
  template <typename U>
  static uint16_t test(...);

 public:
  enum
  {
    value = sizeof(test<T>(0)) == sizeof(uint8_t)
  };
};

template <typename T>
struct has_push_back
{
  template <typename U>
  static uint8_t test(decltype(&U::push_back));
  template <typename U>
  static uint16_t test(...);

 public:
  enum
  {
    value = sizeof(test<T>(0)) == sizeof(uint8_t)
  };
};

template <typename T>
struct has_dump
{
  template <typename U>
  static uint8_t test(decltype(&U::dump));
  template <typename U>
  static uint16_t test(...);

 public:
  enum
  {
    value = sizeof(test<T>(0)) == sizeof(uint8_t)
  };
};

template <typename T>
struct has_load
{
  template <typename U>
  static uint8_t test(decltype(&U::load));
  template <typename U>
  static uint16_t test(...);

 public:
  enum
  {
    value = sizeof(test<T>(0)) == sizeof(uint8_t)
  };
};

};  // end of namespace kmdiff