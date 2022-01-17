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

#include <filesystem>
#include <sstream>
#include <string>
#include <vector>

#if __APPLE__
#include <mach-o/dyld.h>
#include <mach/mach.h>
#endif

#include <fmt/format.h>

#include <kmdiff/exceptions.hpp>

namespace fs = std::filesystem;

namespace kmdiff {

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

  int exec_external_cmd(const std::string& cmd, const std::string& args,
                        const std::string& sout = "", const std::string& serr = "");

  std::string random_dna_seq(size_t size);

  size_t get_current_rss();

  size_t get_peak_rss();

  std::string& str_to_upper(std::string& s);

  template<typename T>
  void check_fstream_good(const std::string& path, const T& stream)
  {
    if (!stream.good())
    {
      if constexpr(std::is_same_v<T, std::ofstream>)
        throw IOError(fmt::format("Unable to write at {}.", path));
      else
        throw IOError(fmt::format("Unable to read at {}.", path));
    }
  }

  template <template<typename> typename Container, typename T>
  void destroy_container(Container<T>& container)
  {
    Container<T>().swap(container);
  }

  template <typename T>
  std::vector<T> slice(const std::vector<T>& v, int i, int j)
  {
    if (j > v.size() - 1) j = v.size() - 1;
    return std::vector<T>(v.cbegin() + i, v.cbegin() + j + 1);
  }

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

  template<typename, typename>
  struct is_same_template : std::false_type {};

  template<template<typename...> typename T,
           typename... A,
           typename... B>
  struct is_same_template<T<A...>, T<B...>> : std::true_type {};

  template<template<typename, std::size_t> typename T,
           typename TA, std::size_t SA,
           typename TB, std::size_t SB>
  struct is_same_template<T<TA, SA>, T<TB, SB>> : std::true_type {};

  template<typename T, typename U>
  constexpr auto is_same_template_v = is_same_template<T, U>::value;

};  // end of namespace kmdiff
