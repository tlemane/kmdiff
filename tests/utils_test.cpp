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

#include <gtest/gtest.h>
#include <spdlog/spdlog.h>

#include <kmdiff/range.hpp>
#include <kmdiff/utils.hpp>
#include <kmdiff/kmer.hpp>
#include <chrono>

using namespace kmdiff;

TEST(utils, commands_exists)
{
  std::string ls_cmd = command_exists(".", "ls");
  EXPECT_EQ(ls_cmd, "ls");

  EXPECT_THROW(command_exists(".", "aqszed"), BinaryNotFound);
}

TEST(utils, get_binary_dir)
{
  std::string dir = get_binary_dir();
  EXPECT_NO_THROW(command_exists(dir, "kmdiff-tests"));
}

TEST(utils, Range)
{
  std::vector<int> v{0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
  std::vector<int> a;

  for (auto& e : Range<int>(v, 2, 7))
    a.push_back(e);

  for (size_t i=0; i<a.size(); i++)
    EXPECT_EQ(a[i], v[i+2]);
}

TEST(utils, slice)
{
  std::vector<int> v{0, 1, 2, 3, 4, 5, 6, 7, 8, 9};

  std::vector<int> from_start = slice(v, 0, 4);
  std::vector<int> mid = slice(v, 2, 6);
  std::vector<int> to_end = slice(v, 5, 9);
  std::vector<int> too_large = slice(v, 0, 15);

  EXPECT_EQ(from_start.size(), 5);
  EXPECT_EQ(mid.size(), 5);
  EXPECT_EQ(to_end.size(), 5);
  EXPECT_EQ(too_large.size(), 10);

  for (size_t i=0; i<from_start.size(); i++)
    EXPECT_EQ(from_start[i], v[i]);
  for (size_t i=0; i<mid.size(); i++)
    EXPECT_EQ(mid[i], v[2+i]);
  for (size_t i=0; i<to_end.size(); i++)
    EXPECT_EQ(to_end[i], v[5+i]);
  for (size_t i=0; i<too_large.size(); i++)
    EXPECT_EQ(too_large[i], v[i]);
}

TEST(utils, verbosity)
{
  EXPECT_EQ(str_to_verbosity_level("debug"), VerbosityLevel::DEBUG);
  EXPECT_EQ(str_to_verbosity_level("info"), VerbosityLevel::INFO);
  EXPECT_EQ(str_to_verbosity_level("warning"), VerbosityLevel::WARNING);
  EXPECT_EQ(str_to_verbosity_level("error"), VerbosityLevel::ERROR);
  EXPECT_EQ(str_to_verbosity_level("def"), VerbosityLevel::WARNING);

  set_verbosity_level("debug");
  EXPECT_EQ(spdlog::get_level(), spdlog::level::debug);
  set_verbosity_level("info");
  EXPECT_EQ(spdlog::get_level(), spdlog::level::info);
  set_verbosity_level("warning");
  EXPECT_EQ(spdlog::get_level(), spdlog::level::warn);
  set_verbosity_level("error");
  EXPECT_EQ(spdlog::get_level(), spdlog::level::err);
}

TEST(utils, rss)
{
  EXPECT_GE(get_peak_rss(), 0);
  EXPECT_GE(get_current_rss(), 0);
}

TEST(utils, ext)
{
  EXPECT_NO_THROW(exec_external_cmd("ls", "", "./tests_tmp/ext_test"));
  EXPECT_THROW(exec_external_cmd("./ret1.sh", ""), ExternalExecFailed);
  EXPECT_GT(get_uname_sr().size(), 0);
}

TEST(utils, sfinae)
{
  EXPECT_FALSE(has_dump<std::vector<int>>::value);
  EXPECT_TRUE(has_dump<KmerSign<32>>::value);
  EXPECT_TRUE(has_load<KmerSign<32>>::value);
  EXPECT_FALSE(has_push_back<KmerSign<32>>::value);
  EXPECT_FALSE(has_insert<KmerSign<32>>::value);
}
