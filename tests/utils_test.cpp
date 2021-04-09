#include <gtest/gtest.h>
#include <kmdiff/utils.hpp>

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
  std::cerr << get_binary_dir() << std::endl;
  EXPECT_NO_THROW(command_exists(dir, "kmtricks-sv-tests"));
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