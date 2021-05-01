#include <gtest/gtest.h>
#include <kmdiff/log_factorial_table.hpp>
#include <kmdiff/linear_model.hpp>

using namespace kmdiff;

TEST(lft, lft)
{
  LogFactorialTable t(50);
  EXPECT_EQ(t[0], 0);
  EXPECT_EQ(t[1], 0);
  EXPECT_EQ(t[10], 15.104412573075514);
  EXPECT_EQ(t[50], 148.47776695177302);
  EXPECT_EQ(t[51], 152.40959258449737);
  EXPECT_EQ(t[100], 363.7393755555635);
}