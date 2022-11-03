#include <iostream>
#include <cmath>

#include <gtest/gtest.h>
#include <kmdiff/corrector.hpp>

using namespace kmdiff;

TEST(corrector, basic_threshold)
{
  basic_threshold c(0.05);
  EXPECT_TRUE(c.apply(0.04));
  EXPECT_FALSE(c.apply(0.06));
}

TEST(corrector, bonferroni)
{
  bonferroni c(0.05, 100);
  EXPECT_TRUE(c.apply(0.0004));
  EXPECT_FALSE(c.apply(0.0006));
}

TEST(corrector, benjamini)
{
  benjamini c(0.25, 25);
  EXPECT_TRUE(c.apply(0.009));
  EXPECT_FALSE(c.apply(0.02));
}

TEST(corrector, sidak)
{
  sidak c(0.05, 100);
  EXPECT_TRUE(c.apply(0.00050));
  EXPECT_FALSE(c.apply(0.00052));
}

TEST(corrector, holm)
{
  holm c(0.05, 100);
  for (std::size_t _ = 0; _ < 90; _++)
    c.apply(0);

  EXPECT_TRUE(c.apply(0.004));
  EXPECT_FALSE(c.apply(0.006));
}
