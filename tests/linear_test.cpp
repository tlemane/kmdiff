#include <gtest/gtest.h>
#include <kmdiff/linear_model.hpp>
#include <cmath>

using namespace kmdiff;

TEST(linear, utils)
{
  vector_t v = {1, 0, 0, 1};
  vector_t v2 = {1, 0, 1, 1};

  EXPECT_TRUE(is_equal_v(v, v));
  EXPECT_FALSE(is_equal_v(v, v2));

  matrix_t m = {
    {0, 1, 1, 0},
    {0, 1, 1, 0},
    {0, 1, 1, 0},
  };

  matrix_t m2 = {
    {0, 1, 1, 0},
    {0, 1, 0, 0},
    {0, 1, 1, 0},
  };

  EXPECT_TRUE(is_equal_d(0.1, 0.2, 0.11));
  EXPECT_FALSE(is_equal_d(0.1, 0.2, 0.09));
  EXPECT_TRUE(is_equal_d(sigmoid(1.0), 0.7310585786300048792512));
  EXPECT_EQ(linear_predictor({1, 2, 3}, {1, 2, 3}), 14);
  EXPECT_TRUE(is_equal_d(predict({1, 2, 3}, {1, 2, 3}), 0.9999991684719723358679));
  EXPECT_TRUE(is_equal_m(m, m));
  EXPECT_FALSE(is_equal_m(m, m2));
  EXPECT_EQ(nrows(m), 3);
  EXPECT_EQ(ncols(m), 4);
}

TEST(linear, transpose)
{
  {
    matrix_t m = {
      {0, 1, 1, 0},
      {0, 1, 1, 0},
      {0, 1, 1, 0},
      {0, 1, 1, 0},
    };

    matrix_t t = {
      {0, 0, 0, 0},
      {1, 1, 1, 1},
      {1, 1, 1, 1},
      {0, 0, 0, 0},
    };

    matrix_t tt = transpose(m);

    EXPECT_TRUE(is_equal_m(t, tt));
  }

  {
    matrix_t m = {
      {0, 1, 1, 0},
      {0, 1, 1, 0},
      {0, 1, 1, 0},
    };

    matrix_t t = {
      {0, 0, 0},
      {1, 1, 1},
      {1, 1, 1},
      {0, 0, 0},
    };

    matrix_t tt = transpose(m);

    EXPECT_TRUE(is_equal_m(t, tt));
  }
}

TEST(linear, multiply)
{
  matrix_t m = {
    {1, 2, 1, 1},
    {1, 1, 6, 1},
    {1, 0, 1, 0},
  };

  matrix_t m2 = {
    {0, 1, 0, 0},
    {1, 2, 1, 1},
    {1, 1, 6, 1},
    {1, 0, 1, 0},
  };

  matrix_t mm = {
    {4, 6, 9, 3},
    {8, 9, 38, 7},
    {1, 2, 6, 1},
  };

  matrix_t mm_ = multiply(m, m2);
  EXPECT_TRUE(is_equal_m(mm, mm_));
}

TEST(linear, lu_decomposition)
{
  matrix_t m = {
    {1, 2, 1, 1},
    {1, 1, 6, 1},
    {1, 0, 1, 0},
    {1, 0, 1, 1},
  };

  matrix_t lower = {
    {1, 0, 0, 0},
    {1, 1, 0, 0},
    {1, 2, 1, 0},
    {1, 2, 1, 1},
  };

  matrix_t upper = {
    {1, 2, 1, 1},
    {0, -1, 5, 0},
    {0, 0, -10, -1},
    {0, 0, 0, 1},
  };

  matrix_t inv = {
    {0.1, -0.2, 1, 0.1},
    {0.5, 0, 0, -0.5},
    {-0.1, 0.2, 0, -0.1},
    {0, 0, -1, 1}
  };

  auto [lower_, upper_] = lu_decomposition(m, nrows(m));

  EXPECT_TRUE(is_equal_m(lower, lower_));
  EXPECT_TRUE(is_equal_m(upper, upper_));

  double epsilon = 1e-15;
  auto [inv_, a, b] = inverse(m, nrows(m));
  for (size_t i=0; i<nrows(inv); i++)
  {
    for (size_t j=0; j<ncols(inv); j++)
    {
      double x = std::signbit(inv[i][j]) ? inv[i][j]+(+0.0) : inv[i][j];
      double y = std::signbit(inv_[i][j]) ? inv_[i][j]+(+0.0) : inv_[i][j];
      EXPECT_TRUE(std::fabs(x-y)<epsilon);
    }
  }
}
