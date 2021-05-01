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
#define private public
#include <kmdiff/model.hpp>

using namespace kmdiff;

TEST(model, convert)
{
  EXPECT_EQ(correction_type_str(CorrectionType::NOTHING), "NOTHING");
  EXPECT_EQ(correction_type_str(CorrectionType::BONFERRONI), "BONFERRONI");
  EXPECT_EQ(correction_type_str(CorrectionType::BENJAMINI), "BENJAMINI");
}

TEST(model, model)
{
  Model<8> model;

  std::vector<uint8_t> v {2, 2, 0, 10, 0, 1, 0, 0, 42, 2};

  EXPECT_FLOAT_EQ(model.compute_mean(v), 5.9);
  std::tuple<double, size_t> t = std::make_tuple(5.9, 6);
  auto [mean, n] = model.compute_mean_e(v);

  EXPECT_FLOAT_EQ(mean, 5.9);
  EXPECT_EQ(n, 6);
}

TEST(model, poisson_likelihood)
{
  size_t nb_control = 30;
  size_t nb_case = 30;

  std::vector<uint32_t> v;
  std::vector<size_t> ct;
  for (size_t i=0; i<nb_control; i++)
  {
    v.push_back(200);
    ct.push_back(1);
  }
  for (size_t i=0; i<nb_case; i++)
    v.push_back(100);

  PoissonLikelihood<100000> p(nb_control, nb_case, ct, ct, 10);

  Range<uint32_t> r1(v, 0, nb_control);
  Range<uint32_t> r2(v, nb_control, nb_case);

  {
    auto [pvalue, sign, mc1, mc2] = p.process(r1, r2);
    EXPECT_EQ(sign, Significance::CONTROL);
  }
  {
    auto [pvalue, sign, mc1, mc2] = p.process(r2, r1);
    EXPECT_EQ(sign, Significance::CASE);
  }

  for (size_t i=0; i<nb_control; i++)
    v[i] = 100;

  {
    auto [pvalue, sign, mc1, mc2] = p.process(r1, r2);
    EXPECT_EQ(sign, Significance::NO);
  }
}

TEST(model, icorrector)
{
  ICorrector c;
  EXPECT_TRUE(c.apply(0.0));
}

TEST(model, bonferroni)
{
  Bonferroni b(0.05, 100);
  EXPECT_TRUE(b.apply(0.0004));
  EXPECT_FALSE(b.apply(0.0006));
}

TEST(model, benjamini)
{
  BenjaminiHochberg b(0.25, 25);
  EXPECT_TRUE(b.apply(0.001));
  EXPECT_FALSE(b.apply(0.021));
}
