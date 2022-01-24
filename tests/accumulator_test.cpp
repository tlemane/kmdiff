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
#include <kmdiff/kmer.hpp>
#include <kmdiff/accumulator.hpp>

using namespace kmdiff;

TEST(accumulator, VectorAccumulator)
{
  acc_t<int> acc = std::make_shared<VectorAccumulator<int>>(100);
  std::vector<int> v {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
  for (auto& e : v)
    acc->push(std::move(e));

  acc->finish();

  EXPECT_EQ(acc->size(), v.size());

  int i = 0;
  while (const std::optional<int>&o = acc->get())
  {
    EXPECT_EQ(*o, v[i]);
    i++;
  }
}

TEST(accumulator, SetAccumulator)
{
  acc_t<int> acc = std::make_shared<SetAccumulator<int>>(100);
  std::vector<int> v {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
  for (auto& e : v)
    acc->push(std::move(e));

  acc->finish();

  EXPECT_EQ(acc->size(), v.size());

  int i = 0;
  while (const std::optional<int>&o = acc->get())
  {
    EXPECT_TRUE(std::find(v.begin(), v.end(), *o) != v.end());
    i++;
  }
}

TEST(accumulator, FileAccumulator)
{
  acc_t<int> acc = std::make_shared<FileAccumulator<int>>("./tests_tmp/acc.txt.lz4");
  std::vector<int> v {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
  for (auto e : v)
    acc->push(std::move(e));

  acc->finish();

  EXPECT_EQ(acc->size(), v.size());

  int i = 0;
  while (const std::optional<int>&o = acc->get())
  {
    EXPECT_EQ(o.value(), v[i]);
    i++;
  }
}

TEST(accumulator, KmerSignVec)
{
  using kmer_sign_t = KmerSign<32>;
  acc_t<kmer_sign_t> acc = std::make_shared<VectorAccumulator<kmer_sign_t>>(100);

  std::vector<kmer_sign_t> v;
  for (size_t i=0; i<10; i++)
  {
    std::string r = random_dna_seq(20);
    km::Kmer<32> kmer(r);
    KmerSign<32> kmer_sign(std::move(kmer), 0.01, Significance::CONTROL);
    v.push_back(kmer_sign);
    EXPECT_EQ(r, kmer_sign.m_kmer.to_string());
  }
  for (auto e : v)
    acc->push(std::move(e));

  acc->finish();
  int i = 0;
  while (const std::optional<kmer_sign_t>&o = acc->get())
  {
    EXPECT_EQ(*o, v[i]);
    i++;
  }
}

