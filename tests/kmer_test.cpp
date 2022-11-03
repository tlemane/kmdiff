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

#include <string>
#include <gtest/gtest.h>
#include <kmdiff/utils.hpp>
#define private public
#define KMTRICKS_PUBLIC
#include <kmdiff/kmer.hpp>
#include <kmtricks/io/lz4_stream.hpp>

using namespace kmdiff;

TEST(kmer, serialize)
{
  std::string r = random_dna_seq(20);
  km::Kmer<32> kmer(r);

  {
    std::shared_ptr<std::ofstream> out =
      std::make_shared<std::ofstream>("tests_tmp/test.out", std::ios::out | std::ios::binary);
    kmer.dump(*out);
  }
  {
    km::Kmer<32> k; k.set_k(20);
    std::shared_ptr<std::ifstream> in =
      std::make_shared<std::ifstream>("tests_tmp/test.out", std::ios::in | std::ios::binary);
    k.load(*in);
    EXPECT_EQ(kmer, k);
  }
}

TEST(kmer, serialize_lz4)
{
  std::string r = random_dna_seq(20);
  std::vector<km::Kmer<32>> v;

  for (size_t i=0; i<10; i++)
    v.push_back(km::Kmer<32>(random_dna_seq(20)));

  {
    std::ofstream out_n("tests_tmp/test.out.lz4", std::ios::out | std::ios::binary);
    std::shared_ptr<std::ostream> out =
      std::make_shared<lz4_stream::basic_ostream<4096>>(out_n);
    for (auto& kmer: v)
      kmer.dump(*out);
  }
  {
    km::Kmer<32> k; k.set_k(20);
    std::shared_ptr<std::istream> in =
      std::make_shared<lz4_stream::basic_istream<4096>>(
         "tests_tmp/test.out.lz4");
    for (auto& kmer : v)
    {
      k.load(*in);
      EXPECT_EQ(kmer, k);
    }
  }
}

TEST(kmerSign, kmerSign)
{
  std::string r = random_dna_seq(20);
  km::Kmer<32> kmer(r);
  KmerSign<32> kmer_sign(std::move(kmer), 0.01, Significance::CONTROL);

  EXPECT_EQ(kmer_sign.m_kmer.to_string(), r);
  EXPECT_EQ(kmer_sign.m_pvalue, 0.01);
  EXPECT_EQ(kmer_sign.m_sign, Significance::CONTROL);
}

TEST(kmerSign, serial)
{
  std::string r = random_dna_seq(20);
  km::Kmer<32> kmer(r);
  KmerSign<32> kmer_sign(std::move(kmer), 0.01, Significance::CONTROL);

  {
    std::shared_ptr<std::ofstream> out =
      std::make_shared<std::ofstream>("./tests_tmp/test2.out", std::ios::out | std::ios::binary);
    kmer_sign.dump(out);
  }
  {
    KmerSign<32> k;
    std::shared_ptr<std::ifstream> in =
      std::make_shared<std::ifstream>("tests_tmp/test2.out", std::ios::in | std::ios::binary);
    k.load(in, 20);
    EXPECT_EQ(kmer_sign, k);
  }
}

TEST(kmerSign, serial_lz4)
{
  std::string r = random_dna_seq(20);
  km::Kmer<32> kmer(r);
  KmerSign<32> kmer_sign(std::move(kmer), 0.01, Significance::CONTROL);

  {
    std::ofstream out_n("./tests_tmp/test2.out.lz4", std::ios::out | std::ios::binary);
    std::shared_ptr<std::ostream> out =
      std::make_shared<lz4_stream::basic_ostream<4096>>(out_n);
    kmer_sign.dump(out);
  }
  {
    KmerSign<32> k;
    std::shared_ptr<std::istream> in =
      std::make_shared<lz4_stream::basic_istream<4096>>("tests_tmp/test2.out.lz4");
    k.load(in, 20);
    EXPECT_EQ(kmer_sign, k);
  }
}
