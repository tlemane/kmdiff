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
#define private public
#include <kmdiff/kmer.hpp>
#define _KM_LIB_INCLUDE_
#include <kmtricks/lz4_stream.hpp>

using namespace kmdiff;

TEST(kemr, set_from_str)
{  
  std::string a = random_dna_seq(20);
  std::string b = random_dna_seq(40);
  std::string c = random_dna_seq(240);
  
  Kmer<32> kmer(a);
  Kmer<64> kmer2(b);
  Kmer<256> kmer3(c);

  EXPECT_EQ(a , kmer.to_string());
  EXPECT_EQ(b , kmer2.to_string());
  EXPECT_EQ(c , kmer3.to_string());
}

TEST(kmer, operator)
{
  const std::string str1 = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAT";
  const std::string str2 = "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTA";

  Kmer<64> kmer1(str1);
  Kmer<64> kmer2(str2);
  Kmer<64> kmer3(str1);

  EXPECT_TRUE(kmer1 < kmer2);
  EXPECT_FALSE(kmer1 > kmer2);
  EXPECT_TRUE(kmer2 > kmer1);
  EXPECT_FALSE(kmer2 < kmer1);
  EXPECT_TRUE(kmer1 == kmer3);
}

TEST(kmer, serialize)
{
  std::string r = random_dna_seq(20);
  Kmer<32> kmer(r);

  {
    std::shared_ptr<std::ofstream> out =
      std::make_shared<std::ofstream>("tests_tmp/test.out", std::ios::out | std::ios::binary);
    kmer.dump(out);
  }
  {
    Kmer<32> k;
    std::shared_ptr<std::ifstream> in =
      std::make_shared<std::ifstream>("tests_tmp/test.out", std::ios::in | std::ios::binary);
    k.load(in, 20);
    EXPECT_EQ(kmer, k);
  }
}

TEST(kmer, serialize_lz4)
{
  std::string r = random_dna_seq(20);
  std::vector<Kmer<32>> v;
  
  for (size_t i=0; i<10; i++)
    v.push_back(Kmer<32>(random_dna_seq(20)));

  {
    std::ofstream out_n("tests_tmp/test.out.lz4", std::ios::out | std::ios::binary);
    std::shared_ptr<lz4_stream::basic_ostream<4096>> out =
      std::make_shared<lz4_stream::basic_ostream<4096>>(out_n);
    for (auto& kmer: v) kmer.dump(out);
  }
  {
    Kmer<32> k;
    std::shared_ptr<lz4_stream::basic_istream<4096>> in =
      std::make_shared<lz4_stream::basic_istream<4096>>(
         "tests_tmp/test.out.lz4");
    for (auto& kmer : v)
    {
      k.load(in, 20);
      EXPECT_EQ(kmer, k);
    }
  }
}

TEST(kmerSign, kmerSign)
{
  std::string r = random_dna_seq(20);
  Kmer<32> kmer(r);
  KmerSign<32> kmer_sign(std::move(kmer), 0.01, Significance::CONTROL);

  EXPECT_EQ(kmer_sign.m_kmer.to_string(), r);
  EXPECT_EQ(kmer_sign.m_pvalue, 0.01);
  EXPECT_EQ(kmer_sign.m_sign, Significance::CONTROL);
}

TEST(kmerSign, serial)
{
  std::string r = random_dna_seq(20);
  Kmer<32> kmer(r);
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
  Kmer<32> kmer(r);
  KmerSign<32> kmer_sign(std::move(kmer), 0.01, Significance::CONTROL);

  {
    std::ofstream out_n("./tests_tmp/test2.out.lz4", std::ios::out | std::ios::binary);
    std::shared_ptr<lz4_stream::basic_ostream<4096>> out =
      std::make_shared<lz4_stream::basic_ostream<4096>>(out_n);
    kmer_sign.dump(out);
  }
  {
    KmerSign<32> k;
    std::shared_ptr<lz4_stream::basic_istream<4096>> in =
      std::make_shared<lz4_stream::basic_istream<4096>>("tests_tmp/test2.out.lz4");
    k.load(in, 20);
    EXPECT_EQ(kmer_sign, k);
  }
}