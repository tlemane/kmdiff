#include <gtest/gtest.h>
#include <kmdiff/utils.hpp>
#include <kmdiff/kmer.hpp>
#include <kmdiff/kff_utils.hpp>

using namespace kmdiff;

TEST(kff, read_write)
{
  size_t kmer_size = 20;
  std::vector<std::string> vs;
  for (size_t i=0; i<100; i++)
    vs.push_back(random_dna_seq(kmer_size));
  {
    KffWriter f("./tests_tmp/test.kff", 20);
    for (auto& s: vs)
    {
      Kmer<32> k(s);
      f.write(k);
    }
    f.close();
  }

  {
    KffReader f("./tests_tmp/test.kff", kmer_size);
    int i = 0;
    while (std::optional<Kmer<32>> k = f.read<32>())
    {
      EXPECT_EQ((*k).to_string(), vs[i]); i++;
    }
  }
}