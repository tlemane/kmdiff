#include <gtest/gtest.h>
#include <kmdiff/sv.hpp>

using namespace kmdiff;

TEST(sv, operator)
{
  SV sv("chr1", "insertion", "ACGT", 100, 101);
  SV sv1("chr1", "insertion", "ACGT", 100, 101);
  SV sv2("chr2", "insertion", "ACGT", 42, 43);

  EXPECT_EQ(sv, sv1);
  EXPECT_LT(sv2, sv1);
}

TEST(sv, bed)
{
  SV sv("chr1", "insertion", "ACGT", 100, 101);

  EXPECT_EQ(sv.to_bed_entry(), "chr1\t100\t101\tinsertion\tACGT\t0");
}

TEST(sv, hash)
{
  SV sv("chr1", "insertion", "ACGT", 100, 101);
  std::string s = "chr1insertionACGT100101";

  std::hash<SV> hasher;
  EXPECT_EQ(hasher(sv), XXH64(s.c_str(), s.size(), 0));
}