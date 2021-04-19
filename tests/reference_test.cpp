#include <filesystem>
#include <gtest/gtest.h>
#include <kmdiff/kmer.hpp>
#include <kmdiff/reference.hpp>
#include <kmdiff/simulator.hpp>

namespace fs = std::filesystem;
using namespace kmdiff;

TEST(Reference, init)
{
  Reference r("./data_test/1.fasta");
}

TEST(Reference, fetch)
{
  Reference r("./data_test/1.fasta");
  std::string seq = r.fetch("chr1", 4, 7);
}

TEST(Reference, sv)
{
  Reference r("./data_test/1.fasta");
  {
    SV sv("chr1", visor_map.at("DEL"), "", 15, 16);
    auto [c1, c2] = r.get_kmer_view(sv, 4);
    EXPECT_EQ(c1, "AACTTCAA");
    EXPECT_EQ(c2, "AACCAA");
  }
  {
    SV sv("chr2", visor_map.at("INS"), "TT", 13, 14);
    auto [c1, c2] = r.get_kmer_view(sv, 4);
    EXPECT_EQ(c1, "AACCAA");
    EXPECT_EQ(c2, "AACTTCAA");
  }
  {
    SV sv("chr3", visor_map.at("INV"), "", 14, 17);
    auto [c1, c2] = r.get_kmer_view(sv, 4);
    EXPECT_EQ(c1, "AAAGTAGAAA");
    EXPECT_EQ(c2, "AAACTACAAA");
  }
  {
    SV sv("chr4", visor_map.at("DUP"), "2", 14, 17);
    auto [c1, c2] = r.get_kmer_view(sv, 4);
    EXPECT_EQ(c1, "AAAGTAGAAA");
    EXPECT_EQ(c2, "AAAGTAGGTAGAAA");
  }
  {
    SV sv("chr5", visor_map.at("IDUP"), "2", 14, 17);
    auto [c1, c2] = r.get_kmer_view(sv, 4);
    EXPECT_EQ(c1, "AAAGTAGAAA");
    EXPECT_EQ(c2, "AAAGTAGCTACAAA");
  }
}