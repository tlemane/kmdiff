#include <gtest/gtest.h>
#include <fmt/format.h>
#include <kmdiff/kmtricks_utils.hpp>

using namespace kmdiff;

const std::string kmtricks_dir = "./data_test/kmtricks-dir";

TEST(kmtricks_utils, config)
{
  kmtricks_config_t config = get_kmtricks_config(kmtricks_dir);
  EXPECT_EQ(config.kmer_size, 20);
  EXPECT_EQ(config.nb_partitions, 4);
}

TEST(kmtricks_utils, fofs)
{
  std::vector<std::string> fofs = get_fofs(kmtricks_dir);
  std::string tmp = "{}/storage/kmers_partitions/partition_{}/{}";
  EXPECT_EQ(fofs[0], fmt::format(tmp, kmtricks_dir, 0, "partition0.fof"));
  EXPECT_EQ(fofs[1], fmt::format(tmp, kmtricks_dir, 1, "partition1.fof"));
  EXPECT_EQ(fofs[2], fmt::format(tmp, kmtricks_dir, 2, "partition2.fof"));
  EXPECT_EQ(fofs[3], fmt::format(tmp, kmtricks_dir, 3, "partition3.fof"));

  int i = 0;
  for (auto& e : fofs)
  {
    std::string line;
    std::ifstream in(e);
    std::getline(in, line);
    EXPECT_EQ(line, fmt::format(tmp, kmtricks_dir, i, "D1.kmer.lz4"));
    std::getline(in, line);
    EXPECT_EQ(line, fmt::format(tmp, kmtricks_dir, i, "D2.kmer.lz4"));
    i++;
  }
}

TEST(kmtricks_utils, total)
{
  std::tuple<std::vector<size_t>, std::vector<size_t>> t = get_total_kmer(
    kmtricks_dir, 1, 1
  );

  EXPECT_EQ(std::get<0>(t).size(), 1);
  EXPECT_EQ(std::get<1>(t).size(), 1);

  EXPECT_EQ(std::get<0>(t)[0], 160);
  EXPECT_EQ(std::get<1>(t)[0], 160);
}