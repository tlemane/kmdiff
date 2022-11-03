#include <iostream>
#include <gtest/gtest.h>

#include <kmdiff/utils.hpp>
#include <kmdiff/merge.hpp>
#include <kmdiff/kmtricks_utils.hpp>
#include <vector>
#include <string>

using namespace kmdiff;

TEST(merge, merge)
{
  const std::string km = "./data_test/km_out_dir";
  const std::string outd = "./tests_tmp/out";

  fs::create_directory(outd);

  kmtricks_config_t config = get_kmtricks_config(km);

  auto part_paths = get_partition_paths(km, config.nb_partitions);

  std::vector<acc_t<KmerSign<32>>> accs(config.nb_partitions);
  for (std::size_t i = 0; i < accs.size(); i++)
    accs[i] = std::make_shared<VectorAccumulator<KmerSign<32>>>(100);

  auto [total_controls, total_cases] = get_total_kmer(km, 1, 1, 1);


  std::shared_ptr<IModel<65536+1>> model {nullptr};
  model = std::make_shared<PoissonLikelihood<65536+1>>(1, 1, total_controls, total_cases, 100);

  std::vector<uint32_t> a_min(1+1, 1);

  global_merge<32, 65536+1> merger(
      part_paths, a_min, model, accs, config.kmer_size, 1, 1, 0.05/10000, 1, nullptr);

  auto T = merger.merge();
  EXPECT_EQ(total_controls[0], 160);
  EXPECT_EQ(total_cases[0], 160);
  EXPECT_EQ(T, 320);

  auto [c1, c2] = merger.signs();
  EXPECT_EQ(c1, 0);
  EXPECT_EQ(c2, 0);
}
