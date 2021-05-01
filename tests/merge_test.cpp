#include <gtest/gtest.h>
#include <kmdiff/merge.hpp>
#include <kmdiff/kmtricks_utils.hpp>
#include <vector>
#include <string>

using namespace kmdiff;

TEST(merge, merge)
{
  const std::string km = "./data_test/kmtricks-dir";
  const std::string outd = "./tests_tmp/out";
  fs::create_directory(outd);

  kmtricks_config_t config = get_kmtricks_config(km);
  std::vector<std::string> fofs = get_fofs(km);

  std::vector<acc_t<KmerSign<32>>> accumulators(config.nb_partitions);
  for (size_t i = 0; i < accumulators.size(); i++)
    accumulators[i] = std::make_shared<VectorAccumulator<KmerSign<32>>>(4096);

  auto [total_controls, total_cases] = get_total_kmer(km, 1, 1);

  std::shared_ptr<Model<8>> model = std::make_shared<PoissonLikelihood<8>>(
      1, 1, total_controls, total_cases, 10);

  std::vector<uint32_t> dummy_a_min = std::vector<uint32_t>(1+1, 1);

  GlobalMerge<32, 8> merger(fofs, dummy_a_min, config.kmer_size, 1, 1, 10.0, outd, 2, model, accumulators, "", "", false, 0.0);

  size_t s = merger.merge();
  EXPECT_EQ(s, 320);
}