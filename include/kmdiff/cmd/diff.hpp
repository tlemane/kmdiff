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

#pragma once

// std
#include <csignal>
#include <memory>
#include <thread>
#include <string>
#include <vector>

// ext
#include <spdlog/spdlog.h>

// int
#include <kmdiff/cmd/cmd_common.hpp>
#include <kmdiff/config.hpp>
#include <kmdiff/kmtricks_utils.hpp>
#include <kmdiff/correction.hpp>
#include <kmdiff/aggregator.hpp>
#include <kmdiff/accumulator.hpp>
#include <kmdiff/merge.hpp>
#include <kmdiff/cmd/diff_opt.hpp>

#define KMTRICKS_PUBLIC
#include <kmtricks/kmdir.hpp>

namespace kmdiff
{

template<std::size_t KSIZE>
void main_diff(kmdiff_options_t options)
{
  diff_options_t opt = std::static_pointer_cast<struct diff_options>(options);
  spdlog::debug(opt->display());

  Timer whole_time;
  kmtricks_config_t config = get_kmtricks_config(opt->kmtricks_dir);

  std::string output_part_dir = fmt::format("{}/partitions", opt->output_directory);
  fs::create_directories(output_part_dir);

  spdlog::info("Process partitions...");
  Timer merge_time;

  km::KmDir::get().init(opt->kmtricks_dir, fmt::format("{}/kmtricks.fof", opt->kmtricks_dir));

  std::vector<std::vector<std::string>> part_paths;
  for (std::size_t i = 0; i < config.nb_partitions; i++)
    part_paths.push_back(km::KmDir::get().get_files_to_merge(i, true, km::KM_FILE::KMER));

  std::vector<acc_t<KmerSign<KSIZE>>> accumulators(config.nb_partitions);
  for (std::size_t i = 0; i < accumulators.size(); i++)
  {
    if (opt->in_memory)
      accumulators[i] = std::make_shared<VectorAccumulator<KmerSign<KSIZE>>>(65536);
    else
      accumulators[i] = std::make_shared<FileAccumulator<KmerSign<KSIZE>>>(
        fmt::format("{}/acc_{}", output_part_dir, i), config.kmer_size
      );
  }

  std::vector<std::uint32_t> ab_mins(opt->nb_controls + opt->nb_cases, 1);

  auto [total_controls, total_cases] = get_total_kmer(opt->kmtricks_dir, opt->nb_controls, opt->nb_cases);

  std::shared_ptr<Model<DMAX_C>> model = std::make_shared<PoissonLikelihood<DMAX_C>>(
    opt->nb_controls, opt->nb_cases, total_controls, total_cases, 500);

  global_merge<KSIZE, DMAX_C> merger(
    part_paths, ab_mins, model, accumulators, config.kmer_size, opt->nb_controls,
    opt->nb_cases, opt->threshold, opt->nb_threads);

  std::size_t total_kmers = merger.merge();

  spdlog::info("Partitions processed ({})", merge_time.formatted());

  spdlog::info("Found {} significant k-mers.", merger.nb_sign());

  std::shared_ptr<ICorrector> corrector {nullptr};
  std::unique_ptr<IAggregator<KSIZE>> aggregator {nullptr};

  spdlog::info("Aggregate and apply correction...");

  if (opt->correction == CorrectionType::BONFERRONI)
  {
    corrector = std::make_shared<Bonferroni>(opt->threshold / 100000, total_kmers);
    aggregator = std::make_unique<BonferonniAggregator<KSIZE>>(accumulators, corrector, opt, config);
  }
  else if (opt->correction == CorrectionType::BENJAMINI)
  {
    corrector = std::make_shared<BenjaminiHochberg>(opt->threshold / 100000, total_kmers);
    aggregator = std::make_unique<BenjaminiAggregator<KSIZE>>(accumulators, corrector, opt, config);
  }
  else
  {
    aggregator = make_uncorrected_aggregator<KSIZE>(accumulators, opt, config, nullptr);
  }

  aggregator->run();

  spdlog::info(
    "Done in {}, Peak RSS -> {} MB.", whole_time.formatted(), static_cast<size_t>(get_peak_rss() * 0.0009765625)
  );

}

};  // namespace kmdiff
