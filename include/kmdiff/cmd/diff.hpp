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
#include <kmdiff/time.hpp>
#include <kmdiff/corrector.hpp>

#ifdef WITH_PLUGIN
  #include <kmdiff/model_manager.hpp>
#endif

#define KMTRICKS_PUBLIC
#include <kmtricks/kmdir.hpp>

namespace kmdiff {

  template<std::size_t KSIZE>
  void main_diff(kmdiff_options_t options)
  {
    diff_options_t opt = std::static_pointer_cast<struct diff_options>(options);
    spdlog::debug(opt->display());

    #ifdef WITH_PLUGIN
      if (!opt->model_lib_path.empty())
        plugin_manager<IModel<DMAX_C>>::get().init(opt->model_lib_path, opt->model_config);
    #endif

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
      {
        accumulators[i] = std::make_shared<VectorAccumulator<KmerSign<KSIZE>>>(65536);
      }
      else
      {
        accumulators[i] = std::make_shared<FileAccumulator<KmerSign<KSIZE>>>(
          fmt::format("{}/acc_{}", output_part_dir, i), config.kmer_size);
      }
    }

    std::vector<std::uint32_t> ab_mins(opt->nb_controls + opt->nb_cases, 1);

    auto [total_controls, total_cases] = get_total_kmer(opt->kmtricks_dir, opt->nb_controls, opt->nb_cases, config.abundance_min);

    spdlog::debug("\nNb k-mers controls: {}\n Nb k-mers cases: {}",
                  str_vector(total_controls),
                  str_vector(total_cases));

    std::shared_ptr<IModel<DMAX_C>> model {nullptr};

    if (opt->model_lib_path.empty())
    {
      model = std::make_shared<PoissonLikelihood<DMAX_C>>(
        opt->nb_controls, opt->nb_cases, total_controls, total_cases, 500);
    }
    #ifdef WITH_PLUGIN
      else
      {
        if (opt->pop_correction)
        {
          spdlog::warn("population stratification correction disabled with custom models.");
        }

        opt->pop_correction = false;
        model = plugin_manager<IModel<DMAX_C>>::get().get_plugin();
      }
    #endif

    std::string pop_dir = fmt::format("{}/popstrat", opt->output_directory);
    fs::create_directory(pop_dir);

    std::string gwas_eigenstratX_geno = fmt::format("{}/gwas_eigenstratX.geno", pop_dir);
    std::string gwas_eigenstratX_snp = fmt::format("{}/gwas_eigenstratX.snp", pop_dir);

    eig_geno_t<DMAX_C> geno {nullptr};
    eig_snp_t snp {nullptr};
    std::shared_ptr<Sampler<DMAX_C>> sampler {nullptr};

    if (opt->pop_correction)
    {
      geno = std::make_shared<EigGenoFile<DMAX_C>>(gwas_eigenstratX_geno);
      snp = std::make_shared<EigSnpFile>(gwas_eigenstratX_snp);
      sampler = std::make_shared<Sampler<DMAX_C>>(geno, snp, opt->kmer_pca);
    }

    global_merge<KSIZE, DMAX_C> merger(
      part_paths, ab_mins, model, accumulators, config.kmer_size, opt->nb_controls,
      opt->nb_cases, opt->threshold/opt->cutoff, opt->nb_threads, sampler);

    std::size_t total_kmers = merger.merge();
    auto [sign_controls, sign_cases] = merger.signs();

    spdlog::info("Partitions processed ({})", merge_time.formatted());

    spdlog::info("{}/{} significant k-mers.", merger.nb_sign(), total_kmers);
    spdlog::info("Before correction: {} (control), {} (case).", sign_controls, sign_cases);

    if (opt->learning_rate) pop_strat_corrector::s_learn_rate = opt->learning_rate;
    if (opt->max_iteration) pop_strat_corrector::s_max_iter = opt->max_iteration;
    if (opt->epsilon) pop_strat_corrector::s_epsilon = opt->epsilon;
    if (opt->stand) pop_strat_corrector::s_stand = opt->stand;
    if (opt->irls) pop_strat_corrector::s_irls = opt->irls;

    if (opt->pop_correction)
    {
      geno->close();
      snp->close();
    }

    #ifdef WITH_POPSTRAT
      std::vector<acc_t<KmerSign<KSIZE>>> pop_accumulators;
      if (opt->pop_correction)
      {
        std::string gwas_eigenstratX_ind = fmt::format("{}/gwas_eigenstratX.ind", pop_dir);
        std::string gwas_eigenstratX_total = fmt::format("{}/gwas_eigenstratX.total", pop_dir);
        std::string pcs_evec = fmt::format("{}/pcs.evec", pop_dir);

        spdlog::info("PCA for population stratification correction...");

        Timer pca_time;

        std::string parfile_path = fmt::format("{}/parfile.txt", pop_dir);
        std::string gwas_info_path = fmt::format("{}/gwas_infos.txt", pop_dir);
        std::string fof = fmt::format("{}/kmtricks.fof", opt->kmtricks_dir);

        write_parfile(parfile_path);
        write_gwas_info(fof, gwas_info_path, opt->nb_controls, opt->nb_cases);
        write_gwas_info(fof, gwas_eigenstratX_ind, opt->nb_controls, opt->nb_cases);
        write_gwas_eigenstrat_total(gwas_eigenstratX_total, total_controls, total_cases);

        std::string log_eigenstrat = "eigenstrat.log";
        run_eigenstrat_smartpca(pop_dir, "parfile.txt", log_eigenstrat, opt->is_diploid);

        spdlog::info("PCA done. ({}).", pca_time.formatted());

        Timer pop_time;

        spdlog::info("Apply population stratification correction...");
        auto pop_corrector = std::make_shared<pop_strat_corrector>(
          opt->nb_controls, opt->nb_cases, total_controls, total_cases, opt->npc);

        pop_corrector->load_Z(pcs_evec);
        pop_corrector->load_Y(gwas_eigenstratX_ind);
        pop_corrector->load_C(opt->covariates);
        pop_corrector->load_ginfo(gwas_info_path);
        pop_corrector->init_global_features();

        pop_accumulators.resize(accumulators.size());

        for (std::size_t p = 0; p < accumulators.size(); p++)
        {
          if (opt->in_memory)
          {
            pop_accumulators[p] = std::make_shared<VectorAccumulator<KmerSign<KSIZE>>>(65536);
          }
          else
          {
            pop_accumulators[p] = std::make_shared<FileAccumulator<KmerSign<KSIZE>>>(
              fmt::format("{}/accp_{}", output_part_dir, p), config.kmer_size);
          }
        }

        pop_corrector->apply(accumulators, pop_accumulators, 1);

        accumulators.swap(pop_accumulators);

        spdlog::info("Population correction done. ({}).", pop_time.formatted());
      }
    #endif

    Timer agg_time;

    std::shared_ptr<ICorrector> corrector {nullptr};
    std::unique_ptr<IAggregator<KSIZE>> agg {nullptr};

    if (opt->correction == CorrectionType::NOTHING)
      spdlog::info("Aggregate partitions...");
    else
      spdlog::info("Aggregate partitions and apply significance correction...");

    indicators::ProgressBar* pb = nullptr;

    if ((spdlog::get_level() != spdlog::level::debug) && isatty_stderr())
    {
      pb = get_progress_bar("progress", config.nb_partitions, 50, indicators::Color::white, false);
      pb->print_progress();
    }

    if (opt->correction != CorrectionType::BENJAMINI && opt->correction != CorrectionType::HOLM)
    {
      if (opt->correction == CorrectionType::BONFERRONI)
      {
        corrector = std::make_shared<bonferroni>(opt->threshold, total_kmers);
      }
      else if (opt->correction == CorrectionType::SIDAK)
      {
        corrector = std::make_shared<sidak>(opt->threshold, total_kmers);
      }
      else
      {
        corrector = std::make_shared<basic_threshold>(opt->threshold);
      }
      agg = std::make_unique<aggregator<KSIZE>>(accumulators,
                                                corrector,
                                                config,
                                                opt->output_directory,
                                                opt->kff,
                                                opt->nb_threads,
                                                pb);

    }
    else
    {
      if (opt->correction == CorrectionType::BENJAMINI)
      {
        corrector = std::make_shared<benjamini>(opt->threshold, total_kmers);
      }
      else if (opt->correction == CorrectionType::HOLM)
      {
        corrector = std::make_shared<holm>(opt->threshold, total_kmers);
      }
      agg = std::make_unique<sorted_aggregator<KSIZE>>(accumulators,
                                                       corrector,
                                                       config,
                                                       opt->output_directory,
                                                       opt->kff,
                                                       opt->nb_threads,
                                                       pb);


    }

    agg->run();

    auto [c_controls, c_cases] = agg->counts();

    delete pb;
    spdlog::info("Partitions aggregated ({})", agg_time.formatted());
    spdlog::info("Significant k-mers: {} (control), {} (case).", c_controls, c_cases);

    spdlog::info(
      "Done in {}, Peak RSS -> {} MB.", whole_time.formatted(), static_cast<size_t>(get_peak_rss() * 0.0009765625)
    );
  }

} // end of namespace kmdiff

