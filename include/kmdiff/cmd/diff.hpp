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

  inline void copy_kdir(const std::string& from, const std::string& to)
  {
    fs::create_directory(to);
    fs::create_directory(to + "/matrices");
    fs::create_directory(to + "/repartition_gatb");
    fs::create_directory(to + "/config_gatb");

    fs::copy(from + "/config_gatb", to + "/config_gatb", fs::copy_options::recursive);
    fs::copy(from + "/repartition_gatb", to + "/repartition_gatb", fs::copy_options::recursive);
    fs::copy(from + "/options.txt", to, fs::copy_options::recursive);
    fs::copy(from + "/kmtricks.fof", to, fs::copy_options::recursive);
  }

  template<std::size_t KSIZE>
  std::size_t do_diff(diff_options_t opt,
               const kmtricks_config_t& config,
               const std::string& output_part_dir,
               std::vector<acc_t<KmerSign<KSIZE>>>& accumulators,
               std::shared_ptr<Sampler<DMAX_C>> sampler)
  {
    Timer merge_time;

    spdlog::info("Process partitions");

    std::vector<std::vector<std::string>> part_paths;
    std::vector<std::string> matrix_paths;

    for (const auto& entry: fs::directory_iterator(km::KmDir::get().m_matrix_storage))
    {
      if (fs::exists(entry.path()))
      {
        matrix_paths.push_back(entry.path().string());
      }
    }

    bool from_matrix = false;

    if (matrix_paths.empty())
    {
      for (std::size_t i = 0; i < config.nb_partitions; i++)
      {
        part_paths.push_back(km::KmDir::get().get_files_to_merge(i, true, km::KM_FILE::KMER));
      }
    }
    else
    {
      from_matrix = true;
      std::ofstream out_opt_c(km::KmDir::get().m_root + "/kmdiff-count.opt");
    }

    for (std::size_t i = 0; i < accumulators.size(); i++)
    {
      accumulators[i] = std::make_shared<FileAccumulator<KmerSign<KSIZE>>>(
        fmt::format("{}/p{}_uncorrected", output_part_dir, i), config.kmer_size, false, !opt->keep_tmp);
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
        opt->nb_controls, opt->nb_cases, total_controls, total_cases, opt->log_size);
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

    std::string sign_matrix_dir = fmt::format("{}/positive_kmer_matrix", opt->output_directory);

    if (opt->save_sk)
    {
      copy_kdir(km::KmDir::get().m_root, sign_matrix_dir);
      sign_matrix_dir += "/matrices";
    }

    global_merge<KSIZE, DMAX_C> merger(
      part_paths, ab_mins, model, accumulators, config.kmer_size, opt->nb_controls,
      opt->nb_cases, opt->threshold/opt->cutoff, opt->nb_threads, sampler, opt->save_sk ? sign_matrix_dir : std::string(""));

    std::size_t total_kmers = 0;

    if (from_matrix)
      total_kmers = merger.merge(matrix_paths);
    else
      total_kmers = merger.merge();

    auto [sign_controls, sign_cases] = merger.signs();

    spdlog::info("Partitions processed ({})", merge_time.formatted());

    spdlog::info("{}/{} significant k-mers.", merger.nb_sign(), total_kmers);
    spdlog::info("Before correction: {} (control), {} (case).", sign_controls, sign_cases);

    return total_kmers;
  }

  #ifdef WITH_POPSTRAT
  template<std::size_t KSIZE>
    void do_pop(std::vector<acc_t<KmerSign<KSIZE>>>& accumulators,
                const std::string& pop_dir,
                const std::string& output_part_dir,
                diff_options_t opt,
                const kmtricks_config_t& config)
    {
      Timer pca_time;

      std::vector<acc_t<KmerSign<KSIZE>>> pop_accumulators;

      std::string gwas_eigenstratX_ind = fmt::format("{}/gwas_eigenstratX.ind", pop_dir);
      std::string gwas_eigenstratX_total = fmt::format("{}/gwas_eigenstratX.total", pop_dir);
      std::string pcs_evec = fmt::format("{}/pcs.evec", pop_dir);

      std::string parfile_path = fmt::format("{}/parfile.txt", pop_dir);
      std::string gwas_info_path = fmt::format("{}/gwas_infos.txt", pop_dir);
      std::string fof = fmt::format("{}/kmtricks.fof", opt->kmtricks_dir);

      auto [total_controls, total_cases] = get_total_kmer(opt->kmtricks_dir, opt->nb_controls, opt->nb_cases, config.abundance_min);

      write_parfile(parfile_path);
      write_gwas_info(fof, gwas_info_path, opt->nb_controls, opt->nb_cases, opt->gender);
      write_gwas_info(fof, gwas_eigenstratX_ind, opt->nb_controls, opt->nb_cases, opt->gender);
      write_gwas_eigenstrat_total(gwas_eigenstratX_total, total_controls, total_cases);

      std::string log_eigenstrat = "eigenstrat.log";
      run_eigenstrat_smartpca(
        pop_dir, "parfile.txt", log_eigenstrat, opt->is_diploid, opt->nb_controls + opt->nb_cases);

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
        pop_accumulators[p] = std::make_shared<FileAccumulator<KmerSign<KSIZE>>>(
          fmt::format("{}/p{}_popstrat_uncorrected", output_part_dir, p), config.kmer_size, false, !opt->keep_tmp);
      }

      pop_corrector->apply(accumulators, pop_accumulators, opt->nb_threads);

      accumulators.swap(pop_accumulators);

      spdlog::info("Population correction done. ({}).", pop_time.formatted());
    }
  #endif

  template<std::size_t KSIZE>
  void do_correction(std::vector<acc_t<KmerSign<KSIZE>>>& accumulators,
                     diff_options_t opt,
                     const kmtricks_config_t& config,
                     std::size_t total_kmers)
  {
    Timer agg_time;

    if (opt->correction == CorrectionType::NOTHING)
      spdlog::info("Aggregate partitions...");
    else
      spdlog::info("Aggregate partitions and apply significance correction...");

    indicators::ProgressBar* pb = nullptr;

    if ((spdlog::get_level() != spdlog::level::debug) && isatty_stderr())
    {
      pb = get_progress_bar("progress", config.nb_partitions, 50, indicators::Color::white, false);
      pb->set_progress(0);
      pb->print_progress();
    }

    auto corrector = make_corrector(opt->correction, opt->threshold, total_kmers);
    auto agg = make_aggregator<KSIZE>(
        accumulators, corrector, config, opt->output_directory, opt->kff, opt->nb_threads, pb);

    agg->run();

    auto [c_controls, c_cases] = agg->counts();

    delete pb;
    spdlog::info("Partitions aggregated ({})", agg_time.formatted());
    spdlog::info("Significant k-mers: {} (control), {} (case).", c_controls, c_cases);
  }

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
    km::KmDir::get().init(opt->kmtricks_dir, fmt::format("{}/kmtricks.fof", opt->kmtricks_dir));
    kmtricks_config_t config = get_kmtricks_config(opt->kmtricks_dir);
    km::Kmer<KSIZE>::m_kmer_size = config.kmer_size;

    bool prev_run = fs::exists(fmt::format("{}/options.bin", opt->output_directory));
    diff_options_t prev_opt = nullptr;

    bool prev_1 = false;
    bool prev_2 = false;
    bool prev_f = false;
    unsigned action = 0;

    std::string output_part_dir = fmt::format("{}/partitions", opt->output_directory);
    fs::create_directories(output_part_dir);

    if (prev_run)
    {
      prev_opt = load_opt(fmt::format("{}/options.bin", opt->output_directory));
      spdlog::debug(fmt::format("Previous {}", prev_opt->display()));
      action = compare_opt(opt, prev_opt);
      prev_1 = partitions_exist("{}/p{}_uncorrected", config.nb_partitions, output_part_dir);
      prev_2 = partitions_exist("{}/p{}_popstrat_uncorrected", config.nb_partitions, output_part_dir);
      prev_f = fs::exists(fmt::format("{}/control_kmers.fasta", opt->output_directory)) &&
               fs::exists(fmt::format("{}/case_kmers.fasta", opt->output_directory));

      spdlog::debug("prev1 -> {}", prev_1);
      spdlog::debug("prev2 -> {}", prev_2);
      spdlog::debug("prevf -> {}", prev_f);
      spdlog::debug("action -> {}", action);
    }

    eig_geno_t<DMAX_C> geno {nullptr};
    eig_snp_t snp {nullptr};
    std::shared_ptr<Sampler<DMAX_C>> sampler {nullptr};

    std::string pop_dir;

    if (opt->pop_correction)
    {
      pop_dir = fmt::format("{}/popstrat", opt->output_directory);
      fs::create_directory(pop_dir);
    }

    bool redo_c = false;
    std::vector<acc_t<KmerSign<KSIZE>>> accumulators(config.nb_partitions);

    if (!prev_1 || (action & 0b1))
    {
      std::string gwas_eigenstratX_geno = fmt::format("{}/gwas_eigenstratX.geno", pop_dir);
      std::string gwas_eigenstratX_snp = fmt::format("{}/gwas_eigenstratX.snp", pop_dir);

      if (opt->pop_correction)
      {
        geno = std::make_shared<EigGenoFile<DMAX_C>>(gwas_eigenstratX_geno);
        snp = std::make_shared<EigSnpFile>(gwas_eigenstratX_snp);
        sampler = std::make_shared<Sampler<DMAX_C>>(geno, snp, opt->kmer_pca, opt->seed);
      }

      opt->total_kmers = do_diff<KSIZE>(opt, config, output_part_dir, accumulators, sampler);
      redo_c = true;

      if (opt->pop_correction)
      {
        geno->close();
        snp->close();
      }
    }
    else
    {
      opt->total_kmers = prev_opt->total_kmers;
      for (std::size_t i = 0; i < accumulators.size(); ++i)
      {
        accumulators[i] = std::make_shared<FileAccumulator<KmerSign<KSIZE>>>(
          fmt::format("{}/p{}_uncorrected", output_part_dir, i), config.kmer_size, true, !opt->keep_tmp);
      }
    }

    dump_opt(opt, fmt::format("{}/options.bin", opt->output_directory));

    #ifdef WITH_POPSTRAT
      pop_strat_corrector::set_params(opt->max_iteration,
                                      opt->learning_rate,
                                      opt->epsilon,
                                      opt->stand,
                                      opt->irls);

      if (opt->pop_correction && ((!prev_2 || (action & 0b10)) || ((action & 0b1) || !prev_1)))
      {
        do_pop<KSIZE>(accumulators, pop_dir, output_part_dir, opt, config);
        redo_c = true;
      }
    #endif

    if ((!prev_f || (action > 0)) || redo_c)
    {
      do_correction<KSIZE>(accumulators, opt, config, opt->total_kmers);
    }

    spdlog::info(
      "Done in {}, Peak RSS -> {} MB.",
      whole_time.formatted(),
      static_cast<size_t>(get_peak_rss() * 0.0009765625)
    );
  }

} // end of namespace kmdiff

