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

//#include <kmdiff/aggregator.hpp>
#include <kmdiff/cmd.hpp>
//#include <kmdiff/io/bam.hpp>
//#include <kmdiff/popstrat.hpp>
//#include <kmdiff/validator.hpp>

namespace kmdiff
{
void main_count(kmdiff_options_t options)
{
  count_options_t opt = std::static_pointer_cast<struct count_options>(options);

  std::string cmd_path = command_exists(get_binary_dir(), "kmtricks");
  spdlog::debug(opt->display());

  std::string kmtricks_args =
      "pipeline --file {} --run-dir {} --kmer-size {} --hard-min {} "
      "--threads {} --minimizer-type {} --mode kmer:count:bin "
      "--minimizer-size {} --repartition-type {} --nb-partitions {} --cpr --until count --hist";

  std::string fmt_args = fmt::format(
      kmtricks_args, opt->file, opt->dir, opt->kmer_size, opt->abundance_min,
      opt->nb_threads, opt->minimizer_type, opt->minimizer_size,
      opt->repartition_type, opt->nb_partitions);

  exec_external_cmd(cmd_path, fmt_args);

  std::ofstream out_opt(fmt::format("{}/{}-count.opt", opt->dir, KMD_PROJECT_NAME));
  out_opt << opt->display() << std::endl;
}

//void main_diff(kmdiff_options_t options)
//{
//  diff_options_t opt = std::static_pointer_cast<struct diff_options>(options);
//  spdlog::debug(opt->display());
//
//  Timer whole;
//  kmtricks_config_t config = get_kmtricks_config(opt->kmtricks_dir);
//
//  std::vector<std::string> fofs = get_fofs(opt->kmtricks_dir);
//  std::vector<uint32_t> dummy_a_min = std::vector<uint32_t>(opt->nb_controls + opt->nb_cases, 1);
//
//#ifdef KMDIFF_DEV_MODE
//  if (options->signal) std::raise(options->signal);
//#endif
//
//  spdlog::info("Process partitions...");
//  Timer merge_timer;
//
//#ifdef KMDIFF_DEV_MODE
//  PopStratCorrector::s_learn_rate = opt->learning_rate;
//  PopStratCorrector::s_max_iteration = opt->max_iteration;
//  PopStratCorrector::s_epsilon = opt->epsilon;
//  PopStratCorrector::s_stand = opt->stand;
//#endif
//
//  auto [total_controls, total_cases] =
//      get_total_kmer(opt->kmtricks_dir, opt->nb_controls, opt->nb_cases);
//
//  std::shared_ptr<Model<DEF_MAX_COUNT>> model = std::make_shared<PoissonLikelihood<DEF_MAX_COUNT>>(
//      opt->nb_controls, opt->nb_cases, total_controls, total_cases, 50000);
//
//  std::string output_part_dir = fmt::format("{}/partitions", opt->output_directory);
//  fs::create_directories(output_part_dir);
//
//  {
//    std::ofstream out_opt(fmt::format("{}/{}-diff.opt", opt->output_directory, PROJECT_NAME));
//    out_opt << opt->display() << std::endl;
//  }
//
//  std::vector<acc_t<KmerSign<DEF_MAX_KMER>>> accumulators(config.nb_partitions);
//  for (size_t i = 0; i < accumulators.size(); i++)
//  {
//    if (opt->in_memory)
//      accumulators[i] = std::make_shared<VectorAccumulator<KmerSign<DEF_MAX_KMER>>>(65536);
//    else
//      accumulators[i] = std::make_shared<FileAccumulator<KmerSign<DEF_MAX_KMER>>>(
//          fmt::format("{}/acc_{}", output_part_dir, i), config.kmer_size);
//  }
//
//#ifdef WITH_POPSTRAT
//  std::string popstrat_dir = fmt::format("{}/popstrat", opt->output_directory);
//  std::string gwas_eigenstratX_geno = fmt::format("{}/gwas_eigenstratX.geno", popstrat_dir);
//  std::string gwas_eigenstratX_snp = fmt::format("{}/gwas_eigenstratX.snp", popstrat_dir);
//  std::string gwas_eigenstratX_ind = fmt::format("{}/gwas_eigenstratX.ind", popstrat_dir);
//  std::string gwas_eigenstratX_total = fmt::format("{}/gwas_eigenstratX.total", popstrat_dir);
//  std::string pcs_evec = fmt::format("{}/pcs.evec", popstrat_dir);
//
//  GlobalMerge<DEF_MAX_KMER, DEF_MAX_COUNT> merger(
//      fofs, dummy_a_min, config.kmer_size, opt->nb_controls, opt->nb_cases, opt->threshold,
//      output_part_dir, opt->nb_threads, model, accumulators, gwas_eigenstratX_geno,
//      gwas_eigenstratX_snp, opt->pop_correction, opt->kmer_pca);
//#else
//  GlobalMerge<DEF_MAX_KMER, DEF_MAX_COUNT> merger(
//      fofs, dummy_a_min, config.kmer_size, opt->nb_controls, opt->nb_cases, opt->threshold,
//      output_part_dir, opt->nb_threads, model, accumulators, "", "", false, 0);
//#endif
//
//  size_t total_kmers = merger.merge();
//
//  spdlog::info(
//      "Partitions processed ({} seconds).", merge_timer.elapsed<std::chrono::seconds>().count());
//
//  spdlog::info("Found {} significant k-mers.", merger.nb_sign());
//
//#ifdef WITH_POPSTRAT
//  std::shared_ptr<PopStratCorrector> pop_corrector = nullptr;
//  if (opt->pop_correction)
//  {
//    Timer pop_time;
//    spdlog::info("Population stratification correction...");
//
//    Timer pca_time;
//    spdlog::info("Run eigenstrat/smartpca...");
//    fs::create_directory(popstrat_dir);
//    std::string parfile_path = "parfile.txt";
//    std::string gwas_info_path = fmt::format("{}/gwas_infos.txt", popstrat_dir);
//    std::string kmtricks_fof = fmt::format("{}/storage/fof.txt", opt->kmtricks_dir);
//    write_parfile(parfile_path);
//    write_gwas_info(kmtricks_fof, gwas_info_path, opt->nb_controls, opt->nb_cases);
//    write_gwas_info(kmtricks_fof, gwas_eigenstratX_ind, opt->nb_controls, opt->nb_cases);
//    write_gwas_eigenstrat_total(opt->kmtricks_dir, gwas_eigenstratX_total);
//
//    std::string log_eigenstrat = "eigenstrat.log";
//    run_eigenstrat_smartpca(popstrat_dir, parfile_path, log_eigenstrat, opt->is_diploid);
//    spdlog::info("smartpca done. ({} seconds).", pca_time.elapsed<std::chrono::seconds>().count());
//
//    spdlog::info(
//        "stratification correction done. ({} seconds).",
//        pop_time.elapsed<std::chrono::seconds>().count());
//
//    pop_corrector = std::make_shared<PopStratCorrector>(
//        opt->nb_controls, opt->nb_cases, total_controls, total_cases, opt->npc);
//    pop_corrector->load_Z(pcs_evec);
//    pop_corrector->load_Y(gwas_eigenstratX_ind);
//    pop_corrector->load_C(opt->covariates);
//    pop_corrector->load_ginfo(gwas_info_path);
//    pop_corrector->init_global_features();
//  }
//#endif
//
//  std::shared_ptr<ICorrector> corrector = nullptr;
//  std::unique_ptr<IAggregator<DEF_MAX_KMER>> aggregator;
//
//  spdlog::info("Aggregate and apply correction...", correction_type_str(opt->correction));
//
//  if (opt->correction == CorrectionType::BONFERRONI)
//  {
//    corrector = std::make_shared<Bonferroni>(opt->threshold / 100000, total_kmers);
//#ifdef WITH_POPSTRAT
//    aggregator =
//        std::make_unique<BonferonniAggregator<DEF_MAX_KMER>>(accumulators, corrector, opt, config);
//    aggregator->add_pop_corrector(pop_corrector);
//#else
//    aggregator =
//        std::make_unique<BonferonniAggregator<DEF_MAX_KMER>>(accumulators, corrector, opt, config);
//#endif
//  }
//  else if (opt->correction == CorrectionType::BENJAMINI)
//  {
//    corrector = std::make_shared<BenjaminiHochberg>(opt->threshold / 100000, total_kmers);
//#ifdef WITH_POPSTRAT
//    aggregator =
//        std::make_unique<BenjaminiAggregator<DEF_MAX_KMER>>(accumulators, corrector, opt, config);
//    aggregator->add_pop_corrector(pop_corrector);
//#else
//    aggregator =
//        std::make_unique<BenjaminiAggregator<DEF_MAX_KMER>>(accumulators, corrector, opt, config);
//#endif
//  }
//  else
//  {
//#ifdef WITH_POPSTRAT
//    aggregator = make_uncorrected_aggregator<DEF_MAX_KMER>(accumulators,
//                                                           opt,
//                                                           config,
//                                                           pop_corrector);
//#else
//    aggregator = make_uncorrected_aggregator<DEF_MAX_KMER>(accumulators,
//                                                           opt,
//                                                           config,
//                                                           nullptr);
//#endif
//  }
//
//  aggregator->run();
//
//  if ((!opt->seq_control.empty() || !opt->seq_case.empty()) && !opt->kff)
//  {
//    std::string control_out = fmt::format("{}/control_kmers{}", opt->output_directory, ".fasta");
//    std::string case_out = fmt::format("{}/case_kmers{}", opt->output_directory, ".fasta");
//    std::string out_sam_control = fmt::format("{}/control_align.sam", opt->output_directory);
//    std::string out_sam_case = fmt::format("{}/case_align.sam", opt->output_directory);
//
//    Validator control_validator(opt->seq_control, control_out, out_sam_control, config.kmer_size);
//    Validator case_validator(opt->seq_case, case_out, out_sam_case, config.kmer_size);
//
//    size_t nb_target_control, nb_covered_control;
//    size_t nb_target_case, nb_covered_case;
//
//    if (!opt->seq_control.empty())
//    {
//      control_validator.align(config.kmer_size / 4, opt->nb_threads);
//      control_validator.valid(nb_target_control, nb_covered_control);
//    }
//
//    if (!opt->seq_case.empty())
//    {
//      case_validator.align(config.kmer_size / 4, opt->nb_threads);
//      case_validator.valid(nb_target_case, nb_covered_case);
//    }
//
//    if (!opt->seq_control.empty())
//    {
//      double control_ratio =
//          static_cast<double>(nb_covered_control) / static_cast<double>(nb_target_control);
//      spdlog::info(
//          "{}% of control's SVs are covered by at least one k-mer.", control_ratio * 100.0);
//    }
//
//    if (!opt->seq_case.empty())
//    {
//      double case_ratio =
//          static_cast<double>(nb_covered_case) / static_cast<double>(nb_target_case);
//      spdlog::info("{}% of case's SVs are covered by at least one k-mer.", case_ratio * 100.0);
//    }
//  }
//
//  spdlog::info(
//      "Done ({} seconds), Peak RSS -> {} MB.", whole.elapsed<std::chrono::seconds>().count(),
//      static_cast<size_t>(get_peak_rss() * 0.0009765625));
//}
//
//void main_popsim(kmdiff_options_t options)
//{
//  popsim_options_t opt = std::static_pointer_cast<struct popsim_options>(options);
//  spdlog::debug(opt->display());
//
//  Timer sim_timer;
//
//  std::string bin_dir = get_binary_dir();
//  std::string visor_bin = command_exists(bin_dir, "VISOR");
//  std::string random_bin = command_exists(bin_dir, "randomregion.r");
//
//  std::string ref_dir = fmt::format("{}/ref", opt->output_directory);
//  std::string bed_dir = fmt::format("{}/bed", opt->output_directory);
//  fs::create_directory(opt->output_directory);
//  fs::create_directory(ref_dir);
//  fs::create_directory(bed_dir);
//
//  {
//    std::ofstream out_opt(fmt::format("{}/{}-popsim.opt", opt->output_directory, PROJECT_NAME));
//    out_opt << opt->display() << std::endl;
//  }
//
//  std::string dim_file = fmt::format("{}/ref/dim.tsv", opt->output_directory);
//  size_t genome_size = get_dim_file(opt->reference, dim_file);
//
//  std::string control_pool_bed = fmt::format("{}/control_pool.bed", bed_dir);
//  std::string case_pool_bed = fmt::format("{}/case_pool.bed", bed_dir);
//
//  Timer bed_timer;
//  spdlog::info("Generate {}...", control_pool_bed);
//  get_random_region(
//      random_bin, dim_file, convert_to_visor(opt->type_sv_controls), opt->ratio_sv_controls,
//      control_pool_bed, opt->nb_sv_controls, opt->mean_sv_len, opt->sd_sv_len);
//  spdlog::info("Done ({} seconds).", bed_timer.elapsed<std::chrono::seconds>().count());
//
//  bed_timer.reset();
//  spdlog::info("Generate {}...", case_pool_bed);
//  get_random_region(
//      random_bin, dim_file, convert_to_visor(opt->type_sv_cases), opt->ratio_sv_cases,
//      case_pool_bed, opt->nb_sv_cases, opt->mean_sv_len, opt->sd_sv_len);
//  spdlog::info("Done ({} seconds).", bed_timer.elapsed<std::chrono::seconds>().count());
//
//  bed_timer.reset();
//  spdlog::info("Generate individuals...");
//  SVPool control_pool(control_pool_bed);
//  SVPool case_pool(case_pool_bed);
//
//  Reference ref(opt->reference);
//
//  Simulator simulator(
//      ref, control_pool, opt->nb_controls, opt->mean_sv_per_indiv_controls,
//      opt->sd_sv_per_indiv_controls, opt->prob_case, case_pool, opt->nb_cases,
//      opt->mean_sv_per_indiv_cases, opt->sd_sv_per_indiv_cases, opt->prob_control);
//
//  std::string control_dir = fmt::format("{}/controls", opt->output_directory);
//  std::string case_dir = fmt::format("{}/cases", opt->output_directory);
//  fs::create_directory(control_dir);
//  fs::create_directory(case_dir);
//
//  simulator.generate(control_dir, case_dir);
//
//  std::string control_pool_bed_real = fmt::format("{}/real_control.bed", bed_dir);
//  std::string case_pool_bed_real = fmt::format("{}/real_case.bed", bed_dir);
//  std::string shared_bed = fmt::format("{}/shared.bed", bed_dir);
//
//  std::string seq_control_pool_bed_real = fmt::format("{}/seq_real_control.fasta", bed_dir);
//  std::string seq_case_pool_bed_real = fmt::format("{}/seq_real_case.fasta", bed_dir);
//  std::string seq_shared_bed = fmt::format("{}/seq_shared.fasta", bed_dir);
//
//  simulator.dump(
//      control_pool_bed_real, case_pool_bed_real, shared_bed, seq_control_pool_bed_real,
//      seq_case_pool_bed_real, seq_shared_bed, opt->kmer_size);
//
//  spdlog::info(
//      "Individuals generated ({} seconds).", bed_timer.elapsed<std::chrono::seconds>().count());
//
//  bed_timer.reset();
//  spdlog::info("Generate fasta...");
//  simulator.generate_refs(visor_bin, opt->reference, control_dir, case_dir);
//  spdlog::info("Fasta generated ({} seconds).", bed_timer.elapsed<std::chrono::seconds>().count());
//
//  bed_timer.reset();
//  spdlog::info("Generate reads...");
//  std::string km_fof = simulator.generate_reads(
//      control_dir, case_dir, genome_size, opt->read_size, opt->coverage, opt->error_rate,
//      opt->mutation_rate, opt->indel_fraction, opt->extend, opt->nb_threads);
//  spdlog::info("Reads generated ({} seconds).", bed_timer.elapsed<std::chrono::seconds>().count());
//
//  std::string km_fof_path = fmt::format("{}/kmtricks.fof", opt->output_directory);
//  std::ofstream km_fof_out(km_fof_path, std::ios::out);
//  if (!km_fof_out.good()) std::runtime_error(fmt::format("Unable to write at {}.", km_fof_path));
//
//  km_fof_out << km_fof;
//  km_fof_out.flush();
//
//  spdlog::info(
//      "Population simulated ({} seconds).", sim_timer.elapsed<std::chrono::seconds>().count());
//
//  spdlog::info("kmtricks input fof available at {}.", km_fof_path);
//}
//
//void main_call(kmdiff_options_t options)
//{
//  call_options_t opt = std::static_pointer_cast<struct call_options>(options);
//  spdlog::debug(opt->display());
//
//  Timer call_timer;
//  spdlog::info("Map significant k-mers...");
//
//  std::string output_directory = fmt::format("{}/alignments", opt->directory);
//  {
//    std::ofstream out_opt(fmt::format("{}/{}-call.opt", opt->directory, PROJECT_NAME));
//    out_opt << opt->display() << std::endl;
//  }
//  std::string bbmap_bin = command_exists(get_binary_dir(), "bbmap.sh");
//  std::string bbmap_index = "ref={} threads={}";
//  std::string bbmap_align = "in={} out={} k={} threads={} nodisk";
//
//  Timer index_time;
//  spdlog::info("Build reference index...");
//  std::string index_cmd = fmt::format(bbmap_index, opt->reference, opt->nb_threads);
//  exec_external_cmd(bbmap_bin, index_cmd);
//  spdlog::info("Done ({} seconds).", index_time.elapsed<std::chrono::seconds>().count());
//
//  Timer control_time;
//  spdlog::info("Map control k-mers...");
//  std::string control_kmer = fmt::format("{}/control_kmers.fasta", opt->directory);
//  if (!fs::exists(control_kmer)) throw FileNotFound(fmt::format("{} not found.", control_kmer));
//  std::string control_output = fmt::format("{}/control_kmers.sam", output_directory);
//  std::string control_cmd =
//      fmt::format(bbmap_align, control_kmer, control_output, opt->seed_size, opt->nb_threads);
//  exec_external_cmd(bbmap_bin, control_cmd);
//  spdlog::info(
//      "Control k-mers mapped ({} seconds).", control_time.elapsed<std::chrono::seconds>().count());
//
//  Timer case_time;
//  spdlog::info("Map case k-mers...");
//  std::string case_kmer = fmt::format("{}/case_kmers.fasta", opt->directory);
//  if (!fs::exists(case_kmer)) throw FileNotFound(fmt::format("{} not found.", case_kmer));
//  std::string case_output = fmt::format("{}/case_kmers.sam", output_directory);
//  std::string case_cmd =
//      fmt::format(bbmap_align, case_kmer, case_output, opt->seed_size, opt->nb_threads);
//  exec_external_cmd(bbmap_bin, case_cmd);
//  spdlog::info(
//      "Case k-mers mapped ({} seconds).", case_time.elapsed<std::chrono::seconds>().count());
//
//  spdlog::info("Mapping done ({} seconds).", call_timer.elapsed<std::chrono::seconds>().count());
//
//  //std::string abyss_bin = command_exists(get_binary_dir(), "ABYSS");
//  //std::string abyss_args = "-k{} -c0 -e0 {} -o {}";
//
//  //control_time.reset();
//  //spdlog::info("Assemble control k-mers...");
//
//  //std::string control_abyss_cmd = fmt::format(abyss_args);
//  //exec_external_cmd(abyss_bin, control_abyss_cmd);
//
//  //spdlog::info("Controls assembly done ({} seconds).",
//  //              control_time.elapsed<std::chrono::seconds>().count());
//
//  //case_time.reset();
//  //spdlog::info("Assemble case k-mers...");
//  //std::string case_abyss_cmd = fmt::format(abyss_args);
//  //exec_external_cmd(abyss_bin, case_abyss_cmd);
//  //spdlog::info("Cases assembly done ({} seconds).",
//  //             case_time.elapsed<std::chrono::seconds>().count());
//
//  //control_time.reset();
//  //spdlog::info("Map control contigs...");
//
//  //spdlog::info("Mapping of control contigs done ({} seconds).",
//  //              control_time.elapsed<std::chrono::seconds>().count());
//
//
//  //case_time.reset();
//  //spdlog::info("Map case contigs...");
//
//  //spdlog::info("Mapping of case contigs done ({} seconds).",
//  //              case_time.elapsed<std::chrono::seconds>().count());
//
//
//  //spdlog::info("Control k-mers mapping dump at: {}");
//  //spdlog::info("Case k-mers mapping dump at: {}");
//  //spdlog::info("Control contigs mapping dump at: {}");
//  //spdlog::info("Case contigs mapping dump at: {}");
//}

};  // namespace kmdiff
