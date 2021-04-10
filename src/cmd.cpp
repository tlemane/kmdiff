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

#include <kmdiff/cmd.hpp>

namespace kmdiff
{
void main_count(kmdiff_options_t options)
{
  count_options_t opt = std::static_pointer_cast<struct count_options>(options);

  std::string cmd_path = command_exists(get_binary_dir(), "kmtricks.py");
  spdlog::debug(opt->display());

  std::string kmtricks_args =
      "run --file {} --run-dir {} --kmer-size {} --count-abundance-min {} "
      "--max-memory {} --nb-cores {} --minimizer-type {} --mode bin "
      "--minimizer-size {} --repartition-type {} --nb-partitions {} --lz4 --until count";

  std::string fmt_args = fmt::format(
      kmtricks_args, opt->file, opt->dir, opt->kmer_size, opt->abundance_min, opt->memory,
      opt->nb_threads, opt->minimizer_type, opt->minimizer_size, opt->repartition_type,
      opt->nb_partitions);

#ifdef KMDIFF_DEV_MODE
  if (options->signal) std::raise(options->signal);
#endif

  exec_external_cmd(cmd_path, fmt_args);

  std::ofstream out_opt(fmt::format("{}/{}-count.opt", opt->dir, PROJECT_NAME));
  out_opt << opt->display() << std::endl;
}

void main_diff(kmdiff_options_t options)
{
  diff_options_t opt = std::static_pointer_cast<struct diff_options>(options);
  spdlog::debug(opt->display());

  Timer whole;
  kmtricks_config_t config = get_kmtricks_config(opt->kmtricks_dir);

  std::vector<std::string> fofs = get_fofs(opt->kmtricks_dir);
  std::vector<uint32_t> dummy_a_min = std::vector<uint32_t>(opt->nb_controls + opt->nb_cases, 1);

#ifdef KMDIFF_DEV_MODE
  if (options->signal) std::raise(options->signal);
#endif

  spdlog::info("Process partitions...");
  Timer merge_timer;

  auto [total_controls, total_cases] =
      get_total_kmer(opt->kmtricks_dir, opt->nb_controls, opt->nb_cases);

  std::shared_ptr<Model<DEF_MAX_COUNT>> model = std::make_shared<PoissonLikelihood<DEF_MAX_COUNT>>(
      opt->nb_controls, opt->nb_cases, total_controls, total_cases);

  std::string output_part_dir = fmt::format("{}/partitions", opt->output_directory);
  fs::create_directories(output_part_dir);

  {
    std::ofstream out_opt(fmt::format("{}/{}-diff.opt", opt->output_directory, PROJECT_NAME));
    out_opt << opt->display() << std::endl;
  }

  std::vector<acc_t<KmerSign<DEF_MAX_KMER>>> accumulators(config.nb_partitions);
  for (size_t i = 0; i < accumulators.size(); i++)
  {
    if (opt->in_memory)
      accumulators[i] = std::make_shared<VectorAccumulator<KmerSign<DEF_MAX_KMER>>>(4096);
    else
      accumulators[i] = std::make_shared<FileAccumulator<KmerSign<DEF_MAX_KMER>>>(
          fmt::format("{}/acc_{}", output_part_dir, i), config.kmer_size);
  }

  GlobalMerge<DEF_MAX_KMER, DEF_MAX_COUNT> merger(
      fofs, dummy_a_min, config.kmer_size, opt->nb_controls, opt->nb_cases, opt->threshold,
      output_part_dir, opt->nb_threads, model, accumulators);

  size_t total_kmers = merger.merge();

  spdlog::info(
      "Partitions processed ({} seconds).", merge_timer.elapsed<std::chrono::seconds>().count());

  std::shared_ptr<ICorrector> corrector = nullptr;

  spdlog::info("Aggregate and apply {} correction...", correction_type_str(opt->correction));
  if (opt->correction == CorrectionType::BONFERRONI)
    corrector = std::make_shared<Bonferroni>(opt->threshold, total_kmers);

  BlockingQueue<KmerSign<DEF_MAX_KMER>> case_queue(50000, config.nb_partitions);
  BlockingQueue<KmerSign<DEF_MAX_KMER>> control_queue(50000, config.nb_partitions);

  ThreadPool pool(opt->nb_threads < 2 ? 1 : opt->nb_threads);

  for (size_t p = 0; p < config.nb_partitions; p++)
  {
    auto task = [&case_queue, &control_queue, &accumulators, p](int id) {
      while (std::optional<KmerSign<DEF_MAX_KMER>>& o = accumulators[p]->get())
      {
        if ((*o).m_sign == Significance::CONTROL)
          control_queue.push(std::move(*o));
        else if ((*o).m_sign == Significance::CASE)
          case_queue.push(std::move(*o));
      }
      control_queue.end_signal(p);
      case_queue.end_signal(p);
    };
    pool.add_task(task);
  }

  using seq_out_t = std::unique_ptr<klibpp::SeqStreamOut>;

  std::string ext = opt->kff ? ".kff" : ".fasta";
  
  std::string control_out = fmt::format("{}/control_kmers{}", opt->output_directory, ext);
  auto control_consumer = std::thread([&control_queue, &control_out, corrector, &opt, &config]() {
    KmerSign<DEF_MAX_KMER> k;
    klibpp::KSeq record;
    seq_out_t out = nullptr;
    kff_w_t out_kff = nullptr;

    if (opt->kff)
      out_kff = std::make_unique<KffWriter>(control_out, config.kmer_size);
    else
      out = std::make_unique<klibpp::SeqStreamOut>(control_out.c_str());

    size_t i = 0;
    while (control_queue.pop(k))
    {
      bool keep = true;
      if (corrector) keep = corrector->apply(k.m_pvalue);

      if (keep)
      {
        if (!opt->kff)
        {
          record.name = fmt::format("{}_pval={}_control={}_case={}",
                                  i, k.m_pvalue, k.m_mean_control, k.m_mean_case);
          record.seq = k.m_kmer.to_string();
          *out << klibpp::format::fasta << record;
        }
        else
        {
          out_kff->write(k);
        }
      }
      i++;
    }
    if (out_kff) out_kff->close();
    spdlog::info("Over-represented k-mers in controls dumped at {}", control_out);
  });

  std::string case_out = fmt::format("{}/case_kmers{}", opt->output_directory, ext);
  auto case_consumer = std::thread([&case_queue, &case_out, corrector, &opt, &config]() {
    KmerSign<DEF_MAX_KMER> k;
    klibpp::KSeq record;
    seq_out_t out = nullptr;
    kff_w_t out_kff = nullptr;

    if (opt->kff)
      out_kff = std::make_unique<KffWriter>(case_out, config.kmer_size);
    else
      out = std::make_unique<klibpp::SeqStreamOut>(case_out.c_str());
    
    size_t i = 0;
    while (case_queue.pop(k))
    {
      bool keep = true;
      if (corrector) keep = corrector->apply(k.m_pvalue);

      if (keep)
      {
        if (!opt->kff)
        {
          record.name = fmt::format("{}_pval={}_control={}_case={}",
                                  i, k.m_pvalue, k.m_mean_control, k.m_mean_case);
          record.seq = k.m_kmer.to_string();
          *out << klibpp::format::fasta << record;
        }
        else
        {
          out_kff->write(k);
        }
      }
      i++;
    }
    if (out_kff) out_kff->close();
    spdlog::info("Over-represented k-mers in cases dumped at {}", case_out);
  });

  pool.join_all();
  case_consumer.join();
  control_consumer.join();

  spdlog::info("Done ({} seconds), Peak RSS -> {} MB.",
               whole.elapsed<std::chrono::seconds>().count(),
               static_cast<size_t>(get_peak_rss() * 0.0009765625));
}

void main_popsim(kmdiff_options_t options)
{
  popsim_options_t opt = std::static_pointer_cast<struct popsim_options>(options);
  spdlog::debug(opt->display());

  Timer sim_timer;

  std::string bin_dir = get_binary_dir();
  std::string visor_bin = command_exists(bin_dir, "VISOR");
  std::string random_bin = command_exists(bin_dir, "randomregion.r");

  std::string ref_dir = fmt::format("{}/ref", opt->output_directory);
  std::string bed_dir = fmt::format("{}/bed", opt->output_directory);
  fs::create_directory(opt->output_directory);
  fs::create_directory(ref_dir);
  fs::create_directory(bed_dir);

  {
    std::ofstream out_opt(fmt::format("{}/{}-popsim.opt", opt->output_directory, PROJECT_NAME));
    out_opt << opt->display() << std::endl;
  }

  std::string dim_file = fmt::format("{}/ref/dim.tsv", opt->output_directory);
  size_t genome_size = get_dim_file(opt->reference, dim_file);

  std::string control_pool_bed = fmt::format("{}/control_pool.bed", bed_dir);
  std::string case_pool_bed = fmt::format("{}/case_pool.bed", bed_dir);

  Timer bed_timer;
  spdlog::info("Generate {}...", control_pool_bed);
  get_random_region(
      random_bin, dim_file, convert_to_visor(opt->type_sv_controls), opt->ratio_sv_controls,
      control_pool_bed, opt->nb_sv_controls, opt->mean_sv_len, opt->sd_sv_len);
  spdlog::info("Done ({} seconds).", bed_timer.elapsed<std::chrono::seconds>().count());

  bed_timer.reset();
  spdlog::info("Generate {}...", case_pool_bed);
  get_random_region(
      random_bin, dim_file, convert_to_visor(opt->type_sv_cases), opt->ratio_sv_cases,
      case_pool_bed, opt->nb_sv_cases, opt->mean_sv_len, opt->sd_sv_len);
  spdlog::info("Done ({} seconds).", bed_timer.elapsed<std::chrono::seconds>().count());

  bed_timer.reset();
  spdlog::info("Generate individuals...");
  SVPool control_pool(control_pool_bed);
  SVPool case_pool(case_pool_bed);

  Simulator simulator(
      control_pool, opt->nb_controls, opt->mean_sv_per_indiv_controls,
      opt->sd_sv_per_indiv_controls, opt->prob_case, case_pool, opt->nb_cases,
      opt->mean_sv_per_indiv_cases, opt->sd_sv_per_indiv_cases, opt->prob_control);

  std::string control_dir = fmt::format("{}/controls", opt->output_directory);
  std::string case_dir = fmt::format("{}/cases", opt->output_directory);
  fs::create_directory(control_dir);
  fs::create_directory(case_dir);

  simulator.generate(control_dir, case_dir);

  std::string control_pool_bed_real = fmt::format("{}/real_control.bed", bed_dir);
  std::string case_pool_bed_real = fmt::format("{}/real_case.bed", bed_dir);
  std::string shared_bed = fmt::format("{}/shared.bed", bed_dir);

  simulator.dump(control_pool_bed_real, case_pool_bed_real, shared_bed);
  spdlog::info(
      "Individuals generated ({} seconds).", bed_timer.elapsed<std::chrono::seconds>().count());

  bed_timer.reset();
  spdlog::info("Generate fasta...");
  simulator.generate_refs(visor_bin, opt->reference, control_dir, case_dir);
  spdlog::info("Fasta generated ({} seconds).", bed_timer.elapsed<std::chrono::seconds>().count());

  bed_timer.reset();
  spdlog::info("Generate reads...");
  std::string km_fof = simulator.generate_reads(
      control_dir, case_dir, genome_size, opt->read_size, opt->coverage, opt->error_rate,
      opt->mutation_rate, opt->indel_fraction, opt->extend, opt->nb_threads);
  spdlog::info("Reads generated ({} seconds).", bed_timer.elapsed<std::chrono::seconds>().count());

  std::string km_fof_path = fmt::format("{}/kmtricks.fof", opt->output_directory);
  std::ofstream km_fof_out(km_fof_path, std::ios::out);
  if (!km_fof_out.good()) std::runtime_error(fmt::format("Unable to write at {}.", km_fof_path));

  km_fof_out << km_fof;
  km_fof_out.flush();

  spdlog::info(
      "Population simulated ({} seconds).", sim_timer.elapsed<std::chrono::seconds>().count());

  spdlog::info("kmtricks input fof available at {}.", km_fof_path);
}

void main_call(kmdiff_options_t options)
{
  call_options_t opt = std::static_pointer_cast<struct call_options>(options);
  spdlog::debug(opt->display());

  Timer call_timer;
  spdlog::info("Map significant k-mers...");
  
  std::string output_directory = fmt::format("{}/alignments", opt->directory);
  {
    std::ofstream out_opt(fmt::format("{}/{}-call.opt", opt->directory, PROJECT_NAME));
    out_opt << opt->display() << std::endl;
  }
  std::string bbmap_bin = command_exists(get_binary_dir(), "bbmap.sh");
  std::string bbmap_index = "ref={} threads={}";
  std::string bbmap_align = "in={} out={} k={} threads={} nodisk";

  Timer index_time;
  spdlog::info("Build reference index...");
  std::string index_cmd = fmt::format(bbmap_index, opt->reference, opt->nb_threads);
  exec_external_cmd(bbmap_bin, index_cmd);
  spdlog::info("Done ({} seconds).",
                index_time.elapsed<std::chrono::seconds>().count());

  Timer control_time;
  spdlog::info("Map control k-mers...");
  std::string control_kmer = fmt::format("{}/control_kmers.fasta", opt->directory);
  if (!fs::exists(control_kmer)) 
    throw FileNotFound(fmt::format("{} not found.", control_kmer));
  std::string control_output = fmt::format("{}/control_kmers.sam", output_directory);
  std::string control_cmd = fmt::format(bbmap_align,
                                        control_kmer,
                                        control_output,
                                        opt->seed_size,
                                        opt->nb_threads);
  exec_external_cmd(bbmap_bin, control_cmd);
  spdlog::info("Control k-mers mapped ({} seconds).",
                control_time.elapsed<std::chrono::seconds>().count());

  Timer case_time;
  spdlog::info("Map case k-mers...");
  std::string case_kmer = fmt::format("{}/case_kmers.fasta", opt->directory);
  if (!fs::exists(case_kmer)) 
    throw FileNotFound(fmt::format("{} not found.", case_kmer));
  std::string case_output = fmt::format("{}/control_kmers.sam", output_directory);
  std::string case_cmd = fmt::format(bbmap_align,
                                     case_kmer,
                                     case_output,
                                     opt->seed_size,
                                     opt->nb_threads);
  exec_external_cmd(bbmap_bin, case_cmd);
  spdlog::info("Case k-mers mapped ({} seconds).",
                case_time.elapsed<std::chrono::seconds>().count());
  
  spdlog::info("Mapping done ({} seconds).",
                call_timer.elapsed<std::chrono::seconds>().count());
}
};  // namespace kmdiff