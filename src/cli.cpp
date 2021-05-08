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

#include <kmdiff/cli.hpp>

namespace kmdiff
{
kmdiffCli::kmdiffCli(
    const std::string& name,
    const std::string& desc,
    const std::string& version,
    const std::string& authors)
{
  cli = std::make_shared<bc::Parser<1>>(bc::Parser<1>(name, desc, version, authors));
  diff_opt = std::make_shared<struct diff_options>(diff_options{});
  count_opt = std::make_shared<struct count_options>(count_options{});
  call_opt = std::make_shared<struct call_options>(call_options{});
  info_cli(cli);
  count_cli(cli, count_opt);
  diff_cli(cli, diff_opt);
  call_cli(cli, call_opt);
#ifdef WITH_POPSIM
  popsim_opt = std::make_shared<struct popsim_options>(popsim_options{});
  popsim_cli(cli, popsim_opt);
#endif
}

std::tuple<COMMAND, kmdiff_options_t> kmdiffCli::parse(int argc, char* argv[])
{
  try
  {
    (*cli).parse(argc, argv);
  }
  catch (const bc::ex::BCliError& e)
  {
    bc::utils::exit_bcli(e);
    exit(EXIT_FAILURE);
  }

  if (cli->is("diff"))
    return std::make_tuple(COMMAND::DIFF, diff_opt);
  else if (cli->is("count"))
    return std::make_tuple(COMMAND::COUNT, count_opt);
  else if (cli->is("popsim"))
    return std::make_tuple(COMMAND::POPSIM, popsim_opt);
  else if (cli->is("call"))
    return std::make_tuple(COMMAND::CALL, call_opt);
  else
    return std::make_tuple(COMMAND::INFOS, count_opt);
}

void add_common(bc::cmd_t cmd, kmdiff_options_t options)
{
  cmd->add_group("common", "");
  cmd->add_param("-t/--threads", "Number of threads.")
      ->def(std::to_string(std::thread::hardware_concurrency()))
      ->meta("INT")
      ->setter(options->nb_threads)
      ->checker(bc::check::is_number);
  cmd->add_param("-h/--help", "Show this message and exit.")
      ->as_flag()
      ->action(bc::Action::ShowHelp);
  cmd->add_param("--version", "Show version and exit.")->as_flag()->action(bc::Action::ShowVersion);
  cmd->add_param("-v/--verbose", "Verbosity level [DEBUG|INFO|WARNING|ERROR].")
      ->meta("STR")
      ->def("INFO")
      ->checker(bc::check::f::in("DEBUG|INFO|WARNING|ERROR"))
      ->setter(options->verbosity);

#ifdef KMDIFF_DEV_MODE
  cmd->add_group("dev-common", "dev common parameters");
  cmd->add_param("-s/--signal", "Signal to raise.")
      ->checker(bc::check::is_number)
      ->setter(options->signal);
#endif
}

kmdiff_options_t count_cli(std::shared_ptr<bc::Parser<1>> cli, count_options_t options)
{
  bc::cmd_t count_cmd = cli->add_command("count", "Count k-mers with kmtricks.");

  count_cmd->add_param("-f/--file", "fof that contains path of read files")
      ->checker(bc::check::is_file)
      ->meta("FILE")
      ->setter(options->file);

  count_cmd->add_param("-d/--run-dir", "Output directory.")->meta("DIR")->setter(options->dir);

  int max = requiredK<DEF_MAX_KMER>::value / 2;
  int min = max / 2;
  count_cmd->add_param("-k/--kmer-size", fmt::format("size of k-mers [{}, {}]", min, max))
      ->checker(bc::check::is_number)
      ->checker(bc::check::f::range(min+1, max))
      ->setter(options->kmer_size)
      ->def(max == 32 ? "31" : "40")
      ->meta("INT");

  count_cmd->add_param("-c/--count-abundance-min", "min abundance for solid k-mers")
      ->checker(bc::check::is_number)
      ->setter(options->abundance_min)
      ->def("1")
      ->meta("INT");

  count_cmd->add_param("-m/--max-memory", "max memory per core (in mb)")
      ->checker(bc::check::is_number)
      ->setter(options->memory)
      ->def("4000")
      ->meta("INT");

  count_cmd->add_group("advanced performance tweaks", {});

  count_cmd->add_param("--minimizer-type", "minimizer type (0=lexi, 1=freq)")
      ->checker(bc::check::f::range(0, 1))
      ->setter(options->minimizer_type)
      ->def("0")
      ->meta("INT");

  count_cmd->add_param("--minimizer-size", "size of minimizer")
      ->checker(bc::check::is_number)
      ->setter(options->minimizer_size)
      ->def("10")
      ->meta("INT");

  count_cmd->add_param("--repartition-type", "minimizer repartition (0=unordered, 1=ordered)")
      ->checker(bc::check::f::range(0, 1))
      ->setter(options->repartition_type)
      ->def("0")
      ->meta("INT");

  count_cmd->add_param("--nb-partitions", "number of partitions (0=auto)")
      ->checker(bc::check::is_number)
      ->setter(options->nb_partitions)
      ->def("4")
      ->meta("INT");

  add_common(count_cmd, options);

  return options;
}

kmdiff_options_t diff_cli(std::shared_ptr<bc::Parser<1>> cli, diff_options_t options)
{
  bc::cmd_t diff_cmd = cli->add_command("diff", "Differential k-mers analysis.");

  auto is_kmtricks_dir = [](const std::string& p,
                            const std::string& v) -> bc::check::checker_ret_t {
    return std::make_tuple(
        fs::exists(fmt::format("{}/config.log", v)),
        fmt::format("{} {} : Not a kmtricks runtime directory.", p, v));
  };

  diff_cmd->add_param("--km-run", "kmtricks run directory.")
      ->meta("DIR")
      ->checker(bc::check::is_dir)
      ->checker(is_kmtricks_dir)
      ->setter(options->kmtricks_dir);

  diff_cmd->add_param("-o/--output-dir", "output directory.")
      ->meta("DIR")
      ->setter(options->output_directory)
      ->def("./kmdiff_output");

  diff_cmd->add_param("--nb-controls", "Number of controls.")
      ->meta("INT")
      ->checker(bc::check::is_number)
      ->setter(options->nb_controls);

  diff_cmd->add_param("--nb-cases", "Number of cases.")
      ->meta("INT")
      ->checker(bc::check::is_number)
      ->setter(options->nb_cases);

  diff_cmd->add_param("--coverage", "Coverage (running time concern, no impact on results).")
      ->meta("INT")
      ->def("20")
      ->checker(bc::check::is_number)
      ->setter(options->coverage);

  diff_cmd->add_param("--significance", "Significance threshold.")
      ->meta("FLOAT")
      ->checker(bc::check::is_number)
      ->checker(bc::check::f::range(0.0, 1.0))
      ->setter(options->threshold)
      ->def("0.05");

  auto corr_setter = [options](const std::string& v) {
    if (v == "bonferroni")
      options->correction = CorrectionType::BONFERRONI;
    else if (v == "benjamini")
      options->correction = CorrectionType::BENJAMINI;
    else
      options->correction = CorrectionType::NOTHING;
  };

  bc::param_t cp = diff_cmd->add_param("-c/--correction", "Significance correction.")
      ->meta("STR")
      ->def("bonferroni")
      ->checker(bc::check::f::in("bonferroni|benjamini|disable"))
      ->setter_c(corr_setter);

  auto correction_warn = [cp](){
    if (cp->value() == "benjamini")
      spdlog::warn("-c/--correction benjamini: all significants k-mers will live in memory.");
  };

  cp->callback(correction_warn);

  diff_cmd->add_param("--kff-output", "Output significant k-mers in kff format.")
      ->as_flag()
      ->setter(options->kff);

  auto memory_warn = [](){
    spdlog::warn("--in-memory: all significants k-mers will live in memory.");
  };
  diff_cmd->add_param("--in-memory", "Perform correction in memory.")
      ->as_flag()
      ->setter(options->in_memory)
      ->callback(memory_warn);

#ifdef WITH_POPSIM
  diff_cmd->add_param("--seq-control", "Fasta with expected control sv sequences.")
      ->meta("FILE")
      ->def("")
      ->setter(options->seq_control);

  diff_cmd->add_param("--seq-case", "Fasta with expected case sv sequences.")
      ->meta("FILE")
      ->def("")
      ->setter(options->seq_case);
#endif

#ifdef WITH_POPSTRAT
  diff_cmd->add_group("population stratification", "");

  diff_cmd->add_param("--pop-correction", "Apply correction of population stratification.")
      ->as_flag()
      ->setter(options->pop_correction);

  diff_cmd->add_param("--kmer-pca", "Proportion of k-mers used for PCA (in [0.0, 0.05]).")
      ->meta("FLOAT")
      ->def("0.001")
      ->checker(bc::check::f::range(0.0, 0.05))
      ->setter(options->kmer_pca);

  auto ploidy_setter = [options](const std::string& v) {
    options->ploidy = bc::utils::lexical_cast<size_t>(v);
    if (options->ploidy == 2) options->is_diploid = true;
    else options->is_diploid = false;
  };

  diff_cmd->add_param("--ploidy", "Ploidy level.")
      ->meta("INT")
      ->def("2")
      ->checker(bc::check::is_number)
      ->setter_c(ploidy_setter);

  diff_cmd->add_param("--n-pc", "Number of principal components (in [2, 10]).")
      ->meta("INT")
      ->def("2")
      ->checker(bc::check::f::range(2, 10))
      ->setter(options->npc);

  diff_cmd->add_param("--covariates", "Covariate file.")
      ->meta("FILE")
      ->def("")
      ->checker(bc::check::is_file)
      ->setter(options->covariates);
#endif

#ifdef KMDIFF_DEV_MODE
  diff_cmd->add_group("dev", "dev parameters");

  diff_cmd->add_param("--learning-rate", "Learning rate.")
      ->meta("FLOAT")
      ->def("0.1")
      ->checker(bc::check::f::range(0.0, 1.0))
      ->setter(options->learning_rate);

  diff_cmd->add_param("--max-iteration", "Max iteration.")
      ->meta("INT")
      ->def("25")
      ->checker(bc::check::is_number)
      ->setter(options->max_iteration);

  diff_cmd->add_param("--epsilon", "Epsilon.")
      ->meta("FLOAT")
      ->def("1e-15")
      ->setter(options->epsilon);

  diff_cmd->add_param("--stand", "Standardization.")
      ->as_flag()
      ->setter(options->stand);
#endif

  add_common(diff_cmd, options);

  return options;
}

void info_cli(std::shared_ptr<bc::Parser<1>> cli)
{
  bc::cmd_t info_cmd = cli->add_command("infos", "Show build infos.");
}

kmdiff_options_t popsim_cli(std::shared_ptr<bc::Parser<1>> cli, popsim_options_t options)
{
  bc::cmd_t popsim_cmd = cli->add_command("popsim", "Simulate population.");

  popsim_cmd->add_param("-r/--reference", "Reference genome.")
      ->meta("FILE")
      ->checker(bc::check::seems_fastx)
      ->checker(bc::check::is_file)
      ->setter(options->reference);

  popsim_cmd->add_param("-o/--output-dir", "Output directory.")
      ->meta("DIR")
      ->setter(options->output_directory);

  popsim_cmd->add_param("--kmer-size", "Size of k-mers.")
      ->meta("INT")
      ->checker(bc::check::is_number)
      ->setter(options->kmer_size);

  popsim_cmd->add_group("SVs", "");

  popsim_cmd->add_param("--mean-sv-len", "SVs mean length.")
      ->meta("INT")
      ->def("200")
      ->checker(bc::check::is_number)
      ->setter(options->mean_sv_len);

  popsim_cmd->add_param("--sd-sv-len", "SVs sd length.")
      ->meta("INT")
      ->def("10")
      ->checker(bc::check::is_number)
      ->setter(options->sd_sv_len);

  popsim_cmd->add_group("controls", "control parameters");

  popsim_cmd->add_param("--nb-controls", "Number of controls")
      ->meta("INT")
      ->def("10")
      ->checker(bc::check::f::range(1, 100))
      ->setter(options->nb_controls);

  popsim_cmd->add_param("--controls-pool-size", "Number of SVs in control population.")
      ->meta("INT")
      ->def("500")
      ->checker(bc::check::is_number)
      ->setter(options->nb_sv_controls);

  popsim_cmd->add_param("--controls-types", "")
      ->meta("STR")
      ->def("INS:DEL:INV")
      ->setter(options->type_sv_controls);

  popsim_cmd->add_param("--controls-ratio", "")
      ->meta("STR")
      ->def("40:40:20")
      ->setter(options->ratio_sv_controls);

  popsim_cmd->add_param("--controls-sv-per-indiv", "Number of SVs per individual in controls.")
      ->meta("INT")
      ->def("300")
      ->checker(bc::check::is_number)
      ->setter(options->mean_sv_per_indiv_controls);

  popsim_cmd->add_param("--controls-sd-per-indiv", "Standard deviation.")
      ->meta("INT")
      ->def("10")
      ->checker(bc::check::is_number)
      ->setter(options->sd_sv_per_indiv_controls);

  popsim_cmd->add_param("--prob-case", "")
      ->meta("FLOAT")
      ->def("0.01")
      ->checker(bc::check::f::range(0.0, 0.5))
      ->setter(options->prob_case);

  popsim_cmd->add_group("cases", "case parameters");

  popsim_cmd->add_param("--nb-cases", "")
      ->meta("INT")
      ->def("10")
      ->checker(bc::check::f::range(1, 100))
      ->setter(options->nb_cases);

  popsim_cmd->add_param("--cases-pool-size", "Number of SVs in case population.")
      ->meta("INT")
      ->def("500")
      ->checker(bc::check::is_number)
      ->setter(options->nb_sv_cases);

  popsim_cmd->add_param("--cases-types", "")
      ->meta("STR")
      ->def("INS:DEL:INV")
      ->setter(options->type_sv_cases);

  popsim_cmd->add_param("--cases-ratio", "")
      ->meta("STR")
      ->def("40:40:20")
      ->setter(options->ratio_sv_cases);

  popsim_cmd->add_param("--cases-sv-per-indiv", "Number of SVs per individual in cases.")
      ->meta("INT")
      ->def("300")
      ->checker(bc::check::is_number)
      ->setter(options->mean_sv_per_indiv_cases);

  popsim_cmd->add_param("--cases-sd-per-indiv", "Standard deviation.")
      ->meta("INT")
      ->def("10")
      ->checker(bc::check::is_number)
      ->setter(options->sd_sv_per_indiv_cases);

  popsim_cmd->add_param("--prob-control", "")
      ->meta("FLOAT")
      ->def("0.01")
      ->checker(bc::check::f::range(0.0, 0.5))
      ->setter(options->prob_control);

  popsim_cmd->add_group("reads", "wgsim params");

  popsim_cmd->add_param("--read-size", "Size of reads.")
      ->meta("INT")
      ->def("100")
      ->checker(bc::check::is_number)
      ->setter(options->read_size);

  popsim_cmd->add_param("-c/--coverage", "Sequencing coverage.")
      ->meta("INT")
      ->def("10")
      ->checker(bc::check::is_number)
      ->setter(options->coverage);

  popsim_cmd->add_param("-e/--error-rate", "Sequencing error rate.")
      ->meta("FLOAT")
      ->def("0.01")
      ->checker(bc::check::f::range(0.0, 1.0))
      ->setter(options->error_rate);

  popsim_cmd->add_param("-m/--mutation-rate", "Mutation rate.")
      ->meta("FLOAT")
      ->def("0.0")
      ->checker(bc::check::f::range(0.0, 1.0))
      ->setter(options->mutation_rate);

  popsim_cmd->add_param("-i/--indel-fraction", "Fraction of indels.")
      ->meta("FLOAT")
      ->def("0.0")
      ->checker(bc::check::f::range(0.0, 1.0))
      ->setter(options->indel_fraction);

  popsim_cmd->add_param("--extend", "Extend probability.")
      ->meta("FLOAT")
      ->def("0.0")
      ->checker(bc::check::f::range(0.0, 1.0))
      ->setter(options->extend);

  add_common(popsim_cmd, options);

  return options;
}

kmdiff_options_t call_cli(std::shared_ptr<bc::Parser<1>> cli, call_options_t options)
{
  bc::cmd_t call_cmd = cli->add_command("call", "SVs calling from k-mers.");

  call_cmd->add_param("-r/--reference", "Reference genome.")
          ->meta("FILE")
          ->checker(bc::check::seems_fastx)->checker(bc::check::is_file)
          ->setter(options->reference);

  call_cmd->add_param("-d/--diff-dir", fmt::format("Output directory of {} diff.", PROJECT_NAME))
           ->meta("DIR")
           ->checker(bc::check::is_dir)
           ->setter(options->directory);

  int max = requiredK<DEF_MAX_KMER>::value / 2;
  int min = max / 2;
  call_cmd->add_param("-k/--kmer-size", fmt::format("size of k-mers [{}, {}]", min, max))
          ->checker(bc::check::is_number)
          ->checker(bc::check::f::range(min+1, max))
          ->setter(options->kmer_size)
          ->meta("INT");

  add_common(call_cmd, options);

  return options;
}
};  // namespace kmdiff