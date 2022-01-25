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

namespace kmdiff {

  kmdiffCli::kmdiffCli(
      const std::string& name,
      const std::string& desc,
      const std::string& version,
      const std::string& authors)
  {
    cli = std::make_shared<bc::Parser<1>>(bc::Parser<1>(name, desc, version, authors));
    diff_opt = std::make_shared<struct diff_options>(diff_options{});
    count_opt = std::make_shared<struct count_options>(count_options{});
    info_cli(cli);
    count_cli(cli, count_opt);
    diff_cli(cli, diff_opt);
  }

  std::tuple<COMMAND, kmdiff_options_t> kmdiffCli::parse(int argc, char* argv[])
  {
    if (argc > 1 && (std::string(argv[1]) == "--help" || std::string(argv[1]) == "-h"))
    {
      cli->show_help();
      exit(EXIT_FAILURE);
    }

    if (argc > 1 && (std::string(argv[1]) == "--version" || std::string(argv[1]) == "-v"))
    {
      std::cerr << fmt::format("kmdiff {}", KMD_PROJECT_VER) << std::endl;
      exit(EXIT_FAILURE);
    }

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
    if (cli->is("count"))
      return std::make_tuple(COMMAND::COUNT, count_opt);
    else
      return std::make_tuple(COMMAND::INFOS, count_opt);
  }

  void add_common(bc::cmd_t cmd, kmdiff_options_t options)
  {
    cmd->add_group("common", "");
    cmd->add_param("-t/--threads", "number of threads.")
        ->def(std::to_string(std::thread::hardware_concurrency()))
        ->meta("INT")
        ->setter(options->nb_threads)
        ->checker(bc::check::is_number);
    cmd->add_param("-h/--help", "show this message and exit.")
        ->as_flag()
        ->action(bc::Action::ShowHelp);
    cmd->add_param("--version", "show version and exit.")->as_flag()->action(bc::Action::ShowVersion);
    cmd->add_param("-v/--verbose", "Verbosity level [debug|info|warning|error].")
        ->meta("STR")
        ->def("info")
        ->checker(bc::check::f::in("debug|info|warning|error"))
        ->setter(options->verbosity);
  }

  kmdiff_options_t count_cli(std::shared_ptr<bc::Parser<1>> cli, count_options_t options)
  {
    bc::cmd_t count_cmd = cli->add_command("count", "Count k-mers with kmtricks.");

    count_cmd->add_param("-f/--file", "fof that contains path of read files")
        ->checker(bc::check::is_file)
        ->meta("FILE")
        ->setter(options->file);

    count_cmd->add_param("-d/--run-dir", "output directory.")->meta("DIR")->setter(options->dir);

    count_cmd->add_param("-k/--kmer-size", fmt::format("size of k-mers [{}, {}]", 8, KL[KMD_KMERN-1]-1))
        ->def("31")
        ->checker(bc::check::f::range(8, KL[KMD_KMERN-1]-1))
        ->setter(options->kmer_size)
        ->meta("INT");

    count_cmd->add_param("-c/--hard-min", "min abundance for solid k-mers")
        ->checker(bc::check::is_number)
        ->setter(options->abundance_min)
        ->def("1")
        ->meta("INT");

    count_cmd->add_param("-r/--recurrence-min", "min recurrence to keep a k-mer")
        ->def("1")
        ->meta("INT")
        ->checker(bc::check::is_number)
        ->setter(options->recurrence_min);

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
        ->def("0")
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
          fs::exists(fmt::format("{}/kmtricks.fof", v)),
          fmt::format("{} {} : Not a kmtricks runtime directory.", p, v));
    };

    diff_cmd->add_param("-d/--km-run", "kmtricks run directory.")
        ->meta("DIR")
        ->checker(bc::check::is_dir)
        ->checker(is_kmtricks_dir)
        ->setter(options->kmtricks_dir);

    diff_cmd->add_param("-o/--output-dir", "output directory.")
        ->meta("DIR")
        ->setter(options->output_directory)
        ->def("./kmdiff_output");

    diff_cmd->add_param("-1/--nb-controls", "number of controls.")
        ->meta("INT")
        ->checker(bc::check::is_number)
        ->setter(options->nb_controls);

    diff_cmd->add_param("-2/--nb-cases", "number of cases.")
        ->meta("INT")
        ->checker(bc::check::is_number)
        ->setter(options->nb_cases);

    auto sign_check = [](const std::string& p, const std::string& v) -> bc::check::checker_ret_t {
      std::istringstream vs(v); double vd;
      try { vs >> vd; }
      catch (...) { return std::make_tuple(false, fmt::format("{} {}: Not a number!", p, v)); }

      bool in = true;
      if (vd < 0.0 || vd > 0.05)
        in = false;
      return std::make_tuple(in, fmt::format("Not in range [0.0, 0.05]"));
    };

    diff_cmd->add_param("-s/--significance", "significance threshold.")
        ->meta("FLOAT")
        ->checker(sign_check)
        ->setter(options->threshold)
        ->def("0.05");

    diff_cmd->add_param("-u/--cutoff", "cutoff")
      ->meta("INT")
      ->checker(bc::check::is_number)
      ->setter(options->cutoff)
      ->def("100000");

    auto corr_setter = [options](const std::string& v) {
      if (v == "bonferroni")
        options->correction = CorrectionType::BONFERRONI;
      else if (v == "benjamini")
        options->correction = CorrectionType::BENJAMINI;
      else if (v == "sidak")
        options->correction = CorrectionType::SIDAK;
      else if (v == "holm")
        options->correction = CorrectionType::HOLM;
      else
        options->correction = CorrectionType::NOTHING;
    };

    bc::param_t cp = diff_cmd->add_param("-c/--correction",
                                         "significance correction. (bonferroni|benjamini|sidak|holm|disabled)")
        ->meta("STR")
        ->def("bonferroni")
        ->checker(bc::check::f::in("bonferroni|benjamini|sidak|holm|disabled"))
        ->setter_c(corr_setter);

    auto correction_warn = [cp](){
      if (cp->value() == "benjamini")
        spdlog::warn("-c/--correction benjamini: all significants k-mers will live in memory.");
    };

    cp->callback(correction_warn);

    diff_cmd->add_param("-f/--kff-output", "output significant k-mers in kff format.")
        ->as_flag()
        ->setter(options->kff);

    auto memory_warn = [](){
      spdlog::warn("-m/--in-memory: all significants k-mers will live in memory.");
    };
    diff_cmd->add_param("-m/--in-memory", "in-memory correction.")
        ->as_flag()
        ->setter(options->in_memory)
        ->callback(memory_warn);

    diff_cmd->add_param("-r/--cpr", "compress intermediate files.")
        ->as_flag()
        ->setter(options->cpr)
        ->hide();

    #ifdef WITH_PLUGIN
      diff_cmd->add_group("custom model", "");

      diff_cmd->add_param("--cmodel", "path to model shared library.")
        ->meta("STR")
        ->def("")
        ->checker(bc::check::f::ext("so|dylib"))
        ->setter(options->model_lib_path);

      diff_cmd->add_param("--config", "model config")
        ->meta("STR")
        ->def("")
        ->setter(options->model_config);
    #endif

    #ifdef WITH_POPSTRAT
      diff_cmd->add_group("population stratification", "");

      diff_cmd->add_param("--pop-correction", "apply correction for population stratification.")
          ->as_flag()
          ->setter(options->pop_correction);

      diff_cmd->add_param("--gender", "gender file")
          ->meta("FILE")
          ->def("")
          ->checker(bc::check::is_file)
          ->setter(options->gender);

      diff_cmd->add_param("--kmer-pca", "proportion of k-mers used for PCA (in [0.0, 0.05]).")
          ->meta("FLOAT")
          ->def("0.001")
          ->checker(bc::check::f::range(0.0, 0.05))
          ->setter(options->kmer_pca);

      auto ploidy_setter = [options](const std::string& v) {
        options->ploidy = bc::utils::lexical_cast<size_t>(v);
        if (options->ploidy == 2) options->is_diploid = true;
        else options->is_diploid = false;
      };

      diff_cmd->add_param("--ploidy", "ploidy level.")
          ->meta("INT")
          ->def("2")
          ->checker(bc::check::is_number)
          ->setter_c(ploidy_setter);

      diff_cmd->add_param("--n-pc", "number of principal components (in [2, 10]).")
          ->meta("INT")
          ->def("2")
          ->checker(bc::check::f::range(2, 10))
          ->setter(options->npc);

      diff_cmd->add_param("--covariates", "covariates file.")
          ->meta("FILE")
          ->def("")
          ->checker(bc::check::is_file)
          ->setter(options->covariates);
    #endif

    #if KMD_DEV_MODE
      diff_cmd->add_group("dev", "");
    #else
      diff_cmd->add_group("dev", "")->hide();
    #endif

    diff_cmd->add_param("--learning-rate", "learning rate.")
        ->meta("FLOAT")
        ->def("0")
        ->checker(bc::check::f::range(0.0, 1.0))
        ->setter(options->learning_rate);

    diff_cmd->add_param("--max-iteration", "max iteration.")
        ->meta("INT")
        ->def("0")
        ->checker(bc::check::is_number)
        ->setter(options->max_iteration);

    diff_cmd->add_param("--epsilon", "epsilon.")
        ->meta("INT")
        ->def("0")
        ->setter(options->epsilon);

    diff_cmd->add_param("--stand", "standardization.")
        ->as_flag()
        ->setter(options->stand);

    diff_cmd->add_param("--irls", "use irls algorithm.")
        ->as_flag()
        ->setter(options->irls);

    diff_cmd->add_param("--random-seed", "random seed for pca sampling (deterministic only with 1 thread).")
        ->meta("INT")
        ->def("0")
        ->setter(options->seed);

    diff_cmd->add_param("--log-factorial", "size of precomputed table.")
        ->meta("INT")
        ->def("10000")
        ->setter(options->log_size);

    add_common(diff_cmd, options);

    return options;
  }

  void info_cli(std::shared_ptr<bc::Parser<1>> cli)
  {
    bc::cmd_t info_cmd = cli->add_command("infos", "Show build infos.");
  }

} // end of namespace kmdiff
