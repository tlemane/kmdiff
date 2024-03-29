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

#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>

#include <kmdiff/cli.hpp>
#include <kmdiff/cmd.hpp>
#include <kmdiff/config.hpp>
#include <kmdiff/signals.hpp>
#include <kmdiff/utils.hpp>
#include <kmdiff/kmtricks_utils.hpp>

#ifndef KMER_LIST
  #define KMER_LIST KMD_KMER_LIST
#endif

#ifndef KMER_N
  #define KMER_N KMD_KMERN
#endif

#include <kmtricks/loop_executor.hpp>

#define KM_EXEC_HELPER(helper_name, func_name)            \
  template<std::size_t S>                                 \
  struct helper_name                                      \
  {                                                       \
    template<typename... Args>                            \
    auto operator()(Args&&... args)                       \
    {                                                     \
      return func_name<S>(std::forward<Args>(args)...);   \
    }                                                     \
  }

KM_EXEC_HELPER(main_diff_exec, kmdiff::main_diff);

int main(int argc, char* argv[])
{
  using namespace kmdiff;

  SignalHandler::get().init();
  kmdiffCli cli(KMD_PROJECT_NAME, KMD_PROJECT_DESC, KMD_PROJECT_VER, "");

  auto [cmd, options] = cli.parse(argc, argv);

  set_verbosity_level(options->verbosity);
  auto cerr_logger = spdlog::stderr_color_mt("kmdiff");
  cerr_logger->set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%^%l%$] %v");
  spdlog::set_default_logger(cerr_logger);

  try
  {
    if (cmd == COMMAND::COUNT)
    {
      main_count(options);
    }
    else if (cmd == COMMAND::DIFF)
    {
      auto c = get_kmtricks_config(std::static_pointer_cast<struct diff_options>(options)->kmtricks_dir);
      km::const_loop_executor<0, KMER_N>::exec<main_diff_exec>(c.kmer_size, options);
      //main_diff(options);
    }
    else if (cmd == COMMAND::INFOS)
    {
      main_infos();
    }
    else if (cmd == COMMAND::CALL)
    {
      //main_call(options);
    }
#ifdef WITH_POPSIM
    else if (cmd == COMMAND::POPSIM)
    {
      //main_popsim(options);
    }
#endif
  }
  catch (const kmdiff_exception& e)
  {
    spdlog::error(fmt::format("{} - {}", e.get_name(), e.get_msg()));
  }
  catch (const std::exception& e)
  {
    spdlog::error(e.what());
  }

  return 0;
}
