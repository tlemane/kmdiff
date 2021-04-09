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

// ext
#include <spdlog/spdlog.h>

// int
#include <kmdiff/cli.hpp>
#include <kmdiff/cmd.hpp>
#include <kmdiff/config.hpp>
#include <kmdiff/signals.hpp>
#include <kmdiff/utils.hpp>

int main(int argc, char* argv[])
{
  using namespace kmdiff;

  SignalHandler::get().init();
  kmdiffCli cli(PROJECT_NAME, PROJECT_DESC, PROJECT_VER, AUTHOR);

  auto [cmd, options] = cli.parse(argc, argv);

  set_verbosity_level(options->verbosity);

  try
  {
    if (cmd == COMMAND::COUNT)
    {
      main_count(options);
    }
    else if (cmd == COMMAND::DIFF)
    {
      main_diff(options);
    }
    else if (cmd == COMMAND::POPSIM)
    {
      main_popsim(options);
    }
    else if (cmd == COMMAND::INFOS)
    {
      main_infos();
    }
    else if (cmd == COMMAND::CALL)
    {
      main_call(options);
    }
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
