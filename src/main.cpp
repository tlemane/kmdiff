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
