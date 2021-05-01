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
#include <memory>
#include <string>

// ext
#include <bcli/bcli.hpp>

// int
#include <kmdiff/cli/cli_common.hpp>
#include <kmdiff/cli/count.hpp>
#include <kmdiff/cli/diff.hpp>
#include <kmdiff/cli/infos.hpp>
#include <kmdiff/cli/popsim.hpp>
#include <kmdiff/cli/call.hpp>
#include <kmdiff/cmd/cmd_common.hpp>
#include <kmdiff/cmd/infos.hpp>
#include <kmdiff/cmd/call.hpp>

namespace kmdiff
{
class kmdiffCli
{
public:
  kmdiffCli(
      const std::string& name,
      const std::string& desc,
      const std::string& version,
      const std::string& authors);

  std::tuple<COMMAND, kmdiff_options_t> parse(int argc, char* argv[]);

private:
  cli_t cli {nullptr};
  diff_options_t diff_opt {nullptr};
  count_options_t count_opt {nullptr};
  popsim_options_t popsim_opt {nullptr};
  call_options_t call_opt {nullptr};
};

};  // namespace kmdiff