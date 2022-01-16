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

#include <memory>
#include <string>

#include <bcli/bcli.hpp>

#include <kmdiff/config.hpp>
#include <kmdiff/cli/cli_common.hpp>
#include <kmdiff/cli/count.hpp>
#include <kmdiff/cli/diff.hpp>
#include <kmdiff/cli/infos.hpp>

namespace kmdiff {

  constexpr int KL[KMER_N] = {KMER_LIST};

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
  };

}  // end of namespace kmdiff

