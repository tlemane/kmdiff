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

#include <fstream>
#include <spdlog/spdlog.h>
#include <fmt/format.h>

#include <kmdiff/cmd.hpp>

namespace kmdiff {

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

} // end of namespace kmdiff

