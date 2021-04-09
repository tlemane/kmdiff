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
#include <sstream>
#include <string>

// ext
#include <spdlog/spdlog.h>

// int
#include <kmdiff/cmd/cmd_common.hpp>
#include <kmdiff/config.hpp>
#include <kmdiff/utils.hpp>

namespace kmdiff
{
struct count_options : kmdiff_options
{
  std::string file;
  std::string dir;
  int kmer_size;
  int abundance_min;
  int memory;
  int minimizer_type;
  int minimizer_size;
  int repartition_type;
  int nb_partitions;

  std::string display()
  {
    std::stringstream ss;
    ss << this->global_display();
    RECORD(ss, file);
    RECORD(ss, dir);
    RECORD(ss, kmer_size);
    RECORD(ss, abundance_min);
    RECORD(ss, memory);
    RECORD(ss, minimizer_type);
    RECORD(ss, minimizer_size);
    RECORD(ss, repartition_type);
    RECORD(ss, nb_partitions);
    return ss.str();
  }
};

using count_options_t = std::shared_ptr<struct count_options>;

void main_count(kmdiff_options_t options);

};  // namespace kmdiff