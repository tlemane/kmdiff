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
#include <csignal>
#include <memory>
#include <thread>

// ext
#include <spdlog/spdlog.h>
#define _KM_LIB_INCLUDE_
#include <kmtricks/utilities.hpp>

// int
#include <kmdiff/accumulator.hpp>
#include <kmdiff/blocking_queue.hpp>
#include <kmdiff/cmd/cmd_common.hpp>
#include <kmdiff/config.hpp>
#include <kmdiff/kmtricks_utils.hpp>
#include <kmdiff/merge.hpp>
#include <kmdiff/model.hpp>
#include <kmdiff/threadpool.hpp>
#include <kmdiff/utils.hpp>
#include <kmdiff/kff_utils.hpp>

namespace kmdiff
{
struct diff_options : kmdiff_options
{
  std::string kmtricks_dir;
  std::string output_directory;
  size_t nb_controls;
  size_t nb_cases;
  double threshold;
  CorrectionType correction;
  bool in_memory;
  bool kff;
  std::string seq_control;
  std::string seq_case;

  std::string display()
  {
    std::stringstream ss;
    ss << this->global_display();
    RECORD(ss, kmtricks_dir);
    RECORD(ss, output_directory);
    RECORD(ss, nb_controls);
    RECORD(ss, nb_cases);
    RECORD(ss, threshold);
    RECORD(ss, correction_type_str(correction));
    RECORD(ss, in_memory);
    RECORD(ss, kff);
    RECORD(ss, seq_control);
    RECORD(ss, seq_case);
    return ss.str();
  }
};

using diff_options_t = std::shared_ptr<struct diff_options>;

void main_diff(kmdiff_options_t options);

};  // namespace kmdiff