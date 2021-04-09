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

// int
#include <kmdiff/cmd/cmd_common.hpp>
#include <kmdiff/config.hpp>
#include <kmdiff/kmtricks_utils.hpp>
#include <kmdiff/model.hpp>
#include <kmdiff/simulator.hpp>
#include <kmdiff/utils.hpp>

namespace kmdiff
{
struct popsim_options : kmdiff_options
{
  std::string output_directory;
  std::string reference;
  int nb_controls;
  int nb_cases;

  int mean_sv_len;
  int sd_sv_len;

  int nb_sv_controls;
  std::string type_sv_controls;
  std::string ratio_sv_controls;

  int nb_sv_cases;
  std::string type_sv_cases;
  std::string ratio_sv_cases;

  int mean_sv_per_indiv_controls;
  int sd_sv_per_indiv_controls;
  double prob_case;

  int mean_sv_per_indiv_cases;
  int sd_sv_per_indiv_cases;
  double prob_control;

  int read_size;
  int coverage;
  double error_rate;
  double mutation_rate;
  double indel_fraction;
  double extend;

  std::string display()
  {
    std::stringstream ss;
    ss << "\n";
    RECORD(ss, output_directory);
    RECORD(ss, reference);
    RECORD(ss, nb_controls);
    RECORD(ss, nb_cases);
    RECORD(ss, mean_sv_len);
    RECORD(ss, sd_sv_len);
    RECORD(ss, nb_controls);
    RECORD(ss, type_sv_controls);
    RECORD(ss, ratio_sv_controls);
    RECORD(ss, nb_sv_cases);
    RECORD(ss, type_sv_cases);
    RECORD(ss, ratio_sv_cases);
    RECORD(ss, mean_sv_per_indiv_controls);
    RECORD(ss, sd_sv_per_indiv_controls);
    RECORD(ss, prob_case);
    RECORD(ss, mean_sv_per_indiv_cases);
    RECORD(ss, sd_sv_per_indiv_cases);
    RECORD(ss, prob_control);
    RECORD(ss, read_size);
    RECORD(ss, coverage);
    RECORD(ss, error_rate);
    RECORD(ss, mutation_rate);
    RECORD(ss, indel_fraction);
    RECORD(ss, extend);
    ss << this->global_display();
    return ss.str();
  }
};

using popsim_options_t = std::shared_ptr<struct popsim_options>;

void main_popsim(kmdiff_options_t options);

};  // namespace kmdiff