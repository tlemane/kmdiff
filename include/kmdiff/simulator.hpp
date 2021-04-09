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
#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>

// ext
#include <fmt/format.h>
#include <robin_hood.h>
#include <spdlog/spdlog.h>
#include <wgsim.h>

#include <bcli/bcli.hpp>
#include <kseq++/kseq++.hpp>
#include <kseq++/seqio.hpp>

// int
#include <kmdiff/simulator.hpp>
#include <kmdiff/sv.hpp>
#include <kmdiff/threadpool.hpp>
#include <kmdiff/utils.hpp>

namespace kmdiff
{
const std::unordered_map<std::string, std::string> visor_map = {
    {"INS", "insertion"},
    {"DEL", "deletion"},
    {"INV", "inversion"},
    {"TD", "tandem duplication"},
    {"ITD", "inverted tandem duplication"},
    {"TCO", "translocation copy-paste"},
    {"TCU", "translocation cut-paste"},
    {"RT", "reciprocal translocation"}};

std::string convert_to_visor(const std::string& types);

size_t get_dim_file(const std::string& path, const std::string& output_path);

void get_random_region(
    const std::string& script,
    const std::string& dim_file,
    const std::string& svtype,
    const std::string& svratio,
    const std::string& output_file,
    int n,
    int l,
    int s);

void get_mutated(
    const std::string& bin_path,
    const std::string& ref,
    const std::string& bed,
    const std::string& out);

void get_reads(
    const std::string& ref,
    const std::string& out_template,
    size_t reference_size,
    size_t read_size,
    size_t coverage,
    double error_rate,
    double mutation_rate,
    double indel_fraction,
    double extend,
    size_t seed);

class SVPool
{
 private:
  std::string m_bed_path;
  std::vector<SV> m_svs;
  size_t m_size{0};

 public:
  SVPool(const std::string& bed_path);
  std::vector<SV> get(size_t size, std::mt19937 g);

 private:
  void load();
};

class Simulator
{
 private:
  std::random_device m_rd;
  std::mt19937 m_g{m_rd()};
  std::normal_distribution<double> m_norm_dist;
  std::uniform_real_distribution<double> m_dist{0, 1};
  std::uniform_int_distribution<uint32_t> m_dist_int{0, 0xFFFFFFFF};

  SVPool& m_control_pool;
  size_t m_nb_control{0};
  size_t m_mean_control{0};
  size_t m_sd_control{0};
  double m_prob_case{0};

  SVPool& m_case_pool;
  size_t m_nb_case{0};
  size_t m_mean_case{0};
  size_t m_sd_case{0};
  double m_prob_control{0};

  robin_hood::unordered_set<SV> m_real_control_pool;
  robin_hood::unordered_set<SV> m_real_case_pool;
  std::vector<SV> m_shared_pool;

 public:
  Simulator(
      SVPool& control_pool,
      size_t nb_control,
      size_t mean_control,
      size_t sd_control,
      double prob_case,
      SVPool& case_pool,
      size_t nb_case,
      size_t mean_case,
      size_t sd_case,
      double prob_control);

  void generate(const std::string& control_dir, const std::string& case_dir);

  void generate_refs(
      const std::string& visor_bin,
      const std::string& ref,
      const std::string& control_dir,
      const std::string& case_dir);

  std::string generate_reads(
      const std::string& control_dir,
      const std::string& case_dir,
      size_t genome_size,
      size_t read_size,
      size_t coverage,
      double error_rate,
      double mutation_rate,
      double indel_fraction,
      double extend,
      size_t threads);

  std::tuple<size_t, size_t> sample(size_t size, double prob_bad);

  void dump(
      const std::string& real_control, const std::string& real_case, const std::string& shared);
};

};  // end of namespace kmdiff