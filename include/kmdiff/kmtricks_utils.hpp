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
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

// ext
#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include <bcli/bcli.hpp>

// int
#include <kmdiff/exceptions.hpp>

#define KMTRICKS_PUBLIC
#include <kmtricks/io/fof.hpp>

namespace fs = std::filesystem;

namespace kmdiff
{

  struct kmtricks_config
  {
    size_t kmer_size{0};
    size_t nb_partitions{0};
    size_t abundance_min{0};
    std::string to_string() { return fmt::format("k={} n={} a={}", kmer_size, nb_partitions, abundance_min); }
  };

  using part_paths_t = std::vector<std::vector<std::string>>;
  using kmtricks_config_t = struct kmtricks_config;

  kmtricks_config_t get_kmtricks_config(const std::string& run_dir);

  km::Fof get_fofs(const std::string& run_dir);

  std::tuple<std::vector<size_t>, std::vector<size_t>> get_total_kmer(
      const std::string& run_dir, size_t nb_controls, size_t nb_cases, size_t ab_min);

  part_paths_t get_partition_paths(const std::string& kmdir, std::size_t nb_parts);

} // end of namespace kmdiff
