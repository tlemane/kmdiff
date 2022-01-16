#pragma once

// std
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

#include <iostream>

// int
#include <kmdiff/config.hpp>
#include <kmdiff/utils.hpp>
#include <limits>

namespace kmdiff
{

inline void main_infos()
{
  std::cerr << "- HOST -"
            << "\n";
  std::cerr << "build host: " << KMD_HOST_SYSTEM << "\n";
  std::cerr << "run host: " << get_uname_sr() << "\n";
  std::cerr << "- BUILD -"
            << "\n";
  std::cerr << "c compiler: " << KMD_COMPILER_C << "\n";
  std::cerr << "cxx compiler: " << KMD_COMPILER_CXX << "\n";
  std::cerr << "conda: " << KMD_CONDA_BUILD << "\n";
  std::cerr << "static: " << KMD_STATIC_BUILD << "\n";
  std::cerr << "dev: " << KMD_DEV_BUILD << "\n";
  std::cerr << "popsim: " << KMD_POPSIM_BUILD << "\n";
  std::cerr << "popstrat: " << KMD_POPSTRAT_BUILD << "\n";
  std::cerr << "kmer: " << KMD_KMER_LIST_STR << "\n";
  std::cerr << "max_c: " << DMAX_C << "\n";
  std::cerr << "\n";
  std::cerr << "- GIT SHA1 / VERSION -" << "\n";
  std::cerr << "kmdiff: " << KMD_GIT_SHA1 << "\n";
  std::cerr << "kmtricks: " << KMD_KMTRICKS_SHA1 << "\n";
  std::cerr << "bcli: " << KMD_BCLI_SHA1 << "\n";
  std::cerr << "fmt: " << KMD_FMT_SHA1 << "\n";
  std::cerr << "kff: " << KMD_KFF_SHA1 << "\n";
  std::cerr << "lz4: " << KMD_LZ4_SHA1 << "\n";
  std::cerr << "robin-hood: " << KMD_ROBIN_SHA1 << "\n";
  std::cerr << "spdlog: " << KMD_SPDLOG_SHA1 << "\n";
  std::cerr << "xxHash: " << KMD_XXHASH_SHA1 << "\n";
  std::cerr << "wgsim: " << KMD_WGSIM_SHA1 << "\n";
  std::cerr << "htslib: " << KMD_HTSLIB_SHA1 << "\n";
  std::cerr << std::flush;
}

};  // namespace kmdiff
