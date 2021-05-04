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

#define _KM_LIB_INCLUDE_
#include <kmtricks/utilities.hpp>

namespace kmdiff
{
inline void main_infos()
{
  std::cerr << "- HOST -"
            << "\n";
  std::cerr << "build host: " << HOST_SYSTEM << "\n";
  std::cerr << "run host: " << get_uname_sr() << "\n";
  std::cerr << "- BUILD -"
            << "\n";
  std::cerr << "c compiler: " << COMPILER_C << "\n";
  std::cerr << "cxx compiler: " << COMPILER_CXX << "\n";
  std::cerr << "conda: " << CONDA_BUILD << "\n";
  std::cerr << "static: " << STATIC_BUILD << "\n";
  std::cerr << "dev: " << DEV_BUILD << "\n";
  std::cerr << "popsim: " << POPSIM_BUILD << "\n";
  std::cerr << "popstrat: " << POPSTRAT_BUILD << "\n";
  std::cerr << "max_k: " << DEF_MAX_KMER << "\n";
  std::cerr << "max_c: " << std::numeric_limits<selectC<DEF_MAX_COUNT>::type>::max() << "\n";
  std::cerr << "\n";
  std::cerr << "- GIT SHA1 / VERSION -" << "\n";
  std::cerr << "kmdiff: " << GIT_SHA1 << "\n";
  std::cerr << "kmtricks: " << KMTRICKS_SHA1 << "\n";
  std::cerr << "bcli: " << BCLI_SHA1 << "\n";
  std::cerr << "fmt: " << FMT_SHA1 << "\n";
  std::cerr << "kff: " << KFF_SHA1 << "\n";
  std::cerr << "lz4: " << LZ4_SHA1 << "\n";
  std::cerr << "robin-hood: " << ROBIN_SHA1 << "\n";
  std::cerr << "spdlog: " << SPDLOG_SHA1 << "\n";
  std::cerr << "xxHash: " << XXHASH_SHA1 << "\n";
  std::cerr << "wgsim: " << WGSIM_SHA1 << "\n";
  std::cerr << "zlib: " << ZLIB_SHA1 << "\n";
  std::cerr << "htslib: " << HTSLIB_SHA1 << "\n";
  std::cerr << std::flush;
}

};  // namespace kmdiff