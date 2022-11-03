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
#include <sstream>
#include <string>

#define KRECORD(ss, var) ss << #var << "=" << var << ","

namespace kmdiff {

  enum class COMMAND
  {
    DIFF,
    COUNT,
    INFOS,
    POPSIM,
    CALL
  };

  struct kmdiff_options
  {
    std::string verbosity{};
    int nb_threads{1};

  #ifdef KMDIFF_DEV_MODE
    int signal{0};
  #endif
    std::string global_display()
    {
      std::stringstream ss;
      ss << "Options: ";
      KRECORD(ss, verbosity);
      KRECORD(ss, nb_threads);
  #ifdef KMDIFF_DEV_MODE
      KRECORD(ss, signal);
  #endif

      return ss.str();
    }
  };

  using kmdiff_options_t = std::shared_ptr<struct kmdiff_options>;

}; // end of namespace kmdiff
