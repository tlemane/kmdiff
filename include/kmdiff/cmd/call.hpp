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
struct call_options : kmdiff_options
{
  std::string directory;
  std::string reference;
  size_t      seed_size;

  std::string display()
  {
    std::stringstream ss;
    ss << this->global_display();
    RECORD(ss, directory);
    RECORD(ss, reference);
    RECORD(ss, seed_size);
    return ss.str();
  }
};

using call_options_t = std::shared_ptr<struct call_options>;

void main_call(kmdiff_options_t options);

};  // namespace kmdiff