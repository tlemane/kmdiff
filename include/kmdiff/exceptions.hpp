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
#include <stdexcept>
#include <string>

namespace kmdiff
{
class kmdiff_exception : public std::exception
{
 private:
  std::string name{"Base error"};
  std::string msg{"Base error msg should never be printed"};

 public:
  kmdiff_exception(const std::string& name, const std::string& msg) : name(name), msg(msg) {}
  std::string get_name() const { return name; }
  std::string get_msg() const { return msg; }
};

#define kmdiff_EXCEPTION(name)                                     \
  class name : public kmdiff_exception                             \
  {                                                              \
   public:                                                       \
    name(const std::string& msg) : kmdiff_exception(#name, msg) {} \
  }

kmdiff_EXCEPTION(BinaryNotFound);
kmdiff_EXCEPTION(ExternalExecFailed);
kmdiff_EXCEPTION(KmtricksFileNotFound);
kmdiff_EXCEPTION(FileNotFound);
kmdiff_EXCEPTION(ConfigError);
kmdiff_EXCEPTION(IOError);


kmdiff_EXCEPTION(VCFOpenError);
kmdiff_EXCEPTION(VCFHeaderError);

kmdiff_EXCEPTION(BEDOpenError);
kmdiff_EXCEPTION(BEDBadFormat);

kmdiff_EXCEPTION(BAMOpenError);
kmdiff_EXCEPTION(BAMHeaderError);

kmdiff_EXCEPTION(EigenStratError);

kmdiff_EXCEPTION(SingularError);

};  // namespace kmdiff