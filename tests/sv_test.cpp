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

#include <gtest/gtest.h>
#include <kmdiff/sv.hpp>

using namespace kmdiff;

TEST(sv, operator)
{
  SV sv("chr1", "insertion", "ACGT", 100, 101);
  SV sv1("chr1", "insertion", "ACGT", 100, 101);
  SV sv2("chr2", "insertion", "ACGT", 42, 43);

  EXPECT_EQ(sv, sv1);
  EXPECT_LT(sv2, sv1);
}

TEST(sv, bed)
{
  SV sv("chr1", "insertion", "ACGT", 100, 101);

  EXPECT_EQ(sv.to_bed_entry(), "chr1\t100\t101\tinsertion\tACGT\t0");
}

TEST(sv, hash)
{
  SV sv("chr1", "insertion", "ACGT", 100, 101);
  std::string s = "chr1insertionACGT100101";

  std::hash<SV> hasher;
  EXPECT_EQ(hasher(sv), XXH64(s.c_str(), s.size(), 0));
}