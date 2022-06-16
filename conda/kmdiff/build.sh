#!/bin/bash

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#   kmdiff
#   Authors: T. Lemane
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Affero General Public License as
#  published by the Free Software Foundation, either version 3 of the
#  License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Affero General Public License for more details.
#
#  You should have received a copy of the GNU Affero General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

mkdir build-conda
cd build-conda

cmake .. -DCMAKE_PREFIX_PATH=${CONDA_PREFIX} -DWITH_KMTRICKS=OFF -DWITH_POPSTRAT=ON -DCONDA_BUILD=ON

make -j4

mkdir -p $PREFIX/bin

cp -r ./build-conda/bin/kmdiff $PREFIX/bin
cp -r ./build-conda/bin/smartpca $PREFIX/bin
cp -r ./build-conda/bin/evec2pca.perl $PREFIX/bin

