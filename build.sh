#!/bin/bash

# /* SubModST solver, that solves the Cardinality-Constrained Submodular Monotone
#    Subset Maximization problem.
#    Copyright (C) 2024  Henning Woydt
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or any
#    later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.
# ==============================================================================*/

ROOT=${PWD}

function get_num_cores {
  grep -c ^processor /proc/cpuinfo;
}

# update all submodules
git submodule update --init --recursive

# setup 3rd party
cd extern
./setup_3rd_party.sh
cd ${ROOT}

# make directory
mkdir build
cd build

# build
cmake .. -DCMAKE_BUILD_TYPE=Release && cd ${ROOT}
cmake --build build --parallel "$(get_num_cores)" --target submodst
# cmake --build build --parallel "$(get_num_cores)" --target submodst_gtest
# cd ${ROOT}

# create test data - currently disabled
# cd utility && ./generate_all.sh
# cd ${ROOT}

# run tests
# cd build
# ./submodst_gtest
