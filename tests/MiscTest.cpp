/* SubModST solver, that solves the Cardinality-Constrained Submodular Monotone
   Subset Maximization problem.
   Copyright (C) 2024  Henning Woydt

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or any
   later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
==============================================================================*/

#include <vector>
#include <string>

#include <gtest/gtest.h>

#include "../src/algorithms/BFSolver.h"
#include "../src/algorithms/Solver.h"
#include "../src/structures/Graph.h"
#include "../src/structures/Matrix.h"
#include "../src/util/utility.h"

namespace SubModST {

    TEST(DataTest, GroupClosenessCentralityTest) {
        return;
        std::string data_dir_path = "../data/Graph/";
        std::string extension = ".edges";

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            SubModST::Graph graph(file);
            SubModST::Matrix<int> mtx = graph.get_GroupClosenessCentralityMatrix();

            for (size_t i = 0; i < mtx.get_n(); ++i) {
                for (size_t j = 0; j < mtx.get_n(); ++j) {
                    if (i == j) {
                        EXPECT_EQ(mtx.get(i, j), 0);
                    } else {
                        EXPECT_EQ(mtx.get(i, j) < 0, 1);
                    }
                }
            }
        }
    }

}