/* CCSMSM solver, that solves the Cardinality-Constrained Submodular Monotone
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
#include <filesystem>

#include <gtest/gtest.h>

#include "../src/algorithms/BFSolver.h"
#include "../src/algorithms/Solver.h"
#include "../src/structures/Graph.h"
#include "../src/structures/Matrix.h"
#include "../src/util/utility.h"

namespace CCSMSM {

    void test_BFSolverVsSolver_PartialDominatingSet(std::string &data_file_path, std::string &sol_dir_path) {
        auto graph = Graph(data_file_path);
        Matrix<int> mtx = graph.get_PartialDominatingSetMatrix();
        size_t n = mtx.get_n();

        std::string data_file_name = std::filesystem::path(data_file_path).filename();

        for (size_t k = 1; k <= n; ++k) {

            std::vector<uint32_t> bf_solution;

            // check if the solution file exists
            std::string sol_path = sol_dir_path + data_file_name + "-" + std::to_string(k) + ".sol";
            if (file_exists(sol_path)) {
                bf_solution = read_solution(sol_path);
            } else {
                auto bfSolver = BFSolver<Matrix<int>, int>(mtx, k, false);
                bfSolver.search();
                bf_solution = bfSolver.get_solution();
                write_solution(bf_solution, sol_path);
            }

            std::vector<AlgorithmConfiguration> acs = get_all_acs();
            for (auto &ac: acs) {

                auto solver = Solver<Matrix<int>, int>(mtx, k, ac);
                solver.search();

                int bf_res = mtx.evaluate(bf_solution, bf_solution.size());
                int res = solver.get_solution_score();

                EXPECT_EQ(bf_res, res) << data_file_path << " n: " << n << " k: " << k << std::endl;
            }
        }
    }

    TEST(PartialDominatingSetTest, PartialDominatingSet3) {
        std::string data_dir_path = "../data/private/Graph/3/";
        std::string extension = ".edges";
        std::string sol_dir_path = data_dir_path + "PartialDominatingSet/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_PartialDominatingSet(file, sol_dir_path);
        }
    }

    TEST(PartialDominatingSetTest, PartialDominatingSet4) {
        std::string data_dir_path = "../data/private/Graph/4/";
        std::string extension = ".edges";
        std::string sol_dir_path = data_dir_path + "PartialDominatingSet/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_PartialDominatingSet(file, sol_dir_path);
        }
    }

    TEST(PartialDominatingSetTest, PartialDominatingSet5) {
        std::string data_dir_path = "../data/private/Graph/5/";
        std::string extension = ".edges";
        std::string sol_dir_path = data_dir_path + "PartialDominatingSet/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_PartialDominatingSet(file, sol_dir_path);
        }
    }

    TEST(PartialDominatingSetTest, PartialDominatingSet6) {
        std::string data_dir_path = "../data/private/Graph/6/";
        std::string extension = ".edges";
        std::string sol_dir_path = data_dir_path + "PartialDominatingSet/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_PartialDominatingSet(file, sol_dir_path);
        }
    }

    TEST(PartialDominatingSetTest, PartialDominatingSet7) {
        std::string data_dir_path = "../data/private/Graph/7/";
        std::string extension = ".edges";
        std::string sol_dir_path = data_dir_path + "PartialDominatingSet/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_PartialDominatingSet(file, sol_dir_path);
        }
    }

    TEST(PartialDominatingSetTest, PartialDominatingSet8) {
        std::string data_dir_path = "../data/private/Graph/8/";
        std::string extension = ".edges";
        std::string sol_dir_path = data_dir_path + "PartialDominatingSet/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_PartialDominatingSet(file, sol_dir_path);
        }
    }

    TEST(PartialDominatingSetTest, PartialDominatingSet9) {
        std::string data_dir_path = "../data/private/Graph/9/";
        std::string extension = ".edges";
        std::string sol_dir_path = data_dir_path + "PartialDominatingSet/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_PartialDominatingSet(file, sol_dir_path);
        }
    }

    TEST(PartialDominatingSetTest, PartialDominatingSet10) {
        std::string data_dir_path = "../data/private/Graph/10/";
        std::string extension = ".edges";
        std::string sol_dir_path = data_dir_path + "PartialDominatingSet/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_PartialDominatingSet(file, sol_dir_path);
        }
    }

    TEST(PartialDominatingSetTest, PartialDominatingSet11) {
        std::string data_dir_path = "../data/private/Graph/11/";
        std::string extension = ".edges";
        std::string sol_dir_path = data_dir_path + "PartialDominatingSet/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_PartialDominatingSet(file, sol_dir_path);
        }
    }

    TEST(PartialDominatingSetTest, PartialDominatingSet12) {
        std::string data_dir_path = "../data/private/Graph/12/";
        std::string extension = ".edges";
        std::string sol_dir_path = data_dir_path + "PartialDominatingSet/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_PartialDominatingSet(file, sol_dir_path);
        }
    }

    TEST(PartialDominatingSetTest, PartialDominatingSet13) {
        std::string data_dir_path = "../data/private/Graph/13/";
        std::string extension = ".edges";
        std::string sol_dir_path = data_dir_path + "PartialDominatingSet/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_PartialDominatingSet(file, sol_dir_path);
        }
    }

    TEST(PartialDominatingSetTest, PartialDominatingSet14) {
        std::string data_dir_path = "../data/private/Graph/14/";
        std::string extension = ".edges";
        std::string sol_dir_path = data_dir_path + "PartialDominatingSet/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_PartialDominatingSet(file, sol_dir_path);
        }
    }

    TEST(PartialDominatingSetTest, PartialDominatingSet15) {
        std::string data_dir_path = "../data/private/Graph/15/";
        std::string extension = ".edges";
        std::string sol_dir_path = data_dir_path + "PartialDominatingSet/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_PartialDominatingSet(file, sol_dir_path);
        }
    }

    TEST(PartialDominatingSetTest, PartialDominatingSet16) {
        std::string data_dir_path = "../data/private/Graph/16/";
        std::string extension = ".edges";
        std::string sol_dir_path = data_dir_path + "PartialDominatingSet/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_PartialDominatingSet(file, sol_dir_path);
        }
    }

    TEST(PartialDominatingSetTest, PartialDominatingSet17) {
        std::string data_dir_path = "../data/private/Graph/17/";
        std::string extension = ".edges";
        std::string sol_dir_path = data_dir_path + "PartialDominatingSet/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_PartialDominatingSet(file, sol_dir_path);
        }
    }

    TEST(PartialDominatingSetTest, PartialDominatingSet18) {
        std::string data_dir_path = "../data/private/Graph/18/";
        std::string extension = ".edges";
        std::string sol_dir_path = data_dir_path + "PartialDominatingSet/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_PartialDominatingSet(file, sol_dir_path);
        }
    }

    TEST(PartialDominatingSetTest, PartialDominatingSet19) {
        std::string data_dir_path = "../data/private/Graph/19/";
        std::string extension = ".edges";
        std::string sol_dir_path = data_dir_path + "PartialDominatingSet/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_PartialDominatingSet(file, sol_dir_path);
        }
    }

    TEST(PartialDominatingSetTest, PartialDominatingSet20) {
        std::string data_dir_path = "../data/private/Graph/20/";
        std::string extension = ".edges";
        std::string sol_dir_path = data_dir_path + "PartialDominatingSet/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_PartialDominatingSet(file, sol_dir_path);
        }
    }
}