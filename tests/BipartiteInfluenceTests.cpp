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
#include "../src/structures/BipartiteInfluence.h"

namespace SubModST {

    void test_BFSolverVsSolver_BipartiteInfluence([[maybe_unused]] std::string &data_file_path, [[maybe_unused]] std::string &sol_dir_path) {
        auto bipartiteInfluence = BipartiteInfluence<double>(data_file_path);
        size_t n = bipartiteInfluence.get_n();

        std::string data_file_name = std::filesystem::path(data_file_path).filename();

        for (size_t k = 1; k <= n; ++k) {

            std::vector<uint32_t> bf_solution;

            // check if the solution file exists
            std::string sol_path = sol_dir_path + data_file_name + "-" + std::to_string(k) + ".sol";
            if (file_exists(sol_path)) {
                bf_solution = read_solution(sol_path);
            } else {
                auto bfSolver = BFSolver<BipartiteInfluence<double>, double>(bipartiteInfluence, k, false);
                bfSolver.search();
                bf_solution = bfSolver.get_solution();
                write_solution(bf_solution, sol_path);
            }

            std::vector<AlgorithmConfiguration> acs = get_all_acs();
            for (auto &ac: acs) {

                auto solver = Solver<BipartiteInfluence<double>, double>(bipartiteInfluence, k, ac);
                solver.search();

                double bf_res = bipartiteInfluence.evaluate(bf_solution, bf_solution.size());
                double res = solver.get_solution_score();

                EXPECT_NEAR(bf_res, res, 0.000000001) << data_file_path << " n: " << n << " k: " << k << " name: " << ac.nickname << std::endl;
            }
        }
    }

    TEST(BipartiteInfluenceTest, BipartiteInfluence3) {
        std::string data_dir_path = "../data/private/BipartiteInfluence/3/";
        std::string extension = ".binf";
        std::string sol_dir_path = data_dir_path + "BipartiteInfluence/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_BipartiteInfluence(file, sol_dir_path);
        }
    }

    TEST(BipartiteInfluenceTest, BipartiteInfluence4) {
        std::string data_dir_path = "../data/private/BipartiteInfluence/4/";
        std::string extension = ".binf";
        std::string sol_dir_path = data_dir_path + "BipartiteInfluence/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_BipartiteInfluence(file, sol_dir_path);
        }
    }

    TEST(BipartiteInfluenceTest, BipartiteInfluence5) {
        std::string data_dir_path = "../data/private/BipartiteInfluence/5/";
        std::string extension = ".binf";
        std::string sol_dir_path = data_dir_path + "BipartiteInfluence/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_BipartiteInfluence(file, sol_dir_path);
        }
    }

    TEST(BipartiteInfluenceTest, BipartiteInfluence6) {
        std::string data_dir_path = "../data/private/BipartiteInfluence/6/";
        std::string extension = ".binf";
        std::string sol_dir_path = data_dir_path + "BipartiteInfluence/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_BipartiteInfluence(file, sol_dir_path);
        }
    }

    TEST(BipartiteInfluenceTest, BipartiteInfluence7) {
        std::string data_dir_path = "../data/private/BipartiteInfluence/7/";
        std::string extension = ".binf";
        std::string sol_dir_path = data_dir_path + "BipartiteInfluence/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_BipartiteInfluence(file, sol_dir_path);
        }
    }

    TEST(BipartiteInfluenceTest, BipartiteInfluence8) {
        std::string data_dir_path = "../data/private/BipartiteInfluence/8/";
        std::string extension = ".binf";
        std::string sol_dir_path = data_dir_path + "BipartiteInfluence/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_BipartiteInfluence(file, sol_dir_path);
        }
    }

    TEST(BipartiteInfluenceTest, BipartiteInfluence9) {
        std::string data_dir_path = "../data/private/BipartiteInfluence/9/";
        std::string extension = ".binf";
        std::string sol_dir_path = data_dir_path + "BipartiteInfluence/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_BipartiteInfluence(file, sol_dir_path);
        }
    }

    TEST(BipartiteInfluenceTest, BipartiteInfluence10) {
        std::string data_dir_path = "../data/private/BipartiteInfluence/10/";
        std::string extension = ".binf";
        std::string sol_dir_path = data_dir_path + "BipartiteInfluence/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_BipartiteInfluence(file, sol_dir_path);
        }
    }

    TEST(BipartiteInfluenceTest, BipartiteInfluence11) {
        std::string data_dir_path = "../data/private/BipartiteInfluence/11/";
        std::string extension = ".binf";
        std::string sol_dir_path = data_dir_path + "BipartiteInfluence/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_BipartiteInfluence(file, sol_dir_path);
        }
    }

    TEST(BipartiteInfluenceTest, BipartiteInfluence12) {
        std::string data_dir_path = "../data/private/BipartiteInfluence/12/";
        std::string extension = ".binf";
        std::string sol_dir_path = data_dir_path + "BipartiteInfluence/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_BipartiteInfluence(file, sol_dir_path);
        }
    }

    TEST(BipartiteInfluenceTest, BipartiteInfluence13) {
        std::string data_dir_path = "../data/private/BipartiteInfluence/13/";
        std::string extension = ".binf";
        std::string sol_dir_path = data_dir_path + "BipartiteInfluence/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_BipartiteInfluence(file, sol_dir_path);
        }
    }

    TEST(BipartiteInfluenceTest, BipartiteInfluence14) {
        std::string data_dir_path = "../data/private/BipartiteInfluence/14/";
        std::string extension = ".binf";
        std::string sol_dir_path = data_dir_path + "BipartiteInfluence/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_BipartiteInfluence(file, sol_dir_path);
        }
    }

    TEST(BipartiteInfluenceTest, BipartiteInfluence15) {
        std::string data_dir_path = "../data/private/BipartiteInfluence/15/";
        std::string extension = ".binf";
        std::string sol_dir_path = data_dir_path + "BipartiteInfluence/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_BipartiteInfluence(file, sol_dir_path);
        }
    }

    TEST(BipartiteInfluenceTest, BipartiteInfluence16) {
        std::string data_dir_path = "../data/private/BipartiteInfluence/16/";
        std::string extension = ".binf";
        std::string sol_dir_path = data_dir_path + "BipartiteInfluence/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_BipartiteInfluence(file, sol_dir_path);
        }
    }

    TEST(BipartiteInfluenceTest, BipartiteInfluence17) {
        std::string data_dir_path = "../data/private/BipartiteInfluence/17/";
        std::string extension = ".binf";
        std::string sol_dir_path = data_dir_path + "BipartiteInfluence/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_BipartiteInfluence(file, sol_dir_path);
        }
    }

    TEST(BipartiteInfluenceTest, BipartiteInfluence18) {
        std::string data_dir_path = "../data/private/BipartiteInfluence/18/";
        std::string extension = ".binf";
        std::string sol_dir_path = data_dir_path + "BipartiteInfluence/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_BipartiteInfluence(file, sol_dir_path);
        }
    }

    TEST(BipartiteInfluenceTest, BipartiteInfluence19) {
        std::string data_dir_path = "../data/private/BipartiteInfluence/19/";
        std::string extension = ".binf";
        std::string sol_dir_path = data_dir_path + "BipartiteInfluence/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_BipartiteInfluence(file, sol_dir_path);
        }
    }

    TEST(BipartiteInfluenceTest, BipartiteInfluence20) {
        std::string data_dir_path = "../data/private/BipartiteInfluence/20/";
        std::string extension = ".binf";
        std::string sol_dir_path = data_dir_path + "BipartiteInfluence/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_BipartiteInfluence(file, sol_dir_path);
        }
    }
}
