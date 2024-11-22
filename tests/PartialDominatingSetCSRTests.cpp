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
#include <filesystem>

#include <gtest/gtest.h>

#include "../src/algorithms/BFSolver.h"
#include "../src/algorithms/Solver.h"
#include "../src/structures/CSRGraphPDS.h"
#include "../src/util/utility.h"

namespace SubModST {

    void test_BFSolverVsSolver_PartialDominatingSetCSR(std::string &data_file_path, std::string &sol_dir_path) {
        auto graph = CSRGraphPDS<int>(data_file_path);
        size_t n = graph.get_n();

        std::string data_file_name = std::filesystem::path(data_file_path).filename();

        for (size_t k = 1; k <= n; ++k) {

            std::vector<uint32_t> bf_solution;

            // check if the solution file exists
            std::string sol_path = sol_dir_path + data_file_name + "-" + std::to_string(k) + ".sol";
            if (file_exists(sol_path)) {
                bf_solution = read_solution(sol_path);
            } else {
                auto bfSolver = BFSolver<CSRGraphPDS<int>, int>(graph, k, false);
                bfSolver.search();
                bf_solution = bfSolver.get_solution();
                write_solution(bf_solution, sol_path);
            }

            std::vector<AlgorithmConfiguration> acs = get_all_acs();
            for (auto &ac: acs) {

                auto solver = Solver<CSRGraphPDS<int>, int>(graph, k, ac);
                solver.search();

                int bf_res = graph.evaluate(bf_solution, bf_solution.size());
                int res = solver.get_solution_score();

                EXPECT_EQ(bf_res, res) << data_file_path << " n: " << n << " k: " << k << std::endl;
            }
        }
    }

    TEST(PartialDominatingSetCSRTest, PartialDominatingSetCSR3) {
        std::string data_dir_path = "../data/private/CSRGraph/3/";
        std::string extension = ".gr";
        std::string sol_dir_path = data_dir_path + "PartialDominatingSetCSR/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_PartialDominatingSetCSR(file, sol_dir_path);
        }
    }

    TEST(PartialDominatingSetCSRTest, PartialDominatingSetCSR4) {
        std::string data_dir_path = "../data/private/CSRGraph/4/";
        std::string extension = ".gr";
        std::string sol_dir_path = data_dir_path + "PartialDominatingSetCSR/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_PartialDominatingSetCSR(file, sol_dir_path);
        }
    }

    TEST(PartialDominatingSetCSRTest, PartialDominatingSetCSR5) {
        std::string data_dir_path = "../data/private/CSRGraph/5/";
        std::string extension = ".gr";
        std::string sol_dir_path = data_dir_path + "PartialDominatingSetCSR/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_PartialDominatingSetCSR(file, sol_dir_path);
        }
    }

    TEST(PartialDominatingSetCSRTest, PartialDominatingSetCSR6) {
        std::string data_dir_path = "../data/private/CSRGraph/6/";
        std::string extension = ".gr";
        std::string sol_dir_path = data_dir_path + "PartialDominatingSetCSR/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_PartialDominatingSetCSR(file, sol_dir_path);
        }
    }

    TEST(PartialDominatingSetCSRTest, PartialDominatingSetCSR7) {
        std::string data_dir_path = "../data/private/CSRGraph/7/";
        std::string extension = ".gr";
        std::string sol_dir_path = data_dir_path + "PartialDominatingSetCSR/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_PartialDominatingSetCSR(file, sol_dir_path);
        }
    }

    TEST(PartialDominatingSetCSRTest, PartialDominatingSetCSR8) {
        std::string data_dir_path = "../data/private/CSRGraph/8/";
        std::string extension = ".gr";
        std::string sol_dir_path = data_dir_path + "PartialDominatingSetCSR/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_PartialDominatingSetCSR(file, sol_dir_path);
        }
    }

    TEST(PartialDominatingSetCSRTest, PartialDominatingSetCSR9) {
        std::string data_dir_path = "../data/private/CSRGraph/9/";
        std::string extension = ".gr";
        std::string sol_dir_path = data_dir_path + "PartialDominatingSetCSR/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_PartialDominatingSetCSR(file, sol_dir_path);
        }
    }

    TEST(PartialDominatingSetCSRTest, PartialDominatingSetCSR10) {
        std::string data_dir_path = "../data/private/CSRGraph/10/";
        std::string extension = ".gr";
        std::string sol_dir_path = data_dir_path + "PartialDominatingSetCSR/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_PartialDominatingSetCSR(file, sol_dir_path);
        }
    }

    TEST(PartialDominatingSetCSRTest, PartialDominatingSetCSR11) {
        std::string data_dir_path = "../data/private/CSRGraph/11/";
        std::string extension = ".gr";
        std::string sol_dir_path = data_dir_path + "PartialDominatingSetCSR/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_PartialDominatingSetCSR(file, sol_dir_path);
        }
    }

    TEST(PartialDominatingSetCSRTest, PartialDominatingSetCSR12) {
        std::string data_dir_path = "../data/private/CSRGraph/12/";
        std::string extension = ".gr";
        std::string sol_dir_path = data_dir_path + "PartialDominatingSetCSR/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_PartialDominatingSetCSR(file, sol_dir_path);
        }
    }

    TEST(PartialDominatingSetCSRTest, PartialDominatingSetCSR13) {
        std::string data_dir_path = "../data/private/CSRGraph/13/";
        std::string extension = ".gr";
        std::string sol_dir_path = data_dir_path + "PartialDominatingSetCSR/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_PartialDominatingSetCSR(file, sol_dir_path);
        }
    }

    TEST(PartialDominatingSetCSRTest, PartialDominatingSetCSR14) {
        std::string data_dir_path = "../data/private/CSRGraph/14/";
        std::string extension = ".gr";
        std::string sol_dir_path = data_dir_path + "PartialDominatingSetCSR/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_PartialDominatingSetCSR(file, sol_dir_path);
        }
    }

    TEST(PartialDominatingSetCSRTest, PartialDominatingSetCSR15) {
        std::string data_dir_path = "../data/private/CSRGraph/15/";
        std::string extension = ".gr";
        std::string sol_dir_path = data_dir_path + "PartialDominatingSetCSR/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_PartialDominatingSetCSR(file, sol_dir_path);
        }
    }

    TEST(PartialDominatingSetCSRTest, PartialDominatingSetCSR16) {
        std::string data_dir_path = "../data/private/CSRGraph/16/";
        std::string extension = ".gr";
        std::string sol_dir_path = data_dir_path + "PartialDominatingSetCSR/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_PartialDominatingSetCSR(file, sol_dir_path);
        }
    }

    TEST(PartialDominatingSetCSRTest, PartialDominatingSetCSR17) {
        std::string data_dir_path = "../data/private/CSRGraph/17/";
        std::string extension = ".gr";
        std::string sol_dir_path = data_dir_path + "PartialDominatingSetCSR/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_PartialDominatingSetCSR(file, sol_dir_path);
        }
    }

    TEST(PartialDominatingSetCSRTest, PartialDominatingSetCSR18) {
        std::string data_dir_path = "../data/private/CSRGraph/18/";
        std::string extension = ".gr";
        std::string sol_dir_path = data_dir_path + "PartialDominatingSetCSR/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_PartialDominatingSetCSR(file, sol_dir_path);
        }
    }

    TEST(PartialDominatingSetCSRTest, PartialDominatingSetCSR19) {
        std::string data_dir_path = "../data/private/CSRGraph/19/";
        std::string extension = ".gr";
        std::string sol_dir_path = data_dir_path + "PartialDominatingSetCSR/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_PartialDominatingSetCSR(file, sol_dir_path);
        }
    }

    TEST(PartialDominatingSetCSRTest, PartialDominatingSetCSR20) {
        std::string data_dir_path = "../data/private/CSRGraph/20/";
        std::string extension = ".gr";
        std::string sol_dir_path = data_dir_path + "PartialDominatingSetCSR/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_PartialDominatingSetCSR(file, sol_dir_path);
        }
    }
}