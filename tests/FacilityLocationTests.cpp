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
#include "../src/structures/FacilityLocation.h"

namespace SubModST {
    void test_BFSolverVsSolver_FacilityLocation(std::string &data_file_path, std::string &sol_dir_path) {
        auto facilityLocation = FacilityLocation(data_file_path);
        Matrix<double> mtx = facilityLocation.get_FacilityLocationMatrix();
        size_t n = mtx.get_n();

        std::string data_file_name = std::filesystem::path(data_file_path).filename();

        for (size_t k = 1; k <= n; ++k) {

            std::vector<uint32_t> bf_solution;

            // check if the solution file exists
            std::string sol_path = sol_dir_path + data_file_name + "-" + std::to_string(k) + ".sol";
            if (file_exists(sol_path)) {
                bf_solution = read_solution(sol_path);
            } else {
                auto bfSolver = BFSolver<Matrix<double>, double>(mtx, k, false);
                bfSolver.search();
                bf_solution = bfSolver.get_solution();
                write_solution(bf_solution, sol_path);
            }

            std::vector<AlgorithmConfiguration> acs = get_all_acs();
            for (auto &ac: acs) {

                auto solver = Solver<Matrix<double>, double>(mtx, k, ac);
                solver.search();

                double bf_res = mtx.evaluate(bf_solution, bf_solution.size());
                double res = solver.get_solution_score();

                EXPECT_NEAR(bf_res, res, 0.000000001) << data_file_path << " n: " << n << " k: " << k << std::endl;
            }
        }
    }

    TEST(FacilityLocationTest, FacilityLocation3) {
        std::string data_dir_path = "../data/private/FacilityLocation/3/";
        std::string extension = ".csv";
        std::string sol_dir_path = data_dir_path + "FacilityLocation/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_FacilityLocation(file, sol_dir_path);
        }
    }

    TEST(FacilityLocationTest, FacilityLocation4) {
        std::string data_dir_path = "../data/private/FacilityLocation/4/";
        std::string extension = ".csv";
        std::string sol_dir_path = data_dir_path + "FacilityLocation/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_FacilityLocation(file, sol_dir_path);
        }
    }

    TEST(FacilityLocationTest, FacilityLocation5) {
        std::string data_dir_path = "../data/private/FacilityLocation/5/";
        std::string extension = ".csv";
        std::string sol_dir_path = data_dir_path + "FacilityLocation/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_FacilityLocation(file, sol_dir_path);
        }
    }

    TEST(FacilityLocationTest, FacilityLocation6) {
        std::string data_dir_path = "../data/private/FacilityLocation/6/";
        std::string extension = ".csv";
        std::string sol_dir_path = data_dir_path + "FacilityLocation/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_FacilityLocation(file, sol_dir_path);
        }
    }

    TEST(FacilityLocationTest, FacilityLocation7) {
        std::string data_dir_path = "../data/private/FacilityLocation/7/";
        std::string extension = ".csv";
        std::string sol_dir_path = data_dir_path + "FacilityLocation/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_FacilityLocation(file, sol_dir_path);
        }
    }

    TEST(FacilityLocationTest, FacilityLocation8) {
        std::string data_dir_path = "../data/private/FacilityLocation/8/";
        std::string extension = ".csv";
        std::string sol_dir_path = data_dir_path + "FacilityLocation/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_FacilityLocation(file, sol_dir_path);
        }
    }

    TEST(FacilityLocationTest, FacilityLocation9) {
        std::string data_dir_path = "../data/private/FacilityLocation/9/";
        std::string extension = ".csv";
        std::string sol_dir_path = data_dir_path + "FacilityLocation/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_FacilityLocation(file, sol_dir_path);
        }
    }

    TEST(FacilityLocationTest, FacilityLocation10) {
        std::string data_dir_path = "../data/private/FacilityLocation/10/";
        std::string extension = ".csv";
        std::string sol_dir_path = data_dir_path + "FacilityLocation/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_FacilityLocation(file, sol_dir_path);
        }
    }

    TEST(FacilityLocationTest, FacilityLocation11) {
        std::string data_dir_path = "../data/private/FacilityLocation/11/";
        std::string extension = ".csv";
        std::string sol_dir_path = data_dir_path + "FacilityLocation/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_FacilityLocation(file, sol_dir_path);
        }
    }

    TEST(FacilityLocationTest, FacilityLocation12) {
        std::string data_dir_path = "../data/private/FacilityLocation/12/";
        std::string extension = ".csv";
        std::string sol_dir_path = data_dir_path + "FacilityLocation/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_FacilityLocation(file, sol_dir_path);
        }
    }

    TEST(FacilityLocationTest, FacilityLocation13) {
        std::string data_dir_path = "../data/private/FacilityLocation/13/";
        std::string extension = ".csv";
        std::string sol_dir_path = data_dir_path + "FacilityLocation/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_FacilityLocation(file, sol_dir_path);
        }
    }

    TEST(FacilityLocationTest, FacilityLocation14) {
        std::string data_dir_path = "../data/private/FacilityLocation/14/";
        std::string extension = ".csv";
        std::string sol_dir_path = data_dir_path + "FacilityLocation/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_FacilityLocation(file, sol_dir_path);
        }
    }

    TEST(FacilityLocationTest, FacilityLocation15) {
        std::string data_dir_path = "../data/private/FacilityLocation/15/";
        std::string extension = ".csv";
        std::string sol_dir_path = data_dir_path + "FacilityLocation/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_FacilityLocation(file, sol_dir_path);
        }
    }

    TEST(FacilityLocationTest, FacilityLocation16) {
        std::string data_dir_path = "../data/private/FacilityLocation/16/";
        std::string extension = ".csv";
        std::string sol_dir_path = data_dir_path + "FacilityLocation/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_FacilityLocation(file, sol_dir_path);
        }
    }

    TEST(FacilityLocationTest, FacilityLocation17) {
        std::string data_dir_path = "../data/private/FacilityLocation/17/";
        std::string extension = ".csv";
        std::string sol_dir_path = data_dir_path + "FacilityLocation/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_FacilityLocation(file, sol_dir_path);
        }
    }

    TEST(FacilityLocationTest, FacilityLocation18) {
        std::string data_dir_path = "../data/private/FacilityLocation/18/";
        std::string extension = ".csv";
        std::string sol_dir_path = data_dir_path + "FacilityLocation/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_FacilityLocation(file, sol_dir_path);
        }
    }

    TEST(FacilityLocationTest, FacilityLocation19) {
        std::string data_dir_path = "../data/private/FacilityLocation/19/";
        std::string extension = ".csv";
        std::string sol_dir_path = data_dir_path + "FacilityLocation/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_FacilityLocation(file, sol_dir_path);
        }
    }

    TEST(FacilityLocationTest, FacilityLocation20) {
        std::string data_dir_path = "../data/private/FacilityLocation/20/";
        std::string extension = ".csv";
        std::string sol_dir_path = data_dir_path + "FacilityLocation/";

        std::filesystem::create_directories(data_dir_path);
        std::filesystem::create_directories(sol_dir_path);

        std::vector<std::string> files = get_directory_files(data_dir_path, extension);
        for (std::string &file: files) {
            test_BFSolverVsSolver_FacilityLocation(file, sol_dir_path);
        }
    }
}
