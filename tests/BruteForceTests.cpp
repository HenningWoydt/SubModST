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
#include "../src/structures/Clustering.h"
#include "../src/structures/FacilityLocation.h"
#include "../src/structures/WeightedCoverage.h"

namespace SubModST {
    /**
     * Checks the found sets of the algorithm.
     *
     * @tparam Structure Structure to optimize.
     * @tparam SFType Return data type of the score function.
     * @param bfSolver The Brute-Force solver.
     * @param n n.
     * @param k k.
     */
    template<typename Structure, typename SFType>
    void check_found_sets(BFSolver<Structure, SFType> &bfSolver,
                          size_t n,
                          size_t k) {
        std::vector<std::vector<uint32_t>> all_sets = bfSolver.get_all_s();

        // the number of sets found should be equal to (n choose k)
        EXPECT_EQ(n_choose_k(n, k), all_sets.size());

        // each found set must have size k
        for (std::vector<uint32_t> &vec: all_sets) {
            EXPECT_EQ(vec.size(), k);
        }

        // all sets must differ from each other
        for (size_t i = 0; i < all_sets.size(); ++i) {
            for (size_t j = i + 1; j < all_sets.size(); ++j) {
                EXPECT_EQ(false, all_sets[i] == all_sets[j]);
            }
        }
    }

    void test_PartialDominatingSet(std::string &file_path) {
        auto graph = Graph(file_path);
        Matrix<int> mtx = graph.get_PartialDominatingSetMatrix();
        size_t n = mtx.get_n();

        for (size_t k = 1; k <= n; ++k) {
            auto bfSolver = BFSolver<Matrix<int>, int>(mtx, k, true);
            bfSolver.search();

            check_found_sets(bfSolver, n, k);
        }
    }

    void test_PartialVertexSet(std::string &file_path) {
        auto graph = Graph(file_path);
        Matrix<int> mtx = graph.get_PartialVertexSetMatrix();
        size_t n = mtx.get_n();

        for (size_t k = 1; k <= n; ++k) {
            auto bfSolver = BFSolver<Matrix<int>, int>(mtx, k, true);
            bfSolver.search();

            check_found_sets(bfSolver, n, k);
        }
    }

    void test_GroupClosenessCentrality(std::string &file_path) {
        auto graph = Graph(file_path);
        Matrix<int> mtx = graph.get_GroupClosenessCentralityMatrix();
        size_t n = mtx.get_n();

        for (size_t k = 1; k <= n; ++k) {
            auto bfSolver = BFSolver<Matrix<int>, int>(mtx, k, true);
            bfSolver.search();

            check_found_sets(bfSolver, n, k);
        }
    }

    void test_EuclidianClustering(std::string &file_path) {
        auto dataPoints = Clustering(file_path);
        Matrix<double> mtx = dataPoints.get_EuclidianClusteringMatrix();
        size_t n = mtx.get_n();

        for (size_t k = 1; k <= n; ++k) {
            auto bfSolver = BFSolver<Matrix<double>, double>(mtx, k, true);
            bfSolver.search();

            check_found_sets(bfSolver, n, k);
        }
    }

    void test_ManhattanClustering(std::string &file_path) {
        auto dataPoints = Clustering(file_path);
        Matrix<double> mtx = dataPoints.get_ManhattanClusteringMatrix();
        size_t n = mtx.get_n();

        for (size_t k = 1; k <= n; ++k) {
            auto bfSolver = BFSolver<Matrix<double>, double>(mtx, k, true);
            bfSolver.search();

            check_found_sets(bfSolver, n, k);
        }
    }

    void test_FacilityLocation(std::string &file_path) {
        auto facilityLocation = FacilityLocation(file_path);
        Matrix<double> mtx = facilityLocation.get_FacilityLocationMatrix();
        size_t n = mtx.get_n();

        for (size_t k = 1; k <= n; ++k) {
            auto bfSolver = BFSolver<Matrix<double>, double>(mtx, k, true);
            bfSolver.search();

            check_found_sets(bfSolver, n, k);
        }
    }

    void test_WeightedCoverage(std::string &file_path) {
        auto weightedcoverage = WeightedCoverage(file_path);
        Matrix<double> mtx = weightedcoverage.get_WeightedCoverageMatrix();
        size_t n = mtx.get_n();

        for (size_t k = 1; k <= n; ++k) {
            auto bfSolver = BFSolver<Matrix<double>, double>(mtx, k, true);
            bfSolver.search();

            check_found_sets(bfSolver, n, k);
        }
    }

    void test_BipartiteInfluence([[maybe_unused]] std::string &file_path) {
        // TODO: Fill
    }

    TEST(BruteForceTest, PartialDominatingSet3) {
        std::vector<std::string> files = get_directory_files("../data/private/Graph/3/", ".edges");
        for (std::string &file: files) {
            test_PartialDominatingSet(file);
        }
    }

    TEST(BruteForceTest, PartialDominatingSet4) {
        std::vector<std::string> files = get_directory_files("../data/private/Graph/4/", ".edges");
        for (std::string &file: files) {
            test_PartialDominatingSet(file);
        }
    }

    TEST(BruteForceTest, PartialDominatingSet5) {
        std::vector<std::string> files = get_directory_files("../data/private/Graph/5/", ".edges");
        for (std::string &file: files) {
            test_PartialDominatingSet(file);
        }
    }

    TEST(BruteForceTest, PartialDominatingSet6) {
        std::vector<std::string> files = get_directory_files("../data/private/Graph/6/", ".edges");
        for (std::string &file: files) {
            test_PartialDominatingSet(file);
        }
    }

    TEST(BruteForceTest, PartialDominatingSet7) {
        std::vector<std::string> files = get_directory_files("../data/private/Graph/7/", ".edges");
        for (std::string &file: files) {
            test_PartialDominatingSet(file);
        }
    }

    TEST(BruteForceTest, PartialDominatingSet8) {
        std::vector<std::string> files = get_directory_files("../data/private/Graph/8/", ".edges");
        for (std::string &file: files) {
            test_PartialDominatingSet(file);
        }
    }

    TEST(BruteForceTest, PartialDominatingSet9) {
        std::vector<std::string> files = get_directory_files("../data/private/Graph/9/", ".edges");
        for (std::string &file: files) {
            test_PartialDominatingSet(file);
        }
    }

    TEST(BruteForceTest, PartialDominatingSet10) {
        std::vector<std::string> files = get_directory_files("../data/private/Graph/10/", ".edges");
        for (std::string &file: files) {
            test_PartialDominatingSet(file);
        }
    }

    TEST(BruteForceTest, PartialVertexSet3) {
        std::vector<std::string> files = get_directory_files("../data/private/Graph/3/", ".edges");
        for (std::string &file: files) {
            test_PartialVertexSet(file);
        }
    }

    TEST(BruteForceTest, PartialVertexSet4) {
        std::vector<std::string> files = get_directory_files("../data/private/Graph/4/", ".edges");
        for (std::string &file: files) {
            test_PartialVertexSet(file);
        }
    }

    TEST(BruteForceTest, PartialVertexSet5) {
        std::vector<std::string> files = get_directory_files("../data/private/Graph/5/", ".edges");
        for (std::string &file: files) {
            test_PartialVertexSet(file);
        }
    }

    TEST(BruteForceTest, PartialVertexSet6) {
        std::vector<std::string> files = get_directory_files("../data/private/Graph/6/", ".edges");
        for (std::string &file: files) {
            test_PartialVertexSet(file);
        }
    }

    TEST(BruteForceTest, PartialVertexSet7) {
        std::vector<std::string> files = get_directory_files("../data/private/Graph/7/", ".edges");
        for (std::string &file: files) {
            test_PartialVertexSet(file);
        }
    }

    TEST(BruteForceTest, PartialVertexSet8) {
        std::vector<std::string> files = get_directory_files("../data/private/Graph/8/", ".edges");
        for (std::string &file: files) {
            test_PartialVertexSet(file);
        }
    }

    TEST(BruteForceTest, PartialVertexSet9) {
        std::vector<std::string> files = get_directory_files("../data/private/Graph/9/", ".edges");
        for (std::string &file: files) {
            test_PartialVertexSet(file);
        }
    }

    TEST(BruteForceTest, PartialVertexSet10) {
        std::vector<std::string> files = get_directory_files("../data/private/Graph/10/", ".edges");
        for (std::string &file: files) {
            test_PartialVertexSet(file);
        }
    }

    TEST(BruteForceTest, GroupClosenessCentrality3) {
        std::vector<std::string> files = get_directory_files("../data/private/Graph/3/", ".edges");
        for (std::string &file: files) {
            test_GroupClosenessCentrality(file);
        }
    }

    TEST(BruteForceTest, GroupClosenessCentrality4) {
        std::vector<std::string> files = get_directory_files("../data/private/Graph/4/", ".edges");
        for (std::string &file: files) {
            test_GroupClosenessCentrality(file);
        }
    }

    TEST(BruteForceTest, GroupClosenessCentrality5) {
        std::vector<std::string> files = get_directory_files("../data/private/Graph/5/", ".edges");
        for (std::string &file: files) {
            test_GroupClosenessCentrality(file);
        }
    }

    TEST(BruteForceTest, GroupClosenessCentrality6) {
        std::vector<std::string> files = get_directory_files("../data/private/Graph/6/", ".edges");
        for (std::string &file: files) {
            test_GroupClosenessCentrality(file);
        }
    }

    TEST(BruteForceTest, GroupClosenessCentrality7) {
        std::vector<std::string> files = get_directory_files("../data/private/Graph/7/", ".edges");
        for (std::string &file: files) {
            test_GroupClosenessCentrality(file);
        }
    }

    TEST(BruteForceTest, GroupClosenessCentrality8) {
        std::vector<std::string> files = get_directory_files("../data/private/Graph/8/", ".edges");
        for (std::string &file: files) {
            test_GroupClosenessCentrality(file);
        }
    }

    TEST(BruteForceTest, GroupClosenessCentrality9) {
        std::vector<std::string> files = get_directory_files("../data/private/Graph/9/", ".edges");
        for (std::string &file: files) {
            test_GroupClosenessCentrality(file);
        }
    }

    TEST(BruteForceTest, GroupClosenessCentrality10) {
        std::vector<std::string> files = get_directory_files("../data/private/Graph/10/", ".edges");
        for (std::string &file: files) {
            test_GroupClosenessCentrality(file);
        }
    }

    TEST(BruteForceTest, EuclidianClustering3) {
        std::vector<std::string> files = get_directory_files("../data/private/Clustering/3/", ".mtx");
        for (std::string &file: files) {
            test_EuclidianClustering(file);
        }
    }

    TEST(BruteForceTest, EuclidianClustering4) {
        std::vector<std::string> files = get_directory_files("../data/private/Clustering/4/", ".mtx");
        for (std::string &file: files) {
            test_EuclidianClustering(file);
        }
    }

    TEST(BruteForceTest, EuclidianClustering5) {
        std::vector<std::string> files = get_directory_files("../data/private/Clustering/5/", ".mtx");
        for (std::string &file: files) {
            test_EuclidianClustering(file);
        }
    }

    TEST(BruteForceTest, EuclidianClustering6) {
        std::vector<std::string> files = get_directory_files("../data/private/Clustering/6/", ".mtx");
        for (std::string &file: files) {
            test_EuclidianClustering(file);
        }
    }

    TEST(BruteForceTest, EuclidianClustering7) {
        std::vector<std::string> files = get_directory_files("../data/private/Clustering/7/", ".mtx");
        for (std::string &file: files) {
            test_EuclidianClustering(file);
        }
    }

    TEST(BruteForceTest, EuclidianClustering8) {
        std::vector<std::string> files = get_directory_files("../data/private/Clustering/8/", ".mtx");
        for (std::string &file: files) {
            test_EuclidianClustering(file);
        }
    }

    TEST(BruteForceTest, EuclidianClustering9) {
        std::vector<std::string> files = get_directory_files("../data/private/Clustering/9/", ".mtx");
        for (std::string &file: files) {
            test_EuclidianClustering(file);
        }
    }

    TEST(BruteForceTest, EuclidianClustering10) {
        std::vector<std::string> files = get_directory_files("../data/private/Clustering/10/", ".mtx");
        for (std::string &file: files) {
            test_EuclidianClustering(file);
        }
    }

    TEST(BruteForceTest, ManhattanClustering3) {
        std::vector<std::string> files = get_directory_files("../data/private/Clustering/3/", ".mtx");
        for (std::string &file: files) {
            test_ManhattanClustering(file);
        }
    }

    TEST(BruteForceTest, ManhattanClustering4) {
        std::vector<std::string> files = get_directory_files("../data/private/Clustering/4/", ".mtx");
        for (std::string &file: files) {
            test_ManhattanClustering(file);
        }
    }

    TEST(BruteForceTest, ManhattanClustering5) {
        std::vector<std::string> files = get_directory_files("../data/private/Clustering/5/", ".mtx");
        for (std::string &file: files) {
            test_ManhattanClustering(file);
        }
    }

    TEST(BruteForceTest, ManhattanClustering6) {
        std::vector<std::string> files = get_directory_files("../data/private/Clustering/6/", ".mtx");
        for (std::string &file: files) {
            test_ManhattanClustering(file);
        }
    }

    TEST(BruteForceTest, ManhattanClustering7) {
        std::vector<std::string> files = get_directory_files("../data/private/Clustering/7/", ".mtx");
        for (std::string &file: files) {
            test_ManhattanClustering(file);
        }
    }

    TEST(BruteForceTest, ManhattanClustering8) {
        std::vector<std::string> files = get_directory_files("../data/private/Clustering/8/", ".mtx");
        for (std::string &file: files) {
            test_ManhattanClustering(file);
        }
    }

    TEST(BruteForceTest, ManhattanClustering9) {
        std::vector<std::string> files = get_directory_files("../data/private/Clustering/9/", ".mtx");
        for (std::string &file: files) {
            test_ManhattanClustering(file);
        }
    }

    TEST(BruteForceTest, ManhattanClustering10) {
        std::vector<std::string> files = get_directory_files("../data/private/Clustering/10/", ".mtx");
        for (std::string &file: files) {
            test_ManhattanClustering(file);
        }
    }

    TEST(BruteForceTest, FacilityLocation3) {
        std::vector<std::string> files = get_directory_files("../data/private/FacilityLocation/3/", ".csv");
        for (std::string &file: files) {
            test_FacilityLocation(file);
        }
    }

    TEST(BruteForceTest, FacilityLocation4) {
        std::vector<std::string> files = get_directory_files("../data/private/FacilityLocation/4/", ".csv");
        for (std::string &file: files) {
            test_FacilityLocation(file);
        }
    }

    TEST(BruteForceTest, FacilityLocation5) {
        std::vector<std::string> files = get_directory_files("../data/private/FacilityLocation/5/", ".csv");
        for (std::string &file: files) {
            test_FacilityLocation(file);
        }
    }

    TEST(BruteForceTest, FacilityLocation6) {
        std::vector<std::string> files = get_directory_files("../data/private/FacilityLocation/6/", ".csv");
        for (std::string &file: files) {
            test_FacilityLocation(file);
        }
    }

    TEST(BruteForceTest, FacilityLocation7) {
        std::vector<std::string> files = get_directory_files("../data/private/FacilityLocation/7/", ".csv");
        for (std::string &file: files) {
            test_FacilityLocation(file);
        }
    }

    TEST(BruteForceTest, FacilityLocation8) {
        std::vector<std::string> files = get_directory_files("../data/private/FacilityLocation/8/", ".csv");
        for (std::string &file: files) {
            test_FacilityLocation(file);
        }
    }

    TEST(BruteForceTest, FacilityLocation9) {
        std::vector<std::string> files = get_directory_files("../data/private/FacilityLocation/9/", ".csv");
        for (std::string &file: files) {
            test_FacilityLocation(file);
        }
    }

    TEST(BruteForceTest, FacilityLocation10) {
        std::vector<std::string> files = get_directory_files("../data/private/FacilityLocation/10/", ".csv");
        for (std::string &file: files) {
            test_FacilityLocation(file);
        }
    }

    TEST(BruteForceTest, WeightedCoverage3) {
        std::vector<std::string> files = get_directory_files("../data/private/WeightedCoverage/3/", ".csv");
        for (std::string &file: files) {
            test_WeightedCoverage(file);
        }
    }

    TEST(BruteForceTest, WeightedCoverage4) {
        std::vector<std::string> files = get_directory_files("../data/private/WeightedCoverage/4/", ".csv");
        for (std::string &file: files) {
            test_WeightedCoverage(file);
        }
    }

    TEST(BruteForceTest, WeightedCoverage5) {
        std::vector<std::string> files = get_directory_files("../data/private/WeightedCoverage/5/", ".csv");
        for (std::string &file: files) {
            test_WeightedCoverage(file);
        }
    }

    TEST(BruteForceTest, WeightedCoverage6) {
        std::vector<std::string> files = get_directory_files("../data/private/WeightedCoverage/6/", ".csv");
        for (std::string &file: files) {
            test_WeightedCoverage(file);
        }
    }

    TEST(BruteForceTest, WeightedCoverage7) {
        std::vector<std::string> files = get_directory_files("../data/private/WeightedCoverage/7/", ".csv");
        for (std::string &file: files) {
            test_WeightedCoverage(file);
        }
    }

    TEST(BruteForceTest, WeightedCoverage8) {
        std::vector<std::string> files = get_directory_files("../data/private/WeightedCoverage/8/", ".csv");
        for (std::string &file: files) {
            test_WeightedCoverage(file);
        }
    }

    TEST(BruteForceTest, WeightedCoverage9) {
        std::vector<std::string> files = get_directory_files("../data/private/WeightedCoverage/9/", ".csv");
        for (std::string &file: files) {
            test_WeightedCoverage(file);
        }
    }

    TEST(BruteForceTest, WeightedCoverage10) {
        std::vector<std::string> files = get_directory_files("../data/private/WeightedCoverage/10/", ".csv");
        for (std::string &file: files) {
            test_WeightedCoverage(file);
        }
    }

    TEST(BruteForceTest, BipartiteInfluence3) {
        std::vector<std::string> files = get_directory_files("../data/private/BipartiteInfluence/3/", ".csv");
        for (std::string &file: files) {
            test_BipartiteInfluence(file);
        }
    }

    TEST(BruteForceTest, BipartiteInfluence4) {
        std::vector<std::string> files = get_directory_files("../data/private/BipartiteInfluence/4/", ".csv");
        for (std::string &file: files) {
            test_BipartiteInfluence(file);
        }
    }

    TEST(BruteForceTest, BipartiteInfluence5) {
        std::vector<std::string> files = get_directory_files("../data/private/BipartiteInfluence/5/", ".csv");
        for (std::string &file: files) {
            test_BipartiteInfluence(file);
        }
    }

    TEST(BruteForceTest, BipartiteInfluence6) {
        std::vector<std::string> files = get_directory_files("../data/private/BipartiteInfluence/6/", ".csv");
        for (std::string &file: files) {
            test_BipartiteInfluence(file);
        }
    }

    TEST(BruteForceTest, BipartiteInfluence7) {
        std::vector<std::string> files = get_directory_files("../data/private/BipartiteInfluence/7/", ".csv");
        for (std::string &file: files) {
            test_BipartiteInfluence(file);
        }
    }

    TEST(BruteForceTest, BipartiteInfluence8) {
        std::vector<std::string> files = get_directory_files("../data/private/BipartiteInfluence/8/", ".csv");
        for (std::string &file: files) {
            test_BipartiteInfluence(file);
        }
    }

    TEST(BruteForceTest, BipartiteInfluence9) {
        std::vector<std::string> files = get_directory_files("../data/private/BipartiteInfluence/9/", ".csv");
        for (std::string &file: files) {
            test_BipartiteInfluence(file);
        }
    }

    TEST(BruteForceTest, BipartiteInfluence10) {
        std::vector<std::string> files = get_directory_files("../data/private/BipartiteInfluence/10/", ".csv");
        for (std::string &file: files) {
            test_BipartiteInfluence(file);
        }
    }
}
