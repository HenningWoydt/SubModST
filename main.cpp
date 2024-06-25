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

#include <iostream>

#include "src/structures/Clustering.h"
#include "src/structures/Matrix.h"
#include "src/algorithms/Solver.h"
#include "src/structures/Graph.h"
#include "src/structures/FacilityLocation.h"
#include "src/structures/WeightedCoverage.h"
#include "src/structures/BipartiteInfluence.h"

int main(int argc, char *argv[]) {

    SubModST::AlgorithmConfiguration ac(argc, argv);

    if (ac.function == "PartialDominatingSet") {
        auto graph = SubModST::Graph(ac.input_file_path);
        SubModST::Matrix<int> mtx = graph.get_PartialDominatingSetMatrix();
        auto solver = SubModST::Solver<SubModST::Matrix<int>, int>(mtx, ac.k, ac);
        solver.search();
        solver.write_output(ac.output_file_path);

    } else if (ac.function == "GroupClosenessCentrality") {
        auto graph = SubModST::Graph(ac.input_file_path);
        SubModST::Matrix<int> mtx = graph.get_GroupClosenessCentralityMatrix();
        auto solver = SubModST::Solver<SubModST::Matrix<int>, int>(mtx, ac.k, ac);
        solver.search();
        solver.write_output(ac.output_file_path);

    } else if (ac.function == "PartialVertexSet") {
        auto graph = SubModST::Graph(ac.input_file_path);
        SubModST::Matrix<int> mtx = graph.get_PartialVertexSetMatrix();
        auto solver = SubModST::Solver<SubModST::Matrix<int>, int>(mtx, ac.k, ac);
        solver.search();
        solver.write_output(ac.output_file_path);

    } else if (ac.function == "EuclidianClustering") {
        auto dataPoints = SubModST::Clustering(ac.input_file_path);
        SubModST::Matrix<double> mtx = dataPoints.get_EuclidianClusteringMatrix();
        auto solver = SubModST::Solver<SubModST::Matrix<double>, double>(mtx, ac.k, ac);
        solver.search();
        solver.write_output(ac.output_file_path);

    } else if (ac.function == "ManhattanClustering") {
        auto dataPoints = SubModST::Clustering(ac.input_file_path);
        SubModST::Matrix<double> mtx = dataPoints.get_ManhattanClusteringMatrix();
        auto solver = SubModST::Solver<SubModST::Matrix<double>, double>(mtx, ac.k, ac);
        solver.search();
        solver.write_output(ac.output_file_path);

    } else if (ac.function == "FacilityLocation") {
        auto facilityLocation = SubModST::FacilityLocation(ac.input_file_path);
        SubModST::Matrix<double> mtx = facilityLocation.get_FacilityLocationMatrix();
        auto solver = SubModST::Solver<SubModST::Matrix<double>, double>(mtx, ac.k, ac);
        solver.search();
        solver.write_output(ac.output_file_path);

    } else if (ac.function == "WeightedCoverage") {
        auto weightedCoverage = SubModST::WeightedCoverage(ac.input_file_path);
        SubModST::Matrix<double> mtx = weightedCoverage.get_WeightedCoverageMatrix();
        auto solver = SubModST::Solver<SubModST::Matrix<double>, double>(mtx, ac.k, ac);
        solver.search();
        solver.write_output(ac.output_file_path);

    } else if (ac.function == "BipartiteInfluence") {
        auto bipartiteInfluence = SubModST::BipartiteInfluence<double>(ac.input_file_path);
        auto solver = SubModST::Solver<SubModST::BipartiteInfluence<double>, double>(bipartiteInfluence, ac.k, ac);
        solver.search();
        solver.write_output(ac.output_file_path);

    } else {
        std::cerr << "Dont recognize Score function \"" << ac.function << "\"!" << std::endl;
        exit(EXIT_FAILURE);
    }

    return EXIT_SUCCESS;
}
