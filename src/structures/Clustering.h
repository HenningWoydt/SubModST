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

#ifndef CCSMSM_CLUSTERING_H
#define CCSMSM_CLUSTERING_H

#include <cstddef>
#include <vector>
#include <string>
#include <iterator>
#include <cmath>

#include "../util/utility.h"
#include "Matrix.h"

namespace CCSMSM {
    /**
     * Class for defining a clustering instance (set of vectors).
     */
    class Clustering {
    private:
        size_t m_n_vec; // number of vectors
        size_t m_dim; // dimensionality of each point
        std::vector<std::vector<double>> m_vecs;

    public:
        /**
         * Reads a clustering instance from the specified file. File should have
         * the following format:
         * x11 x12 ... x1d\n
         * x21 x22 ... x2d\n
         * ...
         * xn1 xn2 ... xnd
         *
         * Each entry should be a 'double' value.
         * Lines with '%' are ignored.
         *
         * @param file_path Path to the file.
         * @return The set of datapoints.
         */
        explicit Clustering(const std::string &file_path) {
            if (!CCSMSM::file_exists(file_path)) {
                std::cout << "File " << file_path << " was not found!\n";
                exit(EXIT_FAILURE);
            }

            std::ifstream file(file_path);
            std::string line;

            while (std::getline(file, line)) {
                while (line.back() == '\n' || line.back() == '\r') {
                    line.pop_back();
                }

                if (line[0] != '%') {
                    std::vector<std::string> temp = split(line, ' ');
                    std::vector<double> t;
                    for (auto &s: temp) { t.push_back(std::stod(s)); }
                    m_vecs.emplace_back(t);
                }
            }
            file.close();

            m_n_vec = m_vecs.size();
            m_dim = m_vecs[0].size();
        }

        /**
         * Returns the Euclidian distance matrix, but all entries are negative.
         * Optimizing this structure will solve the k-Medoid Clustering problem.
         *
         * @return The matrix.
         */
        inline Matrix<double> get_EuclidianClusteringMatrix() {
            Matrix<double> mtx(m_n_vec, m_n_vec);

            double sum_of_distance = 0.0;
            for (size_t i = 0; i < m_n_vec; ++i) {
                for (size_t j = i; j < m_n_vec; ++j) {

                    double distance = 0.0;
                    for (size_t d = 0; d < m_dim; ++d) {
                        distance += (m_vecs[i][d] - m_vecs[j][d]) * (m_vecs[i][d] - m_vecs[j][d]);
                    }
                    distance = sqrt(distance);
                    sum_of_distance += distance;

                    mtx.set(i, j, -distance); // negative
                    mtx.set(j, i, -distance); // negative
                }
            }

            mtx.set_special_scores(-2 * sum_of_distance, std::numeric_limits<double>::max());
            return mtx;
        }

        /**
         * Returns the Manhattan distance matrix, but all entries are negative.
         * Optimizing this structure will solve the k-Medoid Clustering problem
         * with Manhattan distances.
         *
         * @return The matrix.
         */
        inline Matrix<double> get_ManhattanClusteringMatrix() {
            Matrix<double> mtx(m_n_vec, m_n_vec);

            double sum_of_distance = 0.0;
            for (size_t i = 0; i < m_n_vec; ++i) {
                for (size_t j = i; j < m_n_vec; ++j) {

                    double distance = 0.0;
                    for (size_t d = 0; d < m_dim; ++d) {
                        distance += std::abs(m_vecs[i][d] - m_vecs[j][d]);
                    }
                    sum_of_distance += distance;

                    mtx.set(i, j, -distance); // negative
                    mtx.set(j, i, -distance); // negative
                }
            }

            mtx.set_special_scores(-2 * sum_of_distance, std::numeric_limits<double>::max());
            return mtx;
        }
    };
}

#endif //CCSMSM_CLUSTERING_H
