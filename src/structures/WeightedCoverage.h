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

#ifndef SUBMODST_WEIGHTEDCOVERAGE_H
#define SUBMODST_WEIGHTEDCOVERAGE_H

#include <string>
#include <vector>
#include <iostream>

#include "../util/utility.h"
#include "Matrix.h"

namespace SubModST {
    class WeightedCoverage {
    private:
        size_t m_n_items;
        size_t m_n_sensors;

        std::vector<double> m_weights;
        std::vector<std::vector<int>> m_active;

    public:
        /**
         * Reads a file containing data for weighted coverage. The first line
         * contains n weights for the n items. All other lines represent one
         * sensor each. An entry in each line ist set to 1 if the sensor
         * observes the item and 0 if not. All values are separated by a ','.
         *
         * @param file_path Path to the file.
         */
        explicit WeightedCoverage(std::string &file_path) {
            if (!SubModST::file_exists(file_path)) {
                std::cout << "File " << file_path << " was not found!\n";
                exit(EXIT_FAILURE);
            }

            std::ifstream file(file_path);
            std::string line;
            bool weights_read = false;

            while (std::getline(file, line)) {
                while (line.back() == '\n' || line.back() == '\r') {
                    line.pop_back();
                }

                if (line[0] == '%') {
                    continue;
                }

                if (!weights_read) {
                    std::vector<std::string> weights_str = split(line, ',');
                    for (std::string &s: weights_str) {
                        m_weights.push_back(std::stod(s));
                    }
                    weights_read = true;
                } else {
                    m_active.emplace_back();
                    std::vector<std::string> active_str = split(line, ',');
                    for (std::string &s: active_str) {
                        m_active.back().push_back((int) std::stod(s));
                    }
                }
            }
            file.close();

            m_n_items = m_active.size();
            m_n_sensors = m_active[0].size();
        }


        /**
         * Returns the matrix.
         *
         * @return
         */
        inline Matrix<double> get_WeightedCoverageMatrix() const {
            Matrix<double> mtx(m_n_sensors, m_n_items);

            for (size_t i = 0; i < m_n_sensors; ++i) {
                for (size_t j = 0; j < m_n_items; ++j) {
                    mtx.set(i, j, m_weights[j] * m_active[j][i]);
                }
            }

            mtx.set_special_scores(0, sum(m_weights));
            return mtx;
        }

    };
}

#endif //SUBMODST_WEIGHTEDCOVERAGE_H
