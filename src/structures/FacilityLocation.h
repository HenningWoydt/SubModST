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

#ifndef SUBMODST_FACILITYLOCATION_H
#define SUBMODST_FACILITYLOCATION_H

#include <cstddef>
#include <vector>
#include <string>
#include <iterator>

#include "../util/utility.h"
#include "Matrix.h"

namespace SubModST {
    /**
     * Class for defining the Facility Location problem.
     */
    class FacilityLocation {
    private:
        size_t m_n_facilities; // number of facilities
        size_t m_n_customers; // number of customers
        std::vector<std::vector<double>> m_benefits;

    public:
        /**
         * Reads a set of benefits for locations from the specified file. File
         * should have the following format:
         * x11 x12 ... x1c\n
         * x21 x22 ... x2c\n
         * ...
         * xf1 xf2 ... xfc
         *
         * Each entry should be a non-negative 'double' value.
         * Lines with '%' are ignored.
         *
         * @param file_path Path to the file.
         * @return The set of datapoints.
         */
        explicit FacilityLocation(const std::string &file_path) {
            if (!SubModST::file_exists(file_path)) {
                std::cout << "File " << file_path << " was not found!\n";
                exit(EXIT_FAILURE);
            }

            std::ifstream file(file_path);
            std::string line;
            std::vector<uint32_t> edges;

            while (std::getline(file, line)) {
                while (line.back() == '\n' || line.back() == '\r') {
                    line.pop_back();
                }

                if (line[0] != '%') {
                    std::vector<std::string> temp = split(line, ' ');
                    std::vector<double> t;
                    for (auto &s: temp) { t.push_back(std::stod(s)); }
                    m_benefits.emplace_back(t);
                }
            }
            file.close();

            m_n_facilities = m_benefits[0].size();
            m_n_customers = m_benefits.size();
        }

        /**
         * Returns the matrix.
         *
         * @return The matrix.
         */
        inline Matrix<double> get_FacilityLocationMatrix() {
            Matrix<double> mtx(m_n_facilities, m_n_customers);

            for (size_t i = 0; i < m_n_customers; ++i) {
                for (size_t j = 0; j < m_n_facilities; ++j) {
                    mtx.set(i, j, m_benefits[i][j]);
                }
            }

            mtx.set_special_scores(0.0, std::numeric_limits<double>::max());
            return mtx;
        }

    };
}

#endif //SUBMODST_FACILITYLOCATION_H
