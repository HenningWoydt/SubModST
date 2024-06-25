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

#ifndef SUBMODST_MATRIX_H
#define SUBMODST_MATRIX_H

#include <numeric>

#include "StructureInterface.h"
#include "../util/utility.h"

namespace SubModST {
    /**
     * Matrix structure for general optimization.
     */
    template<typename SFType>
    class Matrix final : public StructureInterface<SFType> {
    private:
        // default matrix vars
        size_t m_n_row; // number of rows
        size_t m_n_col; // number of columns
        std::vector<SFType> m_mtx; // the matrix

        // structures to speed up score function evaluation
        size_t m_depth;
        std::vector<std::vector<SFType>> m_temp_rows;
        std::vector<SFType> m_temp;

    public:
        /**
         * Default constructor.
         *
         * @param n_row Number of columns.
         * @param n_col Number of rows.
         */
        Matrix(size_t n_row, size_t n_col) {
            // default matrix vars
            m_n_row = n_row;
            m_n_col = n_col;
            m_mtx.resize(m_n_row * m_n_col);

            // structures to speed up score function evaluation
            m_depth = 0;
        }

        /**
         * Sets the value at the specified location.
         *
         * @param row The row of the matrix.
         * @param col The column of the matrix.
         * @param val The value to be placed.
         */
        inline void set(size_t row, size_t col, SFType val) {
            m_mtx[row * m_n_col + col] = val;
        }

        /**
         * Returns the value at the specified location.
         *
         * @param row The row of the matrix.
         * @param col The column of the matrix.
         * @return The value at the location.
         */
        inline SFType get(size_t row, size_t col) {
            return m_mtx[row * m_n_col + col];
        }

        /**
         * Returns the number of rows in the matrix.
         *
         * @return Number of rows.
         */
        inline size_t get_n() override {
            return m_n_row;
        }

        /**
         * Prints the matrix to the standard output.
         */
        inline void print() {
            for (size_t j = 0; j < m_n_row; ++j) {
                std::cout << "[";
                for (size_t i = 0; i < m_n_col - 1; ++i) {
                    std::cout << m_mtx[(j * m_n_col) + i] << ", ";
                }
                std::cout << m_mtx[(j * m_n_col) + m_n_col - 1] << "]" << std::endl;
            }
        }

        inline SFType evaluate_empty_set() override {
            return StructureInterface<SFType>::m_empty_set_score;
        };

        inline SFType evaluate_1D(const std::vector<uint32_t> &s, const size_t s_size) override {
            SFType score = 0;
            for (size_t i = 0; i < m_n_col; ++i) {
                score += std::max(m_temp_rows[m_depth][i], m_mtx[s[s_size - 1] * m_n_col + i]);
            }
            return score;
        };

        inline SFType evaluate_2D(const std::vector<uint32_t> &s, const size_t s_size) override {
            SFType score = 0;
            for (size_t i = 0; i < m_n_col; ++i) {
                SFType a = m_temp_rows[m_depth][i];
                SFType b = m_mtx[s[s_size - 2] * m_n_col + i];
                SFType c = m_mtx[s[s_size - 1] * m_n_col + i];
                score += std::max(a, std::max(b, c));
            }
            return score;
        };

        inline SFType evaluate(const std::vector<uint32_t> &s, size_t s_size) override {
            // copy the first row into temp
            for (size_t i = 0; i < m_n_col; ++i) {
                m_temp[i] = m_mtx[s[0] * m_n_col + i];
            }

            // element wise max for all other rows and temp
            for (size_t j = 1; j < s_size; ++j) {
                for (size_t i = 0; i < m_n_col; ++i) {
                    m_temp[i] = std::max(m_temp[i], m_mtx[s[j] * m_n_col + i]);
                }
            }

            // sum temp
            return std::accumulate(m_temp.begin(), m_temp.end(), 0.0);
        };

        inline void initialize_helping_structures(size_t k) override {
            m_depth = 0;
            m_temp_rows.resize(k + 1, std::vector<SFType>(m_n_col, -std::numeric_limits<SFType>::max()));
            m_temp.resize(m_n_col);
        };

        inline void visit_new_depth(const std::vector<uint32_t> &s, size_t s_size) override {
            m_depth += 1;

            // element wise max with the latest row
            std::copy(m_temp_rows[m_depth - 1].begin(), m_temp_rows[m_depth - 1].end(), m_temp_rows[m_depth].begin());
            for (size_t i = 0; i < m_n_col; ++i) {
                m_temp_rows[m_depth][i] = std::max(m_temp_rows[m_depth - 1][i], m_mtx[s[s_size - 1] * m_n_col + i]);
            }
        };

        inline void return_from_last_depth() override {
            m_depth -= 1;
        };
    };
}

#endif //SUBMODST_MATRIX_H
