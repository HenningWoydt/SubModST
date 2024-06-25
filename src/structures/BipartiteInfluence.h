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

#ifndef CCSMSM_BIPARTITEINFLUENCE_H
#define CCSMSM_BIPARTITEINFLUENCE_H

#include <string>
#include <cmath>

namespace CCSMSM {
    template<typename SFType>
    class BipartiteInfluence final : public StructureInterface<SFType> {
    private:
        size_t m_n_items;
        size_t m_n_targets;

        // TODO: The layout of m_targets can further be optimized
        std::vector<std::vector<double>> m_targets;

        // structures to speed up score function evaluation
        size_t m_depth = 0;
        std::vector<std::vector<SFType>> m_temp_rows;
        std::vector<SFType> m_temp;

    public:
        /**
         * Reads an instance for Bipartite Influence.
         *
         * @param file_path Path to the file.
         * @return The set of datapoints.
         */
        explicit BipartiteInfluence(std::string &file_path) {
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
                    std::vector<std::string> temp = split(line, ',');
                    std::vector<double> t;
                    for (auto &s: temp) { t.push_back(std::stod(s)); }
                    m_targets.push_back(t);
                }
            }
            file.close();

            m_n_targets = m_targets.size();
            m_n_items = m_targets[0].size();

            for (size_t i = 0; i < m_n_targets; ++i) {
                for (size_t j = 0; j < m_n_items; ++j) {
                    m_targets[i][j] = std::log(1 - m_targets[i][j]);
                }
            }

            m_temp.resize(m_n_targets);

            StructureInterface<SFType>::set_special_scores(0, std::numeric_limits<SFType>::max());
        }

        void print() {
            for (auto &vec: m_targets) {
                ::CCSMSM::print(vec, vec.size());
            }
        }

        /**
         * Returns the number of rows in the matrix.
         *
         * @return Number of rows.
         */
        inline size_t get_n() override {
            return m_n_items;
        }

        inline SFType evaluate_empty_set() override {
            return StructureInterface<SFType>::m_empty_set_score;
        };

        inline SFType evaluate_1D(const std::vector<uint32_t> &s, const size_t s_size) override {
            SFType sum = 0;
            for (size_t i = 0; i < m_n_targets; ++i) {
                sum += std::exp(m_temp_rows[m_depth][i] + m_targets[i][s[s_size - 1]]);
            }
            return ((SFType) m_n_targets) - sum;
        };

        inline SFType evaluate_2D(const std::vector<uint32_t> &s, const size_t s_size) override {
            SFType sum = 0;
            for (size_t i = 0; i < m_n_targets; ++i) {
                SFType a = m_temp_rows[m_depth][i];
                SFType b = m_targets[i][s[s_size - 2]];
                SFType c = m_targets[i][s[s_size - 1]];
                sum += std::exp(a + b + c);
            }
            return ((SFType) m_n_targets) - sum;
        };

        inline SFType evaluate(const std::vector<uint32_t> &s, size_t s_size) override {
            // fill temp with 0s
            std::fill(m_temp.begin(), m_temp.end(), 0);

            // sum all values
            for (size_t i = 0; i < m_n_targets; ++i) {
                SFType temp_sum = 0;
                for (size_t j = 0; j < s_size; ++j) {
                    temp_sum += m_targets[i][s[j]];
                }
                m_temp[i] = temp_sum;
            }

            // take exponent of all sums
            SFType sum = 0;
            for (size_t i = 0; i < m_n_targets; ++i) {
                sum += std::exp(m_temp[i]);
            }

            return ((SFType) m_n_targets) - sum;
        };

        inline void initialize_helping_structures([[maybe_unused]] size_t k) override {
            m_depth = 0;
            m_temp_rows.resize(k + 1, std::vector<SFType>(m_n_targets, 0));
            m_temp.resize(m_n_targets);
        };

        inline void visit_new_depth([[maybe_unused]] const std::vector<uint32_t> &s, [[maybe_unused]]  size_t s_size) override {
            m_depth += 1;

            // element wise sum with the latest row
            std::copy(m_temp_rows[m_depth - 1].begin(), m_temp_rows[m_depth - 1].end(), m_temp_rows[m_depth].begin());
            for (size_t i = 0; i < m_n_targets; ++i) {
                m_temp_rows[m_depth][i] += m_targets[i][s[s_size - 1]];
            }
        };

        inline void return_from_last_depth() override {
            m_depth -= 1;
        };
    };
}

#endif //CCSMSM_BIPARTITEINFLUENCE_H
