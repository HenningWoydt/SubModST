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

#ifndef CCSMSM_BFSOLVER_H
#define CCSMSM_BFSOLVER_H

#include <cstddef>
#include <vector>
#include <cstdint>

namespace CCSMSM {

    /**
     * Brute-Force solver.
     *
     * @tparam Structure Structure to optimize.
     * @tparam SFType Return data type of the score function.
     */
    template<typename Structure, typename SFType>
    class BFSolver {
    private:
        Structure &m_structure;
        size_t m_n;
        size_t m_k;
        size_t m_sf_evaluated;

        std::vector<uint32_t> m_best_s;
        SFType m_best_score;

        std::vector<uint32_t> m_s;

        bool m_collect_s = false;
        std::vector<std::vector<uint32_t>> m_all_s;

    public:
        /**
         * Constructor for the Brute-Force solver.
         *
         * @param structure The function to maximize.
         * @param k Desired solution size.
         * @param collect_s Whether to collect all found s (used for testing).
         */
        BFSolver(Structure &structure, size_t k, bool collect_s = false) : m_structure(structure) {
            m_n = structure.get_n();
            m_k = k;
            m_sf_evaluated = 0;

            // helping structures
            m_structure.initialize_helping_structures(m_k);

            m_best_s.resize(k);
            m_best_score = -std::numeric_limits<SFType>::max();

            m_s.resize(k);

            m_collect_s = collect_s;
        }

        /**
         * Searches for the set of size k with the greatest value.
         */
        inline void search() {
            // special case at the top of the tree
            for (uint32_t c = 0; c < m_n; ++c) {
                m_s[0] = c;
                recursive_search(1);
            }
        };

        /**
         * Returns the best found solution.
         *
         * @return Vector holding the solution.
         */
        inline std::vector<uint32_t> get_solution() {
            return m_best_s;
        }

        /**
         * Returns the score of the best found solution.
         *
         * @return Score of the best found solution.
         */
        inline SFType get_solution_score() {
            return m_best_score;
        }

        /**
         * Returns all found sets of size k. This function is used for testing
         * purposes.
         *
         * @return
         */
        inline std::vector<std::vector<uint32_t>> get_all_s() {
            return m_all_s;
        }

    private:
        /**
         * Recursively searches for a solution in a tree-like manner.
         *
         * @param depth Current depth of the search algorithm.
         */
        inline void recursive_search(size_t depth) {
            if (depth == m_k) {
                // if the size of S is k then evaluate
                if (m_collect_s) {
                    m_all_s.push_back(m_s);
                }

                m_sf_evaluated += 1;
                SFType score = m_structure.evaluate(m_s, depth);
                if (score > m_best_score) {
                    m_best_score = score;
                    std::copy(m_s.begin(), m_s.end(), m_best_s.begin());
                }
            } else {
                // determine which element to add to the set
                for (uint32_t c = m_s[depth - 1] + 1; c < m_n; ++c) {
                    m_s[depth] = c;
                    recursive_search(depth + 1);
                }
            }
        };
    };
}

#endif //CCSMSM_BFSOLVER_H
