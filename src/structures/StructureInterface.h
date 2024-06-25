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

#ifndef CCSMSM_STRUCTUREINTERFACE_H
#define CCSMSM_STRUCTUREINTERFACE_H

#include <limits>
#include <cstddef>
#include <vector>
#include <cstdint>

namespace CCSMSM {

    /**
     * Interface for structures, so that all score functions can be called in the
     * same way.
     *
     * @tparam SFType Return data type of the score function.
     */
    template<typename SFType>
    class StructureInterface {
    protected:
        SFType m_empty_set_score;
        SFType m_max_score;

    public:
        /**
         * Sets special score value that max be different for each score
         * function. If the max_score is not known, set it to the maximal value
         * for the datatype.
         *
         * @param empty_set_score Score of the empty set.
         * @param max_score Maximum reachable score.
         */
        virtual inline void set_special_scores(SFType empty_set_score, SFType max_score) final {
            m_empty_set_score = empty_set_score;
            m_max_score = max_score;
        }

        /**
         * Returns the maximum reachable score.
         *
         * @return The score.
         */
        virtual inline SFType get_max_score() final {
            return m_max_score;
        }

        /**
         * Returns n.
         *
         * @return n.
         */
        virtual inline size_t get_n() = 0;

        /**
         * Gives the score for the empty set.
         *
         * @return Score of empty set.
         */
        virtual inline SFType evaluate_empty_set() = 0;

        /**
         * Gives the score for the set that has only added one element.
         *
         * @return Score of the set.
         */
        virtual inline SFType evaluate_1D(const std::vector<uint32_t> &s, size_t s_size) = 0;

        /**
         * Gives the score for the set that has only added two elements.
         *
         * @return Score of the set.
         */
        virtual inline SFType evaluate_2D(const std::vector<uint32_t> &s, size_t s_size) = 0;

        /**
         * Evaluates the score function in the most general way. The function
         * should be able to accept any set of any size.
         *
         * @param s The set S.
         * @param s_size The size of set S.
         */
        virtual inline SFType evaluate(const std::vector<uint32_t> &s, size_t s_size) = 0;

        /**
         * Will initialize helping structures, that can help with score function
         * evaluation.
         *
         * @param k Number of depths, the search algorithm can explore.
         */
        virtual inline void initialize_helping_structures(size_t k) = 0;

        /**
         * Signals to the structure, that a new depth will be explored. Use this
         * function to keep the helping structures up to date.
         *
         * @param s The set S, that is present at the new depth.
         * @param s_size Size of the set S.
         */
        virtual inline void visit_new_depth(const std::vector<uint32_t> &s, size_t s_size) = 0;

        /**
         * Signals to the structure, that we returned from the last depth. Use this
         * function to keep the helping structures up to date.
         */
        virtual inline void return_from_last_depth() = 0;
    };
}

#endif //CCSMSM_STRUCTUREINTERFACE_H
