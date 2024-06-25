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

#ifndef SUBMODST_CANDIDATEMANAGER_H
#define SUBMODST_CANDIDATEMANAGER_H

#include <cstdint>
#include <vector>
#include <algorithm>

#include "../util/utility.h"

namespace SubModST {
    /**
     * Represents a candidate with score improvement and accuracy information.
     *
     * @tparam SFType Return data type of the score function.
     */
    template<typename SFType>
    struct Candidate {
        uint32_t c;
        SFType gain;
        bool acc;
    };

    template<typename T>
    std::ostream &operator<<(std::ostream &os,
                             const Candidate<T> &candidate) {
        os << "(" << candidate.c << ", " << candidate.gain << ", " << (int) candidate.acc << ")";
        return os;
    }

    /**
     * Manages a collection of c entries for optimization.
     *
     * @tparam SFType Return data type of the score function.
     */
    template<typename SFType>
    class CandidateManager {
    private:
        std::vector<Candidate<SFType>> m_candidates;
        std::vector<SFType> m_csum_si;

        size_t m_s_idx; // start index
        size_t m_e_idx; // end index

    public:
        /**
         * Constructor.
         *
         * @param n Number of candidates.
         */
        explicit CandidateManager(size_t n) {
            m_candidates.resize(n);
            m_csum_si.resize(n + 1);

            m_s_idx = 0;
            m_e_idx = 0;
        };

        /**
         * Sort the candidates by score improvement in descending order.
         */
        inline void sort() {
            auto first = m_candidates.begin();
            std::advance(first, m_s_idx);

            auto last = m_candidates.begin();
            std::advance(last, m_e_idx);

            std::sort(first, last, [](const Candidate<SFType> &a, const Candidate<SFType> &b) {
                return a.gain > b.gain;
            });
        }

        /**
         * Calculate the cumulative sum of score improvements for candidates.
         */
        inline void calc_csum() {
            SFType cumulative_sum = 0;
            m_csum_si[0] = 0; // Initialize the first element.

            for (size_t i = 0; i < m_e_idx; ++i) {
                cumulative_sum += m_candidates[i].gain;
                m_csum_si[i + 1] = cumulative_sum;
            }
        }

        /**
         * Calculate the sum of score improvements for a range of candidates.
         *
         * @param n The number of candidates to include in the sum.
         * @return The sum of score improvements for the specified range of candidates.
         */
        inline SFType get_partial_sum(size_t n) const { return m_csum_si[m_s_idx + n] - m_csum_si[m_s_idx]; }

        /** vector like functions */
        inline Candidate<SFType> operator[](size_t i) const { return m_candidates[i]; }

        inline Candidate<SFType> &operator[](size_t i) { return m_candidates[i]; }

        inline Candidate<SFType> back() const { return m_candidates[m_e_idx - 1]; }

        inline Candidate<SFType> &back() { return m_candidates[m_e_idx - 1]; }

        inline Candidate<SFType> front() const { return m_candidates[m_s_idx]; }

        inline Candidate<SFType> &front() { return m_candidates[m_s_idx]; }

        inline void push_back(Candidate<SFType> entry) {
            m_candidates[m_e_idx] = entry;
            m_e_idx += 1;
        }

        inline void pop_back() { m_e_idx -= 1; }

        inline void pop_front() { m_s_idx += 1; }

        inline size_t size() const { return m_e_idx - m_s_idx; }

        inline size_t get_start() const { return m_s_idx; }

        inline size_t get_end() const { return m_e_idx; }

        inline void clear() {
            m_s_idx = 0;
            m_e_idx = 0;
        }

        /**
         * Prints the CandidateManager to the standard output stream.
         */
        inline void print() {
            if (m_e_idx == 0) {
                std::cout << "C  : []" << std::endl;
                std::cout << "SI : []" << std::endl;
                std::cout << "ACC: []" << std::endl;
                return;
            }


            std::vector<std::string> c_s;
            std::vector<std::string> si_s;
            std::vector<std::string> acc_s;

            // convert all to strings
            size_t max_str_size = 0;
            for (size_t i = 0; i < m_e_idx; ++i) {
                std::string str_c = std::to_string(m_candidates[i].c);
                std::string str_si = std::to_string(m_candidates[i].gain);
                std::string str_acc = std::to_string(m_candidates[i].acc);

                max_str_size = std::max({max_str_size, str_c.size(), str_si.size(), str_acc.size()});

                c_s.push_back(str_c);
                si_s.push_back(str_si);
                acc_s.push_back(str_acc);
            }

            // pad all strings
            for (size_t i = 0; i < m_e_idx; ++i) {
                while (c_s[i].size() <= max_str_size) c_s[i].push_back(' ');
                while (si_s[i].size() <= max_str_size) si_s[i].push_back(' ');
                while (acc_s[i].size() <= max_str_size) acc_s[i].push_back(' ');
            }

            // print candidates
            std::cout << "C  : [";
            for (size_t i = 0; i < m_e_idx - 1; ++i) {
                if (i >= m_s_idx) {
                    std::cout << "\033[1;32m" << c_s[i] << "\033[0m" << ", ";
                } else {
                    std::cout << c_s[i] << ", ";
                }
            }
            if (m_e_idx - 1 >= m_s_idx) {
                std::cout << "\033[1;32m" << c_s[m_e_idx - 1] << "\033[0m" << "]" << std::endl;
            } else {
                std::cout << c_s[m_e_idx - 1] << "]" << std::endl;
            }

            // print score improvements
            std::cout << "SI : [";
            for (size_t i = 0; i < m_e_idx - 1; ++i) {
                if (i >= m_s_idx) {
                    std::cout << "\033[1;32m" << si_s[i] << "\033[0m" << ", ";
                } else {
                    std::cout << si_s[i] << ", ";
                }
            }
            if (m_e_idx - 1 >= m_s_idx) {
                std::cout << "\033[1;32m" << si_s[m_e_idx - 1] << "\033[0m" << "]" << std::endl;
            } else {
                std::cout << si_s[m_e_idx - 1] << "]" << std::endl;
            }

            // print accuracy
            std::cout << "ACC: [";
            for (size_t i = 0; i < m_e_idx - 1; ++i) {
                if (i >= m_s_idx) {
                    std::cout << "\033[1;32m" << acc_s[i] << "\033[0m" << ", ";
                } else {
                    std::cout << acc_s[i] << ", ";
                }
            }
            if (m_e_idx - 1 >= m_s_idx) {
                std::cout << "\033[1;32m" << acc_s[m_e_idx - 1] << "\033[0m" << "]" << std::endl;
            } else {
                std::cout << acc_s[m_e_idx - 1] << "]" << std::endl;
            }
        }

    };
}

#endif //SUBMODST_CANDIDATEMANAGER_H
