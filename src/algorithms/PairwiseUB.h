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

#ifndef CCSMSM_PAIRWISEUB_H
#define CCSMSM_PAIRWISEUB_H

#include <cstddef>
#include <vector>
#include <iostream>
#include <algorithm>

#include "Cache.h"
#include "../util/utility.h"
#include "../3rd_party/blossom5/PerfectMatching.h"

namespace CCSMSM {

    /**
     * Returns upper bounds based on 2D marginal gains.
     *
     * @tparam SFType Return data type of the score function.
     */
    template<typename SFType>
    class PairwiseUB {
    private:
        std::vector<Candidate<SFType>> m_candidates;
        std::vector<Pair<SFType>> m_pairs;

        bool m_candidates_is_sorted = false;
        bool m_pairs_is_sorted = false;

        // needed for dynamic programming
        std::vector<bool> cache_ready;
        std::vector<SFType> cache;

    public:
        /**
         * Default constructor.
         */
        PairwiseUB() = default;

        /**
         * Clears all stored data.
         */
        inline void clear() {
            m_candidates.clear();
            m_pairs.clear();
            m_candidates_is_sorted = false;
            m_pairs_is_sorted = false;

            cache_ready.clear();
        }

        /**
         * Adds one candidate to the algorithm.
         *
         * @param candidate The candidate.
         */
        inline void add_candidate(Candidate<SFType> &candidate) {
            m_candidates.push_back(candidate);
        }

        /**
         * Adds one pair to the algorithm.
         *
         * @param pair The pair.
         */
        inline void add_pair(Pair<SFType> &pair) {
            m_pairs.push_back(pair);
        }

        /**
         * Returns an upper-bound for k elements.
         *
         * @param k Number of elements.
         * @param alg_type Which algorithm to use.
         * @return The upper bound.
         */
        inline SFType get_upper_bound(size_t k, size_t alg_type) {
            if (k == 0) {
                return 0;
            } else if (k == 1) {
                SFType m = -std::numeric_limits<SFType>::max();
                for (auto &candidate: m_candidates) { m = std::max(m, candidate.gain); }
                return m;
            } else if (k == 2) {
                SFType m = -std::numeric_limits<SFType>::max();
                for (auto &pair: m_pairs) { m = std::max(m, pair.gain); }
                return m;
            }

            if (alg_type == 0) {
                return get_greedy_upper_bound(k);
            } else if (alg_type == 1) {
                return get_matching_upper_bound(k);
            } else if (alg_type == 2) {
                return get_dynamic_upper_bound(k);
            } else {
                std::cout << "Dont recognize algorithm type " << alg_type << " for PUB!" << std::endl;
                exit(EXIT_FAILURE);
            }
        }

    private:
        /**
         * Returns an upper bound using the greedy algorithm.
         *
         * @param k Number of elements.
         * @return The upper bound.
         */
        inline SFType get_greedy_upper_bound(size_t k) {
            if (!m_candidates_is_sorted || !m_pairs_is_sorted) {
                std::sort(m_candidates.begin(), m_candidates.end(), [](const Candidate<SFType> &x, const Candidate<SFType> &y) { return x.gain > y.gain; });
                std::sort(m_pairs.begin(), m_pairs.end(), [](const Pair<SFType> &x, const Pair<SFType> &y) { return x.gain > y.gain; });
                m_candidates_is_sorted = true;
                m_pairs_is_sorted = true;
            }

            size_t n_pairs = (k - (k & 1)) / 2;
            SFType sum = 0;
            for (size_t i = 0; i < n_pairs; ++i) {
                sum += m_pairs[i].gain;
            }

            if (k & 1) {
                // have to add one single candidate
                for (auto &candidate: m_candidates) {
                    // check if the candidate is in the first n_pairs pairs
                    bool found = false;
                    for (size_t i = 0; i < n_pairs; ++i) {
                        found |= m_pairs[i].a == candidate.c;
                        found |= m_pairs[i].b == candidate.c;
                    }

                    if (!found) {
                        sum += candidate.gain;
                        break;
                    }
                }
            }

            return sum;
        }

        /**
         * Returns an upper bound using the matching algorithm.
         *
         * @param k Number of elements.
         * @return The upper bound.
         */
        inline SFType get_matching_upper_bound(size_t k) {
            size_t n_pairs = (k - (k & 1)) / 2;

            size_t n_vertices = m_candidates.size();
            size_t n_edges = m_pairs.size();
            size_t n_aux_vertices = n_vertices - (n_pairs * 2); // auxiliary vertices
            size_t n_aux_edges = n_vertices * n_aux_vertices;

            PerfectMatching pm((int) (n_vertices + n_aux_vertices), (int) (n_edges + n_aux_edges));
            pm.options.fractional_jumpstart = true; // true is better than false
            pm.options.dual_greedy_update_option = 0; // 0, 1 better than 2
            pm.options.dual_LP_threshold = 0.00;
            pm.options.update_duals_before = false;
            pm.options.update_duals_after = false;
            pm.options.single_tree_threshold = 1.00;
            pm.options.verbose = false;

            // add the clique
            double max_weight = -std::numeric_limits<double>::max();
            for (size_t i = 0; i < n_edges; ++i) {
                pm.AddEdge(m_pairs[i].a, m_pairs[i].b, -(double) m_pairs[i].gain);
                max_weight = std::max(max_weight, (double) m_pairs[i].gain);
            }

            // add auxiliary vertices
            for (size_t i = 0; i < n_vertices; ++i) {
                for (size_t j = 0; j < n_aux_vertices; ++j) {
                    pm.AddEdge((int) i, (int) (n_vertices + j), -(max_weight + 1));
                }
            }

            // solve
            pm.Solve();

            // get the solution
            SFType sum = 0;
            for (size_t i = 0; i < n_vertices; ++i) {
                size_t j = pm.GetMatch((int) i);
                if (j >= n_vertices || j < i) { continue; }
                // pair i and j is chosen
                for (auto &pair: m_pairs) {
                    if (pair.a == i && pair.b == j) {
                        sum += pair.gain;
                        break;
                    }
                }
            }

            if (k & 1) {
                // have to add one single candidate
                for (auto &candidate: m_candidates) {
                    // check if the candidate is used
                    size_t j = pm.GetMatch((int) candidate.c);
                    if (j >= n_vertices) {
                        // candidate is not matched
                        sum += candidate.gain;
                        break;
                    }
                }
            }

            return sum;
        }

        /**
         * Initializes the dynamic array, if not already initialized.
         */
        inline void initialize_dynamic_array() {
            if (cache_ready.empty()) {
                size_t n = m_candidates.size();
                size_t c = ceil(n, 2);
                size_t cache_size = pow_2(n);
                cache_ready.resize(c);
                std::fill(cache_ready.begin(), cache_ready.end(), 0);

                cache.resize(cache_size);
                std::fill(cache.begin(), cache.end(), std::numeric_limits<SFType>::max());

                // Populate caches with pair data
                for (size_t i = 0; i < m_pairs.size(); ++i) {
                    size_t bitset = 0 | (1 << m_pairs[i].a) | (1 << m_pairs[i].b);
                    cache[bitset] = m_pairs[i].gain;
                }

                // Mark initialization as complete
                cache_ready[1] = true;
            }
        }

        /**
         * Updates the dynamic array for n_pais pairs.
         *
         * @param max_n_pairs Maximum number of pairs.
         */
        inline void update_dynamic_array(size_t max_n_pairs) {
            // check if already initialized and initialize lower depth
            for (size_t n_pairs = 2; n_pairs <= max_n_pairs; ++n_pairs) {
                if (cache_ready[n_pairs]) {
                    continue;
                }

                // update the array
                size_t n = m_candidates.size();
                size_t set_size = n_pairs * 2;
                size_t bitset = pow_2(set_size) - 1; // first correct bitset
                size_t last_bitset = (pow_2(set_size) - 1) << (n - set_size); // last correct bitset

                while (bitset <= last_bitset) {
                    if ((size_t) __builtin_popcount(bitset) == set_size) {
                        for (size_t i = 0; i < n; ++i) {
                            for (size_t j = i + 1; j < n; ++j) {
                                if (get_bit(bitset, i) && get_bit(bitset, j)) {
                                    // we have the correct number of bits set, and i and j are set to 1

                                    // split bitset into two bitsets
                                    size_t bitset1 = (1ULL << i) | (1ULL << j); // set i and j to true
                                    size_t bitset2 = bitset & ~(1ULL << i) & ~(1ULL << j); // set i and j to false

                                    // look for a new minimum
                                    if (cache[bitset1] + cache[bitset2] < cache[bitset]) {
                                        cache[bitset] = cache[bitset1] + cache[bitset2];
                                    }
                                }
                            }
                        }
                    }
                    bitset += 1;
                }

                cache_ready[n_pairs] = true;
            }
        }

        /**
         * Returns an upper bound using the dynamic algorithm.
         *
         * @param k Number of elements.
         * @return The upper bound.
         */
        inline SFType get_dynamic_upper_bound(size_t k) {
            if (!m_candidates_is_sorted) {
                std::sort(m_candidates.begin(), m_candidates.end(), [](const Candidate<SFType> &x, const Candidate<SFType> &y) { return x.gain > y.gain; });
                m_candidates_is_sorted = true;
            }

            size_t n_pairs = (k - (k & 1)) / 2;
            size_t set_size = n_pairs * 2;
            size_t n = m_candidates.size();

            initialize_dynamic_array();
            update_dynamic_array(n_pairs);

            SFType best_sum = -std::numeric_limits<SFType>::max();
            size_t best_bitset = 0;
            size_t bitset = pow_2(set_size) - 1; // first correct bitset
            size_t last_bitset = (pow_2(set_size) - 1) << (n - set_size); // last correct bitset

            // iterate over all ints with n bits and n_pairs*2 bits set to 1
            while (bitset <= last_bitset) {
                if ((size_t) __builtin_popcount(bitset) == set_size && cache[bitset] > best_sum) {
                    best_sum = cache[bitset];
                    best_bitset = bitset;
                }
                bitset++;
            }

            // have to add one single candidate
            if (k & 1) {
                for (auto &candidate: m_candidates) {
                    // check if the candidate is in the selected pairs
                    if (!get_bit(best_bitset, candidate.c)) {
                        best_sum += candidate.gain;
                        break;
                    }
                }
            }

            return best_sum;
        }
    };
}

#endif //CCSMSM_PAIRWISEUB_H
