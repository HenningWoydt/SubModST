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

#ifndef SUBMODST_SOLVER_H
#define SUBMODST_SOLVER_H

#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>
#include <chrono>
#include <numeric>
#include <cmath>

#include "CandidateManager.h"
#include "Heap.h"
#include "Cache.h"
#include "PairwiseUB.h"
#include "../util/utility.h"
#include "../util/AlgorithmConfiguration.h"
#include "StatisticProfiler.h"

namespace SubModST {

    /**
     * Algorithm to find the set S with k elements that maximizes a score
     * function.
     *
     * @tparam Structure Structure to optimize.
     * @tparam SFType Return data type of the score function.
     */
    template<typename Structure, typename SFType>
    class Solver {
    private:
        // default variables
        Structure &m_structure;
        size_t m_n;
        size_t m_k;

        // config of the algorithm
        AlgorithmConfiguration m_ac;

        // vectors to store s
        std::vector<uint32_t> m_s;
        std::vector<uint32_t> m_initial_s;
        std::vector<uint32_t> m_best_s;
        std::vector<uint32_t> m_best_greedy_s;
        size_t m_s_size;

        // scores
        SFType m_initial_score;
        SFType m_best_score;
        SFType m_best_greedy_score;

        // time points
        std::chrono::steady_clock::time_point m_start_point;
        std::chrono::steady_clock::time_point m_end_point;
        std::chrono::steady_clock::time_point m_exhaustive_start_point;
        std::chrono::steady_clock::time_point m_exhaustive_end_point;
        std::chrono::steady_clock::time_point m_greedy_start_point;
        std::chrono::steady_clock::time_point m_greedy_end_point;

        // information for timeout
        size_t m_tried_calls = 0;
        double m_time_limit = std::numeric_limits<double>::max();
        bool m_time_exceeded = false;

        /* -- internal variables -- */
        // information for each depth
        std::vector<SFType> m_scores;
        std::vector<CandidateManager<SFType>> m_c_managers;
        std::vector<Heap<SFType>> m_heaps;

        // heuristics
        std::vector<Cache<SFType>> m_caches;
        std::vector<PairwiseUB<SFType>> m_PUBs;

        // information for greedy
        std::vector<Candidate<SFType>> m_gains;

        // temporary variables
        std::vector<size_t> help_set;

        /* -- statistic variables -- */
        StatisticProfiler<SFType> stat_profiler;


    public:
        /**
         * Constructor for the Brute-Force solver.
         *
         * @param structure The function to maximize.
         * @param k Desired solution size.
         */
        Solver(Structure &structure,
               size_t k,
               AlgorithmConfiguration ac) : m_structure(structure) {
            // default variables
            m_n = structure.get_n();
            m_k = k;

            // helping structures
            m_structure.initialize_helping_structures(m_k);

            // algorithm configuration
            m_ac = std::move(ac);

            // vectors to store s
            m_s.resize(m_k);
            m_initial_s.resize(m_k);
            m_best_s.resize(m_k);
            m_best_greedy_s.resize(m_k);
            m_s_size = 0;

            if (!m_ac.initial_solution_file_path.empty()) {
                m_initial_s = read_solution(m_ac.initial_solution_file_path);
            } else {
                std::iota(m_initial_s.begin(), m_initial_s.end(), 0);
            }

            // scores
            m_initial_score = -std::numeric_limits<SFType>::max();
            m_best_score = -std::numeric_limits<SFType>::max();
            m_best_greedy_score = -std::numeric_limits<SFType>::max();

            // time points
            m_time_limit = m_ac.time_limit;
            m_exhaustive_start_point = get_time_point();
            m_exhaustive_end_point = m_exhaustive_start_point;
            m_greedy_start_point = m_exhaustive_start_point;
            m_greedy_end_point = m_exhaustive_start_point;

            // information for each depth
            m_scores.resize(m_k);
            m_c_managers.resize(m_k, CandidateManager<SFType>(m_n));
            m_heaps.resize(m_k, Heap<SFType>());

            // heuristics
            m_caches.resize(m_k, Cache<SFType>());
            m_PUBs.resize(m_k, PairwiseUB<SFType>());

            // information for greedy
            m_gains.resize(m_n);

            /* -- statistic variables -- */
            stat_profiler = StatisticProfiler<SFType>(m_n, m_k);
        }

        /**
         * Searches for the set of size k with the greatest value.
         */
        inline void search() {
            m_start_point = get_time_point();
            stat_profiler.set_start_point(m_start_point);

            // initial solution
            m_initial_score = m_structure.evaluate(m_initial_s, m_k);
            update_best(m_initial_s, m_initial_score);

            // special case if k == 1
            if (m_k == 1) {
                for (uint32_t c = 0; c < m_n; ++c) {
                    m_s[0] = c;
                    SFType score = m_structure.evaluate(m_s, 1);
                    if (score > m_best_score) { update_best(m_s, score); }
                }
                m_end_point = get_time_point();
                return;
            }

            // special case if k == n
            if (m_k == m_n) {
                std::iota(m_best_s.begin(), m_best_s.end(), 0);
                m_best_score = m_structure.evaluate(m_best_s, m_k);
                m_end_point = get_time_point();
                return;
            }

            m_greedy_start_point = get_time_point();
            search_greedily();
            m_greedy_end_point = get_time_point();

            m_exhaustive_start_point = get_time_point();
            search_exhaustive_iterative();
            m_exhaustive_end_point = get_time_point();

            m_end_point = get_time_point();
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
         * Returns the initial solution.
         *
         * @return Vector holding the initial solution.
         */
        inline std::vector<uint32_t> get_initial_solution() {
            return m_initial_s;
        }

        /**
         * Returns the score of the initial solution.
         *
         * @return Score of the initial solution.
         */
        inline SFType get_initial_solution_score() {
            return m_best_score;
        }

        /**
         * Returns the greedy solution.
         *
         * @return Vector holding the greedy solution.
         */
        inline std::vector<uint32_t> get_greedy_solution() {
            return m_best_greedy_s;
        }

        /**
         * Returns the score of the greedy solution.
         *
         * @return Score of the greedy solution.
         */
        inline SFType get_greedy_solution_score() {
            return m_best_score;
        }

        /**
         * Returns the needed time to solve the instance.
         *
         * @return The time in seconds.
         */
        inline double get_time_needed() {
            return get_elapsed_seconds(m_start_point, m_end_point);
        }

        /**
         * Writes the output to the specified file in JSON format.
         *
         * @param file_path Path to the file.
         */
        inline void write_output(std::string &file_path) {
            std::stringstream ss;
            ss << "{";

            // default variables
            ss << "\"n\": " << m_n << ",\n";
            ss << "\"k\": " << m_k << ",\n";

            // config of the algorithm
            ss << "\"ac\": " << m_ac.to_JSON() << ",\n";

            // vectors to store s
            ss << "\"initial-s\": " << to_JSON(m_initial_s) << ",\n";
            ss << "\"greedy-s\": " << to_JSON(m_best_greedy_s) << ",\n";
            ss << "\"best-s\": " << to_JSON(m_best_s) << ",\n";

            // scores
            ss << "\"initial-score\": " << m_initial_score << ",\n";
            ss << "\"greedy-score\": " << m_best_greedy_score << ",\n";
            ss << "\"best-score\": " << m_best_score << ",\n";

            // time points
            ss << "\"total-time\": " << get_elapsed_seconds(m_start_point, m_end_point) << ",\n";
            ss << "\"greedy-time\": " << get_elapsed_seconds(m_greedy_start_point, m_greedy_end_point) << ",\n";
            ss << "\"exhaustive-time\": " << get_elapsed_seconds(m_exhaustive_start_point, m_exhaustive_end_point) << ",\n";

            // information for timeout
            ss << "\"time-limit-exceeded\": " << m_time_exceeded << ",\n";

            // statistic information
            ss << "\"statistics\": " << stat_profiler.to_JSON() << "\n";

            ss << "}";
            std::ofstream file(file_path);
            file << ss.str();
        }

    private:
        /**
         * Searches for a solution in a greedy-like manner.
         */
        inline void search_greedily() {
            SFType score = m_structure.evaluate_empty_set();

            uint32_t best_c = 0;
            SFType best_gain = -std::numeric_limits<SFType>::max();
            size_t best_idx = 0;

            // first pass calculate all marginal gains
            for (uint32_t c = 0; c < m_n; ++c) {
                m_s[0] = c;
                SFType gain = m_structure.evaluate(m_s, 1) - score;
                m_gains[c] = {c, gain, 1};

                if (gain > best_gain) {
                    best_c = c;
                    best_gain = gain;
                }
            }
            // insert best, update score, sort candidates
            m_s[0] = best_c;
            score += best_gain;
            std::swap(m_gains[best_c], m_gains.back());
            m_gains.pop_back();
            std::sort(m_gains.begin(), m_gains.end(), [](const Candidate<SFType> &a, const Candidate<SFType> &b) {
                return a.gain > b.gain;
            });

            // all other passes, only update if necessary
            for (size_t k = 1; k < m_k; ++k) {
                best_gain = -std::numeric_limits<SFType>::max();
                for (size_t i = 0; i < m_n - k; ++i) {
                    // break if found score is already greater
                    if (best_gain >= m_gains[i].gain) {
                        break;
                    }

                    // check if updated score has better gain
                    m_s[k] = m_gains[i].c;
                    SFType gain = m_structure.evaluate(m_s, k + 1) - score;
                    m_gains[i].gain = gain;

                    if (gain > best_gain) {
                        best_c = m_gains[i].c;
                        best_gain = gain;
                        best_idx = i;
                    }
                }

                // insert best, update score, sort candidates
                m_s[k] = best_c;
                score += best_gain;
                std::swap(m_gains[best_idx], m_gains.back());
                m_gains.pop_back();
                std::sort(m_gains.begin(), m_gains.end(), [](const Candidate<SFType> &a, const Candidate<SFType> &b) {
                    return a.gain > b.gain;
                });
            }

            // insert best greedy score
            std::copy(m_s.begin(), m_s.end(), m_best_greedy_s.begin());
            m_best_greedy_score = score;

            // check if new best solution
            if (m_best_greedy_score > m_best_score) {
                update_best(m_best_greedy_s, m_best_greedy_score);
            }
        }

#define DEPTH_STAY 0
#define DEPTH_UP 1

        /**
         * Exhaustively searches for the best solution.
         */
        inline void search_exhaustive_iterative() {
            m_scores[0] = m_structure.evaluate_empty_set();
            execute_DCO_depth_0();
            stat_profiler.visited_depth(0);

            m_s_size = 0;
            size_t depth = 0;
            size_t depth_action = DEPTH_STAY;

            while (depth != 0 || depth_action != DEPTH_UP) {
                if (depth_action == DEPTH_STAY) {

                    if (m_best_score >= m_structure.get_max_score()) {
                        // we have found a set with more than the desired score
                        depth_action = DEPTH_UP;
                        continue;
                    }

                    // depth, s_size, and r assumed to have the correct value
                    if (m_c_managers[depth].size() < m_k - m_s_size) {
                        // if not enough candidates go up
                        depth_action = DEPTH_UP;
                        continue;
                    }

                    SFType sub_bound = calculate_SUB(depth);
                    if (sub_bound <= m_best_score) {
                        depth_action = DEPTH_UP;
                        continue;
                    }

                    execute_CR(depth);
                    if (m_c_managers[depth].size() == m_k - m_s_size) {
                        execute_RCR(depth);
                        depth_action = DEPTH_UP;
                        continue;
                    }
                    if (m_c_managers[depth].size() < m_k - m_s_size) {
                        depth_action = DEPTH_UP;
                        continue;
                    }

                    stat_profiler.start_track_PUB();
                    SFType pub_bound = calculate_PUB(depth);
                    stat_profiler.end_track_PUB(pub_bound <= m_best_score);
                    if (pub_bound <= m_best_score) {
                        depth_action = DEPTH_UP;
                        continue;
                    }

                    // no heuristic pruned, so go down the tree
                    m_s[m_s_size] = m_c_managers[depth].front().c;
                    m_scores[depth + 1] = m_c_managers[depth].front().acc ? m_scores[depth] + m_c_managers[depth].front().gain : m_structure.evaluate_1D(m_s, m_s_size + 1);
                    m_c_managers[depth].pop_front();
                    m_structure.visit_new_depth(m_s, m_s_size + 1);

                    // --- DEPTH DOWN ---
                    depth += 1;
                    m_s_size += 1;

                    m_caches[depth].clear();

                    if (m_k - m_s_size == 1) {
                        execute_OCR(depth);
                        depth_action = DEPTH_UP;
                        continue;
                    }

                    if (m_c_managers[depth - 1].size() == m_k - m_s_size) {
                        execute_RCR(depth - 1);
                        depth_action = DEPTH_UP;
                        continue;
                    }

                    if (m_c_managers[depth - 1].size() == (m_k - m_s_size) + 1) {
                        execute_R1CR(depth);
                        depth_action = DEPTH_UP;
                        continue;
                    }

                    if (m_c_managers[depth - 1].size() <= 5 && m_k - m_s_size <= 5) {
                        execute_small_BF(depth);
                        depth_action = DEPTH_UP;
                        continue;
                    }

                    execute_DCO(depth);
                    if (m_heaps[depth].is_abort_early()) {
                        depth_action = DEPTH_UP;
                        continue;
                    }

                    depth_action = DEPTH_STAY;
                    continue;
                } else {
                    // we will go one depth up
                    // depth, s_size and r are assumed to be in wrong state

                    depth -= 1;
                    m_s_size -= 1;
                    m_structure.return_from_last_depth();

                    if (has_time_exceeded()) {
                        return;
                    }

                    depth_action = DEPTH_STAY;
                    continue;
                }
            }
        }

        /**
         * Dynamic Candidate Ordering at depth 0. It will calculate the score
         * improvement for all candidates and sort them descending.
         */
        inline void execute_DCO_depth_0() {
            CandidateManager<SFType> &c_manger = m_c_managers[0];
            c_manger.clear();

            for (uint32_t c = 0; c < m_n; ++c) {
                m_s[0] = c;
                SFType si = m_structure.evaluate_1D(m_s, 1) - m_scores[0];
                c_manger.push_back({c, si, true});
            }

            c_manger.sort();
            c_manger.calc_csum();
        }

        /**
         * Executes Dynamic Candidate Ordering on the given depth.
         *
         * @param depth The depth.
         */
        inline void execute_DCO(size_t depth) {
            const CandidateManager<SFType> &pc_manager = m_c_managers[depth - 1];
            const size_t p_s = pc_manager.get_start();
            const size_t p_e = pc_manager.get_end();

            CandidateManager<SFType> &c_manager = m_c_managers[depth];
            c_manager.clear();

            size_t r = m_k - m_s_size;
            SFType r_score = m_best_score - m_scores[depth];
            double r_score_avg = ((double) r_score / (double) r);

            // initialize the heap
            Heap<SFType> &heap = m_heaps[depth];
            if (m_ac.enabled_DCO_heap) {
                heap.initialize(r);
                for (size_t i = p_s; i < p_s + r; ++i) {
                    uint32_t c = pc_manager[i].c;
                    SFType si = pc_manager[i].gain;
                    heap.push({c, si});
                }
            }

            // default: update all
            const double score_threshold = m_ac.determine_LE_score_threshold(r_score_avg);
            const size_t rank_threshold = m_ac.determine_LE_rank_threshold(m_n, m_k, p_e - p_s, r);
            auto update_scheme = [&](size_t curr_rank, double curr_score) {
                if (m_ac.LE_mode == 0) {
                    return curr_score >= score_threshold || curr_rank <= rank_threshold;
                } else if (m_ac.LE_mode == 1) {
                    return curr_score >= score_threshold && curr_rank <= rank_threshold;
                }

                return true;
            };


            // process candidates
            for (size_t i = p_s; i < p_e; ++i) {
                // get old values
                uint32_t c = pc_manager[i].c;
                SFType si = pc_manager[i].gain;
                bool acc = false;

                // check if we can abort early
                if (m_ac.enabled_DCO_heap && heap.get_heap_sum() <= r_score && heap.get_heap_min() >= si) {
                    heap.mark_abort_early();
                    return;
                }

                // check for update
                bool update = update_scheme(i - p_s, si);
                if (update) {
                    // update candidate
                    m_s[m_s_size] = c;
                    si = m_structure.evaluate_1D(m_s, m_s_size + 1) - m_scores[depth];
                    acc = true;
                }

                c_manager.push_back({c, si, acc});
                if (m_ac.enabled_DCO_heap) {
                    if (i < p_s + r) {
                        heap.update({c, si});
                    } else {
                        heap.push({c, si});
                    }
                }
            }

            // sort and calculate cumulative sum
            c_manager.sort();
            c_manager.calc_csum();
        }

        /**
         * Calculates 'Simple Upper Bound'.
         *
         * @param depth Current depth of the algorithm.
         * @return The upper bound.
         */
        inline SFType calculate_SUB(size_t depth) {
            if (!m_ac.enabled_SUB) {
                return std::numeric_limits<SFType>::max();
            }
            const size_t r = m_k - m_s_size;
            return m_scores[depth] + m_c_managers[depth].get_partial_sum(r);
        }

        /**
         * Calculates 'Pairwise Upper Bound'.
         *
         * @param depth Current depth of the algorithm.
         * @return The upper bound.
         */
        inline SFType calculate_PUB(size_t depth) {
            if (!m_ac.enabled_PUB) {
                return std::numeric_limits<SFType>::max();
            }
            m_PUBs[depth].clear();

            size_t r = m_k - m_s_size;
            size_t n_remaining = m_c_managers[depth].size();

            size_t l = 0;
            if (m_ac.PUB_l_type == 0) {
                l = m_ac.PUB_l;
            } else if (m_ac.PUB_l_type == 1) {
                l = r;
            } else if (m_ac.PUB_l_type == 2) {
                l = (size_t) std::sqrt(n_remaining);
            }
            l = std::clamp(l, (size_t) 2, n_remaining);

            // fill the gains into PUB
            size_t offset = m_c_managers[depth].get_start();
            Candidate<SFType> c = {0, 0, 0};
            for (uint32_t i = 0; i < l; ++i) {
                c.c = i;
                c.gain = m_c_managers[depth][offset + i].gain;
                c.acc = m_c_managers[depth][offset + i].acc;

                m_PUBs[depth].add_candidate(c);
            }

            // fill the 2D gains into PUB
            Pair<SFType> p = {0, 0, 0};
            Pair<SFType> p_cache = {0, 0, 0};
            for (size_t i = 0; i < l; ++i) {
                for (size_t j = i + 1; j < l; ++j) {
                    p_cache.a = m_c_managers[depth][offset + i].c;
                    p_cache.b = m_c_managers[depth][offset + j].c;

                    if (!m_caches[depth].retrieve(p_cache)) {
                        m_s[m_s_size] = p_cache.a;
                        m_s[m_s_size + 1] = p_cache.b;
                        p_cache.gain = m_structure.evaluate_2D(m_s, m_s_size + 2) - m_scores[depth];
                        m_caches[depth].insert(p_cache);
                    }

                    p.a = i;
                    p.b = j;
                    p.gain = p_cache.gain;

                    m_PUBs[depth].add_pair(p);
                }
            }

            // calculate the bounds
            SFType upper_bound = -std::numeric_limits<SFType>::max();

            for (size_t i = 0; i <= r; ++i) {
                size_t n_candidates_p1 = i;
                size_t n_candidates_p2 = r - i;

                if (n_candidates_p1 <= l && n_candidates_p2 <= n_remaining - l) {
                    SFType bound_p1 = m_PUBs[depth].get_upper_bound(n_candidates_p1, m_ac.PUB_algo_type);
                    SFType bound_p2 = m_c_managers[depth].get_partial_sum(n_candidates_p2);
                    SFType bound = m_scores[depth] + bound_p1 + bound_p2;

                    upper_bound = std::max(upper_bound, bound);
                    if (upper_bound > m_best_score) {
                        return upper_bound;
                    }
                }
            }

            return upper_bound;
        }

        /**
         * Executes 'Candidate Reduction'.
         *
         * @param depth Current depth of the algorithm.
         */
        inline void execute_CR(size_t depth) {
            if (!m_ac.enabled_CR) {
                return;
            }

            CandidateManager<SFType> &c_manager = m_c_managers[depth];
            const size_t r = m_k - m_s_size;

            SFType min_req_score = m_best_score - (m_scores[depth] + c_manager.get_partial_sum(r - 1));

            while (c_manager.back().gain <= min_req_score && c_manager.size() > r) {
                c_manager.pop_back();
            }
        }

        /**
         * Executes the special 'One Candidate Remaining' routine. If only one
         * more candidate is required in S, it is more efficient to use this
         * function, instead of searching the tree.
         *
         * @param depth Current depth of the algorithm.
         */
        inline void execute_OCR(const size_t depth) {
            const CandidateManager<SFType> &pc_manager = m_c_managers[depth - 1];
            const size_t p_e = pc_manager.get_end();
            const size_t p_s = pc_manager.get_start();

            SFType remaining_score = m_best_score - m_scores[depth];

            for (size_t i = p_s; i < p_e; ++i) {
                if (pc_manager[i].gain >= remaining_score) {
                    m_s[m_s_size] = pc_manager[i].c;
                    SFType new_score = m_structure.evaluate_1D(m_s, m_s_size + 1);

                    if (new_score > m_best_score) {
                        update_best(m_s, new_score);
                        remaining_score = m_best_score - m_scores[depth];
                    }
                } else {
                    break;
                }
            }
        };

        /**
         * Executes the special 'R Candidates Remaining' routine. If r
         * candidates are available and exactly r candidates are still required
         * in S, then we can just insert them and test the set.
         *
         * @param depth Current depth of the algorithm.
         */
        inline void execute_RCR(const size_t depth) {
            const CandidateManager<SFType> &c_manager = m_c_managers[depth];
            const size_t start = c_manager.get_start();
            const size_t r = m_k - m_s_size;

            for (size_t i = 0; i < r; ++i) {
                m_s[m_s_size + i] = c_manager[start + i].c;
            }

            SFType score = m_structure.evaluate(m_s, m_k);

            if (score > m_best_score) { update_best(m_s, score); }
        }

        /**
         * Executes the special 'R+1 Candidates Remaining' routine. If r+1
         * candidates are available and exactly r candidates are still required
         * in S, then this routine is more efficient then searching the tree.
         *
         * @param depth Current depth of the algorithm.
         */
        inline void execute_R1CR(const size_t depth) {
            const CandidateManager<SFType> &pc_manager = m_c_managers[depth - 1];
            const size_t p_s = pc_manager.get_start();
            const size_t r = m_k - m_s_size;

            for (size_t i = 0; i < r; ++i) {
                m_s[m_s_size + i] = pc_manager[p_s + i].c;
            }
            SFType new_score = m_structure.evaluate(m_s, m_s_size + r);
            if (new_score > m_best_score) { update_best(m_s, new_score); }

            for (size_t i = 0; i < r; ++i) {
                m_s[m_s_size + r - 1 - i] = pc_manager[p_s + r - i].c;
                new_score = m_structure.evaluate(m_s, m_s_size + r);
                if (new_score > m_best_score) { update_best(m_s, new_score); }
            }
        }

        /**
         * Special function, that brute forces all possible combinations. Should
         * be used when number of available candidates and number of required
         * candidates is small.
         *
         * @param depth Current depth of the algorithm.
         */
        inline void execute_small_BF(const size_t depth) {
            const CandidateManager<SFType> &pc_manager = m_c_managers[depth - 1];
            const size_t p_e = pc_manager.get_end();
            const size_t p_s = pc_manager.get_start();

            const size_t size = p_e - p_s;
            const size_t r = m_k - m_s_size;

            help_set.resize(r);
            std::iota(help_set.begin(), help_set.end(), 0);
            help_set[r - 1] -= 1;

            while (next_subset(help_set, r, size)) {
                for (size_t i = 0; i < r; ++i) {
                    m_s[m_s_size + i] = pc_manager[p_s + help_set[i]].c;
                }
                SFType score = m_structure.evaluate(m_s, m_k);

                if (score > m_best_score) { update_best(m_s, score); }
            }
        }

        /**
         * Updates the best found element.
         *
         * @param temp New best set.
         * @param new_score New best score.
         */
        inline void update_best(const std::vector<uint32_t> &temp, const SFType new_score) {
            std::copy(temp.begin(), temp.end(), m_best_s.begin());
            m_best_score = new_score;

            stat_profiler.new_best(new_score);
        }

        /**
         * Check if the execution time has exceeded the specified time limit.
         *
         * @return True if the execution time has exceeded the time limit, else False.
         */
        inline bool has_time_exceeded() {
            if (m_tried_calls > 100) {
                double seconds = get_elapsed_seconds(m_start_point, get_time_point());
                m_time_exceeded = seconds > m_time_limit;
                m_tried_calls = 0;
                return m_time_exceeded;
            }
            m_tried_calls += 1;
            return m_time_exceeded;
        }
    };
}

#endif //SUBMODST_SOLVER_H
