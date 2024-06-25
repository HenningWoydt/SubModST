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

#ifndef CCSMSM_STATISTICPROFILER_H
#define CCSMSM_STATISTICPROFILER_H

#include <cstddef>
#include <vector>
#include <string>
#include <sstream>

#include "../util/utility.h"

#define STAT_PROFILER_ENABLED 1

namespace CCSMSM {

    /**
     * Class to track statistics of the solver.
     */
    template<typename SFType>
    class StatisticProfiler {
    private:
#if STAT_PROFILER_ENABLED
        size_t m_n = 0;
        size_t m_k = 0;

        // tracks the number of visits
        std::vector<size_t> m_recursive_calls_per_depth;

        // start point of the search
        std::chrono::steady_clock::time_point m_start_point;

        // tracks the score improvement over time
        std::vector<double> history_time;
        std::vector<SFType> history_score;

        // track PUB usefulness
        size_t PUB_total = 0;
        size_t PUB_success = 0;
        size_t PUB_fail = 0;
        double PUB_time_total = 0.0;
        double PUB_time_success = 0.0;
        double PUB_time_fail = 0.0;
        std::chrono::steady_clock::time_point PUB_start_point;
#endif

    public:
        StatisticProfiler() = default;

        StatisticProfiler(size_t n, size_t k) {
#if STAT_PROFILER_ENABLED
            m_n = n;
            m_k = k;

            m_recursive_calls_per_depth.resize(m_k, 0);
#endif
        }

        inline void set_start_point(std::chrono::steady_clock::time_point sp) {
#if STAT_PROFILER_ENABLED
            m_start_point = sp;
#endif
        }

        inline void visited_depth(size_t depth) {
#if STAT_PROFILER_ENABLED
            m_recursive_calls_per_depth[depth]++;
#endif
        }

        inline void new_best(SFType score) {
#if STAT_PROFILER_ENABLED
            std::chrono::steady_clock::time_point p = std::chrono::steady_clock::now();
            double t = get_elapsed_seconds(m_start_point, p);
            history_time.push_back(t);
            history_score.push_back(score);
#endif
        }

        inline void start_track_PUB() {
#if STAT_PROFILER_ENABLED
            PUB_start_point = std::chrono::steady_clock::now();
#endif
        }

        inline void end_track_PUB(bool success) {
#if STAT_PROFILER_ENABLED
            double t = get_elapsed_seconds(PUB_start_point, std::chrono::steady_clock::now());
            PUB_total += 1;
            PUB_time_total += t;
            if (success) {
                PUB_success += 1;
                PUB_time_success += t;
            } else {
                PUB_fail += 1;
                PUB_time_fail += t;
            }
#endif
        }

        inline std::string to_JSON() const {
            std::stringstream ss;
            ss << "{";
#if STAT_PROFILER_ENABLED
            ss << "\"recursive-calls-per-depth\": " << CCSMSM::to_JSON(m_recursive_calls_per_depth) << ",\n";

            ss << "\"history-time\": " << CCSMSM::to_JSON(history_time) << ",\n";
            ss << "\"history-score\": " << CCSMSM::to_JSON(history_score) << ",\n";

            ss << "\"PUB-total\": " << PUB_total << ",\n";
            ss << "\"PUB-success\": " << PUB_success << ",\n";
            ss << "\"PUB-fail\": " << PUB_fail << ",\n";
            ss << "\"PUB-time-total\": " << PUB_time_total << ",\n";
            ss << "\"PUB-time-success\": " << PUB_time_success << ",\n";
            ss << "\"PUB-time-fail\": " << PUB_time_fail << "\n";
#endif
            ss << "}";
            return ss.str();
        }

    };

}

#endif //CCSMSM_STATISTICPROFILER_H
