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

#ifndef SUBMODST_ALGORITHMCONFIGURATION_H
#define SUBMODST_ALGORITHMCONFIGURATION_H

#include <string>
#include <numeric>
#include <utility>
#include "CommandLineParser.h"

namespace SubModST {
    /**
     * Class for defining the configuration of the algorithm.
     */
    class AlgorithmConfiguration {
    public:
        // needed variables
        std::string input_file_path;
        std::string function;
        size_t k = 0;
        std::string output_file_path;

        // initial solution
        std::string initial_solution_file_path;

        // timing variables
        double time_limit = std::numeric_limits<double>::max();

        // special variable
        std::string nickname;

        // heuristics and methods
        bool enabled_SUB = false;
        bool enabled_CR = false;
        bool enabled_DCO_heap = false;

        // Lazy evaluation variables
        size_t LE_mode = 0; // 0 - score or rank, 1 - score and rank
        double LE_rank_weight = 2;
        size_t LE_rank_var = 0; // 0 - n, 1 - k, 2 - r_n, 3 - r_k
        double LE_score_weight = 0;

        // Pairwise upper bound
        bool enabled_PUB = false;
        size_t PUB_algo_type = 0; // 0 - greedy, 1 - matching, 2 - dynamic
        size_t PUB_l_type = 0; // 0 - constant, 1 - k', 2 - sqrt(n')
        size_t PUB_l = 0;

    public:
        /**
         * Default constructor.
         */
        AlgorithmConfiguration() = default;

        /**
         * Initializes the configuration, based on the passed arguments.
         *
         * @param argc Number of arguments.
         * @param argv Char arrays to the arguments.
         */
        AlgorithmConfiguration(int argc, char *argv[]) {
            auto cmp = CommandLineParser(argc, argv);

            // read input file
            if (cmp.exists("-i", "--input-file")) {
                std::vector<std::string> args = cmp.get_group("-i", "--input-file");
                input_file_path = args[1];
            } else {
                std::cerr << "Argument -i [  --input-file  ] not specified!" << std::endl;
                std::cerr << "Terminating!" << std::endl;
                exit(EXIT_FAILURE);
            }

            // read score function
            if (cmp.exists("-f", "--function")) {
                std::vector<std::string> args = cmp.get_group("-f", "--function");
                function = args[1];
            } else {
                std::cerr << "Argument -f [  --function  ] not specified!" << std::endl;
                std::cerr << "Terminating!" << std::endl;
                exit(EXIT_FAILURE);
            }

            // read k
            if (cmp.exists("-k", "--k")) {
                std::vector<std::string> args = cmp.get_group("-k", "--k");
                k = std::stoi(args[1]);
            } else {
                std::cerr << "Argument -k [  --k  ] not specified!" << std::endl;
                std::cerr << "Terminating!" << std::endl;
                exit(EXIT_FAILURE);
            }

            // read output file
            if (cmp.exists("-o", "--output-file")) {
                std::vector<std::string> args = cmp.get_group("-o", "--output-file");
                output_file_path = args[1];
            } else {
                std::cerr << "Argument -o [  --output-file  ] not specified!" << std::endl;
                std::cerr << "Terminating!" << std::endl;
                exit(EXIT_FAILURE);
            }

            // read initial solution file
            if (cmp.exists("--initial-solution-file")) {
                std::vector<std::string> args = cmp.get_group("--initial-solution-file");
                initial_solution_file_path = args[1];
            }

            // read time limit
            if (cmp.exists("-t", "--time-limit")) {
                std::vector<std::string> args = cmp.get_group("-t", "--time-limit");
                time_limit = std::stod(args[1]);
                if (time_limit < 0) {
                    std::cerr << "Time limit (" << time_limit << ") is negative!" << std::endl;
                    std::cerr << "Terminating!" << std::endl;
                    exit(EXIT_FAILURE);
                }
                if (time_limit == 0.0) {
                    time_limit = std::numeric_limits<double>::max();
                }
            }

            // read nickname
            if (cmp.exists("-n", "--nickname")) {
                std::vector<std::string> args = cmp.get_group("-n", "--nickname");
                initialize_via_nickname(args[1]);
                return;
            }
        }

        /**
         * Constructs the configuration, based on a nickname.
         *
         * @param nickname Nickname of the configuration.
         */
        explicit AlgorithmConfiguration(std::string &nick_name) {
            initialize_via_nickname(nick_name);
        }

        /**
         * Constructs the configuration, based on a nickname.
         *
         * @param nickname Nickname of the configuration.
         */
        explicit AlgorithmConfiguration(std::string nick_name) {
            initialize_via_nickname(nick_name);
        }

        /**
         * Constructs the configuration, based on a nickname.
         *
         * @param nickname Nickname of the configuration.
         */
        inline void initialize_via_nickname(std::string &nick_name) {
            nickname = std::move(nick_name);

            if (nickname == "Plain") {
                enabled_SUB = false;
                enabled_CR = false;
                enabled_DCO_heap = false;
            } else if (nickname == "Simple") {
                enabled_SUB = true;
                enabled_CR = true;
                enabled_DCO_heap = false;
            } else if (nickname == "Simple+") {
                enabled_SUB = true;
                enabled_CR = true;
                enabled_DCO_heap = true;
            } else if (nickname == "LE-Rank") {
                enabled_SUB = true;
                enabled_CR = true;
                enabled_DCO_heap = true;
                LE_mode = 0;
                LE_rank_weight = 3;
                LE_rank_var = 3;
                LE_score_weight = 0.0;
            } else if (nickname == "LE-Score") {
                enabled_SUB = true;
                enabled_CR = true;
                enabled_DCO_heap = true;
                LE_mode = 0;
                LE_rank_weight = 0;
                LE_rank_var = 0;
                LE_score_weight = 1.0;
            } else if (nickname == "LE-RankOrScore") {
                enabled_SUB = true;
                enabled_CR = true;
                enabled_DCO_heap = true;
                LE_mode = 0;
                LE_rank_weight = 3;
                LE_rank_var = 3;
                LE_score_weight = 1.0;
            } else if (nickname == "LE-RankAndScore") {
                enabled_SUB = true;
                enabled_CR = true;
                enabled_DCO_heap = true;
                LE_mode = 1;
                LE_rank_weight = 3;
                LE_rank_var = 3;
                LE_score_weight = 1.0;
            } else if (nickname == "PWG-k^") {
                enabled_SUB = true;
                enabled_CR = true;
                enabled_DCO_heap = true;
                LE_mode = 0;
                LE_rank_weight = 0;
                LE_rank_var = 0;
                LE_score_weight = 1.0;
                enabled_PUB = true;
                PUB_algo_type = 0;
                PUB_l_type = 1;
                PUB_l = 0;
            } else if (nickname == "PWG-Sqrt-n^") {
                enabled_SUB = true;
                enabled_CR = true;
                enabled_DCO_heap = true;
                LE_mode = 0;
                LE_rank_weight = 0;
                LE_rank_var = 0;
                LE_score_weight = 1.0;
                enabled_PUB = true;
                PUB_algo_type = 0;
                PUB_l_type = 2;
                PUB_l = 0;
            } else if (nickname == "PWM-k^") {
                enabled_SUB = true;
                enabled_CR = true;
                enabled_DCO_heap = true;
                LE_mode = 0;
                LE_rank_weight = 0;
                LE_rank_var = 0;
                LE_score_weight = 1.0;
                enabled_PUB = true;
                PUB_algo_type = 1;
                PUB_l_type = 1;
                PUB_l = 0;
            } else if (nickname == "PWM-Sqrt-n^") {
                enabled_SUB = true;
                enabled_CR = true;
                enabled_DCO_heap = true;
                LE_mode = 0;
                LE_rank_weight = 0;
                LE_rank_var = 0;
                LE_score_weight = 1.0;
                enabled_PUB = true;
                PUB_algo_type = 1;
                PUB_l_type = 2;
                PUB_l = 0;
            } else if (nickname == "PWD-k^") {
                enabled_SUB = true;
                enabled_CR = true;
                enabled_DCO_heap = true;
                LE_mode = 0;
                LE_rank_weight = 0;
                LE_rank_var = 0;
                LE_score_weight = 1.0;
                enabled_PUB = true;
                PUB_algo_type = 2;
                PUB_l_type = 1;
                PUB_l = 0;
            } else if (nickname == "PWD-10") {
                enabled_SUB = true;
                enabled_CR = true;
                enabled_DCO_heap = true;
                LE_mode = 0;
                LE_rank_weight = 0;
                LE_rank_var = 0;
                LE_score_weight = 1.0;
                enabled_PUB = true;
                PUB_algo_type = 2;
                PUB_l_type = 0;
                PUB_l = 10;
            } else {
                std::cout << "Nickname: " << nickname << " not recognized! Switching to Simple+." << std::endl;
                enabled_SUB = true;
                enabled_CR = true;
                enabled_DCO_heap = true;
            }
        }

        /**
         * Determines the score threshold for Lazy Evaluation.
         *
         * @param r_score_avg Average score required.
         * @return The threshold.
         */
        inline double determine_LE_score_threshold(double r_score_avg) const {
            return r_score_avg * LE_score_weight;
        }

        /**
         * Determines the rank threshold for Lazy Evaluation.
         *
         * @param t_n Total number of elements.
         * @param t_k Total solution size.
         * @param t_r_n Currently remaining number of elements.
         * @param t_r_k Currently remaining solution size.
         * @return
         */
        inline size_t determine_LE_rank_threshold(size_t t_n, size_t t_k, size_t t_r_n, size_t t_r_k) const {
            size_t var_value = t_n;
            if (LE_rank_var == 0) { var_value = t_n; }
            if (LE_rank_var == 1) { var_value = t_k; }
            if (LE_rank_var == 2) { var_value = t_r_n; }
            if (LE_rank_var == 3) { var_value = t_r_k; }

            auto threshold = (size_t) (((double) var_value) * LE_rank_weight);

            return std::min(threshold, t_n);
        }

        /**
         * Writes all variable to a string in JSON format.
         *
         * @return The string.
         */
        inline std::string to_JSON() const {
            std::stringstream ss;

            ss << "{";

            // needed variables
            ss << "\"input-file-path\": \"" << input_file_path << "\",\n";
            ss << "\"function\": \"" << function << "\",\n";
            ss << "\"k\": " << k << ",\n";
            ss << "\"output-file-path\": \"" << output_file_path << "\",\n";

            // initial solution
            ss << "\"initial-solution-file-path\": \"" << initial_solution_file_path << "\",\n";

            // timing variables
            ss << "\"time-limit\": " << time_limit << ",\n";

            // special variable
            ss << "\"nickname\": \"" << nickname << "\",\n";

            // heuristics and methods
            ss << "\"enabled-SUB\": " << enabled_SUB << ",\n";
            ss << "\"enabled-CR\": " << enabled_CR << ",\n";
            ss << "\"enabled-DCO-HEAP\": " << enabled_DCO_heap << ",\n";

            // Lazy evaluation variables
            ss << "\"LE-mode\": " << LE_mode << ",\n";
            ss << "\"LE-rank-weight\": " << LE_rank_weight << ",\n";
            ss << "\"LE-rank-var\": " << LE_rank_var << ",\n";
            ss << "\"LE-score-weight\": " << LE_score_weight << ",\n";

            // Pairwise upper bound
            ss << "\"enabled-PUB\": " << enabled_PUB << ",\n";
            ss << "\"PUB-algo-type\": " << PUB_algo_type << ",\n";
            ss << "\"PUB-l\": " << PUB_l << "\n";

            ss << "}";
            return ss.str();
        }
    };
}

#endif //SUBMODST_ALGORITHMCONFIGURATION_H
