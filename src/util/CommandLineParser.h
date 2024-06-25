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

#ifndef SUBMODST_COMMANDLINEPARSER_H
#define SUBMODST_COMMANDLINEPARSER_H

#include <vector>
#include <string>

namespace SubModST {
    /**
     * Class for defining the configuration of the algorithm.
     */
    class CommandLineParser {
    private:
        std::vector<std::vector<std::string>> m_grouped_args;

    public:
        /**
         * Reads in the command line arguments.
         *
         * @param argc Number of arguments.
         * @param argv Arguments.
         */
        CommandLineParser(int argc, char *argv[]) {
            // read in args as string
            std::vector<std::string> args;
            for (int i = 1; i < argc; ++i) {
                args.emplace_back(argv[i]);
            }

            // group all arguments
            for (auto &arg: args) {
                if (arg[0] == '-') {
                    m_grouped_args.emplace_back();
                }
                m_grouped_args.back().emplace_back(arg);
            }
        }

        /**
         * Determines if a group of arguments belong to the specified strings.
         *
         * @param short_str Short version of the string.
         * @param long_str Long version of the string.
         * @return True if the group exists, false else.
         */
        inline bool exists(std::string &short_str, std::string &long_str) {
            for (auto &vec: m_grouped_args) {
                if (vec[0] == short_str || vec[0] == long_str) {
                    return true;
                }
            }
            return false;
        }

        /**
         * Determines if a group of arguments belong to the specified strings.
         *
         * @param short_str Short version of the string.
         * @param long_str Long version of the string.
         * @return True if the group exists, false else.
         */
        inline bool exists(std::string &&short_str, std::string &&long_str) {
            for (auto &vec: m_grouped_args) {
                if (vec[0] == short_str || vec[0] == long_str) {
                    return true;
                }
            }
            return false;
        }

        /**
         * Determines if a group of arguments belong to the specified strings.
         *
         * @param long_str Long version of the string.
         * @return True if the group exists, false else.
         */
        inline bool exists(std::string &long_str) {
            for (auto &vec: m_grouped_args) {
                if (vec[0] == long_str) {
                    return true;
                }
            }
            return false;
        }

        /**
         * Determines if a group of arguments belong to the specified strings.
         *
         * @param long_str Long version of the string.
         * @return True if the group exists, false else.
         */
        inline bool exists(std::string &&long_str) {
            for (auto &vec: m_grouped_args) {
                if (vec[0] == long_str) {
                    return true;
                }
            }
            return false;
        }

        /**
         * Returns the group fitting to the specified string.
         *
         * @param short_str Short version of the string.
         * @param long_str Long version of the string.
         * @return The group of arguments.
         */
        inline std::vector<std::string> get_group(std::string &short_str, std::string &long_str) {
            for (auto &vec: m_grouped_args) {
                if (vec[0] == short_str || vec[0] == long_str) {
                    return vec;
                }
            }
            return {};
        }

        /**
         * Returns the group fitting to the specified string.
         *
         * @param short_str Short version of the string.
         * @param long_str Long version of the string.
         * @return The group of arguments.
         */
        inline std::vector<std::string> get_group(std::string &&short_str, std::string &&long_str) {
            for (auto &vec: m_grouped_args) {
                if (vec[0] == short_str || vec[0] == long_str) {
                    return vec;
                }
            }
            return {};
        }

        /**
         * Returns the group fitting to the specified string.
         *
         * @param long_str Long version of the string.
         * @return The group of arguments.
         */
        inline std::vector<std::string> get_group(std::string &long_str) {
            for (auto &vec: m_grouped_args) {
                if (vec[0] == long_str) {
                    return vec;
                }
            }
            return {};
        }

        /**
         * Returns the group fitting to the specified string.
         *
         * @param long_str Long version of the string.
         * @return The group of arguments.
         */
        inline std::vector<std::string> get_group(std::string &&long_str) {
            for (auto &vec: m_grouped_args) {
                if (vec[0] == long_str) {
                    return vec;
                }
            }
            return {};
        }
    };
}

#endif //SUBMODST_COMMANDLINEPARSER_H
