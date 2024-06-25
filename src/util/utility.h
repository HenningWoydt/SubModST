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

#ifndef SUBMODST_UTILITY_H
#define SUBMODST_UTILITY_H

#include <string>
#include <fstream>
#include <vector>
#include <algorithm>
#include <filesystem>
#include <iostream>
#include "AlgorithmConfiguration.h"

namespace SubModST {

    /**
     * Reads a m_solution from a .sol file.
     *
     * @param file_path Path to the file.
     * @param shift Shifts the read numbers to the left (x -= shift).
     * @return Vector containing the m_solution.
     */
    std::vector<uint32_t> read_solution(const std::string &file_path);

    /**
     * Write a m_solution vector to a file.
     *
     * @param solution Vector holding the solution.
     * @param file_path Path to the file.
     */
    void write_solution(std::vector<uint32_t> &solution, std::string &file_path);

    /**
     * Prints the vector in Numpy style.
     *
     * @tparam T Type of the vector.
     * @param vec The vector.
     */
    template<typename T>
    void print(std::vector<T> &&vec) {
        if (vec.empty()) {
            std::cout << "[]" << std::endl;
            return;
        }

        std::cout << "[";
        for (size_t i = 0; i < vec.size() - 1; ++i) {
            std::cout << vec[i] << ", ";
        }
        std::cout << vec[vec.size() - 1] << "]" << std::endl;
    }

    /**
     * Prints the vector in Numpy style.
     *
     * @tparam T Type of the vector.
     * @param vec The vector.
     */
    template<typename T>
    void print(std::vector<T> &vec) {
        if (vec.empty()) {
            std::cout << "[]" << std::endl;
            return;
        }

        std::cout << "[";
        for (size_t i = 0; i < vec.size() - 1; ++i) {
            std::cout << vec[i] << ", ";
        }
        std::cout << vec[vec.size() - 1] << "]" << std::endl;
    }

    /**
     * Prints the vector in Numpy style.
     *
     * @tparam T Type of the vector.
     * @param vec The vector.
     * @param size The size.
     */
    template<typename T>
    void print(std::vector<T> &vec, size_t size) {
        if (size == 0) {
            std::cout << "[]" << std::endl;
            return;
        }

        std::cout << "[";
        for (size_t i = 0; i < size - 1; ++i) {
            std::cout << vec[i] << ", ";
        }
        std::cout << vec[size - 1] << "]" << std::endl;
    }

    /**
     * Determines if a file exists at a given path.
     *
     * @param file_path Path to the file.
     * @return True if the file exists, False else.
     */
    inline bool file_exists(const std::string &file_path) {
        std::ifstream f(file_path.c_str());
        return f.good();
    }

    /**
     * Gets all files in a specific directory.
     *
     * @param directory_path Path to the directory.
     * @return Vector holding all paths.
     */
    std::vector<std::string> get_directory_files(std::string &directory_path);

    /**
 * Gets all files in a specific directory with the specified extension.
 *
 * @param directory_path Path to the directory.
 * @param extension Allowed extension.
 * @return Vector holding all paths.
 */
    std::vector<std::string> get_directory_files(std::string &directory_path, std::string &extension);

    /**
     * Gets all files in a specific directory with the specified extension.
     *
     * @param directory_path Path to the directory.
     * @param extension Allowed extension.
     * @return Vector holding all paths.
     */
    std::vector<std::string> get_directory_files(std::string &&directory_path, std::string &&extension);

    /**
     * Calculate the binomial coefficient "n choose k."
     *
     * @param n The total number of items in the set.
     * @param k The number of items to choose from the set.
     * @return The binomial coefficient "n choose k."
     */
    size_t n_choose_k(size_t n, size_t k);

    /**
     * Iterates over all possible subsets of {0, ..., n - 1} with size k.
     *
     * @param set The set.
     * @param k Size of the set.
     * @param n Number of elements.
     * @return True if the subset is valid.
     */
    bool next_subset(std::vector<size_t> &set, size_t k, size_t n);

    /**
     * Returns the current time point.
     *
     * @return The time point.
     */
    inline std::chrono::steady_clock::time_point get_time_point() {
        return std::chrono::steady_clock::now();
    }

    /**
     * Get the elapsed seconds between the both time points.
     *
     * @param sp Start point.
     * @param ep End point.
     * @return Elapsed time in seconds
     */
    inline double get_elapsed_seconds(std::chrono::steady_clock::time_point sp,
                                      std::chrono::steady_clock::time_point ep) {
        return (double) std::chrono::duration_cast<std::chrono::nanoseconds>(ep - sp).count() / 1e9;
    }

    /**
     * Returns almost all configurations. Use this function for testing
     * purposes.
     *
     * @return The configurations.
     */
    std::vector<AlgorithmConfiguration> get_all_acs();

    /**
     * This function takes a file path and content as input parameters and writes
     * the content to the specified file. If the file does not exist, it will be
     * created. If the file already exists, its contents will be replaced with the
     * new content.
     *
     * @param file_path The path to the file where content will be written.
     * @param content The content to be written to the file.
     */
    void write_to_file(std::string &file_path, std::string &content);

    /**
     * Splits a string into a vector of substrings based on a delimiter character.
     *
     * @param str The string to split.
     * @param c The delimiter character.
     * @return A vector of substrings.
     */
    std::vector<std::string> split(const std::string &s, char c);

    template<typename T>
    T sum(const std::vector<T> &vec) {
        T sum = 0;
        for (size_t i = 0; i < vec.size(); ++i) {
            sum += vec[i];
        }
        return sum;
    }

    template<typename T>
    std::string to_JSON(const std::vector<T> &vec) {
        if (vec.empty()) {
            return "[]";
        }

        std::stringstream ss;
        ss << "[";
        for (size_t i = 0; i < vec.size() - 1; ++i) {
            ss << vec[i] << ", ";
        }
        ss << vec.back() << "]";
        return ss.str();
    }

    /**
     * Calculate the power of 2 raised to the given exponent.
     *
     * @param n The exponent for the power of 2.
     * @return The result of 2 raised to the power of 'n'.
     */
    inline size_t pow_2(size_t n) {
        return 1 << n;
    }

    /**
     * Calculates ceil of x/y.
     *
     * @param x x.
     * @param y y.
     * @return ceil(c/y)
     */
    inline size_t ceil(size_t x, size_t y) {
        return (x + y - 1) / y;
    }

    /**
     * Returns the bit at position pos of n.
     *
     * @param n Number.
     * @param pos The position.
     * @return True if bit is 1, False else.
     */
    inline bool get_bit(size_t n, size_t pos) {
        return (n & (1 << pos)) != 0;
    }
}

#endif //SUBMODST_UTILITY_H
