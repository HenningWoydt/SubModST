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

#include "utility.h"

std::vector<uint32_t> SubModST::read_solution(const std::string &file_path) {
    std::vector<uint32_t> vec;

    std::ifstream file(file_path);
    if (file.is_open()) {
        for (std::string line; std::getline(file, line);) {
            if (line.back() == '\n' || line.back() == '\r') {
                line.pop_back();
            }
            int vertex = std::stoi(line);
            vec.push_back(vertex);
        }
    } else {
        std::cout << "Could not open file " << file_path << " !" << std::endl;
    }

    return vec;
}

void SubModST::write_solution(std::vector<uint32_t> &solution, std::string &file_path) {
    std::ofstream file(file_path);

    if (file.is_open()) {
        for (uint32_t i: solution) {
            file << i << "\n";
        }
        file.close();
    } else {
        std::cout << "Could not open file " << file_path << " !" << std::endl;
    }
}

std::vector<std::string> SubModST::get_directory_files(std::string &directory_path) {
    std::vector<std::string> paths;
    for (const auto &entry: std::filesystem::directory_iterator(directory_path)) {
        paths.push_back(entry.path().string());
    }
    std::sort(paths.begin(), paths.end());
    return paths;
}

std::vector<std::string> SubModST::get_directory_files(std::string &directory_path, std::string &extension) {
    std::vector<std::string> paths;
    for (const auto &entry: std::filesystem::directory_iterator(directory_path)) {
        if (entry.path().extension() == extension) {
            paths.push_back(entry.path().string());
        }
    }
    std::sort(paths.begin(), paths.end());
    return paths;
}

std::vector<std::string> SubModST::get_directory_files(std::string &&directory_path, std::string &&extension) {
    std::filesystem::create_directories(directory_path);
    std::vector<std::string> paths;
    for (const auto &entry: std::filesystem::directory_iterator(directory_path)) {
        if (entry.path().extension() == extension) {
            paths.push_back(entry.path().string());
        }
    }
    std::sort(paths.begin(), paths.end());
    return paths;
}

size_t SubModST::n_choose_k(size_t n, size_t k) {
    if (k == 0) {
        return 1;
    }
    return (n * n_choose_k(n - 1, k - 1)) / k;
}

bool SubModST::next_subset(std::vector<size_t> &set, size_t set_size, size_t n) {
    set[set_size - 1] += 1;
    if (set[set_size - 1] < n - 1) {
        return true;
    }
    size_t j = set_size - 1;

    for (size_t i = 0; i < set_size - 1; ++i) {
        if (set[set_size - i - 1] == n - i) {
            // this placed has reached its maximum
            set[set_size - i - 2] += 1;
            j = set_size - i - 2;
        } else {
            break;
        }
    }

    for (size_t k = j + 1; k < set_size; ++k) {
        set[k] = set[k - 1] + 1;
    }
    return set[0] <= n - set_size;
}

std::vector<SubModST::AlgorithmConfiguration> SubModST::get_all_acs() {
    std::vector<SubModST::AlgorithmConfiguration> acs;

    acs = {
            // AlgorithmConfiguration("Plain"),
            // AlgorithmConfiguration("Simple"),
            AlgorithmConfiguration("Simple+"),
            // AlgorithmConfiguration("LE-Rank"),
            AlgorithmConfiguration("LE-Score"),
            // AlgorithmConfiguration("LE-RankOrScore"),
            // AlgorithmConfiguration("LE-RankAndScore"),
            // AlgorithmConfiguration("PWG-k^"),
            // AlgorithmConfiguration("PWG-Sqrt-n^"),
            // AlgorithmConfiguration("PWM-k^"),
            // AlgorithmConfiguration("PWM-Sqrt-n^"),
            // AlgorithmConfiguration("PWD-k^"),
            // AlgorithmConfiguration("PWD-10"),
    };

    return acs;
}

void SubModST::write_to_file(std::string &file_path, std::string &content) {
    // Open the file for writing
    std::ofstream file(file_path);

    // Check if the file is opened successfully
    if (file.is_open()) {
        // Write content to the file
        file << content;
        // Close the file
        file.close();
    } else {
        // Print an error message if the file cannot be opened
        std::cerr << "Error: Unable to open file for writing." << std::endl;
    }
}

std::vector<std::string> SubModST::split(const std::string &s, char c) {
    std::vector<std::string> result;
    std::string token;
    std::istringstream iss(s);

    while (std::getline(iss, token, c)) {
        result.push_back(token);
    }

    return result;
}

