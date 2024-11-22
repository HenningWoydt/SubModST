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

#ifndef SUBMODST_CSRGRAPHPDS_H
#define SUBMODST_CSRGRAPHPDS_H

#include <filesystem>
#include <iostream>
#include <string>

#include "Matrix.h"
#include "StructureInterface.h"
#include "../util/utility.h"

namespace SubModST {

    /**
     * Class for defining a graph.
     */
    template<typename SFType>
    class CSRGraphPDS final : public StructureInterface<SFType> {
    private:
        size_t m_n_vertex; // number of nodes
        size_t m_n_edges; // number of edges

        std::vector<uint64_t> m_sizes;
        std::vector<uint32_t> m_neighborhood;

        // structures to speed up score function evaluation
        size_t m_depth;
        std::vector<std::vector<SFType>> m_temp_rows;
        std::vector<SFType> m_temp;

        std::vector<SFType> m_score_cache;

    public:
        /**
         * Reads a graph from the specified file. File should have the following
         * format:
         * p identifier n m
         * idx1 idx2\n
         * idx3 idx4\n
         * ...
         *
         * Vertices have values from 1 to n.
         * Lines with 'c' are ignored.
         *
         * @param file_path Path to the file.
         * @return The Graph in CSR format.
         */
        explicit CSRGraphPDS(const std::string &file_path) {
            if (!SubModST::file_exists(file_path)) {
                std::cout << "File " << file_path << " was not found!\n";
                exit(EXIT_FAILURE);
            }

            std::ifstream file(file_path);
            std::string line;
            std::vector<std::vector<uint32_t>> adj_list;
            bool header_processed = false;

            while (std::getline(file, line)) {
                if (line[0] != 'c') {
                    if(!header_processed) {
                        std::istringstream iss(line);
                        std::string p, ds, n_vertex_str, n_edges_str;
                        iss >> p >> ds >> n_vertex_str >> n_edges_str;

                        m_n_vertex = std::stoi(n_vertex_str);
                        m_n_edges = std::stoi(n_edges_str);

                        m_sizes.resize(m_n_vertex+1);
                        m_neighborhood.resize(m_n_edges*2);

                        adj_list.resize(m_n_vertex);

                        header_processed = true;
                    } else {
                        std::istringstream iss(line);
                        uint32_t a, b;
                        iss >> a >> b;

                        adj_list[a-1].push_back(b-1);
                        adj_list[b-1].push_back(a-1);
                    }
                }
            }
            file.close();

            m_sizes[0] = 0;

            for(size_t i = 0; i < m_n_vertex; ++i) {
                m_sizes[i + 1] = m_sizes[i] + adj_list[i].size();
                for(size_t j = 0; j < adj_list[i].size(); ++j) {
                    m_neighborhood[m_sizes[i] + j] = adj_list[i][j];
                }
            }

            // structures to speed up score function evaluation
            m_depth = 0;

            StructureInterface<SFType>::set_special_scores(0, m_n_vertex);
        }

        /**
         * Returns the number of rows in the matrix.
         *
         * @return Number of rows.
         */
        inline size_t get_n() override {
            return m_n_vertex;
        }

        inline SFType evaluate_empty_set() override {
            return StructureInterface<SFType>::m_empty_set_score;
        };

        inline SFType evaluate_1D(const std::vector<uint32_t> &s, const size_t s_size) override {
            SFType score = m_score_cache[m_depth];
            uint32_t v = s[s_size - 1];
            score += m_temp_rows[m_depth][v] == 0;
            for(size_t i = 0; i < m_sizes[v + 1] - m_sizes[v]; ++i) {
                uint32_t u = m_neighborhood[m_sizes[v] + i];
                score += m_temp_rows[m_depth][u] == 0;
            }
            return score;
        };

        inline SFType evaluate_2D(const std::vector<uint32_t> &s, const size_t s_size) override {
            return evaluate(s, s_size);
            /*
            SFType score = 0;
            for (size_t i = 0; i < m_n_col; ++i) {
                SFType a = m_temp_rows[m_depth][i];
                SFType b = m_mtx[s[s_size - 2] * m_n_col + i];
                SFType c = m_mtx[s[s_size - 1] * m_n_col + i];
                score += std::max(a, std::max(b, c));
            }
            return score;
            */
        };

        inline SFType evaluate(const std::vector<uint32_t> &s, size_t s_size) override {
            // reset m_temp
            std::fill(m_temp.begin(), m_temp.end(), 0);

            SFType score = 0;

            // mark each dominated vertex and count how many
            for (size_t j = 0; j < s_size; ++j) {
                uint32_t v = s[j];
                score += m_temp[v] == 0;
                m_temp[v] = 1;
                for(size_t i = 0; i < m_sizes[v + 1] - m_sizes[v]; ++i) {
                    uint32_t u = m_neighborhood[m_sizes[v] + i];
                    score += m_temp[u] == 0;
                    m_temp[u] = 1;
                }
            }

            return score;
        };

        inline void initialize_helping_structures(size_t k) override {
            m_depth = 0;
            m_temp_rows.resize(k + 1, std::vector<SFType>(m_n_vertex, 0));
            m_temp.resize(m_n_vertex);

            m_score_cache.resize(k, 0);
        };

        inline void visit_new_depth(const std::vector<uint32_t> &s, size_t s_size) override {
            m_depth += 1;

            // copy last row
            std::copy(m_temp_rows[m_depth - 1].begin(), m_temp_rows[m_depth - 1].end(), m_temp_rows[m_depth].begin());
            m_score_cache[m_depth] = m_score_cache[m_depth - 1];

            // insert newly dominated vertices
            uint32_t v = s[s_size - 1];
            m_score_cache[m_depth] += m_temp_rows[m_depth][v] == 0;
            m_temp_rows[m_depth][v] = 1;

            for(size_t i = 0; i < m_sizes[v + 1] - m_sizes[v]; ++i) {
                uint32_t u = m_neighborhood[m_sizes[v] + i];
                m_score_cache[m_depth] += m_temp_rows[m_depth][u] == 0;
                m_temp_rows[m_depth][u] = 1;
            }
        };

        inline void return_from_last_depth() override {
            m_depth -= 1;
        };
    };
}

#endif //SUBMODST_CSRGRAPHPDS_H
