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

#ifndef SUBMODST_GRAPH_H
#define SUBMODST_GRAPH_H

#include <string>
#include <iostream>
#include <filesystem>

#include "StructureInterface.h"
#include "../util/utility.h"
#include "Matrix.h"

namespace SubModST {

    /**
     * Class for defining a graph.
     */
    class Graph {
    private:
        size_t m_n_vertex; // number of nodes
        size_t m_n_edges; // number of edges
        std::vector<std::vector<uint32_t>> m_adj_list;

    public:
        /**
         * Reads a graph from the specified file. File should have the following
         * format:
         * idx1 idx2\n
         * idx3 idx4\n
         * ...
         *
         * Vertices have values from 0 to n-1.
         * Lines with '%' are ignored.
         *
         * @param file_path Path to the file.
         * @return The Graph.
         */
        explicit Graph(const std::string &file_path) {
            if (!SubModST::file_exists(file_path)) {
                std::cout << "File " << file_path << " was not found!\n";
                exit(EXIT_FAILURE);
            }

            std::ifstream file(file_path);
            std::string line;
            std::vector<uint32_t> edges;

            while (std::getline(file, line)) {
                if (line[0] != '%') {
                    std::istringstream iss(line);
                    uint32_t a, b;
                    iss >> a >> b;

                    edges.push_back(a);
                    edges.push_back(b);
                }
            }
            file.close();

            m_n_vertex = *std::max_element(edges.begin(), edges.end()) + 1;
            m_n_edges = edges.size() / 2;

            m_adj_list.resize(m_n_vertex, {});
            for (size_t i = 0; i < m_n_edges; ++i) {
                m_adj_list[edges[i * 2]].push_back(edges[i * 2 + 1]);
                m_adj_list[edges[i * 2 + 1]].push_back(edges[i * 2]);
            }
        };

        /**
         * Returns the matrix used to optimize Partial Dominating Set. The
         * matrix has nxn entries. If vertices i and j are connected in the
         * graph then m_ij = 1, else m_ij = 0.
         *
         * @return The matrix.
         */
        inline Matrix<int> get_PartialDominatingSetMatrix() const {
            Matrix<int> mtx(m_n_vertex, m_n_vertex);
            mtx.set_special_scores(0, (int) m_n_vertex);

            // set all entries to 0
            for (size_t i = 0; i < m_n_vertex; ++i) {
                for (size_t j = 0; j < m_n_vertex; ++j) {
                    mtx.set(i, j, 0);
                }
            }

            // set edges to 1
            for (size_t i = 0; i < m_n_vertex; ++i) {
                for (auto &j: m_adj_list[i]) {
                    mtx.set(i, j, 1);
                }
            }

            return mtx;
        }

        /**
         * Returns the matrix used to optimize Negative Group Farness. The
         * matrix has nxn entries. The entry m_ij is the distance of the
         * shortest path between i an j, but negative. The matrix represents the
         * negative distance matrix of the graph. The distance matrix has to be
         * negative, because the solver expects a maximization problem.
         * The graph must be connected, else it results in undefined behaviour.
         *
         * @return The matrix.
         */
        inline Matrix<int> get_GroupClosenessCentralityMatrix() const {
            Matrix<int> mtx(m_n_vertex, m_n_vertex);
            mtx.set_special_scores(-(int) (m_n_vertex * m_n_vertex), std::numeric_limits<int>::max());

            // initialize arrays to help
            std::vector<uint32_t> stack1(m_n_vertex);
            std::vector<uint32_t> stack2(m_n_vertex);
            std::vector<uint8_t> bool_arr(m_n_vertex);

            for (uint32_t i = 0; i < m_n_vertex; ++i) {
                // for every vertex, calculate distance to all other vertices
                int curr_distance = 0;
                size_t stack1_size = 0;
                size_t stack2_size = 0;
                std::fill(bool_arr.begin(), bool_arr.end(), 0);

                stack1[0] = i;
                stack1_size += 1;

                while (stack1_size != 0) {
                    // set distance for current stack
                    for (size_t j = 0; j < stack1_size; ++j) {
                        mtx.set(i, stack1[j], -curr_distance); // negative
                        bool_arr[stack1[j]] = 1;
                    }

                    // get the next stack
                    for (size_t j = 0; j < stack1_size; ++j) {
                        uint32_t vertex = stack1[j];

                        for (uint32_t neighbour: m_adj_list[vertex]) {
                            if (!bool_arr[neighbour]) {
                                bool_arr[neighbour] = 1;
                                stack2[stack2_size] = neighbour;
                                stack2_size += 1;
                            }
                        }
                    }
                    curr_distance += 1;

                    // swap stacks
                    stack1.swap(stack2);
                    std::swap(stack1_size, stack2_size);
                    stack2_size = 0;
                }
            }
            return mtx;
        }

        /**
         * Returns the matrix used to optimize Partial Vertex Set. The
         * matrix has nxm entries. If vertex i and edge j are connected in the
         * graph then m_ij = 1, else m_ij = 0.
         *
         * @return The matrix.
         */
        inline Matrix<int> get_PartialVertexSetMatrix() const {
            Matrix<int> mtx(m_n_vertex, m_n_edges);
            mtx.set_special_scores(0, (int) m_n_edges);

            // set all entries to 0
            for (size_t i = 0; i < m_n_vertex; ++i) {
                for (size_t j = 0; j < m_n_edges; ++j) {
                    mtx.set(i, j, 0);
                }
            }

            struct Edge {
                uint32_t v1;
                uint32_t v2;
            };

            // get all edges
            std::vector<Edge> edges;
            for (size_t i = 0; i < m_n_vertex; ++i) {
                for (auto &j: m_adj_list[i]) {
                    uint32_t v1 = std::min((uint32_t) i, j);
                    uint32_t v2 = std::max((uint32_t) i, j);
                    edges.push_back({v1, v2});
                }
            }

            // remove duplicate edges
            std::sort(edges.begin(), edges.end(), [](const Edge &a, const Edge &b) {
                return a.v1 < b.v1 || (a.v1 == b.v1 && a.v2 < b.v2);
            });
            auto it = std::unique(edges.begin(), edges.end(), [](const Edge &a, const Edge &b) {
                return a.v1 == b.v1 && a.v2 == b.v2;
            });
            edges.erase(it, edges.end());
            edges.shrink_to_fit();

            // set edges to 1
            for (size_t i = 0; i < edges.size(); ++i) {
                uint32_t v1 = edges[i].v1;
                uint32_t v2 = edges[i].v2;
                mtx.set(v1, i, 1);
                mtx.set(v2, i, 1);
            }

            return mtx;
        }
    };
}

#endif //SUBMODST_GRAPH_H
