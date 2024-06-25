"""
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
"""
import os

import numpy as np
import networkx as nx


def generate_graph(n: int, file_path: str):
    """
    Generates a connected graph and save it to a file.

    :param n: Number of nodes.
    :param file_path: Path to store the file.
    :return: None.
    """
    connected = False
    graph = None

    while not connected:
        m = np.random.randint(n, n * n + 1)
        graph = nx.gnm_random_graph(n, m)
        connected = nx.is_connected(graph)

    nx.write_edgelist(graph, file_path, data=False)


def main() -> None:
    """
    Main function.

    :return: None
    """
    graph_folder_path = f'../data/private/Graph/'
    if not os.path.exists(graph_folder_path):
        os.makedirs(graph_folder_path, exist_ok=True)

    n_nodes = [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
    for n in n_nodes:
        graph_n_folder_path = graph_folder_path + f'{n}/'
        if not os.path.exists(graph_n_folder_path):
            os.makedirs(graph_n_folder_path, exist_ok=True)

        for i in range(50):
            file_path = graph_n_folder_path + f'{i}.edges'
            if not os.path.exists(file_path):
                generate_graph(n, file_path)


if __name__ == '__main__':
    main()
