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


def generate_clustering_instance(n: int, d: int, file_path: str):
    """
    Generates a connected graph and save it to a file.

    :param n: Number of points.
    :param d: Dimensionality of each point.
    :param file_path: Path to store the file.
    :return: None.
    """
    vecs = np.random.rand(n, d)
    vecs = np.matrix(vecs)
    with open(file_path, 'wb') as f:
        for line in vecs:
            np.savetxt(f, line, fmt='%.10f')


def main() -> None:
    """
    Main function.

    :return: None
    """
    clustering_folder_path = f'../data/private/Clustering/'
    if not os.path.exists(clustering_folder_path):
        os.makedirs(clustering_folder_path, exist_ok=True)

    n_vecs = [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
    for n in n_vecs:
        clustering_n_folder_path = clustering_folder_path + f'{n}/'
        if not os.path.exists(clustering_n_folder_path):
            os.makedirs(clustering_n_folder_path, exist_ok=True)

        for i in range(10):
            file_path = clustering_n_folder_path + f'{i}.mtx'
            if not os.path.exists(file_path):
                generate_clustering_instance(n, i + 3, file_path)


if __name__ == '__main__':
    main()
