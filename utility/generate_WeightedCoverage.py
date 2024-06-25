"""
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
"""

import os

import numpy as np


def generate_weighted_coverage_instance(n_s: int, n_i: int, file_path: str):
    """
    Generates a connected graph and save it to a file.

    :param n_s: Number of sensors.
    :param n_i: Number of items.
    :param file_path: Path to store the file.
    :return: None.
    """
    weights = np.matrix(np.random.rand(n_i))
    mask = np.random.choice([0, 1], (n_s, n_i), p=[0.85, 0.15])
    mask = np.matrix(mask).transpose()
    with open(file_path, 'wb') as f:
        np.savetxt(f, weights, fmt='%.10f', delimiter=',')
        for line in mask:
            np.savetxt(f, line, fmt='%.10f', delimiter=',')


if __name__ == '__main__':
    weightedcoverage_folder_path = f'../data/private/WeightedCoverage/'
    if not os.path.exists(weightedcoverage_folder_path):
        os.makedirs(weightedcoverage_folder_path, exist_ok=True)

    n_sensors = [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
    for n in n_sensors:
        weightedcoverage_n_folder_path = weightedcoverage_folder_path + f'{n}/'
        if not os.path.exists(weightedcoverage_n_folder_path):
            os.makedirs(weightedcoverage_n_folder_path, exist_ok=True)

        for i in range(1000):
            file_path = weightedcoverage_n_folder_path + f'{i}.csv'
            if not os.path.exists(file_path):
                generate_weighted_coverage_instance(n, 1 + (i // 50), file_path)
