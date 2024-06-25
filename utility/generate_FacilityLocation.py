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


def generate_facility_location_instance(n_f: int, n_c: int, file_path: str):
    """
    Generates a matrix and saves it to a file.

    :param n_f: Number of facilities.
    :param n_c: Number of customers.
    :param file_path: Path to store the file.
    :return: None.
    """
    mtx = np.random.rand(n_f, n_c)
    mtx = np.matrix(mtx).transpose()
    with open(file_path, 'wb') as f:
        for line in mtx:
            np.savetxt(f, line, fmt='%.10f', delimiter=',')


def main() -> None:
    """
    Main function.

    :return: None
    """
    facilitylocation_folder_path = f'../data/private/FacilityLocation/'
    if not os.path.exists(facilitylocation_folder_path):
        os.makedirs(facilitylocation_folder_path, exist_ok=True)

    n_facility = [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
    for n in n_facility:
        facilitylocation_n_folder_path = facilitylocation_folder_path + f'{n}/'
        if not os.path.exists(facilitylocation_n_folder_path):
            os.makedirs(facilitylocation_n_folder_path, exist_ok=True)

        for i in range(10):
            file_path = facilitylocation_n_folder_path + f'{i}.csv'
            if not os.path.exists(file_path):
                generate_facility_location_instance(n, i+3, file_path)


if __name__ == '__main__':
    main()
