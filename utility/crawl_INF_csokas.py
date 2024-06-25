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
import glob
import os

import numpy as np

RAW_FOLDER_PATH = "../data/raw/INF_new_data_ALL"
INF_FOLDER_PATH = "../data/BipartiteInfluence"


def reformat(file_path: str):
    """
    Formats the content of the file to be readable like all other INF instances.

    :param file_path: Path to the file.
    :return: Formatted content of the file
    """
    lines = []
    f = open(file_path, 'r')

    while True:
        l = f.readline()
        if not l:
            break
        lines.append(l.strip())
    f.close()

    # read set N
    set_N = []
    assert lines[0] == 'set N :='
    lines.pop(0)

    while lines[0] != ';':
        set_N.append(int(lines[0]))
        lines.pop(0)
    lines.pop(0)

    # read set M
    set_M = []
    assert lines[0] == 'set M :='
    lines.pop(0)

    while lines[0] != ';':
        set_M.append(int(lines[0]))
        lines.pop(0)
    lines.pop(0)

    # read weights
    weights = {}
    assert lines[0] == 'param: M_N: g :='
    lines.pop(0)

    while lines[0] != ';':
        temp = lines[0].split(' ')
        m, n, w = int(temp[0]), int(temp[1]), float(temp[2]),
        weights[(n, m)] = w
        lines.pop(0)
    lines.pop(0)

    # make matrix
    n, m = len(set_N), len(set_M)
    mtx = np.full((m, n), 0.0, dtype=float)

    # fill matrix
    n_offset, m_offset = min(set_N), min(set_M)
    for (n, m), w in weights.items():
        mtx[m - m_offset][n - n_offset] = w

    mtx = np.matrix(mtx)

    filename = os.path.splitext(os.path.basename(file_path))[0]
    with open(f'{INF_FOLDER_PATH}/{filename}.csv', 'wb') as f:
        for line in mtx:
            np.savetxt(f, line, fmt='%.16f', delimiter=',')


if __name__ == '__main__':
    # create directories
    if not os.path.exists(INF_FOLDER_PATH):
        os.mkdir(INF_FOLDER_PATH)

    # move the data into the right folders
    for file in glob.glob(f'{RAW_FOLDER_PATH}/*.dat'):
        if '(' in file:
            continue

        reformat(file)
