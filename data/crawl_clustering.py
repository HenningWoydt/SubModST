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
import os.path
import urllib.request

import numpy as np

URLS = [
    "https://cs.uef.fi/sipu/datasets/s1.txt",
    "https://cs.uef.fi/sipu/datasets/s2.txt",
    "https://cs.uef.fi/sipu/datasets/s3.txt",
    "https://cs.uef.fi/sipu/datasets/s4.txt",
    "https://cs.uef.fi/sipu/datasets/a1.txt",
    "https://cs.uef.fi/sipu/datasets/a2.txt",
    "https://cs.uef.fi/sipu/datasets/a3.txt",
    "https://cs.uef.fi/sipu/datasets/dim032.txt",
    "https://cs.uef.fi/sipu/datasets/dim064.txt",
    "https://cs.uef.fi/sipu/datasets/dim128.txt",
    "https://cs.uef.fi/sipu/datasets/dim256.txt",
    "https://cs.uef.fi/sipu/datasets/dim512.txt",
    "https://cs.uef.fi/sipu/datasets/dim1024.txt",
    "https://cs.uef.fi/sipu/datasets/unbalance2.txt",
    "https://cs.uef.fi/sipu/datasets/asymmetric.txt",
    "https://cs.uef.fi/sipu/datasets/overlap.txt",
    "https://cs.uef.fi/sipu/datasets/skewed.txt"
]

RAW_FOLDER_PATH = "raw/cs.uef.fi/"
FOLDER_PATH = "Clustering/"


def crawl_data():
    """
    Crawls the data from https://cs.uef.fi/sipu/datasets/.

    :return: None
    """
    if not os.path.exists(RAW_FOLDER_PATH):
        os.makedirs(RAW_FOLDER_PATH, exist_ok=True)

    for url in URLS:
        file_name = os.path.basename(url)
        file_path = f"{RAW_FOLDER_PATH}{file_name}"
        if os.path.exists(file_path):
            continue

        data = urllib.request.urlopen(url).read().decode('utf-8')
        f = open(file_path, 'w')
        f.write(data)
        f.close()


def clean_data():
    """
    Cleans the downloaded data.

    :return: None
    """
    if not os.path.exists(FOLDER_PATH):
        os.makedirs(FOLDER_PATH, exist_ok=True)

    for url in URLS:
        file_name = os.path.basename(url)
        raw_file_path = f"{RAW_FOLDER_PATH}{file_name}"
        file_path = f"{FOLDER_PATH}{file_name}"

        x = np.loadtxt(raw_file_path)
        x_norm = (x - x.min(0)) / x.ptp(0)
        np.savetxt(file_path, x_norm)


if __name__ == '__main__':
    crawl_data()
    clean_data()
