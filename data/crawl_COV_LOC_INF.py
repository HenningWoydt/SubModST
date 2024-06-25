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
import glob
import io
import os
import shutil
import zipfile

import requests

RAW_FOLDER_PATH = "raw/COV-INF-LOC"
COV_FOLDER_PATH = "WeightedCoverage"
LOC_FOLDER_PATH = "FacilityLocation"
INF_FOLDER_PATH = "BipartiteInfluence"

if __name__ == '__main__':
    # create directories
    if not os.path.exists(COV_FOLDER_PATH):
        os.mkdir(COV_FOLDER_PATH)
    if not os.path.exists(INF_FOLDER_PATH):
        os.mkdir(INF_FOLDER_PATH)
    if not os.path.exists(LOC_FOLDER_PATH):
        os.mkdir(LOC_FOLDER_PATH)

    # download the data
    if not os.path.exists(RAW_FOLDER_PATH):
        r = requests.get('https://drive.usercontent.google.com/u/0/uc?id=1Ra3ctxgtLbZeNzmVEYggHh3K1tlTUXjs&export=download')
        z = zipfile.ZipFile(io.BytesIO(r.content))
        z.extractall(RAW_FOLDER_PATH)

    # move the data into the right folders
    for file in glob.glob(f'{RAW_FOLDER_PATH}/instance/COV/*.csv'):
        shutil.copy2(file, COV_FOLDER_PATH)

    for file in glob.glob(f'{RAW_FOLDER_PATH}/instance/INF/*.csv'):
        shutil.copy2(file, INF_FOLDER_PATH)

    for file in glob.glob(f'{RAW_FOLDER_PATH}/instance/LOC/*.csv'):
        shutil.copy2(file, LOC_FOLDER_PATH)
