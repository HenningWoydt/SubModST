: '
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
'

wget https://pub.ista.ac.at/~vnk/software/blossom5-v2.05.src.tar.gz -O blossom5.tar.gz
tar -xf blossom5.tar.gz
mv blossom5-v2.05.src blossom5

sed -i '40s/.*/#define PERFECT_MATCHING_DOUBLE/' blossom5/PerfectMatching.h
