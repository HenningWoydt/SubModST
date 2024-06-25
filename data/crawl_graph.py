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

import io
import os
import shutil
import zipfile

import networkx as nx
import requests

# each file has a name, whether it has a header and which char splits the
# vertices
FILES = [
    ("ca/ca-netscience", False, ' '),
    ("soc/soc-wiki-Vote", True, ' '),
    ("bio/bio-yeast", True, ' '),
    ("econ/econ-orani678", True, ' '),
    ("soc/soc-advogato", False, ' '),
    ("bio/bio-dmela", True, ' '),
    ("ia/ia-escorts-dynamic", False, ' '),
    ("soc/soc-anybeat", False, ' '),
    ("ca/ca-AstroPh", True, ' '),
    ("soc/fb-pages-media", False, ','),
    # ("bio/bio-diseasome", True, ' '),
    # ("ca/ca-cit-HepTh", False, ' '),
    # ("ca/ca-CondMat", True, ' '),
    # ("ca/ca-CSphd", False, ' '),
    # ("ca/ca-Erdos992", False, ' '),
    # ("ca/ca-GrQc", True, ' '),
    # ("ca/ca-HepPh", False, ' '),
    # ("econ/econ-beause", True, ' '),
    # ("econ/econ-mahindas", True, ' '),
    # ("econ/econ-poli-large", True, ' '),
    # ("inf/inf-openflights", False, ' '),
    # ("inf/inf-power", True, ' '),
    # ("soc/fb-pages-artist", False, ','),
    # ("soc/soc-gemsec-RO", False, ','),
    # ("soc/soc-gemsec-HU", False, ',')
]

RAW_FOLDER_PATH = "raw/network-repository"
FOLDER_PATH = "Graph"


def crawl_data():
    """
    Crawls data from https://networkrepository.com/.

    :return: None
    """
    if not os.path.exists(RAW_FOLDER_PATH):
        os.makedirs(RAW_FOLDER_PATH, exist_ok=True)

    for file in FILES:
        file_name = file[0]
        name = file[0].split('/')[1]

        if os.path.exists(f'{RAW_FOLDER_PATH}/{name}.edges'):
            continue

        r = requests.get(f'https://nrvis.com/download/data/{file_name}.zip')
        z = zipfile.ZipFile(io.BytesIO(r.content))
        z.extractall(f'{RAW_FOLDER_PATH}/{name}')

        if os.path.exists(f'{RAW_FOLDER_PATH}/{name}/{name}.edges'):
            os.rename(f'{RAW_FOLDER_PATH}/{name}/{name}.edges',
                      f'{RAW_FOLDER_PATH}/{name}.edges')

        if os.path.exists(f'{RAW_FOLDER_PATH}/{name}/{name}.mtx'):
            os.rename(f'{RAW_FOLDER_PATH}/{name}/{name}.mtx',
                      f'{RAW_FOLDER_PATH}/{name}.edges')

        shutil.rmtree(f'{RAW_FOLDER_PATH}/{name}/', ignore_errors=False, onerror=None)


def clean_data():
    """
    Cleans the data crawled from https://networkrepository.com/ and converts it
    into one unified format.
    - Keep only the largest component of the graph (size is measured in number of nodes)
    - make the graph undirected (if necessary)
    - delete duplicate edges
    - delete self-loops (edges that connect a node to itself)
    - remove all weights (if present)
    - Rename the remaining nodes, so they are labeled from 0 to n-1
    - Save the resulting graph in the data/network-repository/ directory

    :return: None
    """
    if not os.path.exists(FOLDER_PATH):
        os.mkdir(FOLDER_PATH)

    for file in FILES:
        name = file[0].split('/')[1]
        has_header = file[1]
        split_char = file[2]
        file_path = f'{RAW_FOLDER_PATH}/{name}.edges'
        graph = nx.Graph()

        file = open(file_path, 'r')
        lines = file.readlines()
        file.close()

        # remove all lines that start with "%"
        for i in reversed(range(len(lines))):
            if lines[i].startswith('%'):
                lines.pop(i)

        # pop the header
        if has_header:
            lines.pop(0)

        # read all edges and build graph
        for line in lines:
            line = line.strip()
            temp = line.split(split_char)

            n1, n2 = int(temp[0]), int(temp[1])
            if n1 != n2:
                graph.add_edge(n1, n2)

        # filter out largest component
        largest_cc = max(nx.connected_components(graph), key=len)
        subgraph = graph.subgraph(largest_cc)
        subgraph = nx.convert_node_labels_to_integers(subgraph, first_label=0, ordering='default', label_attribute=None)

        print(name)
        assert nx.is_connected(subgraph)
        mtx = dict(nx.all_pairs_shortest_path_length(subgraph))
        for i in range(subgraph.number_of_nodes()):
            for j in range(subgraph.number_of_nodes()):
                if i == j:
                    assert mtx[i][j] == 0
                else:
                    assert mtx[i][j] > 0

        # get all edges
        edges = list(subgraph.edges())

        # write edges to string
        s = ''
        for edge in edges:
            s += f'{edge[0]} {edge[1]}\n'

        # write to file
        file_path = f'{FOLDER_PATH}/{name}.edges'
        file = open(file_path, 'w')
        file.write(s)
        file.close()


if __name__ == '__main__':
    crawl_data()
    clean_data()
