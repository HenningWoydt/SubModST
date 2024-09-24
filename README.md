
# Submodular Subset Maximization

This repository contains code to exactly solve the **Cardinality-Constrained Submodular Monotone Subset Maximization** problem.

Given a universe $\mathcal{U}$ consisting of $n$ arbitrary items. Let $f : 2^{\mathcal{U}} \to \mathbb{R}$ be a set function.      
Let $\Delta(e \mid A) \coloneqq f(A \cup \{e\}) - f(A)$ be the **marginal gain** of $e$.          
A set function $f$ is **submodular** if for every $A \subseteq B \subseteq \mathcal{U}$ and $e \in \mathcal{U} \setminus B$ it holds that $\Delta(e \mid A) \geq \Delta(e \mid B)$.        
A set function $f$ is **monotone (increasing)** if for every $A \subseteq \mathcal{U}$ and every $e \in \mathcal{U}$ it holds that $\Delta(e \mid A) \geq 0$.

#### Cardinality-Constrained Submodular Monotone Subset Maximization

Given a submodular and monotone set function $f$ and an integer $k$ determine a set $S$ with $|S| = k$ and        
$$S = \text{arg}\max_{S' \subseteq \mathcal{U}, |S'| = k} f(S').$$

## Available Functions

Currently, six functions can be optimized.

#### Group Closeness Centrality
Given a graph $G=(V, E)$ determine a set $S \subseteq V$ of size $k$, such that $$f(S) \coloneqq \sum_{v \in V \setminus S} \min_{u \in S} \text{dist}(u, v)$$ is minimal. Since this originally is a minimization problem, we will use Negative Group Farness $-f(S)$ as the function to maximize.

#### Partial Dominating Set
Given a graph $G=(V, E)$ determine a set $S \subseteq V$ of size $k$, such that $$f(S) \coloneqq |\bigcup_{v \in S} N[v]|$$ is maximized. $N[v] = \{v\} \cup \{u \mid \{u, v\} \in E\}$ is the closed neighborhood of $v$.

#### $k$-Medoid Clustering
Given a set $X$ consisting of $n$ datapoints, each of dimensionality $d$, determine a set $S \subseteq X$ of size $k$ such that $$f(S) \coloneqq \sum_{i=1}^{n} \min_{x_j \in S} d(x_i, x_j)$$ is minimized. $$d(x, y) = \sqrt{ \sum_{i = 1}^{d} (x_i - y_i)^{2} }$$ is the euclidian distance. Like Group Farness, we will maximize $-f(S)$.

#### Facility Location
Given a set of $n$ locations $N$ and a set of $m$ customers $M$. By $g_{ij} \geq 0$ we denote the benefit for customer $j$ when building a facility at location $i$. Determine a set $S \subseteq N$ of size $k$, such that $$f(S) = \sum_{j \in M} \max_{i \in S} g_{ij}$$ is maximized.

#### Weighted Coverage
Given a collection $N = \{s_1, \ldots, s_n\}$ of subsets of an item set $M = \{1, \ldots, m\}$ and a weight function $\omega: M \to \mathbb{R}$, choose $S \subseteq N$ of size $k$ such that the summed weights are maximized, e.g. $$f(S) = \sum_{i \in \bigcup S} \omega(i). $$

#### Bipartite Influence
Let $N = \{1, \ldots, n\}$ be a set of sources and $M = \{1, \ldots, m\}$ a set of targets. Let $G=(N \cup M, E)$ be a bipartite graph with $E \subseteq N \times M$. Let $0 \leq p_{ij} \leq 1$ be the activation probability of edge $(j, i)$ for target $i \in M$ and source $j \in N$. If the edge does not exist in $G$, then $p_{ij} = 0$. A target $i \in M$ is activated by a set $S \subseteq N$ of source with probability $1 - \prod_{j \in S}(1 - p_{ij})$. Determine a set $S \subseteq N$ of size $k$ that maximizes $$f(S) = \sum_{i \in M} \Big(1 - \prod_{j \in S}(1 - p_{ij})\Big). $$

## Installation

### Prerequisites

Komogorv's [Blossom V](https://pub.ista.ac.at/~vnk/software.html#BLOSSOM5) implementation is required.

> Vladimir Kolmogorov. "Blossom V: A new implementation of a minimum  
> cost perfect matching algorithm." In Mathematical Programming  
> Computation (MPC), July 2009, 1(1):43-67.

Its redistribution is prohibited, so we offer a script that will automatically download and extract it.  
Move into directory `src/3rd_party/` and execute `setup_3rd_party.sh`.  
This script will automatically download and extract the files to the correct folder.  
At the end the folder `blossom5` should contain the downloaded code.

### Build
To build the binary use `./build.sh`. The binary `submodst` will be in the folder `build`.

## Usage

Use the following options:

- `-i [ --input-file  ]` Path to the file.
- `-f [ --function  ]` Which function to optimize.
- `-k [ --k  ]` The solution budget.
- `-o [ --output-file  ]` Path to the output file.
- `--initial-solution-file` (Optional) Path to a file holding an initial solution.
- `-t [ --time-limit  ]` (Optional) Maximum amount of running time (preprocessing not included).
- `-n [ --nickname  ]` Name of the algorithm. Currently `Plain`, `Simple`, `Simple+`, `LE-Rank`, `LE-Score`, `LE-RankOrScore`, `LE-RankAndScore`, `PWG-k^`, `PWG-Sqrt-n^`, `PWM-k^`, `PWM-Sqrt-n^`, `PWD-k^`, `PWD-10` are available. We recommend `LE-Score` as the fastest.

## Data

#### Graph

The file format should be

```  
e11 e12  
e21 e22  
...  
en1 en2  
```
with `ei1 ei2` denoting edge $i$ of the graph.

- Each $e_{ij}$ should be a positive integer $\geq 0$.
- The graph should be fully connected.
- The algorithm assumes that the vertex with the smallest ID is 0.
- Lines starting with a `%` are comments and will be ignored.

#### Datapoints

The file format should be

```
x11 x12 ... x1d
x21 x22 ... x2d
...
xn1 xn2 ... xnd
```   
each line denotes one point of the dataset.

- Each $x_{ij}$ should be a double value.
- Lines starting with a `%` are comments and will be ignored.

#### Facility Location

The file format should be
```
g11 g12 ... g1m
g21 g22 ... g2m
...
gn1 gn2 ... gnm
```
- Each entry $g_{ij}$ should be a non-negative double value for customer $i$ and location $j$.
- Lines starting with a `%` are comments and will be ignored.

#### Weighted Coverage

The file format should be
```
w1 w2 ... wm
b11 b12 ... b1m
...
bn1 bn2 ... bnm
```
- Each $w_j$ should be a non-negative double value, that denotes the weight of item $j$.
- Each entry $b_{ij}$ should be either 1.0 (item $j$ included in subset $i$) or 0.0 (not included).
- Lines starting with a `%` are comments and will be ignored.


#### Bipartite Influence

The file format should be
```
p11 p12 ... p1m
p21 p22 ... p2m
...
pn1 pn2 ... pnm
```
- Each entry $p_{ij}$ is the edge activation from source $i$ to target $j$. They should be double values in the range $[0, 1]$.
- Lines starting with a `%` are comments and will be ignored.


## License

GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007

## Project status

Active development.

## Bugs, Questions, Comments and Ideas

If any bugs arise, questions occur, comments want to be shared, or ideas discussed, please do not hesitate to contact the current repository owner (henning.woydt@informatik.uni-heidelberg.de) or leave a GitHub Issue or Discussion. Thanks!

## Reference

If you use this work in any academic work, please cite
```
@inproceedings{DBLP:conf/esa/WoydtKS24,
  author       = {Henning Martin Woydt and
                  Christian Komusiewicz and
                  Frank Sommer},
  editor       = {Timothy Chan and
                  Johannes Fischer and
                  John Iacono and
                  Grzegorz Herman},
  title        = {SubModST: {A} Fast Generic Solver for Submodular Maximization with
                  Size Constraints},
  booktitle    = {32nd Annual European Symposium on Algorithms, {ESA} 2024, September
                  2-4, 2024, Royal Holloway, London, United Kingdom},
  series       = {LIPIcs},
  volume       = {308},
  pages        = {102:1--102:18},
  publisher    = {Schloss Dagstuhl - Leibniz-Zentrum f{\"{u}}r Informatik},
  year         = {2024},
  url          = {https://doi.org/10.4230/LIPIcs.ESA.2024.102},
  doi          = {10.4230/LIPICS.ESA.2024.102},
  timestamp    = {Mon, 23 Sep 2024 12:27:20 +0200},
  biburl       = {https://dblp.org/rec/conf/esa/WoydtKS24.bib},
  bibsource    = {dblp computer science bibliography, https://dblp.org}
}
```
