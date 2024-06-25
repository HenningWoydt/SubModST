# Submodular Subset Maximization

This repository contains code to exactly solve the **Cardinality-Constrained Submodular Monotone Subset Maximization** problem.

Given a universe $\mathcal{U}$ consisting of $n$ arbitrary items. Let $f : 2^{\mathcal{U}} \to \mathbb{R}$ be a set function.    
Let $\Delta(e \mid A) \coloneqq f(A \cup \{e\}) - f(A)$ be the **marginal gain** of $e$.        
Set function $f$ is **submodular** if for every $A \subseteq B \subseteq \mathcal{U}$ and $e \in \mathcal{U} \setminus B$ it holds that $\Delta(e \mid A) \geq \Delta(e \mid B)$.      
Set function $f$ is **monotone (increasing)** if for every $A \subseteq \mathcal{U}$ and every $e \in \mathcal{U}$ it holds that $\Delta(e \mid A) \geq 0$.

#### Cardinality-Constrained Submodular Monotone Subset Maximization

Given a submodular and monotone set function $f$ and an integer $k$ determine a set $S$ with $|S| = k$ and      
$$S = \text{arg}\max_{S' \subseteq \mathcal{U}, |S'| = k} f(S').$$

## Available Functions

Currently, six functions can be optimized.

- Group Closeness Centrality: Given a graph $G=(V, E)$ determine a set $S \subseteq V$ of size $k$, such that $$f(S) \coloneqq \sum_{v \in V \setminus S} \min_{u \in S} \text{dist}(u, v)$$ is minimal. Since this originally is a minimization problem, we will use Negative Group Farness $-f(S)$ as the function to maximize.

- Partial Dominating Set: Given a graph $G=(V, E)$ determine a set $S \subseteq V$ of size $k$, such that $$f(S) \coloneqq |\bigcup_{v \in S} N[v]|$$ is maximized. $N[v] = \{v\} \cup \{u \mid \{u, v\} \in E\}$ is the closed neighborhood of $v$.

- $k$-Medoid Clustering: Given a set $X$ consisting of $n$ datapoints, each of dimensionality $d$, determine a set $S \subseteq X$ of size $k$ such that $$f(S) \coloneqq \sum_{i=1}^{n} \min_{x_j \in S} d(x_i, x_j)$$ is minimized. $$d(x, y) = \sqrt{ \sum_{i = 1}^{d} (x_i - y_i)^{2} }$$ is the euclidian distance. Like Group Farness, we will maximize $-f(S)$.

- Facility Location: Given a set of $n$ locations $N$ and a set of $m$ customers $M$. By $g_{ij} \geq 0$ we denote the benefit for customer $j$ when building a facility at location $i$. Determine a set $S \subseteq N$ of size $k$, such that $$f(S) = \sum_{j \in M} \max_{i \in S} g_{ij}$$ is maximized.

- Weighted Coverage: Description to be added

- Bipartite Influence: Description to be added

## Installation

### Prerequisites

Komogorv's [Blossom V](https://pub.ista.ac.at/~vnk/software.html#BLOSSOM5) implementation is required.

> Vladimir Kolmogorov. "Blossom V: A new implementation of a minimum
> cost perfect matching algorithm." In Mathematical Programming
> Computation (MPC), July 2009, 1(1):43-67.

Sadly its redistribution is prohibited, so we offer a script that will automatically download and extract it.
Move into directory `src/3rd_party/` and execute `setup_3rd_party.sh`.
This script will automatically download and extract the files to the correct folder.
At the end the folder `blossom5` should contain the downloaded code.

### Build
To build the binary use

```
cmake -G 'Unix Makefiles' -DCMAKE_BUILD_TYPE=Release -S . -B ./build
cmake --build build --target CCSMSM
```

The binary `CCSMSM` will be in the folder `build`.

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

The folder `utility` contains scripts to generate test data and download most of the benchmark data.

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

- Description To be added

#### Weighted Coverage

- Description To be added

#### Bipartite Influence

- Description to be added

## License

GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007

## Project status

Active development.
