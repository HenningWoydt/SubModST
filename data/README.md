# How to obtain data

## User Data
Use the script `crawl_all.sh` in to crawl most of the benchmarking data.
The data from Eszter Julianna Csókás and Tamás Vinkó is not publicly available.

- [networkrepository.com](https://networkrepository.com/) for the graphs. The resulting files will be present in `data/Graph/`. The graphs were used in
> Luca Pascal Staus, Christian Komusiewicz, Nils Morawietz, and Frank Sommer.
> Exact algorithms for group closeness centrality.
> In SIAM Conference on Applied and Computational Discrete Algorithms (ACDA ’23), pages 1–12. SIAM, 2023.

- [sites.google.com/view/umepon/benchmark](https://sites.google.com/view/umepon/benchmark) for LOC, COV and INF instances.  The resulting files will be present in `data/FacilityLocation/`, `data/WeightedCoverage/` and `data/BipartiteInfluence/`. The files were used in
> Naoya Uematsu, Shunji Umetani, and Yoshinobu Kawahara.
> An efficient branch-and-cut algorithm for submodular function maximization.
> Journal of the Operations Research Society of Japan, 63(2):41–59, 2020.

- [cs.uef.fi/sipu/datasets/](https://cs.uef.fi/sipu/datasets/) for clustering instances. The resulting files will be present in `data/clustering/`. The files were used in
> Pasi Fränti and Sami Sieranoja.
> K-means properties on six clustering benchmark datasets.
> Applied Intelligence, 48(12):4743–4759, 2018.

> Mohammad Rezaei and Pasi Fränti.
> Can the number of clusters be determined by external indices?
> IEEE Access, 8:89239–89257, 2020.

- Ask Csókás and Vinkó for more INF instances. Place the folder `INF_new_data_all` into `data/raw/` and then run the `python crawl_INF_csokas` in `utility`.
> Eszter Julianna Csókás, Tamás Vinkó.
> Constraint generation approaches for submodular function maximization leveraging graph properties.
> J. Glob. Optim. 88(2): 377-394 (2024)

## Developer Data
Use the script `generate_all.sh`  in `utility` to create random test instances for each function.
The resulting files will be present in
- `data/private/Graph` for graphs
- `data/private/Clustering` for clustering instances
- `data/private/FacilityLocation` for LOC instances
- `data/private/WeightedCoverage` for COV instances
- `data/private/BipartiteInfluence` for INF instances