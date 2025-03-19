# colSBM 0.4.2

* Reworked clustering to have it de-recursified
* Added a clustering function based on graphon distance
* Added a `$cluster` vector to clustering output
* Removing `extract_best_partition()` and `extract_clustering_dendrogram()`

# colSBM 0.4.1

# colSBM 0.4.0

* Adding `colBiSBM`, supporting bipartite collections
* Changing parallelization support, using a `colsbm_lapply` generic function that the user can customize with their own backend (parallel, future)
* Adding `dorebipartite` a new dataset containing a bipartite collection dataset

# colSBM 0.3.0

* Improved generic `plot()` functions

# colSBM 0.2.1

* BREAKING CHANGE: emission distribution arguments are now named "distribution"
instead of "model"
* Fixed an error in `extract_best_partition()`
* Fixed a bug when using "spectal_init" as a `global_opts` arguement

# colSBM 0.2.0

* Added sparse matrix with the `Matrix` package leading to less intensive
  memory usage
* Added `extract_best_partition()` function
* Improved documentation, references and vignette

# colSBM 0.1.1

# colSBM 0.1.0

## Breaking changes

* Changed the name of `EstimateBMPOP()` to `clusterize_networks()`
* Most internal functions are not exported anymore

## Cosmetic changes

* Cleaned documentation and user functions
* Added a `NEWS.md` file to track changes to the package.
