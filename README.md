
<!-- README.md is generated from README.Rmd. Please edit that file -->

# colSBM

<!-- badges: start -->

[![](https://img.shields.io/badge/devel%20version-0.4.6-green.svg)](https://github.com/GrossSBM/colSBM)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![](https://img.shields.io/github/last-commit/GrossSBM/colSBM.svg)](https://github.com/GrossSBM/colSBM/commits/main)
[![](https://badgen.net/badge/DOI/10.1214%2F23-AOAS1831/yellow)](https://doi.org/10.1214/23-AOAS1831)
[![](https://www.r-pkg.org/badges/version/colSBM)](https://CRAN.R-project.org/package=colSBM)
[![Codecov test
coverage](https://codecov.io/gh/GrossSBM/colSBM/graph/badge.svg)](https://app.codecov.io/gh/GrossSBM/colSBM)
<!-- badges: end -->

colSBM is an R package which implements methods for clustering and
inferring the mesoscale structure for collection of networks. In
particular, it allows to:

- Find a common mesoscale structure in a collection of networks by using
  a Stochastic Block Model (SBM) extension for the joint modeling of
  collection of networks.

- Provide a clustering of the nodes of each networks. Classifies
  networks in groups of networks with similar mesoscale structures.

Mathematical details of the methods as well as simulation studies and
applications can be find in Chabert-Liddell, Barbillon, and Donnet
(2024) .

## Installation

You can install the development version of colSBM from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("GrossSBM/colSBM")
```

# References

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-collection" class="csl-entry">

Chabert-Liddell, Saint-Clair, Pierre Barbillon, and Sophie Donnet. 2024.
“Learning Common Structures in a Collection of Networks. An Application
to Food Webs.” *The Annals of Applied Statistics* 18 (2): 1213–35.
<https://doi.org/10.1214/23-AOAS1831>.

</div>

</div>
