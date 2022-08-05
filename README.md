
<!-- README.md is generated from README.Rmd. Please edit that file -->

# colSBM

<!-- badges: start -->

[![](https://img.shields.io/badge/devel%20version-0.1.0-green.svg)](https://github.com/colSBM)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![](https://img.shields.io/github/last-commit/Chabert-Liddell/colSBM.svg)](https://github.com/Chabert-Liddell/colSBM/commits/main)
[![](https://img.shields.io/badge/doi-10.48550/arXiv.2206.00560-blue.svg)](https://doi.org/10.48550/arXiv.2206.00560)
<!-- badges: end -->

colSBM is an R package which implements methods for clustering and
inferring the mesoscale structure for collection of networks. In
particular, it allows to:

-   Find a common mesoscale structure in a collection of networks by
    using a Stochastic Block Model (SBM) extension for the joint
    modeling of collection of networks.

-   Provide a clustering of the nodes of each networks. Classifies
    networks in groups of networks with similar mesoscale structures.

Mathematical details of the methods as well as simulation studies and
applications can be find in Chabert-Liddell, Barbillon, and Donnet
(2022) .

## Installation

You can install the development version of colSBM from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Chabert-Liddell/colSBM")
```

# References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-collection" class="csl-entry">

Chabert-Liddell, Saint-Clair, Pierre Barbillon, and Sophie Donnet. 2022.
“Learning Common Structures in a Collection of Networks. An Application
to Food Webs.” *arXiv Preprint arXiv:2206.00560*.

</div>

</div>
