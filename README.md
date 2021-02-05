[![](https://www.r-pkg.org/badges/version/backbone?color=orange)](https://cran.r-project.org/package=backbone)
[![](http://cranlogs.r-pkg.org/badges/grand-total/backbone?color=blue)](https://cran.r-project.org/package=backbone)
[![](http://cranlogs.r-pkg.org/badges/last-month/backbone?color=green)](https://cran.r-project.org/package=backbone)
<!-- README.md is generated from README.Rmd. Please edit that file -->

# backbone

<!-- badges: start -->

<!-- badges: end -->

## Welcome

Thank you for your interest in the backbone package\! Here you will find
a short example of how to use this package to extract the backbone of a
bipartite projection. For more details on these functions and methods,
please see our latest manuscript on backbone here:

“Domagalski R, Neal ZP, Sagan B (2021) Backbone: an R package for
extracting the backbone of bipartite projections. PLoS ONE 16(1):
e0244363.” <https://doi.org/10.1371/journal.pone.0244363>

For additional resources on how to use the backbone package, please see
<https://www.zacharyneal.com/backbone>

## The Backbone of a Graph

The `backbone` package provides methods for extracting from a weighted
graph a binary or signed backbone that retains only the significant
edges. The user may input a weighted graph, or a bipartite graph from
which a weighted graph is first constructed via projection. Backbone
extraction methods include the stochastic degree sequence model [(Neal,
Z. P. (2014))](https://doi.org/10.1016/j.socnet.2014.06.001),
hypergeometric model [(Neal, Z.
(2013))](https://doi.org/10.1007/s13278-013-0107-y), the fixed degree
sequence model [(Zweig, K. A., and Kaufmann, M.
(2011))](https://doi.org/10.1007/s13278-011-0021-0), as well as a
universal threshold method.

In a graph `G`, edges are either present (i.e. `G_{ij}=1`) or absent
(i.e. `G_{ij}=0`). However in a weighted or valued graph, edges can take
a range of values that may capture such properties as the strength or
capacity of the edge. Although weighted graphs contain a large amount of
information, there are some cases (e.g. visualization, application of
statistical models not developed for weighted graphs) where it is useful
to reduce this information by focusing on an unweighted subgraph that
contains only the most important edges. We call this subgraph the
backbone of `G`, which we denote as `G’`. Extracting `G’` from `G`
requires deciding which edges to preserve. This usually involves
selecting a threshold `T_{ij}` such that edges are preserved if they are
above the threshold (i.e. `G_{ij}’=1` if `G_{ij} > T_{ij}`), and omitted
if they are below the threshold (i.e. `G_{ij}’=0` if `G_{ij} < T_{ij}`).
It is also possible to extract a signed backbone by selecting upper
`T_{ij}` and lower `T’_{ij}` thresholds such that `G_{ij}’=1` if
`G_{ij}>T_{ij}`, `G_{ij}’=-1` if `G_{ij} < T’_{ij}`, and `G_{ij}’=0` if
`G_{ij} > T’_{ij}` and `G_{ij} < T_{ij}`. The key to all backbone
extraction methods lies in the selection of `T`. The `backbone` package
provides several different methods for selecting `T` and thus extracting
`G’` from `G`.

## Installation

You can install the released version of backbone from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("backbone")
```

You can install from GitHub with:

``` r
library(devtools)
install_github("domagal9/backbone", build_vignettes = TRUE)
```

## Example

This is a basic example which shows you how to solve a common problem
using the included Davis’ Southern Women dataset. Using the
`bipartite.null()` function, one can decide if the null model should be
constrained based on the row and/or column sums.

``` r
library(backbone)
#>  ____
#> |  _ \   backbone v1.3.0
#> |#|_) |  Cite: Domagalski, R., Neal, Z. P., & Sagan, B. (2021).
#> |# _ <   Backbone: An R package for extracting the backbone of bipartite projections.
#> |#|_) |  PLoS ONE. https://doi.org/10.1371/journal.pone.0244363
#> |____/   For help: type vignette("backbone"); email zpneal@msu.edu; github domagal9/backbone
data(davis)
null_model_probabilities <- bipartite.null(davis, rows = TRUE, cols = FALSE)
#> Finding the distribution using hypergeometric distribution
backbone <- backbone.extract(null_model_probabilities, signed = TRUE, alpha = 0.05)
```

For more detailed examples and background on the topic, see
`vignette("backbone_introduction", package = "backbone")` or our
manuscript on the backbone package:
<https://doi.org/10.1371/journal.pone.0244363>
