# backbone <img src='man/figures/logo.png' align="right" height="139" />

<!-- badges: start -->

[![](https://www.r-pkg.org/badges/version/backbone?color=orange)](https://cran.r-project.org/package=backbone)
[![](http://cranlogs.r-pkg.org/badges/grand-total/backbone?color=blue)](https://cran.r-project.org/package=backbone)
[![](http://cranlogs.r-pkg.org/badges/last-month/backbone?color=green)](https://cran.r-project.org/package=backbone)
<!-- badges: end -->

## Welcome
Thank you for your interest in the backbone package\! The backbone package implements methods to extract the *backbone* of a network, which is a sparse and unweighted subgraph that contains only the most ‘important’ or ‘significant’ edges. A backbone can be useful when the original network is too dense, when edge weights are not needed, or when edge weights are difficult to interpret. Methods are available for:

* Weighted bipartite projections
* Non-projection weighted networks
* Unweighted networks

In addition, the package implements some other utility functions to:

* Generate random bipartite networks
* Randomize matrices while preserving the row and column sums
* Estimate the Bipartite Configuration Model (BiCM)

For more details on these functions and methods, please see:

* Domagalski R., Neal, Z.P., and Sagan, B. (2021). Backbone: an R package for extracting the backbone of bipartite projections. *PLoS ONE, 16*, e0244363. <https://doi.org/10.1371/journal.pone.0244363>
* Neal, Z.P., Domagalski, R., and Sagan, B. (2021). Comparing Alternatives to the Fixed Degree Sequence Model for Extracting the Backbone of Bipartite Projections. *Scientific Reports*. <https://arxiv.org/abs/2105.13396>
* [www.rbackbone.net](https://www.zacharyneal.com/backbone)

## Installation
The /master branch contains the current CRAN release of the backbone package. You can install it from [CRAN](https://CRAN.R-project.org) with:
``` r
install.packages("backbone")
```

The /develop branch contains the working beta version of the next release of the backbone package. All the functions are documented and have undergone various levels of preliminary debugging, so they should mostly work, but there are no guarantees. Feel free to use the develop version (with caution), and let us know if you run into any problems. You can install it You can install from GitHub with:
``` r
library(devtools)
install_github("zpneal/backbone", ref = "develop", build_vignettes = TRUE)
```

*NOTE: Prior to backbone v2.0.0, the code was maintained on a different GitHub repo: https://github.com/domagal9/backbone. That repo is still available, but is no longer updated. All future updates to backbone will be maintained here.*
