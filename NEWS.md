---
title: "NEWS"
output: html_document
---
## backbone 2.0.2

* fixed bug in `fastball()` so it will work with R < 4.1.0

## backbone 2.0.1

* minor bug fixes
* faster implementation of `fastball()` algorithm
* set `alpha = 0.05` as default in all statistical models
* renamed `fwer` (familywise error rate) parameter as `mtc` (multiple test correction)

## backbone 2.0.0

* remove `davis` example data; add examples using synthetic data
* add support for unweighted graphs: `sparsify()`
* add support for weighted bipartite graphs: `osdsm()`
* add support for non-projection weighted graphs: `disparity()`
* new vignette illustrating all functions
* add implementation of `fastball()` algorithm for marginal-preserving matrix randomization
* re-add `testthat` tests
* allow backbone functions to directly output a backbone, eliminating the need for the `backbone.extract()` function
* add support for any `p.adjust()` method of correcting for familywise error rates
* Minor bug fixes

## backbone 1.5.1

* removed `testthat` tests due to unknown MKL error; will be restored in a future version

## backbone 1.5.0

* add four functions to generate random bipartite graphs: bipartite.from.probability(), bipartite.from.sequence(), bipartite.from.distribution(), and bipartite.add.blocks()
* set diagonal in `positive` and `negative` backbone object matrices to NA
* corrected p-value computation in fixedfill()
* remove running time from backbone object summary dataframe
* update documentation, readme, vignette

## backbone 1.4.0

* add fixedcol() function - null model where column degrees are fixed and row sums are allowed to vary
* add fixedfill() function - null model where the number of 1's in the matrix (number of edges in the graph) are fixed
* replace class.convert() with tomatrix() and frommatrix()
* use updated Poisson binomial calculations (more accurate approximation)
* hyperg() now called fixedrow()
* remove bipartite.null function
* update documentation, readme, vignette
* include logo

## backbone 1.3.1

* speedups to sdsm

## backbone 1.3.0

* update sdsm to use the bicm model - a new, fast, approximation of the probabilities
* remove all other models from sdsm
* if an older model is called in sdsm, show warning that model has changed
* add new function bipartite.null which lets the user pick if they want rows/cols to be fixed or vary
* update fwer m parameter

## backbone 1.2.2

* fix fdsm to accept all graph inputs
* rename sdsm "chi2" model to "rcn"
* universal function can now return weighted projection
* universal function now has a narrative parameter
* class.convert now drops (with warning) rows and columns with zero sum before sending output to universal, sdsm, fdsm, or hyperg.
* update citations

## backbone 1.2.1

* add narrative parameter to backbone.extract for suggested manuscript text
* add scobit model to sdsm
* add time unit to runtime calculation
* minor spelling and comment fixes

## backbone 1.2.0

* add support for sparse matrix, igraph, network, and edgelist objects (see 'class.convert')
* add family-wise error rate test corrections (see 'backbone.extract')
* sdsm: add multiple methods for computing initial probabilities (see 'sdsm' details) one of which uses convex optimization (see 'polytope')
* sdsm: update poisson binomial computation method to increase speed (see 'sdsm' and 'rna')
* add more descriptives to summary dataframe output of backbone object
* update documentation of functions
* update vignette to reflect package changes
* bug fixes for R 4.0.0

## backbone 1.1.0

* add support for sparse matrices
* add support for speedglm in sdsm
* add poisson binomial approx. in sdsm
* add summary output to sdsm, fdsm, hyperg, universal
* update vignette to reflect package changes
* bug fixes

## backbone 1.0.0

* initial release
