---
title: "NEWS"
output: html_document
---

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
