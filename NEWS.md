---
title: "NEWS"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## backbone 1.2.0

* add support for igraph, network, and edgelist objects (see 'class.convert')
* add family-wise error rate test corrections in backbone.extract (see 'holm.bonferroni')
* sdsm: add multiple methods for computing initial probabilities (see 'sdsm' details and 'polytope'), update poisson binomial computation method (removed poibin dependency, see 'rna')
* add more descriptives to summary output
* update documentation

## backbone 1.1.0

* add support for sparse matrices
* add support for speedglm in sdsm
* add poisson binomial approx. in sdsm
* add summary output to sdsm, fdsm, hyperg, universal
* update vignette to reflect package changes
* bug fixes

## backbone 1.0.0

* initial release
