test_that("frommatrix classes match", {
  library(Matrix)
  library(igraph)
  library(network)

  #NOTE: All matrices passed to frommatrix() will be binary or signed unipartite

  ### matrix (signed unipartite) --> Matrix
  M <- matrix(sample(c(-1,0,1),5*5,replace=TRUE),5,5)
  test <- frommatrix(M,"Matrix")
  expect_is(test, "Matrix")

  ### matrix (signed unipartite) --> Sparse Matrix
  M <- matrix(sample(c(-1,0,1),5*5,replace=TRUE),5,5)
  test <- frommatrix(M,"sparseMatrix")
  expect_is(test, "sparseMatrix")

  ### matrix (signed unipartite) --> network
  M <- matrix(sample(c(-1,0,1),5*5,replace=TRUE),5,5)
  test <- frommatrix(M,"network") #OK
  expect_is(test, "network")

  ### matrix (signed unipartite) --> igraph
  M <- matrix(sample(c(-1,0,1),5*5,replace=TRUE),5,5)
  test <- frommatrix(M,"igraph") #OK
  expect_is(test, "igraph")

  ### matrix (signed unipartite) --> edgelist
  M <- matrix(sample(c(-1,0,1),5*5,replace=TRUE),5,5)
  test <- frommatrix(M,"edgelist") #OK
  expect_is(test, "data.frame")

})
