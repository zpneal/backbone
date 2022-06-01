library(backbone)

test_that("BiCM", {
  M <- rbind(c(0,0,1),c(0,1,0),c(1,0,1))
  test <- round(bicm(M),3)
  expect_equal(test, rbind(c(.216,.216,.568),c(.216,.216,.568),c(.568,.568,.863)))
})

test_that("Fastball", {
  M <- matrix(rbinom(100*1000,1,0.5),100,1000)
  test <- fastball(M)
  expect_equal(rowSums(test), rowSums(M))
  expect_equal(colSums(test), colSums(M))
})

test_that("SDSM output", {
  M <- rbind(c(1,0,1,1),c(0,1,0,0),c(1,0,0,1))
  test <- sdsm(M, alpha = NULL)
  expect_equal(test$G, M%*%t(M))  #Weighted projection
  expect_equal(round(test$Pupper,3), rbind(c(.432,.977,.549),c(.977,.304,.972),c(.549,.972,.329)))
  expect_equal(round(test$Plower,3), rbind(c(.909,.354,.844),c(.354,.954,.476),c(.844,.476,.942)))
  expect_equal(test$model, "sdsm")
})

test_that("FDSM output", {
  set.seed(1)
  M <- rbind(c(1,0,1,1),c(0,1,0,0),c(1,0,0,1))
  test <- fdsm(M, trials = 1000, alpha = NULL)
  expect_equal(test$G, M%*%t(M))  #Weighted projection
  expect_equal(round(test$Pupper,3), rbind(c(1,1,.265),c(1,1,1),c(.265,1,1)))
  expect_equal(round(test$Plower,3), rbind(c(1,.510,1),c(.51,1,.755),c(1,.755,1)))
  expect_equal(test$model, "fdsm")
})

test_that("FixedFill output", {
  M <- rbind(c(1,0,1,1),c(0,1,0,0),c(1,0,0,1))
  test <- fixedfill(M, alpha = NULL)
  expect_equal(test$G, M%*%t(M))  #Weighted projection
  expect_equal(round(test$Pupper,3), rbind(c(.004,1,.173),c(1,.732,1),c(.173,1,.173)))
  expect_equal(round(test$Plower,3), rbind(c(1,.268,.996),c(.268,.827,.268),c(.996,.268,.996)))
  expect_equal(test$model, "fixedfill")
})

test_that("FixedRow output", {
  M <- rbind(c(1,0,1,1),c(0,1,0,0),c(1,0,0,1))
  test <- fixedrow(M, alpha = NULL)
  expect_equal(test$G, M%*%t(M))  #Weighted projection
  expect_equal(round(test$Pupper,3), rbind(c(.25,1,.5),c(1,.25,1),c(.5,1,.167)))
  expect_equal(round(test$Plower,3), rbind(c(1,.25,1),c(.25,1,.5),c(1,.5,1)))
  expect_equal(test$model, "fixedrow")
})

test_that("FixedColumn output", {
  M <- rbind(c(1,0,1,1),c(0,1,0,0),c(1,0,0,1))
  test <- fixedcol(M, alpha = NULL)
  expect_equal(test$G, M%*%t(M))  #Weighted projection
  expect_equal(round(test$Pupper,3), rbind(c(.008,.975,.114),c(.975,.568,.975),c(.114,.975,.114)))
  expect_equal(round(test$Plower,3), rbind(c(1,.432,.992),c(.432,.886,.432),c(.992,.432,.992)))
  expect_equal(test$model, "fixedcol")
})

test_that("Disparity Filter output", {
  M <- rbind(c(1,0,1,1),c(0,1,0,0),c(1,0,0,1),c(1,1,0,1))
  test <- disparity(M, alpha = NULL)
  expect_equal(test$G, M)  #Weighted projection
  expect_equal(round(test$Pupper,3), rbind(c(.444,1,.444,.444),c(1,.5,1,1),c(.444,1,1,.444),c(.444,.444,1,.444)))
  expect_equal(round(test$Plower,3), rbind(c(.556,1,.556,.556),c(1,.5,1,1),c(.556,1,1,.556),c(.556,.556,1,.556)))
  expect_equal(test$model, "disparity")
})

test_that("frommatrix classes match", {
  #NOTE: All matrices passed to frommatrix() will be binary or signed unipartite
  # matrix (signed unipartite) --> Matrix
  M <- matrix(sample(c(-1,0,1),5*5,replace=TRUE),5,5)
  test <- frommatrix(M,convert="Matrix")
  expect_s4_class(test, "Matrix")

  # matrix (signed unipartite) --> igraph
  M <- matrix(sample(c(-1,0,1),5*5,replace=TRUE),5,5)
  test <- frommatrix(M,convert="igraph") #OK
  expect_s3_class(test, "igraph")

  # matrix (signed unipartite) --> edgelist
  M <- matrix(sample(c(-1,0,1),5*5,replace=TRUE),5,5)
  test <- frommatrix(M,convert="edgelist") #OK
  expect_s3_class(test, "data.frame")
})
