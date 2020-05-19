test_that("class.convert classes match", {
  library(Matrix)
  library(igraph)
  library(network)

  #normal matrix
  class(davis) #"matrix"
  G <- davis%*%t(davis)
  class(G)

  #sparse matrix
  davis.sp <- as(davis, "sparseMatrix")
  class(davis.sp) #"Matrix"
  methods::is(davis.sp, "sparseMatrix") #TRUE
  G.sp <- as(G, "sparseMatrix")
  class(G.sp)

  #igraph object
  davis.graph <- graph.incidence(davis)
  class(davis.graph) #"igraph
  igraph::is.bipartite(davis.graph) #TRUE
  G.graph <- graph.adjacency(G)

  #edge list
  davis.el <- igraph::get.edgelist(davis.graph)
  class(davis.el)
  methods::is(davis.el, "edge list")
  dim(davis.el)
  G.el <- igraph::get.edgelist(G.graph)

  #network object
  davis.nw <- network::network(davis, ignore.eval = FALSE, names.eval = "weight", loops = TRUE)
  class(davis.nw)
  G.nw <- network::network(G, ignore.eval = FALSE, names.eval = "weight", directed = FALSE, loops = TRUE)
  davis.sn <- network::network(davis, bipartite = TRUE)

  ## test file

  a1<-class.convert(davis)
  a2<-class.convert(davis.sp)
  a3<-class.convert(davis.graph)
  a4<-class.convert(davis.nw)
  a5<-class.convert(davis.el)
  a6<-class.convert(davis.sn)

  b1<-a1[[2]]
  b2<-a2[[2]]
  b3<-a3[[2]]
  b4<-a4[[2]]
  b5<-a5[[2]]
  b6<-a6[[2]]

  expect_equal(davis,b1)
  expect_equal(davis.sp, b2)
  expect_equal(davis, b3)
  expect_equal(davis, b4)
  expect_equal(davis[,order(colnames(davis))], b5[, order(colnames(b5))])
  expect_equal(davis, b6)
  expect_equal(b4,b6)

  c1<-class.convert(G)
  c2<-class.convert(G.sp)
  c3<-class.convert(G.graph)
  c4<-class.convert(G.nw)

  d1<-c1[[2]]
  d2<-c2[[2]]
  d3<-c3[[2]]
  d4<-c4[[2]]

  expect_equal(G,d1)
  expect_equal(G.sp,d2)
  expect_equal(G.sp,d3)
  expect_equal(G, d4)

  ######## 2nd Function #################
  bb.sdsm <- sdsm(davis)
  bb <- backbone.extract(bb.sdsm, signed = TRUE, alpha = .2)
  bb.pos <- backbone.extract(bb.sdsm, signed = FALSE, alpha = .1)
  bb.sp <- as(bb, "sparseMatrix")
  bb.pos.sp <- as(bb.pos, "sparseMatrix")
  class(davis)
  class(bb)

  e1 <- class.convert(bb, "matrix")
  e2 <- class.convert(bb, "sparseMatrix")
  e3 <- class.convert(bb, "igraph")
  e4 <- class.convert(bb, "network")

  f1<-class.convert(e1[[2]])
  f2<-class.convert(e2[[2]])
  f3<-class.convert(e3[[2]])
  f4<-class.convert(e4[[2]])

  expect_equal(bb,f1[[2]])
  expect_equal(bb.sp,f2[[2]])
  expect_equal(bb.sp,f3[[2]])
  expect_equal(bb,f4[[2]])

  ##################################
  bb.sdsm <- sdsm(davis.sp)
  bb <- backbone.extract(bb.sdsm, signed = TRUE, alpha = .2)
  bb.sp <- as(bb, "sparseMatrix")
  expect_true(methods::is(bb, "sparseMatrix"))
  expect_equal(bb, bb.sp)

  ##################################
  bb.sdsm <- sdsm(davis.graph)
  bb <- backbone.extract(bb.sdsm, signed = TRUE, alpha = .2)
  expect_equal(class(davis.graph), class(bb))

  ##################################
  bb.sdsm <- sdsm(davis.nw)
  bb <- backbone.extract(bb.sdsm, signed = TRUE, alpha = .2)
  expect_equal(class(davis.nw), class(bb))
})
