#### Exported Utility Functions ####
## BICM
M <- rbind(c(0,0,1),c(0,1,0),c(1,0,1))
test <- round(bicm(M),3)
expect_equal(test, rbind(c(.216,.216,.568),c(.216,.216,.568),c(.568,.568,.863)), info = "bicm")

## FASTBALL
M <- matrix(rbinom(100*1000,1,0.5),100,1000)
test <- fastball(M)
expect_equal(rowSums(test), rowSums(M), info = "fastball rows")
expect_equal(colSums(test), colSums(M), info = "fastball columns")

## PB (poisson-binomial)
test <- pb(5, c(.123,.234,.345,.456,.567,.678,.789,.801,.911))
expect_equal(test[1], 0.6773302, info = "pb lower")
expect_equal(test[2], 0.6268476, info = "pb upper")

#### Internal Utility Functions ####
## TOMATRIX: convert from matrix
M <- rbind(c(0,0,1),c(0,1,0),c(1,0,1))
test <- backbone:::tomatrix(M)
expect_equal(test$summary$class,"matrix", info = "unipartite matrix")
expect_false(test$summary$bipartite, info = "unipartite matrix")
expect_true(test$summary$symmetric, info = "unipartite matrix")
expect_false(test$summary$weighted, info = "unipartite matrix")
expect_equal(test$G,M, info = "unipartite matrix")

M <- rbind(c(0,0,1),c(1,1,0),c(1,0,1),c(2,5,7))
test <- backbone:::tomatrix(M)
expect_equal(test$summary$class,"matrix", info = "bipartite matrix")
expect_true(test$summary$bipartite, info = "bipartite matrix")
expect_false(test$summary$symmetric, info = "bipartite matrix")
expect_true(test$summary$weighted, info = "bipartite matrix")
expect_equal(test$G,M, info = "bipartite matrix")

M <- rbind(c(0,0,1),c(1,1,0),c(1,0,1))
rownames(M) <- c("A", "B", "C")
colnames(M) <- c("D", "E", "F")
test <- backbone:::tomatrix(M)
expect_equal(test$summary$class,"matrix", info = "square bipartite matrix")
expect_true(test$summary$bipartite, info = "square bipartite matrix")
expect_false(test$summary$symmetric, info = "square bipartite matrix")
expect_false(test$summary$weighted, info = "square bipartite matrix")
expect_equal(test$G,M, info = "square bipartite matrix")

## TOMATRIX: convert from Matrix
M <- rbind(c(0,0,1),c(0,1,0),c(1,0,1))
M <- Matrix::Matrix(M)
test <- backbone:::tomatrix(M)
expect_equal(test$summary$class,"Matrix", info = "unipartite Matrix")
expect_false(test$summary$bipartite, info = "unipartite Matrix")
expect_true(test$summary$symmetric, info = "unipartite Matrix")
expect_false(test$summary$weighted, info = "unipartite Matrix")
expect_equal(test$G,as.matrix(M), info = "unipartite Matrix")

M <- rbind(c(0,0,1),c(1,1,0),c(1,0,1),c(2,5,7))
M <- Matrix::Matrix(M)
test <- backbone:::tomatrix(M)
expect_equal(test$summary$class,"Matrix", info = "bipartite Matrix")
expect_true(test$summary$bipartite, info = "bipartite Matrix")
expect_false(test$summary$symmetric, info = "bipartite Matrix")
expect_true(test$summary$weighted, info = "bipartite Matrix")
expect_equal(test$G,as.matrix(M), info = "bipartite Matrix")

## TOMATRIX: convert from igraph
M <- igraph::erdos.renyi.game(10, 0.5, type="gnp", directed=FALSE)
igraph::V(M)$names <- LETTERS[1:10]
test <- backbone:::tomatrix(M)
expect_equal(test$summary$class,"igraph", info = "unipartite igraph")
expect_false(test$summary$bipartite, info = "unipartite igraph")
expect_true(test$summary$symmetric, info = "unipartite igraph")
expect_false(test$summary$weighted, info = "unipartite igraph")
expect_equal(test$G,igraph::as_adjacency_matrix(M, sparse = F), info = "unipartite igraph")
expect_equal(as.vector(test$attribs$names),LETTERS[1:10], info = "unipartite igraph")

M <- igraph::bipartite.random.game(n1 = 10, n2 = 100, type = "gnm", m = 250, directed=FALSE)
igraph::V(M)$names <- LETTERS[1:10]
igraph::E(M)$weight <- runif(250)
test <- backbone:::tomatrix(M)
expect_equal(test$summary$class,"igraph", info = "bipartite igraph")
expect_true(test$summary$bipartite, info = "bipartite igraph")
expect_false(test$summary$symmetric, info = "bipartite igraph")
expect_true(test$summary$weighted, info = "bipartite igraph")
expect_equal(test$G,igraph::as_incidence_matrix(M, sparse = F, attr = 'weight'), info = "bipartite igraph")
expect_equal(as.vector(test$attribs$names),LETTERS[1:10], info = "bipartite igraph")

## TOMATRIX: convert from edgelist
M <- data.frame(i = c("A", "B", "C"), j = c("B", "A", "A"), z = c(1,5,9))
test <- backbone:::tomatrix(M)
expect_equal(test$summary$class,"edgelist", info = "unipartite edgelist")
expect_false(test$summary$bipartite, info = "unipartite edgelist")
expect_false(test$summary$symmetric, info = "unipartite edgelist")
expect_true(test$summary$weighted, info = "unipartite edgelist")
expect_equal(test$G,
             matrix(c(0,5,9,1,0,0,0,0,0),3,3,dimnames=list(c("A","B","C"),c("A","B","C"))),
             info = "unipartite edgelist")

M <- data.frame(i = c("A", "B", "C"), j = c("D", "E", "F"))
test <- backbone:::tomatrix(M)
expect_equal(test$summary$class,"edgelist", info = "bipartite edgelist")
expect_true(test$summary$bipartite, info = "bipartite edgelist")
expect_false(test$summary$symmetric, info = "bipartite edgelist")
expect_false(test$summary$weighted, info = "bipartite edgelist")
expect_equal(test$G,
             matrix(c(1,0,0,0,1,0,0,0,1),3,3,dimnames=list(c("A","B","C"),c("D","E","F"))),
             info = "bipartite edgelist")

## FROMMATRIX: Convert to matrix
M <- rbind(c(-1,-1,1),c(0,0,0),c(1,0,0))
test <- backbone:::frommatrix(M, convert = "matrix")
expect_equal(test,M, info = "convert to matrix")

## FROMMATRIX: Convert to Matrix
M <- rbind(c(-1,-1,1),c(0,0,0),c(1,0,0))
test <- backbone:::frommatrix(M, convert = "Matrix")
expect_equal(test,Matrix::Matrix(M), info = "convert to Matrix")

## FROMMATRIX: Convert to Edgelist
M <- rbind(c(-1,-1,1),c(0,0,0),c(1,0,0))
test <- backbone:::frommatrix(M, convert = "edgelist")
df <- data.frame(from = c(1,1,1,3), to = c(1,2,3,1), weight = c(-1,-1,1,1))
expect_equal(test,df, info = "convert to edgelist")

## FROMMATRIX: Convert to igraph
M <- rbind(c(-1,-1,1),c(0,0,0),c(1,0,0))
attrib <- data.frame(name = c("A", "B", "C"), gender = c("M", "F", "F"))
test <- backbone:::frommatrix(M, attribs = attrib, convert = "igraph")
expect_true(methods::is(test,"igraph"), info = "convert to igraph")
expect_equal(igraph::as_adjacency_matrix(test, sparse = F, attr = 'weight'),
             matrix(c(-1,0,1,-1,0,0,1,0,0),3,3,dimnames=list(c("A","B","C"),c("A","B","C"))),
             info = "convert to igraph")
expect_equal(igraph::V(test)$name, c("A", "B", "C"))
expect_equal(igraph::V(test)$gender, c("M", "F", "F"))
expect_equal(igraph::E(test)$weight, c(-1,-1,1,1))
expect_equal(igraph::E(test)$sign, c(-1,-1,1,1))

#### Bipartite Bipartite Models ####
## SDSM
M <- rbind(c(1,0,1,1),c(0,1,0,0),c(1,0,0,1))
test <- sdsm(M, alpha = NULL, signed = TRUE)
expect_equal(test$G, M%*%t(M), info = "sdsm")
expect_equal(round(test$Pupper,3), rbind(c(.432,.977,.549),c(.977,.304,.972),c(.549,.972,.329)), info = "sdsm")
expect_equal(round(test$Plower,3), rbind(c(.909,.354,.844),c(.354,.954,.476),c(.844,.476,.942)), info = "sdsm")
expect_equal(test$model, "sdsm", info = "sdsm")

## FDSM
set.seed(1)
M <- rbind(c(1,0,1,1),c(0,1,0,0),c(1,0,0,1))
test <- fdsm(M, trials = 1000, alpha = NULL, signed = TRUE)
expect_equal(test$G, M%*%t(M), info = "fdsm")
expect_equal(round(test$Pupper,3), rbind(c(1,1,.265),c(1,1,1),c(.265,1,1)), info = "fdsm")
expect_equal(round(test$Plower,3), rbind(c(1,.510,1),c(.51,1,.755),c(1,.755,1)), info = "fdsm")
expect_equal(test$model, "fdsm", info = "fdsm")

## FIXEDFILL
M <- rbind(c(1,0,1,1),c(0,1,0,0),c(1,0,0,1))
test <- fixedfill(M, alpha = NULL, signed = TRUE)
expect_equal(test$G, M%*%t(M), info = "fixedfill")
expect_equal(round(test$Pupper,3), rbind(c(.004,1,.173),c(1,.732,1),c(.173,1,.173)), info = "fixedfill")
expect_equal(round(test$Plower,3), rbind(c(1,.268,.996),c(.268,.827,.268),c(.996,.268,.996)), info = "fixedfill")
expect_equal(test$model, "fixedfill", info = "fixedfill")

## FIXEDROW
M <- rbind(c(1,0,1,1),c(0,1,0,0),c(1,0,0,1))
test <- fixedrow(M, alpha = NULL, signed = TRUE)
expect_equal(test$G, M%*%t(M), info = "fixedrow")
expect_equal(round(test$Pupper,3), rbind(c(.25,1,.5),c(1,.25,1),c(.5,1,.167)), info = "fixedrow")
expect_equal(round(test$Plower,3), rbind(c(1,.25,1),c(.25,1,.5),c(1,.5,1)), info = "fixedrow")
expect_equal(test$model, "fixedrow", info = "fixedrow")

## FIXEDCOL
M <- rbind(c(1,0,1,1),c(0,1,0,0),c(1,0,0,1))
test <- fixedcol(M, alpha = NULL, signed = TRUE)
expect_equal(test$G, M%*%t(M), info = "fixedcol")
expect_equal(round(test$Pupper,3), rbind(c(.008,.975,.114),c(.975,.568,.975),c(.114,.975,.114)), info = "fixedcol")
expect_equal(round(test$Plower,3), rbind(c(1,.432,.992),c(.432,.886,.432),c(.992,.432,.992)), info = "fixedcol")
expect_equal(test$model, "fixedcol", info = "fixedcol")

#### Weighted Unipartite Models ####
## DISPARITY
M <- rbind(c(1,0,1,1),c(0,1,0,0),c(1,0,0,1),c(1,1,0,1))
test <- disparity(M, alpha = NULL, signed = TRUE)
expect_equal(test$G, M, info = "disparity")
expect_equal(round(test$Pupper,3), rbind(c(.444,1,.444,.444),c(1,.5,1,1),c(.444,1,1,.444),c(.444,.444,1,.444)), info = "disparity")
expect_equal(round(test$Plower,3), rbind(c(.556,1,.556,.556),c(1,.5,1,1),c(.556,1,1,.556),c(.556,.556,1,.556)), info = "disparity")
expect_equal(test$model, "disparity", info = "disparity")

M <- rbind(c(0,2,3,4),c(1,0,5,2),c(8,2,0,2),c(1,7,2,0))
test <- global(M, upper = 5, lower = 2)
expect_equal(test, rbind(c(0,0,0,0), c(-1,0,0,0), c(1,0,0,0), c(-1,1,0,0)), info = "global threshold")

M <- rbind(c(0,2,3,4),c(1,0,5,2),c(8,2,0,2),c(1,7,2,0))
test <- global(M, upper = function(x)mean(x))
expect_equal(test, rbind(c(0,0,1,1), c(0,0,1,0), c(1,0,0,0), c(0,1,0,0)), info = "global threshold")

#### Unweighted Unipartite Models (sparsify) ####
## Unit tests awaiting further revisions to sparisfy functions




