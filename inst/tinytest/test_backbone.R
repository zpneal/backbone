#### Utility Functions ####
## bicm()
M <- rbind(c(0,0,1),c(0,1,0),c(1,0,1))
test <- round(bicm(M),3)
expect_equal(test, rbind(c(.216,.216,.568),c(.216,.216,.568),c(.568,.568,.863)))  #BiCM probabilities

## fastball()
M <- matrix(rbinom(100*1000,1,0.5),100,1000)
test <- fastball(M)
expect_equal(rowSums(test), rowSums(M))  #Row sums match
expect_equal(colSums(test), colSums(M))  #Column sums match

## .retain()
upper <- rbind(c(.01,.02,.03),  #Unsigned
               c(.05,.06,.07),
               c(0.5,0.6,0.7))
p <- list(upper = upper)
test <- backbone:::.retain(p, alpha = 0.05, mtc = "none")
expect_equal(test, rbind(c(0,1,1),
                         c(0,0,0),
                         c(0,0,0)))

upper <- rbind(c(.01,.02,.03),  #Signed
               c(.05,.06,.07),
               c(0.5,0.6,0.7))
lower <- rbind(c(0.5,0.6,0.7),
               c(.05,.06,.07),
               c(.01,.02,.03))
p <- list(lower = lower, upper = upper)
test <- backbone:::.retain(p, alpha = 0.1, mtc = "none")  #Higher alpha because this is a two-tailed test
expect_equal(test, rbind(c(0,1,1),
                         c(0,0,0),
                         c(-1,-1,0)))

#### Null Model Functions ####
## .sdsm()
M <- rbind(c(0,0,1),c(0,1,0),c(1,0,1))
test <- backbone:::.sdsm(M, signed = TRUE, missing_as_zero = TRUE)
test$upper <- round(test$upper,3)
test$lower <- round(test$lower,3)
expect_true(is(test, "list") & length(test)==2)  #Output is a two-item list
expect_true(all(is.na(diag(test$upper))))  #Upper-tail diagonal is missing
expect_true(all(is.na(diag(test$lower))))  #Lower-tail diagonal is missing
expect_true(isSymmetric(test$upper))  #Upper-tail is symmetric
expect_true(isSymmetric(test$lower))  #Lower-tail is symmetric
expect_true(all(test$upper[upper.tri(test$upper)]>=0 & test$upper[upper.tri(test$upper)]<=1))  #Upper-tail p-values between 0 and 1
expect_true(all(test$lower[upper.tri(test$lower)]>=0 & test$lower[upper.tri(test$lower)]<=1))  #Lower-tail p-values between 0 and 1

## .sdsm_ec()
M <- rbind(c(10,0,1),c(0,11,0),c(1,0,1))
test <- backbone:::.sdsm_ec(M, signed = TRUE, missing_as_zero = TRUE)
test$upper <- round(test$upper,3)
test$lower <- round(test$lower,3)
expect_true(is(test, "list") & length(test)==2)  #Output is a two-item list
expect_true(all(is.na(diag(test$upper))))  #Upper-tail diagonal is missing
expect_true(all(is.na(diag(test$lower))))  #Lower-tail diagonal is missing
expect_true(isSymmetric(test$upper))  #Upper-tail is symmetric
expect_true(isSymmetric(test$lower))  #Lower-tail is symmetric
expect_true(all(test$upper[upper.tri(test$upper)]>=0 & test$upper[upper.tri(test$upper)]<=1))  #Upper-tail p-values between 0 and 1
expect_true(all(test$lower[upper.tri(test$lower)]>=0 & test$lower[upper.tri(test$lower)]<=1))  #Lower-tail p-values between 0 and 1

## .fixedrow()
M <- rbind(c(0,0,1),c(0,1,0),c(1,0,1))
test <- backbone:::.fixedrow(M, signed = TRUE, missing_as_zero = TRUE)
test$upper <- round(test$upper,3)
test$lower <- round(test$lower,3)
expect_true(is(test, "list") & length(test)==2)  #Output is a two-item list
expect_true(all(is.na(diag(test$upper))))  #Upper-tail diagonal is missing
expect_true(all(is.na(diag(test$lower))))  #Lower-tail diagonal is missing
expect_true(isSymmetric(test$upper))  #Upper-tail is symmetric
expect_true(isSymmetric(test$lower))  #Lower-tail is symmetric
expect_true(all(test$upper[upper.tri(test$upper)]>=0 & test$upper[upper.tri(test$upper)]<=1))  #Upper-tail p-values between 0 and 1
expect_true(all(test$lower[upper.tri(test$lower)]>=0 & test$lower[upper.tri(test$lower)]<=1))  #Lower-tail p-values between 0 and 1

## .fixedcol()
M <- rbind(c(0,0,1),c(0,1,0),c(1,0,1))
test <- backbone:::.fixedcol(M, signed = TRUE, missing_as_zero = TRUE)
test$upper <- round(test$upper,3)
test$lower <- round(test$lower,3)
expect_true(is(test, "list") & length(test)==2)  #Output is a two-item list
expect_true(all(is.na(diag(test$upper))))  #Upper-tail diagonal is missing
expect_true(all(is.na(diag(test$lower))))  #Lower-tail diagonal is missing
expect_true(isSymmetric(test$upper))  #Upper-tail is symmetric
expect_true(isSymmetric(test$lower))  #Lower-tail is symmetric
expect_true(all(test$upper[upper.tri(test$upper)]>=0 & test$upper[upper.tri(test$upper)]<=1))  #Upper-tail p-values between 0 and 1
expect_true(all(test$lower[upper.tri(test$lower)]>=0 & test$lower[upper.tri(test$lower)]<=1))  #Lower-tail p-values between 0 and 1

## .fixedfill()
M <- rbind(c(0,0,1),c(0,1,0),c(1,0,1))
test <- backbone:::.fixedfill(M, signed = TRUE, missing_as_zero = TRUE)
test$upper <- round(test$upper,3)
test$lower <- round(test$lower,3)
expect_true(is(test, "list") & length(test)==2)  #Output is a two-item list
expect_true(all(is.na(diag(test$upper))))  #Upper-tail diagonal is missing
expect_true(all(is.na(diag(test$lower))))  #Lower-tail diagonal is missing
expect_true(isSymmetric(test$upper))  #Upper-tail is symmetric
expect_true(isSymmetric(test$lower))  #Lower-tail is symmetric
expect_true(all(test$upper[upper.tri(test$upper)]>=0 & test$upper[upper.tri(test$upper)]<=1))  #Upper-tail p-values between 0 and 1
expect_true(all(test$lower[upper.tri(test$lower)]>=0 & test$lower[upper.tri(test$lower)]<=1))  #Lower-tail p-values between 0 and 1

## .fdsm()
M <- rbind(c(0,0,1),c(0,1,0),c(1,0,1))
test <- backbone:::.fdsm(M, signed = TRUE, missing_as_zero = TRUE, alpha = 0.05, mtc = "none", trials = 1000)
test$upper <- round(test$upper,3)
test$lower <- round(test$lower,3)
expect_true(is(test, "list") & length(test)==2)  #Output is a two-item list
expect_true(all(is.na(diag(test$upper))))  #Upper-tail diagonal is missing
expect_true(all(is.na(diag(test$lower))))  #Lower-tail diagonal is missing
expect_true(isSymmetric(test$upper))  #Upper-tail is symmetric
expect_true(isSymmetric(test$lower))  #Lower-tail is symmetric
expect_true(all(test$upper[upper.tri(test$upper)]>=0 & test$upper[upper.tri(test$upper)]<=1))  #Upper-tail p-values between 0 and 1
expect_true(all(test$lower[upper.tri(test$lower)]>=0 & test$lower[upper.tri(test$lower)]<=1))  #Lower-tail p-values between 0 and 1

## .disparity()
M <- rbind(c(0,0,1),c(0,1,0),c(1,0,1))
test <- backbone:::.disparity(M, signed = TRUE, missing_as_zero = TRUE)
test$upper <- round(test$upper,3)
test$lower <- round(test$lower,3)
expect_true(is(test, "list") & length(test)==2)  #Output is a two-item list
expect_true(all(test$upper[upper.tri(test$upper)]>=0 & test$upper[upper.tri(test$upper)]<=1))  #Upper-tail p-values between 0 and 1
expect_true(all(test$lower[upper.tri(test$lower)]>=0 & test$lower[upper.tri(test$lower)]<=1))  #Lower-tail p-values between 0 and 1

## .lans()
M <- rbind(c(0,0,1),c(0,1,0),c(1,0,1))
test <- backbone:::.lans(M, signed = TRUE, missing_as_zero = TRUE)
test$upper <- round(test$upper,3)
test$lower <- round(test$lower,3)
expect_true(is(test, "list") & length(test)==2)  #Output is a two-item list
expect_true(all(test$upper[upper.tri(test$upper)]>=0 & test$upper[upper.tri(test$upper)]<=1))  #Upper-tail p-values between 0 and 1
expect_true(all(test$lower[upper.tri(test$lower)]>=0 & test$lower[upper.tri(test$lower)]<=1))  #Lower-tail p-values between 0 and 1

## .mlf()
M <- rbind(c(0,0,1),c(0,1,0),c(1,0,1))
test <- backbone:::.mlf(M, signed = TRUE, missing_as_zero = TRUE)
test$upper <- round(test$upper,3)
test$lower <- round(test$lower,3)
expect_true(is(test, "list") & length(test)==2)  #Output is a two-item list
expect_true(all(test$upper[upper.tri(test$upper)]>=0 & test$upper[upper.tri(test$upper)]<=1))  #Upper-tail p-values between 0 and 1
expect_true(all(test$lower[upper.tri(test$lower)]>=0 & test$lower[upper.tri(test$lower)]<=1))  #Lower-tail p-values between 0 and 1

#### Bipartite Backbone Functions ####
## Define function to compute triangle index
trace <- function(x){sum(diag(x))}
matcube <- function(x){x%*%x%*%x}
triangle_index <- function(x){(trace(matcube(x)) + trace(matcube(abs(x))))/(2 * trace(matcube(abs(x))))}

## Bipartite from matrix
B <- rbind(cbind(matrix(rbinom(250,1,.8),10),   #An example block incidence matrix
                 matrix(rbinom(250,1,.2),10),
                 matrix(rbinom(250,1,.2),10)),
           cbind(matrix(rbinom(250,1,.2),10),
                 matrix(rbinom(250,1,.8),10),
                 matrix(rbinom(250,1,.2),10)),
           cbind(matrix(rbinom(250,1,.2),10),
                 matrix(rbinom(250,1,.2),10),
                 matrix(rbinom(250,1,.8),10)))

bb <- backbone_from_bipartite(B, model = "sdsm", signed = TRUE)  #Extract SDSM matrix as signed
expect_true(is(bb,"matrix"))         #Returns as matrix
expect_true(all(bb %in% c(-1,0,1)))  #Contains only -1, 0, 1
expect_true(any(bb %in% c(-1)))      #Contains some negative edges
expect_true(any(bb %in% c(0)))       #Contains some missing edges
expect_true(any(bb %in% c(1)))       #Contains some positive edges
expect_true(triangle_index(bb)>.8)   #Is nearly balanced

bb <- backbone_from_bipartite(B, model = "fdsm", signed = TRUE, trials = 250)  #Extract FDSM matrix as signed
expect_true(is(bb,"matrix"))         #Returns as matrix
expect_true(all(bb %in% c(-1,0,1)))  #Contains only -1, 0, 1
expect_true(any(bb %in% c(-1)))      #Contains some negative edges
expect_true(any(bb %in% c(0)))       #Contains some missing edges
expect_true(any(bb %in% c(1)))       #Contains some positive edges
expect_true(triangle_index(bb)>.8)   #Is nearly balanced

bb <- backbone_from_bipartite(B, model = "fixedrow", signed = TRUE)  #Extract fixedrow matrix as signed
expect_true(is(bb,"matrix"))         #Returns as matrix
expect_true(all(bb %in% c(-1,0,1)))  #Contains only -1, 0, 1
expect_true(any(bb %in% c(-1)))      #Contains some negative edges
expect_true(any(bb %in% c(0)))       #Contains some missing edges
expect_true(any(bb %in% c(1)))       #Contains some positive edges
expect_true(triangle_index(bb)>.8)   #Is nearly balanced

bb <- backbone_from_bipartite(B, model = "fixedcol", signed = TRUE)  #Extract fixedcol matrix as signed
expect_true(is(bb,"matrix"))         #Returns as matrix
expect_true(all(bb %in% c(-1,0,1)))  #Contains only -1, 0, 1
expect_true(any(bb %in% c(-1)))      #Contains some negative edges
expect_true(any(bb %in% c(0)))       #Contains some missing edges
expect_true(any(bb %in% c(1)))       #Contains some positive edges
expect_true(triangle_index(bb)>.8)   #Is nearly balanced

bb <- backbone_from_bipartite(B, model = "fixedfill", signed = TRUE)  #Extract fixedfill matrix as signed
expect_true(is(bb,"matrix"))         #Returns as matrix
expect_true(all(bb %in% c(-1,0,1)))  #Contains only -1, 0, 1
expect_true(any(bb %in% c(-1)))      #Contains some negative edges
expect_true(any(bb %in% c(0)))       #Contains some missing edges
expect_true(any(bb %in% c(1)))       #Contains some positive edges
expect_true(triangle_index(bb)>.8)   #Is nearly balanced

B <- as.vector(B)
make_prohibited <- sample(which(B==0), 5, replace = FALSE)  #Pick some missing edges to prohibit
B[make_prohibited] <- 10
make_required <- sample(which(B==1), 5, replace = FALSE)  #Pick some present edges to require
B[make_required] <- 11
B <- matrix(B, 30, 75)  #Reassemble as matrix
bb <- backbone_from_bipartite(B, model = "sdsm", signed = TRUE)  #Extract SDSM matrix as signed, considering structural values
expect_true(is(bb,"matrix"))         #Returns as matrix
expect_true(all(bb %in% c(-1,0,1)))  #Contains only -1, 0, 1
expect_true(any(bb %in% c(-1)))      #Contains some negative edges
expect_true(any(bb %in% c(0)))       #Contains some missing edges
expect_true(any(bb %in% c(1)))       #Contains some positive edges
expect_true(triangle_index(bb)>.8)   #Is nearly balanced

## Bipartite from igraph
B <- rbind(cbind(matrix(rbinom(250,1,.8),10),   #An example block incidence matrix
                 matrix(rbinom(250,1,.2),10),
                 matrix(rbinom(250,1,.2),10)),
           cbind(matrix(rbinom(250,1,.2),10),
                 matrix(rbinom(250,1,.8),10),
                 matrix(rbinom(250,1,.2),10)),
           cbind(matrix(rbinom(250,1,.2),10),
                 matrix(rbinom(250,1,.2),10),
                 matrix(rbinom(250,1,.8),10)))
B <- igraph::graph_from_biadjacency_matrix(B)          #Convert to igraph
igraph::V(B)$agent_attrib <- c(c(1:30),rep(NA,75))     #Add agent attribute
igraph::V(B)$artifact_attrib <- c(rep(NA,30),c(1:75))  #Add artifact attribute

bb <- backbone_from_bipartite(B, model = "sdsm")                              #Extract SDSM igraph with defaults
expect_true(is(bb,"igraph"))                                                  #Returns as igraph
expect_identical(igraph::vertex_attr_names(bb), c("agent_attrib"))            #Contains correct vertex attributes
expect_identical(igraph::edge_attr_names(bb), c("oldweight"))                 #Contains correct edge attributes
expect_true(igraph::modularity(bb, c(rep(1,10), rep(2,10), rep(3,10))) > .5)  #Backbone has high modularity

bb <- backbone_from_bipartite(B, model = "fdsm", trials = 250)                #Extract FDSM igraph with defaults
expect_true(is(bb,"igraph"))                                                  #Returns as igraph
expect_identical(igraph::vertex_attr_names(bb), c("agent_attrib"))            #Contains correct vertex attributes
expect_identical(igraph::edge_attr_names(bb), c("oldweight"))                 #Contains correct edge attributes
expect_true(igraph::modularity(bb, c(rep(1,10), rep(2,10), rep(3,10))) > .5)  #Backbone has high modularity

bb <- backbone_from_bipartite(B, model = "fixedrow")                          #Extract fixedrow igraph with defaults
expect_true(is(bb,"igraph"))                                                  #Returns as igraph
expect_identical(igraph::vertex_attr_names(bb), c("agent_attrib"))            #Contains correct vertex attributes
expect_identical(igraph::edge_attr_names(bb), c("oldweight"))                 #Contains correct edge attributes
expect_true(igraph::modularity(bb, c(rep(1,10), rep(2,10), rep(3,10))) > .5)  #Backbone has high modularity

bb <- backbone_from_bipartite(B, model = "fixedrow")                          #Extract fixedcol igraph with defaults
expect_true(is(bb,"igraph"))                                                  #Returns as igraph
expect_identical(igraph::vertex_attr_names(bb), c("agent_attrib"))            #Contains correct vertex attributes
expect_identical(igraph::edge_attr_names(bb), c("oldweight"))                 #Contains correct edge attributes
expect_true(igraph::modularity(bb, c(rep(1,10), rep(2,10), rep(3,10))) > .5)  #Backbone has high modularity

bb <- backbone_from_bipartite(B, model = "fixedfill")                         #Extract fixedcol igraph with defaults
expect_true(is(bb,"igraph"))                                                  #Returns as igraph
expect_identical(igraph::vertex_attr_names(bb), c("agent_attrib"))            #Contains correct vertex attributes
expect_identical(igraph::edge_attr_names(bb), c("oldweight"))                 #Contains correct edge attributes
expect_true(igraph::modularity(bb, c(rep(1,10), rep(2,10), rep(3,10))) > .5)  #Backbone has high modularity

igraph::E(B)$weight <- NA
igraph::E(B)$weight <- sample(c(1,11), length(igraph::E(B)$weight), replace = TRUE, prob = c(.9,.1))
bb <- backbone_from_bipartite(B, model = "sdsm")                              #Extract SDSM igraph with defaults, considering structural values
expect_true(is(bb,"igraph"))                                                  #Returns as igraph
expect_identical(igraph::vertex_attr_names(bb), c("agent_attrib"))            #Contains correct vertex attributes
expect_identical(igraph::edge_attr_names(bb), c("oldweight"))                 #Contains correct edge attributes
expect_true(igraph::modularity(bb, c(rep(1,10), rep(2,10), rep(3,10))) > .5)  #Backbone has high modularity

#### Weighted Backbone Functions ####
## Multiscale weighted matrix
W <- matrix(c(0,10,10,10,10,75,0,0,0,0,
              10,0,1,1,1,0,0,0,0,0,
              10,1,0,1,1,0,0,0,0,0,
              10,1,1,0,1,0,0,0,0,0,
              10,1,1,1,0,0,0,0,0,0,
              75,0,0,0,0,0,100,100,100,100,
              0,0,0,0,0,100,0,10,10,10,
              0,0,0,0,0,100,10,0,10,10,
              0,0,0,0,0,100,10,10,0,10,
              0,0,0,0,0,100,10,10,10,0),10)

bb <- backbone_from_weighted(W, model = "disparity")  #Extract disparity backbone
expect_true(is(bb,"matrix"))                          #Returns as matrix
bb <- igraph::graph_from_adjacency_matrix(bb, mode = "undirected")
expect_true(igraph::is_tree(bb))                      #Backbone is a tree

bb <- backbone_from_weighted(W, model = "lans")       #Extract lans backbone
expect_true(is(bb,"matrix"))                          #Returns as matrix
bb <- igraph::graph_from_adjacency_matrix(bb, mode = "undirected")
expect_true(igraph::is_tree(bb))                      #Backbone is a tree

bb <- backbone_from_weighted(W, model = "mlf")        #Extract mlf backbone
expect_true(is(bb,"matrix"))                          #Returns as matrix
bb <- igraph::graph_from_adjacency_matrix(bb, mode = "undirected")
expect_true(igraph::is_tree(bb))                      #Backbone is a tree

bb <- backbone_from_weighted(W, model = "global")     #Extract global backbone (unsigned)
expect_true(is(bb,"matrix"))                          #Returns as matrix
expect_true(table(bb)[1]==58 & table(bb)[2]==42)      #Contains 58 0s and 42 1s
bb <- backbone_from_weighted(W, model = "global", parameter = c(10,74))     #Extract global backbone (signed)
expect_true(table(bb)[1]==12 & table(bb)[2]==78 & table(bb)[3]==10)      #Contains 12 -1s, 78 0s, and 10 1s

## Multiscale weighted igraph
W <- igraph::graph_from_adjacency_matrix(W, mode = "undirected", weighted = TRUE)

bb <- backbone_from_weighted(W, model = "disparity")  #Extract disparity backbone
expect_true(is(bb,"igraph"))                          #Returns as igraph
expect_true(igraph::is_tree(bb))                      #Backbone is a tree

bb <- backbone_from_weighted(W, model = "lans")       #Extract lans backbone
expect_true(is(bb,"igraph"))                          #Returns as igraph
expect_true(igraph::is_tree(bb))                      #Backbone is a tree

bb <- backbone_from_weighted(W, model = "mlf")        #Extract mlf backbone
expect_true(is(bb,"igraph"))                          #Returns as igraph
expect_true(igraph::is_tree(bb))                      #Backbone is a tree

bb <- backbone_from_weighted(W, model = "global")     #Extract global backbone (unsigned)
expect_true(is(bb,"igraph"))                          #Returns as matrix
bb <- igraph::as_adjacency_matrix(bb, sparse = FALSE)   #Get matrix
expect_true(table(bb)[1]==58 & table(bb)[2]==42)      #Contains 58 0s and 42 1s
bb <- backbone_from_weighted(W, model = "global", parameter = c(10,74))     #Extract global backbone (signed)
bb <- igraph::as_adjacency_matrix(bb, sparse = FALSE, attr = "sign")   #Get matrix
expect_true(table(bb)[1]==12 & table(bb)[2]==78 & table(bb)[3]==10)      #Contains 12 -1s, 78 0s, and 10 1s

## Projection of bipartite matrix
W <- rbind(cbind(matrix(rbinom(250,1,.8),10),
                 matrix(rbinom(250,1,.2),10),
                 matrix(rbinom(250,1,.2),10)),
           cbind(matrix(rbinom(250,1,.2),10),
                 matrix(rbinom(250,1,.8),10),
                 matrix(rbinom(250,1,.2),10)),
           cbind(matrix(rbinom(250,1,.2),10),
                 matrix(rbinom(250,1,.2),10),
                 matrix(rbinom(250,1,.8),10)))
W <- W%*%t(W)
diag(W) <- 0

bb <- backbone_from_weighted(W, model = "disparity", signed = TRUE, alpha = 0.5)  #Extract signed disparity matrix
expect_true(is(bb,"matrix"))         #Returns as matrix
expect_true(all(bb %in% c(-1,0,1)))  #Contains only -1, 0, 1
expect_true(any(bb %in% c(-1)))      #Contains some negative edges
expect_true(any(bb %in% c(0)))       #Contains some missing edges
expect_true(any(bb %in% c(1)))       #Contains some positive edges
expect_true(triangle_index(bb)>.8)   #Is nearly balanced

bb <- backbone_from_weighted(W, model = "lans", signed = TRUE, alpha = 0.5)  #Extract signed lans matrix
expect_true(is(bb,"matrix"))         #Returns as matrix
expect_true(all(bb %in% c(-1,0,1)))  #Contains only -1, 0, 1
expect_true(any(bb %in% c(-1)))      #Contains some negative edges
expect_true(any(bb %in% c(0)))       #Contains some missing edges
expect_true(any(bb %in% c(1)))       #Contains some positive edges
expect_true(triangle_index(bb)>.8)   #Is nearly balanced

bb <- backbone_from_weighted(W, model = "mlf", signed = TRUE, alpha = 0.5)  #Extract signed mlf matrix
expect_true(is(bb,"matrix"))         #Returns as matrix
expect_true(all(bb %in% c(-1,0,1)))  #Contains only -1, 0, 1
expect_true(any(bb %in% c(-1)))      #Contains some negative edges
expect_true(any(bb %in% c(0)))       #Contains some missing edges
expect_true(any(bb %in% c(1)))       #Contains some positive edges
expect_true(triangle_index(bb)>.8)   #Is nearly balanced

upper <- mean(W) + sd(W)             #Use mean + sd as positive edge threshold
lower <- mean(W) - sd(W)             #Use mean - sd as negative edge threshold
bb <- backbone_from_weighted(W, model = "global", parameter = c(lower, upper))  #Extract signed global matrix
expect_true(is(bb,"matrix"))         #Returns as matrix
expect_true(all(bb %in% c(-1,0,1)))  #Contains only -1, 0, 1
expect_true(any(bb %in% c(-1)))      #Contains some negative edges
expect_true(any(bb %in% c(0)))       #Contains some missing edges
expect_true(any(bb %in% c(1)))       #Contains some positive edges
triangle_index(bb)
expect_true(triangle_index(bb)>.8)   #Is nearly balanced

## Projection of bipartite igraph
W <- rbind(cbind(matrix(rbinom(250,1,.8),10),
                 matrix(rbinom(250,1,.2),10),
                 matrix(rbinom(250,1,.2),10)),
           cbind(matrix(rbinom(250,1,.2),10),
                 matrix(rbinom(250,1,.8),10),
                 matrix(rbinom(250,1,.2),10)),
           cbind(matrix(rbinom(250,1,.2),10),
                 matrix(rbinom(250,1,.2),10),
                 matrix(rbinom(250,1,.8),10)))
W <- igraph::graph_from_biadjacency_matrix(W)
W <- igraph::bipartite_projection(W, which = "false")
igraph::V(W)$agent_attrib <- c(c(1:30))     #Add agent attribute

bb <- backbone_from_weighted(W, model = "disparity", alpha = 0.25)            #Extract unweighted disparity igraph
expect_true(is(bb,"igraph"))                                                  #Returns as igraph
expect_identical(igraph::vertex_attr_names(bb), c("agent_attrib"))            #Contains correct vertex attributes
expect_identical(igraph::edge_attr_names(bb), c("oldweight"))                 #Contains correct edge attributes
expect_true(igraph::modularity(bb, c(rep(1,10), rep(2,10), rep(3,10))) > .5)  #Backbone has high modularity

bb <- backbone_from_weighted(W, model = "lans", alpha = 0.25)                 #Extract unweighted lans igraph
expect_true(is(bb,"igraph"))                                                  #Returns as igraph
expect_identical(igraph::vertex_attr_names(bb), c("agent_attrib"))            #Contains correct vertex attributes
expect_identical(igraph::edge_attr_names(bb), c("oldweight"))                 #Contains correct edge attributes
expect_true(igraph::modularity(bb, c(rep(1,10), rep(2,10), rep(3,10))) > .5)  #Backbone has high modularity

bb <- backbone_from_weighted(W, model = "mlf", alpha = 0.25)                  #Extract unweighted mlf igraph
expect_true(is(bb,"igraph"))                                                  #Returns as igraph
expect_identical(igraph::vertex_attr_names(bb), c("agent_attrib"))            #Contains correct vertex attributes
expect_identical(igraph::edge_attr_names(bb), c("oldweight"))                 #Contains correct edge attributes
expect_true(igraph::modularity(bb, c(rep(1,10), rep(2,10), rep(3,10))) > .5)  #Backbone has high modularity

threshold <- mean(igraph::E(W)$weight) + sd(igraph::E(W)$weight)              #Use mean + sd as edge threshold
bb <- backbone_from_weighted(W, model = "global", parameter = mean_edge+sd_edge)      #Extract unweighted global igraph
expect_true(is(bb,"igraph"))                                                  #Returns as igraph
expect_identical(igraph::vertex_attr_names(bb), c("agent_attrib"))            #Contains correct vertex attributes
expect_identical(igraph::edge_attr_names(bb), c("oldweight"))                 #Contains correct edge attributes
expect_true(igraph::modularity(bb, c(rep(1,10), rep(2,10), rep(3,10))) > .5)  #Backbone has high modularity
