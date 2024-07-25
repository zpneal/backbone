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
test <- backbone:::.sdsm(M, signed = TRUE, missing_as_zero = FALSE)
test$upper <- round(test$upper,3)
test$lower <- round(test$lower,3)
expect_true(is(test, "list") & length(test)==2)  #Output is a two-item list
expect_equal(test$upper, rbind(c(.380,NA,.606),  #Upper-tail p-values
                               c(NA,.38,NA),
                               c(.606,NA,.437)))  
expect_equal(test$lower, rbind(c(.949,NA,.864),  #Lower-tail p-values
                               c(NA,.949,NA),
                               c(.864,NA,.916)))  

#### Statistical Backbone Functions ####
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

bb <- backbone_from_bipartite(B, signed = TRUE)  #Extract from matrix as signed
expect_true(is(bb,"matrix"))         #Returns as matrix
expect_true(all(bb %in% c(-1,0,1)))  #Contains only -1, 0, 1
expect_true(any(bb %in% c(-1)))      #Contains some negative edges
expect_true(any(bb %in% c(0)))       #Contains some missing edges
expect_true(any(bb %in% c(1)))       #Contains some positive edges

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

bb <- backbone_from_bipartite(B)                                              #Extract from igraph with defaults
expect_true(is(bb,"igraph"))                                                  #Returns as igraph
expect_identical(igraph::vertex_attr_names(bb), c("agent_attrib"))            #Contains correct vertex attributes
expect_identical(igraph::edge_attr_names(bb), c("oldweight"))                 #Contains correct edge attributes
expect_true(igraph::modularity(bb, c(rep(1,10), rep(2,10), rep(3,10))) > .5)  #Backbone has high modularity