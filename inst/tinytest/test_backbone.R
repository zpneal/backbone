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
#ADD HERE

#### Null Model Functions ####
## .sdsm()
M <- rbind(c(0,0,1),c(0,1,0),c(1,0,1))
test <- backbone:::.sdsm(M, signed = TRUE, missing_as_zero = FALSE)
test$upper <- round(test$upper,3)
test$lower <- round(test$lower,3)
expect_true(is(test, "list") & length(test)==2)  #Output is a two-item list
expect_equal(test$upper, rbind(c(.380,NA,.606),c(NA,.38,NA),c(.606,NA,.437)))  #Upper-tail p-values
expect_equal(test$lower, rbind(c(.949,NA,.864),c(NA,.949,NA),c(.864,NA,.916)))  #Lower-tail p-values

#### Statistical Backbone Functions ####
#ADD HERE

