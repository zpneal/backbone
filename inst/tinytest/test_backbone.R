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