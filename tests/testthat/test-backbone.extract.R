test_that("backbone.extract fwer = none", {
  sdsm <- sdsm(davis)
  fdsm <- fdsm(davis)
  hyperg <- hyperg(davis)

  sdsm.signed.none <- backbone.extract(sdsm, signed = TRUE, alpha = 0.4, fwer = "none")
  sdsm.signed.none[sdsm.signed.none < 0] <- 0
  sdsm.positive.none <- backbone.extract(sdsm, signed = FALSE, alpha = 0.4, fwer = "none")
  expect_equal(sdsm.signed.none,sdsm.positive.none)

  fdsm.signed.none <- backbone.extract(fdsm, signed = TRUE, alpha = 0.4, fwer = "none")
  fdsm.signed.none[fdsm.signed.none < 0] <- 0
  fdsm.positive.none <- backbone.extract(fdsm, signed = FALSE, alpha = 0.4, fwer = "none")
  expect_equal(fdsm.signed.none,fdsm.positive.none)

  hyperg.signed.none <- backbone.extract(hyperg, signed = TRUE, alpha = 0.4, fwer = "none")
  hyperg.signed.none[hyperg.signed.none < 0] <- 0
  hyperg.positive.none <- backbone.extract(hyperg, signed = FALSE, alpha = 0.4, fwer = "none")
  expect_equal(hyperg.signed.none,hyperg.positive.none)
})

test_that("backbone.extract fwer = holm", {
  sdsm <- sdsm(davis)
  fdsm <- fdsm(davis)
  hyperg <- hyperg(davis)

  sdsm.signed.none <- backbone.extract(sdsm, signed = TRUE, alpha = 0.4, fwer = "holm")
  sdsm.signed.none[sdsm.signed.none < 0] <- 0
  sdsm.positive.none <- backbone.extract(sdsm, signed = FALSE, alpha = 0.4, fwer = "holm")
  expect_equal(sdsm.signed.none,sdsm.positive.none)

  fdsm.signed.none <- suppressWarnings(backbone.extract(fdsm, signed = TRUE, alpha = 0.4, fwer = "holm"))
  fdsm.signed.none[fdsm.signed.none < 0] <- 0
  fdsm.positive.none <- suppressWarnings(backbone.extract(fdsm, signed = FALSE, alpha = 0.4, fwer = "holm"))
  expect_equal(fdsm.signed.none,fdsm.positive.none)

  hyperg.signed.none <- backbone.extract(hyperg, signed = TRUE, alpha = 0.4, fwer = "holm")
  hyperg.signed.none[hyperg.signed.none < 0] <- 0
  hyperg.positive.none <- backbone.extract(hyperg, signed = FALSE, alpha = 0.4, fwer = "holm")
  expect_equal(hyperg.signed.none,hyperg.positive.none)
})


test_that("backbone.extract fwer = bonferroni", {
  sdsm <- sdsm(davis)
  fdsm <- fdsm(davis)
  hyperg <- hyperg(davis)

  sdsm.signed.none <- backbone.extract(sdsm, signed = TRUE, alpha = 0.4, fwer = "bonferroni")
  sdsm.signed.none[sdsm.signed.none < 0] <- 0
  sdsm.positive.none <- backbone.extract(sdsm, signed = FALSE, alpha = 0.4, fwer = "bonferroni")
  expect_equal(sdsm.signed.none,sdsm.positive.none)

  fdsm.signed.none <- suppressWarnings(backbone.extract(fdsm, signed = TRUE, alpha = 0.4, fwer = "bonferroni"))
  fdsm.signed.none[fdsm.signed.none < 0] <- 0
  fdsm.positive.none <- suppressWarnings(backbone.extract(fdsm, signed = FALSE, alpha = 0.4, fwer = "bonferroni"))
  expect_equal(fdsm.signed.none,fdsm.positive.none)

  hyperg.signed.none <- backbone.extract(hyperg, signed = TRUE, alpha = 0.4, fwer = "bonferroni")
  hyperg.signed.none[hyperg.signed.none < 0] <- 0
  hyperg.positive.none <- backbone.extract(hyperg, signed = FALSE, alpha = 0.4, fwer = "bonferroni")
  expect_equal(hyperg.signed.none,hyperg.positive.none)
})

test_that("backbone.extract class outputs", {
  probs <- sdsm(davis)
  bb <- backbone.extract(probs, alpha = .2, signed = TRUE, fwer = "none")
  bb.matrix <- backbone.extract(probs, alpha = .2, signed = TRUE, fwer = "none", class = "matrix")
  bb.spmatrix <- backbone.extract(probs, alpha = .2, signed = TRUE, fwer = "none", class = "sparseMatrix")
  bb.igraph <- backbone.extract(probs, alpha = .2, signed = TRUE, fwer = "none", class = "igraph")
  bb.network <- backbone.extract(probs, alpha = .2, signed = TRUE, fwer = "none", class = "network")
  bb.edgelist <- backbone.extract(probs, alpha = .2, signed = TRUE, fwer = "none", class = "edgelist")

  expect_equal(bb, bb.matrix)
  expect_equal(class.convert(bb, "sparseMatrix")[[2]], bb.spmatrix)
  expect_true(igraph::identical_graphs(class.convert(bb, "igraph")[[2]], bb.igraph))
  expect_equal(class.convert(bb, "network")[[2]], bb.network)
  expect_equal(class.convert(bb, "edgelist")[[2]], bb.edgelist)
})


















