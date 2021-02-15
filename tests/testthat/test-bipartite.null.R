test_that("bipartite.null matches", {
  B <- davis
  unwrapped <- sdsm(B)
  wrapped <- bipartite.null(B, rows = TRUE, cols = TRUE)
  expect_equal(unwrapped$positive, wrapped$positive)
  expect_equal(unwrapped$negative, wrapped$negative)

  # unwrapped <- fdsm(B, trials = 100)
  # wrapped <- bipartite.null(B, rows = TRUE, cols = TRUE, trials = 100)
  # expect_equal(unwrapped$positive, wrapped$positive)
  # expect_equal(unwrapped$negative, wrapped$negative)

  unwrapped <- hyperg(B)
  wrapped <- bipartite.null(B, rows = TRUE, cols = FALSE)
  expect_equal(unwrapped$positive, wrapped$positive) #OK
  expect_equal(unwrapped$negative, wrapped$negative) #OK
})
