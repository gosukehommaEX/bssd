test_that("group_patterns counts subjects with the same data", {
  g <- group_patterns(c(1, 0, 1, 0, 1), c(2, 3, 2, 3, 1))
  expect_identical(g$y, c(0, 1, 1))
  expect_identical(g$t, c(3, 1, 2))
  expect_identical(g$w, c(2, 1, 2))
  g1 <- group_patterns(1, 5)
  expect_identical(g1$w, 1)
})
