test_that("post1d gives the distribution function and the density", {
  lp <- function(nodes) dgamma(nodes[[1]], 3, 2, log = TRUE)
  g <- post1d(lp, list(lo = 0, hi = 15), 81L)
  expect_equal(sum(g$wq), 1, tolerance = 1e-10)
  expect_equal(g$cdf(c(1, 2)), pgamma(c(1, 2), 3, 2), tolerance = 1e-8)
  expect_equal(g$dens(c(1, 2)), dgamma(c(1, 2), 3, 2), tolerance = 1e-8)
})
