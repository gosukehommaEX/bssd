test_that("post2d integrates a product of Gaussian densities", {
  # u ~ N(1, 0.5^2) and v ~ N(-2, 0.3^2), independent
  lp <- function(nodes) {
    outer(dnorm(nodes[[1]], 1, 0.5, log = TRUE),
          dnorm(nodes[[2]], -2, 0.3, log = TRUE), "+")
  }
  g <- post2d(lp, list(lo = c(-5, -5), hi = c(7, 1)), 48L, 61L)
  expect_equal(sum(g$wq), 1, tolerance = 1e-12)
  expect_equal(g$cdf_v(-1.8), pnorm(-1.8, -2, 0.3), tolerance = 1e-9)
  expect_equal(g$cdf_v(c(-2.2, -1.8)), pnorm(c(-2.2, -1.8), -2, 0.3),
               tolerance = 1e-9)
  # a threshold that depends on u: P(v <= u - 3)
  ref <- pnorm(0, mean = -2 - (1 - 3), sd = sqrt(0.3^2 + 0.5^2))
  expect_equal(g$cdf_rows(g$u - 3), ref, tolerance = 1e-8)
  expect_equal(sum(g$wq * g$u), 1, tolerance = 1e-10)
})
