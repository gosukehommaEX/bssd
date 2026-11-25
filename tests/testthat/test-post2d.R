test_that("post2d integrates a product of Gaussian densities", {
  # u ~ N(1, 0.5^2) and v ~ N(-2, 0.3^2), independent
  g0 <- gl_nodes(48L)
  u <- 6 * g0$x + 1
  gw <- 6 * g0$w
  lo <- rep(-5, length(u))
  hi <- rep(1, length(u))
  v <- row_nodes(lo, hi, 61L)
  L <- dnorm(u, 1, 0.5, log = TRUE) + dnorm(v, -2, 0.3, log = TRUE)
  g <- post2d(L, u, gw, lo, hi)
  expect_equal(sum(g$wq), 1, tolerance = 1e-12)
  expect_equal(g$v, v)
  expect_equal(g$cdf_rows(-1.8), pnorm(-1.8, -2, 0.3), tolerance = 1e-9)
  # a threshold that depends on u: P(v <= u - 3)
  ref <- pnorm(0, mean = -2 - (1 - 3), sd = sqrt(0.3^2 + 0.5^2))
  expect_equal(g$cdf_rows(g$u - 3), ref, tolerance = 1e-8)
  expect_equal(sum(g$wq * g$u), 1, tolerance = 1e-10)
})

test_that("post2d uses a separate interval in each row", {
  # u ~ N(0, 1) and v | u ~ N(3 u, 0.1^2): a narrow ridge, each row with
  # the interval 3 u +- 1
  g0 <- gl_nodes(48L)
  u <- 7 * g0$x
  gw <- 7 * g0$w
  lo <- 3 * u - 1
  hi <- 3 * u + 1
  L <- dnorm(u, log = TRUE) +
    dnorm(row_nodes(lo, hi, 49L), 3 * u, 0.1, log = TRUE)
  g <- post2d(L, u, gw, lo, hi)
  expect_equal(g$cdf_rows(3 * u + 0.05), pnorm(0.5), tolerance = 1e-6)
  expect_equal(sum(g$wq * g$v^2), 9.01, tolerance = 1e-8)
  # thresholds below and above the intervals give 0 and 1
  expect_equal(g$cdf_rows(-Inf), 0, tolerance = 1e-12)
  expect_equal(g$cdf_rows(Inf), 1, tolerance = 1e-12)
})
