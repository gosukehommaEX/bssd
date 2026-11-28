test_that("beta_coord uses the logit when no component has a heavy tail", {
  mix <- data.frame(w = c(0.7, 0.3), shape1 = c(20, 2), shape2 = c(180, 18))
  cd <- beta_coord(mix, 1e-8)
  expect_identical(cd$coordinate, "logit")
  u <- c(-4, -2.2, -1)
  expect_identical(cd$c_of(u), plogis(u))
  expect_identical(cd$lprior_u(u), lmix_dens(u, mix))
})

test_that("beta_coord uses a normal quantile for a shape parameter below 1", {
  cd <- beta_coord(data.frame(w = 1, shape1 = 0.5, shape2 = 5), 1e-8)
  expect_identical(cd$coordinate, "normal")
  u <- c(-3, -0.5, 0, 2)
  expect_equal(cd$c_of(u), qbeta(pnorm(u), 0.5, 5), tolerance = 1e-10)
  expect_equal(cd$lprior_u(u), dnorm(u, log = TRUE), tolerance = 1e-12)
  expect_equal(cd$u_of(cd$c_of(u)), u, tolerance = 1e-9)
  # MAP prior from few historical events: the prior density of u
  # integrates to 1
  cd <- beta_coord(map_beta(events = c(1, 1, 0), n = c(40, 40, 40)), 1e-8)
  expect_identical(cd$coordinate, "normal")
  val <- stats::integrate(function(u) exp(cd$lprior_u(u)), -9, 9,
                          rel.tol = 1e-10, subdivisions = 500L)$value
  expect_equal(val, 1, tolerance = 1e-7)
})
