test_that("gamma_coord uses the log for a shape parameter of at least 1", {
  cd <- gamma_coord(30, 1500)
  expect_identical(cd$coordinate, "log")
  u <- log(c(0.01, 0.02, 0.03))
  expect_equal(cd$lprior_u(u), dgamma(exp(u), 30, 1500, log = TRUE) + u,
               tolerance = 1e-12)
  expect_identical(cd$c_of(u), exp(u))
})

test_that("gamma_coord uses a normal quantile for a shape parameter below 1", {
  cd <- gamma_coord(0.5, 2)
  expect_identical(cd$coordinate, "normal")
  u <- c(-6, -1, 0, 1.5)
  expect_equal(cd$c_of(u), qgamma(pnorm(u), 0.5, 2), tolerance = 1e-10)
  expect_equal(cd$u_of(cd$c_of(u)), u, tolerance = 1e-9)
  expect_identical(cd$lprior_u(u), dnorm(u, log = TRUE))
})
