test_that("lmix_dens is the log density of the logit", {
  mix <- data.frame(w = c(0.6, 0.4), shape1 = c(2, 0.3), shape2 = c(5, 4))
  x <- c(-20, -2, 0, 1.5, 5)
  p <- plogis(x)
  ref <- log((0.6 * dbeta(p, 2, 5) + 0.4 * dbeta(p, 0.3, 4)) * p * (1 - p))
  expect_equal(lmix_dens(x, mix), ref, tolerance = 1e-12)
  expect_failure(expect_equal(lmix_dens(x, mix), ref + 1e-8,
                              tolerance = 1e-12))
  # the density integrates to the difference of the distribution function
  val <- stats::integrate(function(x) exp(lmix_dens(x, mix)), -200, 3,
                          rel.tol = 1e-11, subdivisions = 500L)$value
  expect_equal(val, exp(lmix_cdf(3, mix)) - exp(lmix_cdf(-200, mix)),
               tolerance = 1e-9)
  expect_true(is.finite(lmix_dens(-1e10, mix)))
})
