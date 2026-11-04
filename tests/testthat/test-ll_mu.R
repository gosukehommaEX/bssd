test_that("ll_mu inverts ll_cdf", {
  # reference value computed with Python: log(104) - 0.8 * logit(0.035)
  expect_equal(ll_mu(0.035, 104, 0.8), 7.29781493102103, tolerance = 1e-12)
  mu <- ll_mu(c(0.05, 0.2, 0.6), c(52, 104, 30), c(2, 0.8, 1.2))
  expect_equal(ll_cdf(c(52, 104, 30), mu, c(2, 0.8, 1.2)), c(0.05, 0.2, 0.6),
               tolerance = 1e-12)
  expect_error(ll_mu(1, 104, 0.8), "out of range")
})
