# Reference values: inst/validation/python/ref_multibin.py (nested adaptive
# quadrature with scipy.integrate.quad).

test_that("multibin_post matches the reference values for four studies", {
  fit <- multibin_post(y = c(6, 10, 25, 47), n = c(200, 400, 800, 1600),
                       k = 1, prior_c = c(2, 98))
  expect_s3_class(fit, "bssd_post")
  expect_identical(fit$model, "MULTIBIN")
  expect_equal(fit$prob, 0.7692121710498092, tolerance = 1e-6)
  expect_failure(expect_equal(fit$prob, 0.7692121710498092 * (1 + 1e-4),
                              tolerance = 1e-6))
  expect_equal(unname(fit$rr[c("lower", "median", "upper")]),
               c(0.172947747423795, 1.6745101788866592, 5.239173113461459),
               tolerance = 1e-5)
  expect_true(is.na(fit$rr[["mean"]]))
  expect_equal(unname(fit$estimates["log_delta"]), 0.369402353108627,
               tolerance = 1e-5)
  expect_lt(fit$integration$error, 1e-5)
})

test_that("multibin_post matches the reference values with study priors", {
  fit <- multibin_post(y = c(3, 9), n = c(150, 300), k = c(1, 2),
                       prior_c = rbind(c(2, 98), c(4, 196)),
                       lrr_prior = c(0, 2), rr0 = 1.5)
  expect_equal(fit$prob, 0.4312016377483243, tolerance = 1e-6)
  expect_equal(unname(fit$rr[c("lower", "median", "upper")]),
               c(0.16691008930051165, 1.3252939461633837,
                 4.826140568787256), tolerance = 1e-5)
  expect_equal(unname(fit$estimates["log_delta"]), 0.18961941043329184,
               tolerance = 1e-5)
})

test_that("multibin_post checks its input", {
  expect_error(multibin_post(y = c(3, 9), n = c(150), prior_c = c(2, 98)),
               "'y'")
  expect_error(multibin_post(y = c(3, 400), n = c(150, 300),
                             prior_c = c(2, 98)), "exceed")
  expect_error(multibin_post(y = c(3, 9), n = c(150, 300),
                             prior_c = matrix(1, 3, 2)), "one row")
  expect_error(multibin_post(y = c(3, 9), n = c(150, 300),
                             prior_c = c(2, 98), lrr_prior = c(0, 0)),
               "standard deviation")
})
