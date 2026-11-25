# Reference values: inst/validation/python/ref_exptime.py (composite
# Gauss-Legendre rules on (log lambda_c, log lambda_d), split exactly at the
# thresholds), without the box refinement and the Chebyshev interpolation of
# the package.

test_that("exptime_post matches the reference values with rate0", {
  d <- data_e1()
  fit <- exptime_post(d$y, d$t, k = 2, prior_c = c(30, 1500),
                      prior_d = c(2, 70), rate0 = 0.03, rr0 = 1.2)
  expect_s3_class(fit, "bssd_post")
  expect_identical(fit$model, "EXPTIME")
  expect_equal(fit$prob, 0.3912496214049057, tolerance = 1e-6)
  # the comparison is not vacuous: a shifted reference value fails
  expect_failure(expect_equal(fit$prob, 0.3912496214049057 * (1 + 1e-4),
                              tolerance = 1e-6))
  expect_equal(unname(fit$probs), c(0.3912496214049057, 0.6384825953389326),
               tolerance = 1e-6)
  expect_equal(unname(fit$estimates[c("lambda_c", "lambda_d", "r")]),
               c(0.020358086360643606, 0.028341655149505953,
                 1.4488414415439201), tolerance = 1e-6)
  expect_equal(unname(fit$rate[c("lower", "median", "upper")]),
               c(0.013681907314679885, 0.027626614108717192,
                 0.04707284676474393), tolerance = 1e-5)
  expect_equal(unname(fit$rr[c("lower", "median", "upper")]),
               c(0.6016135660796988, 1.3731330635099064,
                 2.7305536454663413), tolerance = 1e-5)
  expect_identical(fit$signal, fit$prob > 0.95)
  expect_lt(fit$integration$error, 1e-5)
})

test_that("exptime_post gives unnamed results for named priors", {
  d <- data_e1()
  f1 <- exptime_post(d$y, d$t, k = 2, prior_c = c(30, 1500),
                     prior_d = c(2, 70), rate0 = 0.03, details = FALSE)
  f2 <- exptime_post(d$y, d$t, k = 2,
                     prior_c = c(shape = 30, rate = 1500),
                     prior_d = c(shape = 2, rate = 70), rate0 = 0.03,
                     details = FALSE)
  expect_null(names(f2$prob))
  expect_identical(f2$prob, f1$prob)
})

test_that("exptime_post matches the reference values for the rate ratio", {
  d <- data_e2()
  fit <- exptime_post(d$y, d$t, k = 1, prior_c = c(50, 300),
                      prior_d = c(1, 5))
  expect_identical(fit$rr_label, "r")
  expect_equal(fit$prob, 0.3977976747936025, tolerance = 1e-6)
  expect_true(is.na(fit$probs[["drug"]]))
  expect_equal(unname(fit$estimates[c("lambda_c", "lambda_d", "r")]),
               c(0.16828719791981767, 0.1528404804903628,
                 0.9443465170571876), tolerance = 1e-6)
  expect_equal(unname(fit$rr[c("lower", "median", "upper")]),
               c(0.16317195087743558, 0.8693380934443319,
                 2.1576335854057658), tolerance = 1e-5)
})

test_that("exptime_post with details = FALSE gives the same probability", {
  d <- data_e1()
  f1 <- exptime_post(d$y, d$t, k = 2, prior_c = c(30, 1500),
                     prior_d = c(2, 70), rate0 = 0.03)
  f2 <- exptime_post(d$y, d$t, k = 2, prior_c = c(30, 1500),
                     prior_d = c(2, 70), rate0 = 0.03, details = FALSE)
  expect_equal(f2$prob, f1$prob, tolerance = 1e-12)
  expect_true(all(is.na(f2$rr[c("median", "lower", "upper")])))
})

test_that("the posterior probability of exptime_post increases with events", {
  d <- data_e1()
  f0 <- exptime_post(d$y, d$t, k = 2, prior_c = c(30, 1500),
                     prior_d = c(2, 70), rate0 = 0.03, details = FALSE)
  y2 <- d$y
  y2[1:5] <- y2[1:5] + 1
  f1 <- exptime_post(y2, d$t, k = 2, prior_c = c(30, 1500),
                     prior_d = c(2, 70), rate0 = 0.03, details = FALSE)
  expect_gt(f1$prob, f0$prob)
})

test_that("exptime_post checks its input", {
  d <- data_e1()
  expect_error(exptime_post(d$y, d$t[-1], prior_c = c(1, 1),
                            prior_d = c(1, 1)), "'t'")
  expect_error(exptime_post(d$y - 1, d$t, prior_c = c(1, 1),
                            prior_d = c(1, 1)), "'y'")
  expect_error(exptime_post(d$y, d$t, prior_c = c(1, -1),
                            prior_d = c(1, 1)), "'prior_c'")
  expect_error(exptime_post(d$y, d$t, prior_c = c(30, 1500),
                            prior_d = c(1, 1), rate0 = -1), "'rate0'")
  expect_error(exptime_post(d$y, d$t, prior_c = c(30, 1500),
                            prior_d = c(1, 1), control = list(bad = 1)),
               "Unknown")
})
