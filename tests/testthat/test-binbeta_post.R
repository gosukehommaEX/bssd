# Reference values: inst/validation/python/ref_binbeta.py (exact finite
# mixture of products of Beta distributions from the binomial expansion of
# the likelihood; cases B1 to B5).

test_that("binbeta_post matches the exact values with p0", {
  fit <- binbeta_post(y = 40, n = 300, k = 1, prior_c = c(31, 271),
                      prior_d = c(1, 1), p0 = 0.15)
  expect_s3_class(fit, "bssd_post")
  expect_identical(fit$model, "BINBETA")
  expect_equal(fit$prob, 0.6614733754153749, tolerance = 1e-7)
  expect_failure(expect_equal(fit$prob, 0.6614733754153749 * (1 + 1e-5),
                              tolerance = 1e-7))
  expect_equal(unname(fit$probs), c(0.6614733754153749, 0.8969458860909176),
               tolerance = 1e-7)
  expect_equal(unname(fit$estimates[c("theta_c", "theta_d", "rr")]),
               c(0.10264862816559746, 0.1688753774488444,
                 1.7242946661122427), tolerance = 1e-6)
  expect_equal(unname(fit$rate[c("lower", "median", "upper")]),
               c(0.0879347404014438, 0.16767113692809782,
                 0.2566240990607299), tolerance = 1e-6)
  expect_equal(unname(fit$rr[c("lower", "median", "upper")]),
               c(0.7135180420982211, 1.6516881211086338,
                 3.1514671174895152), tolerance = 1e-6)
})

test_that("binbeta_post matches the exact values for the ratio", {
  fit <- binbeta_post(y = 12, n = 80, k = 2, prior_c = c(5, 45),
                      prior_d = c(0.5, 0.5), rr0 = 1.5)
  expect_identical(fit$rr_label, "rr")
  expect_equal(fit$prob, 0.6135251833617793, tolerance = 1e-7)
  expect_equal(unname(fit$estimates[c("theta_c", "theta_d", "rr")]),
               c(0.10271122477026959, 0.17791744455013891,
                 2.2482454769923024), tolerance = 1e-5)
  expect_equal(unname(fit$rr[c("lower", "median", "upper")]),
               c(0.36994837263017366, 1.8189007890849602,
                 6.678978679651721), tolerance = 1e-6)
})

test_that("binbeta_post is accurate when the data conflict with the priors", {
  # priors with mean 0.02 for both proportions and 12 of 100 subjects with
  # the event; P(theta_d > theta_c) = 0.5 by symmetry
  f3 <- binbeta_post(y = 12, n = 100, k = 1, prior_c = c(1, 49),
                     prior_d = c(1, 49), p0 = 0.05)
  expect_equal(unname(f3$probs), c(0.6276917768448429, 0.5),
               tolerance = 1e-6)
  # all subjects with the event: the posterior of theta_d lies far in the
  # upper tail of its prior Beta(1, 99)
  f4 <- binbeta_post(y = 100, n = 100, k = 1, prior_c = c(1, 99),
                     prior_d = c(1, 99), p0 = 0.3, rr0 = 2)
  expect_equal(f4$prob, 0.7584271492339151, tolerance = 1e-6)
  expect_equal(f4$probs[["ratio"]], 0.010088886231998995, tolerance = 1e-5)
  expect_equal(unname(f4$rate[c("lower", "median", "upper")]),
               c(0.2292367577406037, 0.3385467317214596,
                 0.4420681123754414), tolerance = 1e-5)
  # few events: the posterior of theta_d lies in the lower tail of its prior
  # Beta(30, 270)
  f5 <- binbeta_post(y = 5, n = 1000, k = 1, prior_c = c(5, 495),
                     prior_d = c(30, 270), p0 = 0.04, rr0 = 8)
  expect_equal(unname(f5$probs), c(0.6669027639364123, 0.5333887099059609),
               tolerance = 1e-6)
  expect_equal(unname(f5$rr[c("lower", "median", "upper")]),
               c(3.4809226654317693, 8.335877831033217, 25.13355630126387),
               tolerance = 1e-5)
})

test_that("binbeta_post is accurate with a control prior shape below 1", {
  # case BS1 of inst/validation/python/ref_shape_below1.py (exact finite
  # mixture of products of Beta distributions)
  fit <- binbeta_post(y = 12, n = 100, k = 1, prior_c = c(0.5, 5),
                      prior_d = c(1, 1), p0 = 0.15, details = FALSE)
  expect_identical(fit$integration$coordinate, "normal")
  expect_equal(unname(fit$probs), c(0.7117667148324871, 0.8234809625922909),
               tolerance = 1e-6)
  expect_failure(expect_equal(fit$probs[["ratio"]],
                              0.8234809625922909 * (1 + 1e-4),
                              tolerance = 1e-6))
  expect_equal(unname(fit$estimates[c("theta_c", "theta_d")]),
               c(0.06158654927598209, 0.19698035646229917), tolerance = 1e-6)
  # the posterior density of theta_c near 0 is proportional to the prior
  # density, so the posterior mean of the ratio is infinite
  expect_identical(fit$rr[["mean"]], Inf)
  # the logit is kept for shape parameters of at least 1
  f2 <- binbeta_post(y = 12, n = 100, k = 1, prior_c = c(1, 49),
                     prior_d = c(1, 49), details = FALSE)
  expect_identical(f2$integration$coordinate, "logit")
})

test_that("binbeta_post handles no events and all events", {
  f0 <- binbeta_post(y = 0, n = 50, prior_c = c(5, 45), prior_d = c(1, 1),
                     p0 = 0.2, details = FALSE)
  f1 <- binbeta_post(y = 50, n = 50, prior_c = c(5, 45), prior_d = c(1, 1),
                     p0 = 0.2, details = FALSE)
  expect_lt(f0$prob, 0.05)
  expect_gt(f1$prob, 0.95)
  expect_error(binbeta_post(y = 51, n = 50, prior_c = c(5, 45),
                            prior_d = c(1, 1)), "'y'")
  expect_error(binbeta_post(y = 5, n = 50, prior_c = c(5, 45),
                            prior_d = c(1, 1), p0 = 1), "'p0'")
})
