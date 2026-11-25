# Reference values: inst/validation/python/ref_binbeta.py (exact finite
# mixture of products of Beta distributions from the binomial expansion of
# the likelihood).

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
