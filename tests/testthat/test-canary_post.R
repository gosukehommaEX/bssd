# Reference values: inst/validation/python/ref_canary.py (binomial: composite
# Gauss-Legendre rules on (logit p_c, log M) or (logit p_c, logit p_d),
# split exactly at the threshold and, for the ratio, at the truncation point
# of each row; Poisson: adaptive quadrature).

test_that("canary_post matches the reference values for the ratio", {
  fit <- canary_post(y = 12, n = 100, k = 1, prior_c = c(1, 49))
  expect_s3_class(fit, "bssd_post")
  expect_identical(fit$model, "CANARY")
  expect_identical(fit$rr_label, "R")
  expect_equal(fit$prob, 0.8416853521605938, tolerance = 1e-6)
  expect_failure(expect_equal(fit$prob, 0.8416853521605938 * (1 + 1e-4),
                              tolerance = 1e-6))
  expect_equal(unname(fit$rr[c("lower", "median", "upper")]),
               c(0.35212880033048716, 2.3890803857663077,
                 10.534896851962829), tolerance = 1e-5)
  expect_equal(unname(fit$rr["mean"]), 3.1557892713153413, tolerance = 1e-5)
  expect_equal(unname(fit$estimates["p_c"]), 0.06121356291115411,
               tolerance = 1e-6)
  expect_lt(fit$integration$error, 1e-5)
})

test_that("canary_post matches the reference values for the difference", {
  fit <- canary_post(y = 12, n = 100, k = 1, prior_c = c(1, 49),
                     metric = "difference", crit = 0.02)
  expect_identical(fit$rr_label, "D")
  expect_equal(fit$prob, 0.9924243970067982, tolerance = 1e-6)
  expect_equal(unname(fit$rr[c("lower", "median", "upper")]),
               c(0.06536359178201151, 0.21129425233471963,
                 0.35038996637360437), tolerance = 1e-5)
  expect_equal(unname(fit$rr["mean"]), 0.21079059416469043,
               tolerance = 1e-5)
})

test_that("canary_post matches the reference values for the odds ratio", {
  fit <- canary_post(y = 20, n = 120, k = 2, prior_c = c(5, 45),
                     metric = "odds", crit = 1.2)
  expect_identical(fit$rr_label, "Q")
  expect_equal(fit$prob, 0.6197084082546314, tolerance = 1e-6)
  expect_equal(unname(fit$rr[c("lower", "median", "upper")]),
               c(0.3932327219238604, 1.4472301845829256,
                 4.597715602699213), tolerance = 1e-5)
  expect_equal(unname(fit$estimates["p_c"]), 0.12701718711647839,
               tolerance = 1e-6)
})

test_that("canary_post matches the reference values with Poisson counts", {
  fit <- canary_post(y = 250, n = 20000, k = 1, prior_c = c(75.71, 8500),
                     likelihood = "poisson", crit = 2)
  expect_equal(fit$prob, 0.22133664341858217, tolerance = 1e-6)
  expect_equal(unname(fit$rr[c("lower", "median", "upper")]),
               c(1.1159461907040649, 1.7135203692142589,
                 2.5190360486264805), tolerance = 1e-5)
  expect_equal(unname(fit$estimates[c("lambda_c", "metric")]),
               c(0.009160493109283355, 1.7405131865276489),
               tolerance = 1e-6)
})

test_that("canary_post is accurate when all subjects had the event", {
  # the posterior of p_c lies far above its prior, and the density of log(R)
  # is largest at its limit -log(p_c)
  fit <- canary_post(y = 100, n = 100, k = 1, prior_c = c(1, 49), crit = 2)
  expect_equal(fit$prob, 0.30070256225719694, tolerance = 1e-6)
  expect_equal(unname(fit$rr["mean"]), 1.9155563213386304, tolerance = 1e-6)
  expect_equal(unname(fit$estimates["p_c"]), 0.5163759207909249,
               tolerance = 1e-6)
  expect_lt(fit$integration$error, 1e-6)
})

test_that("canary_post uses the lognormal prior of Gould and Wang by default", {
  f1 <- canary_post(y = 12, n = 100, prior_c = c(1, 49), details = FALSE)
  f2 <- canary_post(y = 12, n = 100, prior_c = c(1, 49), details = FALSE,
                    metric_prior = c(-0.5475867818334186, 1.0465054054647005))
  expect_null(names(f1$prob))
  expect_equal(f1$prob, f2$prob, tolerance = 1e-9)
  expect_equal(f1$prior$metric_prior,
               c(meanlog = -0.5475867818334186, sdlog = 1.0465054054647005),
               tolerance = 1e-9)
})

test_that("canary_post checks its input", {
  expect_error(canary_post(y = 12, n = 100, prior_c = c(1, 49),
                           metric = "odds", likelihood = "poisson"),
               "Poisson")
  expect_error(canary_post(y = 120, n = 100, prior_c = c(1, 49)), "'y'")
  expect_error(canary_post(y = 12, n = 100, prior_c = c(1, 49),
                           metric_prior = c(0, -1)), "sdlog")
  expect_error(canary_post(y = 12, n = 100, prior_c = c(1, 49),
                           metric = "difference", crit = 2), "'crit'")
})
