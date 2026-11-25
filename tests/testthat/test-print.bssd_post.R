test_that("print.bssd_post prints the main results and returns its input", {
  fit <- bdribs_post(y = 5, exposure = 500, k = 2, lambda_hp = 18 / 4000,
                     pt_hp = 4000, t_star = 1)
  out <- capture.output(res <- print(fit))
  expect_identical(res, fit)
  expect_true(any(grepl("Model: BDRIBS", out)))
  expect_true(any(grepl("P(r > 1 | data) = 0.8199", out, fixed = TRUE)))
  expect_true(any(grepl("Signal (posterior probability > 0.95): no", out,
                        fixed = TRUE)))
  expect_true(any(grepl("Cumulative incidence at t_star = 1", out,
                        fixed = TRUE)))
})

test_that("print.bssd_post prints the tested quantity of the new models", {
  fit <- binbeta_post(y = 40, n = 300, k = 1, prior_c = c(31, 271),
                      prior_d = c(1, 1), p0 = 0.15)
  out <- capture.output(res <- print(fit))
  expect_identical(res, fit)
  expect_true(any(grepl("Model: BINBETA", out, fixed = TRUE)))
  expect_true(any(grepl("Tested quantity: proportion theta_d of study drug",
                        out, fixed = TRUE)))
  expect_true(any(grepl("P(theta_d > 0.15 | data) =", out, fixed = TRUE)))
  expect_true(any(grepl("Posterior of rr: mean", out, fixed = TRUE)))
  expect_true(any(grepl("Posterior of theta_d: mean", out, fixed = TRUE)))
  fit2 <- multibin_post(y = c(3, 9), n = c(150, 300), prior_c = c(2, 98))
  out2 <- capture.output(print(fit2))
  expect_true(any(grepl("P(delta > 1 | data) =", out2, fixed = TRUE)))
})

test_that("print.bssd_post does not show the summaries that are NA", {
  # multibin_post() does not report the posterior mean of delta
  fit <- multibin_post(y = c(3, 9), n = c(150, 300), prior_c = c(2, 98))
  out <- capture.output(print(fit))
  expect_true(any(grepl("Posterior of delta: median", out, fixed = TRUE)))
  expect_false(any(grepl("NA", out, fixed = TRUE)))
  # without the quantiles only the mean is shown
  fit2 <- binbeta_post(y = 40, n = 300, prior_c = c(31, 271),
                       prior_d = c(1, 1), details = FALSE)
  out2 <- capture.output(print(fit2))
  expect_true(any(grepl("Posterior of rr: mean [0-9.]+$", out2)))
  expect_false(any(grepl("NA", out2, fixed = TRUE)))
})
