# Reference values were computed independently in Python: all four
# parameters on their prior distribution-function scale with Gauss-Legendre
# nodes (48 per coordinate for the means and probabilities, 32 for the
# quantiles), the indicator of RR(t_star) > q handled by splitting the mu_d
# coordinate at the threshold, and no box refinement or Chebyshev
# interpolation. The same values were obtained within 3e-4 by a separate Python
# implementation of the algorithm of the package.

make_data <- function() {
  list(
    t = c(66.5, 93.7, 81.6, 26.5, 34, 91.4, 4.5, 86.1, 83.7, 50.8, 34.3,
      31.8, 29.5, 48.5, 54.5, 59.3, 103.6, 83.3, 66.2, 102.9, 25.5, 20, 65.3,
      8.4, 7.6, 55.5, 50.6, 95.7, 66.9, 55.4, 53.7, 28.8, 5.2, 23.2, 73.2,
      24.1, 41, 4.4, 87, 19.4, 30.8, 92, 55, 88.7, 68, 78.2, 13.1, 58.1,
      54.8, 91.1, 40.1, 63.8, 9.9, 42.8, 36.3, 19, 85.6, 41.9, 101.9, 63),
    y = c(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0)
  )
}

test_that("llbern_post matches the reference values", {
  d <- make_data()
  fit <- llbern_post(d$y, d$t, k = 1, mu_hp = 6, sigma_hp = 0.8, pt_hp = 200,
                     df_p = 1, t_star = 104)
  expect_s3_class(fit, "bssd_post")
  expect_equal(fit$prob, 0.6964229, tolerance = 1e-3)
  # the comparison is not vacuous: a shifted reference value fails
  expect_failure(expect_equal(fit$prob, 0.6964229 * (1 + 1e-2),
                              tolerance = 1e-3))
  expect_equal(unname(fit$rr[c("median", "upper")]),
               c(1.5201310901064538, 4.273130301999516), tolerance = 2e-3)
  expect_true(fit$rr["lower"] > 0 && fit$rr["lower"] < fit$rr["median"])
  expect_equal(unname(fit$rr["mean"]), 1.6384462, tolerance = 2e-3)
  expect_equal(unname(fit$estimates),
               c(5.9799629, 0.8079559, 6.4949205, 1.3470213),
               tolerance = 2e-3)
  expect_equal(unname(fit$incidence), c(0.162272, 0.2430457, 0.0807737),
               tolerance = 2e-3)
  expect_lt(fit$integration$error, 1e-3)
})

test_that("llbern_post matches the reference values for other settings", {
  d <- make_data()
  f2 <- llbern_post(d$y, d$t, k = 1, mu_hp = 6, sigma_hp = 0.8, pt_hp = 200,
                    df_p = 1, t_star = 104, rr0 = 2, details = FALSE,
                    control = list(check = FALSE))
  expect_equal(f2$prob, 0.33036, tolerance = 2e-3)
  expect_true(all(is.na(f2$rr[c("median", "lower", "upper")])))
  f3 <- llbern_post(d$y, d$t, k = 2, mu_hp = 6, sigma_hp = 0.8, pt_hp = 200,
                    df_p = 2, t_star = 52, rr0 = 1.5,
                    control = list(check = FALSE))
  expect_equal(f3$prob, 0.6092996462522098, tolerance = 1e-3)
})

test_that("llbern_post reports twice the difference from the coarser grid", {
  d <- make_data()
  args <- list(y = d$y, t = d$t, k = 1, mu_hp = 6, sigma_hp = 0.8,
               pt_hp = 200, df_p = 1, t_star = 104, details = FALSE)
  fit <- do.call(llbern_post, args)
  f0 <- do.call(llbern_post, c(args, list(control = list(check = FALSE))))
  f1 <- do.call(llbern_post, c(args, list(control = list(
    n_outer = 20, n_inner = 181, check = FALSE))))
  expect_equal(fit$prob, f1$prob)
  expect_equal(fit$integration$error, 2 * abs(f1$prob - f0$prob))
  # the comparison is not vacuous: the difference itself is not accepted
  expect_failure(expect_equal(fit$integration$error,
                              abs(f1$prob - f0$prob)))
})

test_that("llbern_post checks its input", {
  d <- make_data()
  expect_error(llbern_post(d$y, d$t, mu_hp = 6, sigma_hp = 0.8, pt_hp = 200),
               "t_star")
  expect_error(llbern_post(d$y, d$t, mu_hp = 6, sigma_hp = 0, pt_hp = 200,
                           t_star = 104), "out of range")
  expect_error(llbern_post(d$y, d$t, mu_hp = 6, sigma_hp = 0.8, pt_hp = 200,
                           t_star = 104, control = list(n = 3)), "Unknown")
})
