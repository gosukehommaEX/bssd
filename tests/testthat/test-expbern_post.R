# Reference values were computed independently in Python on a dense fixed
# Gauss-Legendre grid over (log lambda, logit p), split at logit(p0),
# without the box refinement and Chebyshev interpolation of the package.

make_data <- function() {
  list(
    t = c(38.3, 50.3, 31.8, 45.1, 38.6, 8.5, 19.3, 42.3, 38.7, 40.7, 25.5, 46.3,
      10.7, 5, 22.9, 46.9, 40, 30.7, 40.9, 9.1, 16.2, 50.6, 31.6, 33.1, 40.6,
      20.6, 13.8, 25, 42.5, 30.1, 16.4, 15, 27.2, 47.9, 23.8, 9.7, 7.9, 6.8,
      39.6, 14, 43.4, 29.5, 24.7, 38, 39.2, 25, 36, 25.8, 27.4, 48.6, 44.2,
      6, 6.6, 33.1, 43.9, 32.4, 49.9, 44.5, 26.3, 47.7),
    y = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1,
      0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
      0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  )
}

test_that("expbern_post matches the reference values with Beta(0.5, 0.5) prior", {
  d <- make_data()
  fit <- expbern_post(d$y, d$t, k = 1, lambda_hp = 0.004, pt_hp = 2000,
                      df = 2, rr0 = 1, t_star = 52)
  expect_s3_class(fit, "bssd_post")
  expect_equal(fit$prob, 0.6313681141356096, tolerance = 1e-6)
  # the comparison is not vacuous: a shifted reference value fails
  expect_failure(expect_equal(fit$prob, 0.6313681141356096 * (1 + 1e-4),
                              tolerance = 1e-6))
  expect_equal(unname(fit$rr[c("median", "lower", "upper")]),
               c(1.4773060987634894, 0.01768901635499989,
                 9.242458364590444), tolerance = 1e-5)
  expect_equal(unname(fit$estimates[c("p", "lambda_p")]),
               c(0.5488989493123925, 0.0047056533207496705),
               tolerance = 1e-6)
  expect_equal(unname(fit$incidence),
               c(0.21224098863287888, 0.2857866446967653,
                 0.07354565606388641), tolerance = 1e-6)
  # the posterior means of r and lambda_d are infinite when b <= 1
  expect_identical(unname(fit$estimates[c("r", "lambda_d")]), c(Inf, Inf))
  expect_lt(fit$integration$error, 1e-4)
})

test_that("expbern_post matches the reference values with k = 2 and rr0 = 1.5", {
  d <- make_data()
  fit <- expbern_post(d$y, d$t, k = 2, lambda_hp = 0.005, pt_hp = 3000,
                      df = 1, p_prior = c(1, 1), rr0 = 1.5, t_star = 26)
  expect_s3_class(fit, "bssd_post")
  expect_equal(fit$prob, 0.26752577330727506, tolerance = 1e-6)
  expect_equal(unname(fit$rr[c("median", "lower", "upper")]),
               c(1.0411213957847538, 0.1920458923174792,
                 3.033535502130235), tolerance = 1e-5)
  expect_equal(unname(fit$estimates[c("p", "lambda_p")]),
               c(0.6480464328256317, 0.0052915273529549915),
               tolerance = 1e-6)
  expect_equal(unname(fit$incidence),
               c(0.12802702563099594, 0.136190599702356,
                 0.008163574071360064), tolerance = 1e-6)
  # the posterior means of r and lambda_d are infinite when b <= 1
  expect_identical(unname(fit$estimates[c("r", "lambda_d")]), c(Inf, Inf))
  expect_lt(fit$integration$error, 1e-4)
})

test_that("expbern_post matches the reference values with Beta(1, 3) prior", {
  d <- make_data()
  fit <- expbern_post(d$y, d$t, k = 1, lambda_hp = 0.004, pt_hp = 1000,
                      df = 1, p_prior = c(1, 3), rr0 = 1, t_star = 52)
  expect_s3_class(fit, "bssd_post")
  expect_equal(fit$prob, 0.29560487790378887, tolerance = 1e-4)
  expect_equal(unname(fit$rr[c("median", "lower", "upper")]),
               c(0.5934364836828607, 0.026881672076536658,
                 3.2104777171292267), tolerance = 1e-4)
  expect_equal(unname(fit$estimates[c("p", "lambda_p")]),
               c(0.3752890244665631, 0.005920649289616942),
               tolerance = 1e-4)
  expect_equal(unname(fit$incidence),
               c(0.2608934271794771, 0.17993734650825244,
                 -0.08095608067122467), tolerance = 1e-4)
  expect_equal(unname(fit$estimates[c("r", "lambda_d")]),
               c(0.8530456928091014, 0.0040438138780992785),
               tolerance = 1e-3)
  expect_lt(fit$integration$error, 1e-4)
})

test_that("expbern_post works without the convergence check", {
  d <- make_data()
  f1 <- expbern_post(d$y, d$t, lambda_hp = 0.004, pt_hp = 2000, df = 2,
                     control = list(check = FALSE))
  expect_true(is.na(f1$integration$error))
  expect_null(f1$incidence)
  expect_equal(f1$prob, 0.6313681141, tolerance = 1e-6)
  f2 <- expbern_post(d$y, d$t, lambda_hp = 0.004, pt_hp = 2000, df = 2,
                     details = FALSE, control = list(check = FALSE))
  expect_equal(f2$prob, f1$prob)
  expect_true(all(is.na(f2$rr[c("median", "lower", "upper")])))
})

test_that("expbern_post checks its input", {
  d <- make_data()
  expect_error(expbern_post(d$y, d$t[-1], lambda_hp = 0.004, pt_hp = 2000),
               "length")
  expect_error(expbern_post(d$y, d$t, lambda_hp = 0.004, pt_hp = 2000,
                            control = list(n_z = 10)), "Unknown")
})
