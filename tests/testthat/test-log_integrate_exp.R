test_that("log_integrate_exp integrates a narrow Gaussian", {
  lf <- function(z) -0.5 * ((z - 30.1234) / 0.01)^2
  expect_equal(log_integrate_exp(lf), log(sqrt(2 * pi) * 0.01),
               tolerance = 1e-9)
  expect_equal(log_integrate_exp(lf, 30.1234, Inf),
               log(sqrt(2 * pi) * 0.01 / 2), tolerance = 1e-9)
  # far tail: log of the normal tail probability times sqrt(2 * pi) * 0.01
  ref <- log(sqrt(2 * pi) * 0.01) +
    pnorm((31 - 30.1234) / 0.01, lower.tail = FALSE, log.p = TRUE)
  expect_equal(log_integrate_exp(lf, 31, 32), ref, tolerance = 1e-9)
  expect_identical(log_integrate_exp(lf, 2, 1), -Inf)
})

test_that("log_integrate_exp handles heavy exponential tails", {
  # integral of plogis(z)^a * plogis(-z)^b over the real line is beta(a, b)
  lf <- function(z) {
    0.5 * plogis(z, log.p = TRUE) + 0.5 * plogis(-z, log.p = TRUE)
  }
  expect_equal(log_integrate_exp(lf), log(pi), tolerance = 1e-9)
  lf2 <- function(z) 2 * plogis(z, log.p = TRUE) + 3 * plogis(-z, log.p = TRUE)
  expect_equal(log_integrate_exp(lf2), lbeta(2, 3), tolerance = 1e-9)
})
