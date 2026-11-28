test_that("mix_logit_quantile inverts the distribution function", {
  q1 <- mix_logit_quantile(data.frame(w = 1, shape1 = 2, shape2 = 5))
  z <- c(-5, -1, 0, 0.7, 3)
  expect_equal(q1(z), qlogis(qbeta(pnorm(z), 2, 5)), tolerance = 1e-10)
  # envelope of a MAP prior from no historical events, with a heavy left
  # tail of the logit
  env <- mix_envelope(map_beta(events = c(0, 0, 0), n = c(40, 40, 40)))
  q2 <- mix_logit_quantile(env)
  z <- seq(-8.4, 8.4, by = 0.6)
  x <- q2(z)
  low <- z <= 0
  expect_equal(lmix_cdf(x[low], env), pnorm(z[low], log.p = TRUE),
               tolerance = 1e-10)
  expect_equal(lmix_cdf(x[!low], env, lower = FALSE),
               pnorm(-z[!low], log.p = TRUE), tolerance = 1e-10)
  expect_lt(min(x), -1e6)
  expect_true(all(diff(x) > 0))
})

test_that("mix_logit_quantile works for MAP priors from common events", {
  # envelopes of such priors have lower tails whose logarithm is close to 0
  # with rounding errors, which must not stop findInterval()
  env <- mix_envelope(map_beta(events = c(8, 12, 8), n = c(100, 100, 100),
                               control = list(check = FALSE)))
  q <- mix_logit_quantile(env)
  z <- seq(-8, 8, by = 0.5)
  x <- q(z)
  expect_true(all(is.finite(x)))
  expect_true(all(diff(x) > 0))
  low <- z <= 0
  expect_equal(lmix_cdf(x[low], env), pnorm(z[low], log.p = TRUE),
               tolerance = 1e-10)
})
