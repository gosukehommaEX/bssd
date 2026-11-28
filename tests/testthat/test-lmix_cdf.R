test_that("lmix_cdf matches pbeta", {
  mix <- data.frame(w = c(0.6, 0.4), shape1 = c(2, 0.3), shape2 = c(5, 4))
  # each tail where pbeta() computes it accurately
  x <- c(-30, -3, 0, 1.5, 3)
  p <- plogis(x)
  lo <- log(0.6 * pbeta(p, 2, 5) + 0.4 * pbeta(p, 0.3, 4))
  expect_equal(lmix_cdf(x, mix), lo, tolerance = 1e-12)
  expect_failure(expect_equal(lmix_cdf(x, mix), lo * (1 + 1e-8),
                              tolerance = 1e-12))
  x <- c(-3, 0, 1.5, 3, 6)
  p <- plogis(x)
  up <- log(0.6 * pbeta(p, 2, 5, lower.tail = FALSE) +
              0.4 * pbeta(p, 0.3, 4, lower.tail = FALSE))
  expect_equal(lmix_cdf(x, mix, lower = FALSE), up, tolerance = 1e-12)
})

test_that("lmix_cdf is continuous where the asymptotic tails start", {
  s <- data.frame(w = 1, shape1 = 0.3, shape2 = 4)
  expect_equal(lmix_cdf(-700, s), lmix_cdf(-700 + 1e-10, s),
               tolerance = 1e-12)
  expect_equal(lmix_cdf(700, s, lower = FALSE),
               lmix_cdf(700 - 1e-10, s, lower = FALSE), tolerance = 1e-12)
  expect_true(all(is.finite(lmix_cdf(c(-1e20, 1e20), s))))
  expect_true(all(is.finite(lmix_cdf(c(-1e20, 1e20), s, lower = FALSE))))
})

test_that("lmix_cdf is finite where pbeta() underflows", {
  # tail probabilities far below the smallest double with large shape
  # parameters: the leading term of the series is used
  s <- data.frame(w = 1, shape1 = 30, shape2 = 300)
  x <- seq(-720, -640, by = 10)
  lo <- lmix_cdf(x, s)
  expect_true(all(is.finite(lo)))
  expect_true(all(diff(lo) > 0))
  up <- lmix_cdf(-x, s, lower = FALSE)
  expect_true(all(is.finite(up)))
  expect_true(all(diff(up) > 0))
})
