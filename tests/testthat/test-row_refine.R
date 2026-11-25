test_that("row_refine finds the interval of each row", {
  # v | u ~ N(3 u, 0.1^2) for u = -1, 0 and 2
  u <- c(-1, 0, 2)
  lp <- function(u, v) dnorm(v, 3 * u, 0.1, log = TRUE)
  b <- row_refine(lp, u, rep(-50, 3), rep(50, 3), 41, 1e-8, 12,
                  rep(-Inf, 3), rep(Inf, 3))
  expect_true(all(b$ok))
  # each interval covers at least 5.5 and at most 15 standard deviations on
  # each side of the mean of its row
  expect_true(all(b$lo < 3 * u - 0.55 & b$lo > 3 * u - 1.5))
  expect_true(all(b$hi > 3 * u + 0.55 & b$hi < 3 * u + 1.5))
})

test_that("row_refine respects the limits of each row", {
  # a density that increases up to a limit that differs between the rows
  u <- c(0, 1)
  lp <- function(u, v) 5 * v
  lim_hi <- c(0.5, 2)
  b <- row_refine(lp, u, c(-10, -10), c(10, 10), 41, 1e-8, 12,
                  c(-Inf, -Inf), lim_hi)
  expect_identical(b$hi, lim_hi)
  # the mass below the lower limit, exp(-5 * (hi - lo)), is below 1e-8
  expect_true(all(b$lo < lim_hi - 3.6 & b$lo > lim_hi - 6))
  # an initial interval above the limit is moved below it
  b2 <- row_refine(lp, u, c(5, 5), c(10, 10), 41, 1e-8, 12, c(-Inf, -Inf),
                   lim_hi)
  expect_identical(b2$hi, lim_hi)
  expect_true(all(b2$lo < lim_hi - 3.6 & b2$lo > lim_hi - 6))
})
