# Reference values: inst/validation/python/ref_misc.py (posterior
# probabilities of canary_post() at counts 5 and 6, and 16 and 17, by the
# dense rule of ref_canary.py).

test_that("signal_count finds the critical count of Gould and Wang", {
  # Table III: prior Beta(20, 980), R > 1 with probability 0.8, N = 100
  res <- signal_count(n = 100, method = "canary", prior_c = c(20, 980),
                      cutoff = 0.8)
  expect_identical(res$y_crit, 6)
  expect_equal(res$prob, 0.842132592597768, tolerance = 1e-6)
  expect_equal(res$prob_below, 0.7167683667445133, tolerance = 1e-6)
})

test_that("signal_count stops at the first count with a signal", {
  # prior Beta(1, 49), R > 2 with probability 0.8, N = 100: the probability
  # exceeds 0.8 from the count 17, and is below 0.8 again at the count 100
  res <- signal_count(n = 100, method = "canary", prior_c = c(1, 49),
                      crit = 2, cutoff = 0.8)
  expect_identical(res$y_crit, 17)
  expect_equal(res$prob, 0.8143444542663405, tolerance = 1e-6)
  expect_equal(res$prob_below, 0.7781208918194892, tolerance = 1e-6)
})

test_that("signal_count works for binbeta_post and for several sizes", {
  res <- signal_count(n = c(50, 100), method = "binbeta",
                      prior_c = c(31, 271), prior_d = c(1, 1), p0 = 0.15,
                      cutoff = 0.8)
  expect_identical(nrow(res), 2L)
  for (j in 1:2) {
    yc <- res$y_crit[j]
    p1 <- binbeta_post(yc, res$n[j], prior_c = c(31, 271), prior_d = c(1, 1),
                       p0 = 0.15, details = FALSE)$prob
    p0 <- binbeta_post(yc - 1, res$n[j], prior_c = c(31, 271),
                       prior_d = c(1, 1), p0 = 0.15, details = FALSE)$prob
    expect_true(p1 > 0.8 && p0 <= 0.8)
  }
  expect_error(signal_count(100, "canary", prior_c = c(1, 49), y = 3),
               "set by")
  expect_error(signal_count(100, "canary", prior_c = c(1, 49),
                            likelihood = "poisson"), "y_max")
})
