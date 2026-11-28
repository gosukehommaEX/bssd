test_that("mix_envelope keeps small mixtures unchanged", {
  mix <- data.frame(w = c(0.7, 0.3), shape1 = c(20, 2), shape2 = c(180, 18))
  expect_identical(mix_envelope(mix), mix)
})

test_that("mix_envelope bounds the ratio of the densities", {
  prior <- map_beta(events = c(0, 0, 0), n = c(40, 40, 40))
  env <- mix_envelope(prior)
  expect_lt(nrow(env), nrow(prior))
  expect_equal(sum(env$w), 1, tolerance = 1e-12)
  x <- c(-10^seq(6, -2, length.out = 300), seq(-5, 5, length.out = 101),
         10^seq(-2, 3, length.out = 100))
  expect_lt(max(lmix_dens(x, prior) - lmix_dens(x, env)), log(10))
  # a single Beta distribution with lighter tails does not bound the ratio
  one <- data.frame(w = 1, shape1 = 1, shape2 = 50)
  expect_gt(max(lmix_dens(x, prior) - lmix_dens(x, one)), log(10))
})
