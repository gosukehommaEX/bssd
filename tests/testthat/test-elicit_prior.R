# Reference values: inst/validation/python/ref_misc.py (two-dimensional root
# finding with scipy.optimize.fsolve).

test_that("elicit_prior reproduces the Gamma prior of Schnell and Ball", {
  p <- elicit_prior(c(0.5, 2), c(0.5, 0.75), "gamma")
  expect_identical(names(p), c("shape", "rate"))
  expect_equal(unname(p), c(0.346004956477767, 0.20898323220514092),
               tolerance = 1e-8)
  expect_equal(pgamma(c(0.5, 2), p[["shape"]], p[["rate"]]), c(0.5, 0.75),
               tolerance = 1e-10)
  # Gould and Wang (2017), Table I: median 1, ratio 2, zeta = 0.8
  p2 <- elicit_prior(c(1, 2), c(0.5, 0.8), "gamma")
  expect_equal(unname(p2), c(1.4233581235293404, 1.1073865262540488),
               tolerance = 1e-8)
})

test_that("elicit_prior finds Beta and lognormal priors", {
  p <- elicit_prior(c(0.02, 0.04), c(0.5, 0.95), "beta")
  expect_identical(names(p), c("shape1", "shape2"))
  expect_equal(unname(p), c(4.549595083088185, 207.16232292428424),
               tolerance = 1e-7)
  q <- elicit_prior(c(0.5, 4), c(0.1, 0.9), "lnorm")
  expect_equal(unname(q), c(0.3465735902799726, 0.8112984282439579),
               tolerance = 1e-10)
  m <- elicit_prior(c(0.1, 10), 0.95, "lnorm", mean = 1)
  expect_equal(unname(m), c(-0.5475867818334186, 1.0465054054647005),
               tolerance = 1e-9)
  m8 <- elicit_prior(c(0.1, 8), 0.95, "lnorm", mean = 1)
  expect_equal(unname(m8), c(-0.5369581373824451, 1.036299317168978),
               tolerance = 1e-9)
})

test_that("elicit_prior checks its input", {
  expect_error(elicit_prior(c(2, 1), c(0.5, 0.75)), "increasing")
  expect_error(elicit_prior(c(0.1, 10), 0.95, "gamma", mean = 1), "lnorm")
  expect_error(elicit_prior(c(0.1, 0.5), 0.95, "lnorm", mean = 1),
               "contains")
  expect_error(elicit_prior(c(0.2, 2), c(0.5, 0.75), "beta"), "'q'")
})
