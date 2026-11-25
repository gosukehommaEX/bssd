# Reference values: inst/validation/python/ref_mapbern.py (dense composite
# Gauss-Legendre grid over (logit a, logit b); for subject-level data, the
# arm likelihoods are integrated with scipy.integrate.quad).

test_that("map_beta matches the reference MAP prior for summary data", {
  prior <- map_beta(events = c(9, 11, 10), n = c(100, 100, 100))
  expect_s3_class(prior, "data.frame")
  expect_identical(names(prior), c("w", "shape1", "shape2"))
  expect_equal(sum(prior$w), 1, tolerance = 1e-12)
  expect_equal(unname(attr(prior, "moments")),
               c(0.10334757285043465, 0.03403531432783694), tolerance = 1e-6)
  expect_failure(expect_equal(unname(attr(prior, "moments"))[1],
                              0.10334757285043465 * (1 + 1e-4),
                              tolerance = 1e-6))
  expect_equal(pmixbeta(c(0.05, 0.08, 0.1, 0.12, 0.15), prior),
               c(0.03943187847234971, 0.25753992102464807,
                 0.49397605128441335, 0.712316272633485, 0.9082824788540658),
               tolerance = 1e-5)
  expect_equal(unname(attr(prior, "hyper")["b"]), 0.5199936162577368,
               tolerance = 1e-5)
  expect_identical(attr(prior, "n_total"), 300)
  expect_lt(attr(prior, "error"), 1e-5)
})

test_that("map_beta matches the reference for unequal arms", {
  prior <- map_beta(events = c(3, 9, 1, 6), n = c(40, 60, 30, 80))
  expect_equal(unname(attr(prior, "moments")),
               c(0.09143255411408796, 0.038375835893030445),
               tolerance = 1e-6)
})

test_that("map_beta with subject-level data agrees with summary data", {
  ev <- c(9, 11, 10)
  y <- unlist(lapply(ev, function(e) c(rep(1, e), rep(0, 100 - e))))
  study <- rep(1:3, each = 100)
  p1 <- map_beta(events = ev, n = c(100, 100, 100))
  p2 <- map_beta(y = y, t = rep(1, 300), study = study)
  expect_equal(attr(p2, "moments"), attr(p1, "moments"), tolerance = 1e-6)
})

test_that("map_beta matches the reference for subject-level data", {
  # the data of ref_mapbern.py (map_subject_data)
  ref <- c(mean = 0.1403765506737122, sd = 0.055414999422045655)
  dat <- map_subject_data()
  prior <- map_beta(y = dat$y, t = dat$t, study = dat$study)
  expect_equal(attr(prior, "moments"), ref, tolerance = 1e-5)
})

test_that("map_beta checks its input", {
  expect_error(map_beta(events = c(5, 200), n = c(100, 100)), "exceed")
  expect_error(map_beta(events = c(5, 2), n = c(100, 100), y = 1), "either")
  expect_error(map_beta(y = c(0, 1), t = c(1, 1)), "study")
})
