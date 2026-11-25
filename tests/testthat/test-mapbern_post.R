# Reference values: inst/validation/python/ref_mapbern.py (composite
# Gauss-Legendre rules on (logit pi_c, logit pi_d), split exactly at
# pi_d = rr0 * pi_c).

test_that("mapbern_post matches the reference values", {
  d <- data_m1()
  fit <- mapbern_post(d$y, d$t, k = 3, prior_c = mix_m1())
  expect_s3_class(fit, "bssd_post")
  expect_identical(fit$model, "MAPBERN")
  expect_equal(fit$prob, 0.956972148630543, tolerance = 2e-6)
  expect_failure(expect_equal(fit$prob, 0.956972148630543 * (1 + 1e-4),
                              tolerance = 2e-6))
  expect_equal(unname(fit$estimates[c("pi_c", "pi_d")]),
               c(0.10141842634504168, 0.2635633972184154), tolerance = 1e-6)
  expect_equal(unname(fit$rate[c("lower", "median", "upper")]),
               c(0.11685489263341933, 0.2588641424427206,
                 0.4366714510236672), tolerance = 1e-5)
  expect_equal(unname(fit$rr[c("lower", "median", "upper")]),
               c(0.8211048286636914, 2.660555390969741, 11.177990376362308),
               tolerance = 1e-5)
})

test_that("mapbern_post matches the reference value with rr0 = 1.5", {
  d <- data_m1()
  fit <- mapbern_post(d$y, d$t, k = 3, prior_c = mix_m1(), rr0 = 1.5,
                      details = FALSE)
  expect_equal(fit$prob, 0.8670458703679775, tolerance = 2e-6)
})

test_that("mapbern_post matches the reference with an informative prior_d", {
  # prior of pi_d with mean 0.02, which conflicts with the data
  d <- data_m1()
  fit <- mapbern_post(d$y, d$t, k = 3, prior_c = mix_m1(), prior_d = c(2, 98),
                      rr0 = 0.25)
  expect_equal(fit$prob, 0.4532473341175618, tolerance = 2e-6)
  expect_equal(unname(fit$estimates[c("pi_c", "pi_d")]),
               c(0.21092135853218738, 0.045348003438555634), tolerance = 1e-6)
  expect_equal(unname(fit$rate[c("lower", "median", "upper")]),
               c(0.006748974569721206, 0.04194418468409068,
                 0.1035246193589002), tolerance = 1e-5)
})

test_that("mapbern_post accepts a single Beta prior and a map_beta() prior", {
  d <- data_m1()
  f1 <- mapbern_post(d$y, d$t, k = 3, prior_c = c(20, 180), details = FALSE)
  f2 <- mapbern_post(d$y, d$t, k = 3,
                     prior_c = data.frame(w = 1, shape1 = 20, shape2 = 180),
                     details = FALSE)
  expect_equal(f1$prob, f2$prob, tolerance = 1e-12)
  prior <- map_beta(events = c(9, 11, 10), n = c(100, 100, 100))
  f3 <- mapbern_post(d$y, d$t, k = 3, prior_c = prior, details = FALSE)
  expect_true(f3$prob > 0 && f3$prob < 1)
})

test_that("mapbern_post checks its input", {
  d <- data_m1()
  expect_error(mapbern_post(d$y, d$t, prior_c = list(1, 2)), "prior_c")
  expect_error(mapbern_post(d$y + 1, d$t, prior_c = c(1, 9)), "'y'")
  expect_error(mapbern_post(d$y, d$t, prior_c = c(1, 9), rate0 = 2),
               "'rate0'")
})
