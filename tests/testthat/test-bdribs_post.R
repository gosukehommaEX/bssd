# Reference values were computed independently in Python by adaptive
# quadrature over p (scipy.integrate.quad); the incidences of case A were
# also confirmed by two-dimensional integration over (delta, p).

test_that("bdribs_post matches the reference values for case A (k = 2)", {
  fit <- bdribs_post(y = 5, exposure = 500, k = 2, lambda_hp = 18 / 4000,
                     pt_hp = 4000, df = 1, rr0 = 1, t_star = 1)
  expect_s3_class(fit, "bssd_post")
  expect_equal(fit$prob, 0.8198879688642029, tolerance = 1e-7)
  # the comparison is not vacuous: a shifted reference value fails
  expect_failure(expect_equal(fit$prob, 0.8198879688642029 * (1 + 1e-5),
                              tolerance = 1e-7))
  expect_equal(unname(fit$rr), c(2.4506025985748776, 2.123375405874829,
                              0.19502828213481802, 6.588037134695716),
               tolerance = 1e-6)
  expect_equal(unname(fit$estimates),
               c(0.763240530730424, 2.4506025985748776,
                 0.004667289090498797, 0.010658886368765089),
               tolerance = 1e-7)
  expect_equal(unname(fit$incidence),
               c(0.004655827293755355, 0.01058191207375836,
                 0.005926084780003005),
               tolerance = 1e-7)
})

test_that("bdribs_post matches the reference values for case B (rr0 = 1.2, df = 2)", {
  fit <- bdribs_post(y = 30, exposure = 8000, k = 1, lambda_hp = 0.002,
                     pt_hp = 10000, df = 2, rr0 = 1.2, t_star = 52)
  expect_s3_class(fit, "bssd_post")
  expect_equal(fit$prob, 0.8989448489448683, tolerance = 1e-7)
  expect_equal(unname(fit$rr), c(2.7703125402513282, 2.480943072102248,
                              0.7265536698551257, 6.499772815937226),
               tolerance = 1e-6)
  expect_equal(unname(fit$estimates),
               c(0.6949784588244851, 2.7703125402513282,
                 0.002149513854898876, 0.005163593826477517),
               tolerance = 1e-7)
  expect_equal(unname(fit$incidence),
               c(0.10521650657678538, 0.23305854503242118,
                 0.1278420384556358),
               tolerance = 1e-7)
})

test_that("bdribs_post matches the reference values for case C (no events, uniform prior)", {
  fit <- bdribs_post(y = 0, exposure = 300, k = 1, lambda_hp = 0.01,
                     pt_hp = 500, df = 1, p_prior = c(1, 1), rr0 = 1,
                     t_star = 10)
  expect_s3_class(fit, "bssd_post")
  expect_equal(fit$prob, 0.13823124712087123, tolerance = 1e-7)
  expect_equal(unname(fit$rr), c(0.5120391977307874, 0.3006361296511069,
                              0.010103644553785414, 2.2616702781822595),
               tolerance = 1e-6)
  expect_equal(unname(fit$estimates),
               c(0.2659706017019098, 0.5120391977307874,
                 0.006972217236005878, 0.0031203919773078737),
               tolerance = 1e-7)
  expect_equal(unname(fit$incidence),
               c(0.06687032234964418, 0.030100291817997005,
                 -0.03677003053164718),
               tolerance = 1e-7)
})

test_that("bdribs_post behaves sensibly", {
  p <- vapply(c(2, 5, 10), function(y) {
    bdribs_post(y = y, exposure = 500, k = 1, lambda_hp = 0.005,
                pt_hp = 2000)$prob
  }, numeric(1))
  expect_true(all(diff(p) > 0))
  fit <- bdribs_post(y = 5, exposure = 500, k = 1, lambda_hp = 0.005,
                     pt_hp = 2000, cutoff = 0.5)
  expect_identical(fit$signal, fit$prob > 0.5)
  expect_null(fit$incidence)
  expect_true(fit$rr["lower"] < fit$rr["median"] &&
                fit$rr["median"] < fit$rr["upper"])
  # P(r > rr0) decreases with rr0
  f2 <- bdribs_post(y = 5, exposure = 500, k = 1, lambda_hp = 0.005,
                    pt_hp = 2000, rr0 = 2)
  expect_lt(f2$prob, fit$prob)
})

test_that("bdribs_post checks its input", {
  expect_error(bdribs_post(y = -1, exposure = 500, lambda_hp = 0.01,
                           pt_hp = 100), "out of range")
  expect_error(bdribs_post(y = 1.5, exposure = 500, lambda_hp = 0.01,
                           pt_hp = 100), "integer")
  expect_error(bdribs_post(y = 1, exposure = 0, lambda_hp = 0.01,
                           pt_hp = 100), "out of range")
  expect_error(bdribs_post(y = 1, exposure = 10, lambda_hp = 0.01,
                           pt_hp = 100, cutoff = 1), "out of range")
})

# Reference values for a known background rate (pt_hp = Inf) were computed in
# Python by adaptive quadrature over p (inst/validation/python/
# ref_bdribs_fixed.py). Cases F1 and F2 are hypothetical example 2 of
# Mukhopadhyay et al. (2018) with nine events.

test_that("bdribs_post with a known background rate matches case F1", {
  fit <- bdribs_post(y = 9, exposure = 540, k = 2, lambda_hp = 0.011,
                     pt_hp = Inf, p_prior = c(1, 1), rr0 = 1, t_star = 0.5,
                     level = 0.9)
  expect_equal(fit$prob, 0.7529735703820455, tolerance = 1e-7)
  expect_equal(unname(fit$rr), c(1.521837857620912, 1.4375784755419578,
                              0.509441171140254, 2.8207258819915),
               tolerance = 1e-6)
  expect_equal(unname(fit$estimates),
               c(0.7180682737398301, 1.521837857620912,
                 0.011, 0.01674021643383003),
               tolerance = 1e-7)
  expect_equal(unname(fit$incidence),
               c(0.005484902691080928, 0.008327562899862584,
                 0.0028426602087816555),
               tolerance = 1e-7)
})

test_that("bdribs_post with a known background rate matches case F2", {
  fit <- bdribs_post(y = 9, exposure = 540, k = 2, lambda_hp = 0.011,
                     pt_hp = Inf, p_prior = c(1, 1), rr0 = 1.3, t_star = 0.5,
                     level = 0.9)
  expect_equal(fit$prob, 0.5802463151874379, tolerance = 1e-7)
  expect_equal(unname(fit$rr), c(1.521837857620912, 1.4375784755419578,
                              0.509441171140254, 2.820725881991503),
               tolerance = 1e-6)
  expect_equal(unname(fit$estimates),
               c(0.7180682737398303, 1.521837857620912,
                 0.011, 0.01674021643383003),
               tolerance = 1e-7)
  expect_equal(unname(fit$incidence),
               c(0.005484902691080928, 0.008327562899862362,
                 0.0028426602087814334),
               tolerance = 1e-7)
})

test_that("bdribs_post with a known background rate matches case F3", {
  fit <- bdribs_post(y = 3, exposure = 900, k = 1, lambda_hp = 0.004,
                     pt_hp = Inf, rr0 = 1.2, t_star = 1)
  expect_equal(fit$prob, 0.15757098956025947, tolerance = 1e-7)
  expect_equal(unname(fit$rr), c(0.5944156591825007, 0.34667188674653143,
                              0.0008707915521784481, 2.4930188814046663),
               tolerance = 1e-6)
  expect_equal(unname(fit$estimates),
               c(0.28497409326424966, 0.5944156591825007,
                 0.004, 0.002377662636730003),
               tolerance = 1e-7)
  expect_equal(unname(fit$incidence),
               c(0.003992010656008516, 0.0023709679336300082,
                 -0.0016210427223785073),
               tolerance = 1e-7)
})

test_that("bdribs_post with details = FALSE returns the probability only", {
  f1 <- bdribs_post(y = 5, exposure = 500, k = 2, lambda_hp = 18 / 4000,
                    pt_hp = 4000, t_star = 1)
  f2 <- bdribs_post(y = 5, exposure = 500, k = 2, lambda_hp = 18 / 4000,
                    pt_hp = 4000, t_star = 1, details = FALSE)
  expect_equal(f2$prob, f1$prob)
  expect_equal(f2$incidence, f1$incidence)
  expect_true(all(is.na(f2$rr)) && all(is.na(f2$estimates)))
  expect_error(bdribs_post(y = 5, exposure = 500, lambda_hp = 0.01,
                           pt_hp = 100, details = NA), "TRUE or FALSE")
})
