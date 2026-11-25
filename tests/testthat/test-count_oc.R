# Reference values: inst/validation/python/ref_misc.py (convolution of the
# binomial probability functions, and binomial and Poisson tails).

test_that("count_oc matches the reference values", {
  expect_equal(count_oc(6, 100, rate_c = 0.02, rate_d = c(0.02, 0.04, 0.08)),
               c(0.015483640641779373, 0.08048694600912454,
                 0.38392899931647617), tolerance = 1e-10)
  expect_equal(count_oc(6, 100, rate_c = 0.02, rate_d = c(0.02, 0.04, 0.08),
                        k = 2, allocation = "random"),
               c(0.015483640641779385, 0.11768375216471988,
                 0.5593072758258119), tolerance = 1e-10)
  expect_equal(count_oc(215, 20000, rate_c = 0.007,
                        rate_d = 0.007 * c(1, 2, 4), likelihood = "poisson"),
               c(2.4954832886901208e-09, 0.3741442355892674,
                 0.9999999999999969), tolerance = 1e-8)
})

test_that("count_oc gives 0 for a missing critical value", {
  expect_identical(count_oc(NA, 100, 0.02, 0.04), 0)
  expect_error(count_oc(6, 100, 0.02, 1.5), "'rate_d'")
})
