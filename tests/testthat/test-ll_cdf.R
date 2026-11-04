test_that("ll_cdf reproduces incidences of Table 1 of Duncan (2026)", {
  # reference values computed with scipy.special.expit
  expect_equal(ll_cdf(104, 7.3, 0.8), 0.03490786619283047, tolerance = 1e-12)
  expect_equal(ll_cdf(104, 5.07, 0.7), 0.3525125923317835, tolerance = 1e-12)
  expect_equal(ll_cdf(52, 6.17, 1.6), 0.19993147844301123, tolerance = 1e-12)
  expect_equal(ll_cdf(104, 6.845, 0.8), 0.06004366491216468, tolerance = 1e-12)
  # scenario 4: 5.0% is the week 52 value, not the week 104 value
  expect_equal(ll_cdf(c(52, 104), 9.85, 2),
               c(0.049765910595029036, 0.06895814220916327), tolerance = 1e-12)
})

test_that("ll_cdf is a distribution function and checks its input", {
  expect_equal(ll_cdf(exp(6), 6, 0.5), 0.5)
  expect_true(all(diff(ll_cdf(1:200, 5, 0.7)) > 0))
  expect_error(ll_cdf(0, 6, 1), "out of range")
  expect_error(ll_cdf(1, 6, -1), "out of range")
})
