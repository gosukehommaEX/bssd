test_that("pooled_ratio gives the ratio of Waterhouse et al. (2022)", {
  # Example 3: 250 patient-years at 1:1 and 80 patient-years at 2:1
  expect_equal(pooled_ratio(c(250, 80), c(1, 2)), 1.1758241758241759,
               tolerance = 1e-12)
  expect_equal(round(pooled_ratio(c(250, 80), c(1, 2)), 2), 1.18)
  # one study gives its own ratio
  expect_equal(pooled_ratio(100, 3), 3)
  expect_error(pooled_ratio(c(1, 2), 1), "'k'")
})
