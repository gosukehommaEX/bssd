test_that("refine_box finds the region of a two-dimensional Gaussian", {
  fun <- function(nodes) {
    outer(-0.5 * ((nodes[[1]] - 3) / 0.1)^2, -0.5 * ((nodes[[2]] + 2) / 2)^2,
          "+")
  }
  for (start in list(c(-50, 50), c(10, 20))) {
    b <- refine_box(fun, c(start[1], -50), c(start[2], 50), 41, 1e-8, 12,
                    c(-Inf, -Inf), c(Inf, Inf))
    # the box covers at least 5.5 and at most 15 standard deviations each side
    expect_true(b$lo[1] < 3 - 5.5 * 0.1 && b$lo[1] > 3 - 15 * 0.1)
    expect_true(b$hi[1] > 3 + 5.5 * 0.1 && b$hi[1] < 3 + 15 * 0.1)
    expect_true(b$lo[2] < -2 - 5.5 * 2 && b$hi[2] > -2 + 5.5 * 2)
  }
})

test_that("refine_box respects the limits", {
  b <- refine_box(function(nodes) array(-5 * nodes[[1]], 41), 0, 1, 41, 1e-8,
                  12, 0, 1)
  expect_equal(b$lo, 0)
  expect_equal(b$hi, 1)
})
