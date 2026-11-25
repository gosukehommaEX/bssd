test_that("row_nodes maps the Chebyshev points to the interval of each row", {
  v <- row_nodes(c(0, -2), c(1, 2), 5L)
  expect_identical(dim(v), c(2L, 5L))
  expect_equal(v[1, ], (cheb_nodes(5L) + 1) / 2)
  expect_equal(v[2, ], 2 * cheb_nodes(5L))
  expect_equal(v[, 1], c(0, -2))
  expect_equal(v[, 5], c(1, 2))
})
