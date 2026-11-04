test_that("cheb_coef_matrix recovers Chebyshev coefficients", {
  n <- 9
  x <- cheb_nodes(n)
  # f = 2 T_0 - T_3 + 0.5 T_8
  f <- 2 - cos(3 * acos(x)) + 0.5 * cos(8 * acos(x))
  a <- drop(cheb_coef_matrix(n) %*% f)
  expect_equal(a, c(2, 0, 0, -1, 0, 0, 0, 0, 0.5), tolerance = 1e-12)
})
