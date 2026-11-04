test_that("cheb_antideriv gives the antiderivative of a polynomial", {
  n <- 6
  x <- cheb_nodes(n)
  # f(x) = 3 x^2, antiderivative from -1 is x^3 + 1
  anti <- cheb_antideriv(cheb_coef_matrix(n) %*% (3 * x^2))
  expect_equal(nrow(anti$B), n + 1L)
  xq <- 0.3
  val <- sum(anti$B[, 1] * cos((0:n) * acos(xq))) - anti$base
  expect_equal(val, xq^3 + 1, tolerance = 1e-12)
})
