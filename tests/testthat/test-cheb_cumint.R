test_that("cheb_cumint matches reference integrals", {
  # reference values computed with scipy.integrate.quad
  f <- function(x) exp(-3 * (x - 0.2)^2) * (1 + 0.3 * sin(4 * x))
  n <- 33
  anti <- cheb_antideriv(cheb_coef_matrix(n) %*% f(cheb_nodes(n)))
  ref <- c(0, 0.03645541789030094, 0.36053011138890756, 0.9911433572179201,
           1.0602334784388836)
  xq <- c(-1, -0.5, 0.1, 0.73, 1)
  got <- vapply(xq, function(v) cheb_cumint(anti, v), numeric(1))
  expect_equal(got, ref, tolerance = 1e-12)
})

test_that("cheb_cumint handles several columns and clamps the limits", {
  n <- 17
  x <- cheb_nodes(n)
  vals <- cbind(rep(1, n), x + 1)
  anti <- cheb_antideriv(cheb_coef_matrix(n) %*% vals)
  got <- cheb_cumint(anti, c(0, 2))
  # integral of 1 from -1 to 0 is 1; integral of x + 1 from -1 to 1 is 2
  expect_equal(got, c(1, 2), tolerance = 1e-12)
})
