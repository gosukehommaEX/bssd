test_that("cc_weights integrate polynomials and smooth functions", {
  n <- 11
  x <- cheb_nodes(n)
  w <- cc_weights(n)
  expect_equal(sum(w), 2, tolerance = 1e-13)
  expect_equal(sum(w * x^4), 2 / 5, tolerance = 1e-13)
  expect_equal(sum(w * x^9), 0, tolerance = 1e-13)
  w21 <- cc_weights(21)
  expect_equal(sum(w21 * exp(cheb_nodes(21))), exp(1) - exp(-1),
               tolerance = 1e-13)
})
