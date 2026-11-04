test_that("cheb_nodes returns increasing Lobatto points on [-1, 1]", {
  x <- cheb_nodes(5)
  expect_equal(x, c(-1, -sqrt(2) / 2, 0, sqrt(2) / 2, 1), tolerance = 1e-14)
  expect_true(all(diff(cheb_nodes(41)) > 0))
})
