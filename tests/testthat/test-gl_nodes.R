test_that("gl_nodes reproduces the three-point rule", {
  g <- gl_nodes(3)
  expect_equal(g$x, c(-sqrt(3 / 5), 0, sqrt(3 / 5)), tolerance = 1e-13)
  expect_equal(g$w, c(5, 8, 5) / 9, tolerance = 1e-13)
})

test_that("gl_nodes integrates polynomials of degree 2n - 1 exactly", {
  g <- gl_nodes(10)
  expect_true(all(diff(g$x) > 0))
  expect_equal(sum(g$w), 2, tolerance = 1e-13)
  expect_equal(sum(g$w * g$x^18), 2 / 19, tolerance = 1e-12)
  expect_equal(sum(g$w * (g$x^19 + g$x^7)), 0, tolerance = 1e-12)
})
