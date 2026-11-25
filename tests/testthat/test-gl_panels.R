test_that("gl_panels without split points is the Gauss-Legendre rule", {
  g <- gl_panels(-1, 1, 20L)
  expect_equal(g$x, gl_nodes(20L)$x)
  expect_equal(g$w, gl_nodes(20L)$w)
})

test_that("gl_panels integrates a function with a kink at the split point", {
  # integral of |x - 0.3| over (-1, 1) is (1.3^2 + 0.7^2) / 2 = 1.09
  g1 <- gl_panels(-1, 1, 20L, 0.3)
  expect_equal(sum(g1$w * abs(g1$x - 0.3)), 1.09, tolerance = 1e-12)
  # without the split, the error is not negligible
  g0 <- gl_panels(-1, 1, 20L)
  expect_gt(abs(sum(g0$w * abs(g0$x - 0.3)) - 1.09), 1e-6)
})

test_that("gl_panels gives each panel at least half of the nodes", {
  # panels (0, 9.5) and (9.5, 10): 19 and 10 nodes; points outside the
  # interval are ignored
  g <- gl_panels(0, 10, 20L, c(9.5, -1, 20))
  expect_length(g$x, 29L)
  expect_true(all(diff(g$x) > 0))
  expect_equal(sum(g$w), 10, tolerance = 1e-12)
})
