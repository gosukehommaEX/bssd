test_that("row_tail_share is small for smooth rows and large for rough rows", {
  x <- cheb_nodes(33)
  L <- rbind(rep(0, 33), -x^2, log(abs(x) + 1e-3), rep(-Inf, 33))
  s <- row_tail_share(L)
  # constant and Gaussian densities: only rounding errors in the highest
  # coefficients (about 1e-14 in Python)
  expect_lt(s[1], 1e-12)
  expect_lt(s[2], 1e-12)
  # a kink at 0: the coefficients decay slowly (share 0.0116 in Python)
  expect_gt(s[3], 1e-3)
  expect_lt(s[3], 0.05)
  # a row whose density is zero everywhere
  expect_identical(s[4], Inf)
})

test_that("row_tail_share treats NA as a zero density", {
  x <- cheb_nodes(33)
  L1 <- rbind(-x^2, -x^2)
  L2 <- L1
  L1[2L, 6L] <- NA
  L2[2L, 6L] <- -Inf
  expect_identical(row_tail_share(L1), row_tail_share(L2))
  expect_gt(row_tail_share(L1)[2L], 0.01)
})
