test_that("expbern_loglik_cpp matches a direct R computation", {
  set.seed(11)
  n <- 50
  t <- runif(n, 0.5, 104)
  y <- rbinom(n, 1, 0.3)
  lam <- exp(c(-11, -6, -3, -1))
  r <- exp(c(-4, -1, 0, 2, 4))
  wp <- 1 / 3
  got <- expbern_loglik_cpp(t, as.integer(y), wp, lam, r)
  ref <- outer(lam, r, Vectorize(function(l, rr) {
    h1 <- l * t
    h2 <- rr * l * t
    m1 <- wp * (-expm1(-h1)) + (1 - wp) * (-expm1(-h2))
    a <- log(wp) - h1
    b <- log(1 - wp) - h2
    l0 <- pmax(a, b) + log1p(exp(-abs(a - b)))
    sum(ifelse(y == 1, log(m1), l0))
  }))
  expect_equal(dim(got), c(4L, 5L))
  expect_equal(got, ref, tolerance = 1e-12)
})
