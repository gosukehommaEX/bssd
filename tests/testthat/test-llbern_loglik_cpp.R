test_that("llbern_loglik_cpp matches a direct R computation", {
  set.seed(12)
  n <- 50
  t <- runif(n, 0.5, 104)
  y <- rbinom(n, 1, 0.3)
  mp <- c(3, 6, 9)
  sp <- c(0.05, 0.8, 2.5)
  md <- c(2, 5, 7, 10)
  sd <- c(0.005, 0.3, 1, 4)
  wp <- 0.4
  got <- llbern_loglik_cpp(t, y, rep(1, n), wp, mp, sp, md, sd)
  lse <- function(a, b) pmax(a, b) + log1p(exp(-abs(a - b)))
  ref <- matrix(NA_real_, 3, 4)
  for (i in 1:3) {
    for (j in 1:4) {
      zp <- (log(t) - mp[i]) / sp[i]
      zd <- (log(t) - md[j]) / sd[j]
      l1 <- lse(log(wp) + plogis(zp, log.p = TRUE),
                log(1 - wp) + plogis(zd, log.p = TRUE))
      l0 <- lse(log(wp) + plogis(-zp, log.p = TRUE),
                log(1 - wp) + plogis(-zd, log.p = TRUE))
      ref[i, j] <- sum(ifelse(y == 1, l1, l0))
    }
  }
  expect_equal(got, ref, tolerance = 1e-12)
})

test_that("llbern_loglik_cpp stays accurate when the products are small", {
  # all subjects with the event and small incidences: the partial products
  # of the C++ code fall below 1e-70 and are moved to the sum of logarithms
  t <- seq(1, 50, length.out = 403)
  y <- rep(1, length(t))
  mp <- c(14, 15)
  sp <- c(0.8, 1)
  md <- c(13, 16, 20)
  sd <- c(1, 0.5, 2)
  wp <- 0.5
  got <- llbern_loglik_cpp(t, y, rep(1, length(t)), wp, mp, sp, md, sd)
  lse <- function(a, b) pmax(a, b) + log1p(exp(-abs(a - b)))
  ref <- matrix(NA_real_, 2, 3)
  for (i in 1:2) {
    for (j in 1:3) {
      zp <- (log(t) - mp[i]) / sp[i]
      zd <- (log(t) - md[j]) / sd[j]
      ref[i, j] <- sum(lse(log(wp) + plogis(zp, log.p = TRUE),
                           log(1 - wp) + plogis(zd, log.p = TRUE)))
    }
  }
  expect_lt(max(ref), -2000)
  expect_equal(got, ref, tolerance = 1e-12)
})

test_that("llbern_loglik_cpp gives the same values for grouped subjects", {
  set.seed(14)
  t <- c(rep(c(26, 52), c(30, 60)), runif(13, 1, 52))
  y <- rbinom(length(t), 1, 0.25)
  pat <- group_patterns(y, t)
  mp <- c(4, 6)
  sp <- c(0.5, 1.5)
  md <- c(3, 5, 8)
  sd <- c(0.2, 1, 3)
  expect_lt(length(pat$t), length(t) / 2)
  expect_equal(llbern_loglik_cpp(pat$t, pat$y, pat$w, 0.5, mp, sp, md, sd),
               llbern_loglik_cpp(t, y, rep(1, length(t)), 0.5, mp, sp, md,
                                 sd),
               tolerance = 1e-13)
})
