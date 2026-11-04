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
  got <- llbern_loglik_cpp(t, as.integer(y), wp, mp, sp, md, sd)
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
