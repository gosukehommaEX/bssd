test_that("multibin_inner equals the Beta function ratio at delta = 1", {
  # with delta = 1, the integral is B(a + y, b + n - y) / B(a, b)
  for (cs in list(c(47, 1600, 2, 98), c(6, 200, 2, 98), c(0, 50, 1, 1),
                  c(1, 10, 0.5, 0.5))) {
    v <- multibin_inner(0, cs[1], cs[2], 0.5, cs[3], cs[4], 12L)
    expect_equal(v, lbeta(cs[3] + cs[1], cs[4] + cs[2] - cs[1]) -
                   lbeta(cs[3], cs[4]), tolerance = 1e-9)
  }
})

test_that("multibin_inner matches quadrature with integrate", {
  for (w in c(-3, 0.5, 2)) {
    cc <- 1 - 0.5 + 0.5 * exp(w)
    f <- function(th) {
      exp((2 + 6 - 1) * log(th) + (98 - 1) * log1p(-th) +
            (200 - 6) * log1p(-pmin(cc * th, 1)) - lbeta(2, 98) + 30)
    }
    up <- min(1, exp(-w))
    ref <- log(integrate(f, 0, 0.1, rel.tol = 1e-12,
                         subdivisions = 1000L)$value +
                 integrate(f, 0.1, up, rel.tol = 1e-12,
                           subdivisions = 1000L)$value) - 30
    expect_equal(multibin_inner(w, 6, 200, 0.5, 2, 98, 12L), ref,
                 tolerance = 1e-7)
  }
})
