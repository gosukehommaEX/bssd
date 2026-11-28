test_that("fast_exp_cpp agrees with exp() and expm1()", {
  set.seed(15)
  x <- c(stats::runif(2e4, -745, 709), stats::runif(2e4, -1, 1),
         10^stats::runif(2e3, -300, 0), -10^stats::runif(2e3, -300, 0),
         -708 + c(-1e-9, 0, 1e-9), 709 + c(-1e-9, 0, 1e-9), 0, NaN, Inf,
         -Inf)
  for (m1 in c(FALSE, TRUE)) {
    got <- fast_exp_cpp(x, m1)
    ref <- if (m1) expm1(x) else exp(x)
    ok <- is.finite(ref) & abs(ref) >= .Machine$double.xmin
    # about four units in the last place, including the error of ref
    expect_lt(max(abs(got[ok] / ref[ok] - 1)), 1e-15)
    # NaN, infinities, zeros and results below the smallest normal number
    expect_identical(is.na(got[!ok]), is.na(ref[!ok]))
    inf <- !ok & is.infinite(ref)
    expect_identical(got[inf], ref[inf])
    small <- !ok & is.finite(ref)
    expect_lt(max(abs(got[small] - ref[small])), 1e-300)
  }
})
