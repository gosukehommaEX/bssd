test_that("print.bssd_post prints the main results and returns its input", {
  fit <- bdribs_post(y = 5, exposure = 500, k = 2, lambda_hp = 18 / 4000,
                     pt_hp = 4000, t_star = 1)
  out <- capture.output(res <- print(fit))
  expect_identical(res, fit)
  expect_true(any(grepl("Model: BDRIBS", out)))
  expect_true(any(grepl("P(r > 1 | data) = 0.8199", out, fixed = TRUE)))
  expect_true(any(grepl("Signal (posterior probability > 0.95): no", out,
                        fixed = TRUE)))
  expect_true(any(grepl("Cumulative incidence at t_star = 1", out,
                        fixed = TRUE)))
})
