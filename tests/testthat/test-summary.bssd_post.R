test_that("summary.bssd_post returns a one-row data frame", {
  fit <- bdribs_post(y = 5, exposure = 500, k = 2, lambda_hp = 18 / 4000,
                     pt_hp = 4000)
  s <- summary(fit)
  expect_s3_class(s, "data.frame")
  expect_identical(nrow(s), 1L)
  expect_identical(s$model, "BDRIBS")
  expect_equal(s$prob, fit$prob)
  expect_equal(s$rr_median, fit$rr[["median"]])
  expect_true(is.na(s$ip_placebo))
  fit2 <- bdribs_post(y = 5, exposure = 500, k = 2, lambda_hp = 18 / 4000,
                      pt_hp = 4000, t_star = 1)
  expect_equal(summary(fit2)$ip_drug, fit2$incidence[["drug"]])
})
