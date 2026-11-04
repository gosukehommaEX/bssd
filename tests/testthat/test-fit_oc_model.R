test_that("fit_oc_model passes the data and priors to the posterior functions", {
  set.seed(4)
  d <- sim_blinded(n = 50, accrual = 0, follow_up = 52, par_p = 0.004,
                   par_d = 0.008)
  pr <- list(bdribs = list(lambda_hp = 0.004, pt_hp = 1000),
             expbern = list(lambda_hp = 0.004, pt_hp = 1000))
  f1 <- fit_oc_model("BDRIBS1", d, 1, pr, 1, 0.95, 52, list())
  expect_equal(f1$prob, bdribs_post(sum(d$event), sum(d$follow_up),
                                    lambda_hp = 0.004, pt_hp = 1000)$prob)
  f2 <- fit_oc_model("BDRIBS2", d, 1, pr, 1, 0.95, 52, list())
  expect_equal(f2$prob, bdribs_post(sum(d$event), sum(d$time_at_risk),
                                    lambda_hp = 0.004, pt_hp = 1000)$prob)
  f3 <- fit_oc_model("EXP_BERN", d, 1, pr, 1, 0.95, 52, list())
  expect_equal(f3$prob, expbern_post(d$event, d$follow_up, lambda_hp = 0.004,
                                     pt_hp = 1000)$prob)
  expect_true(is.na(f3$rr[["median"]]))
})
