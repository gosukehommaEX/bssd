test_that("sim_oc returns summaries for each analysis and model", {
  design <- list(n = 60, k = 1, accrual = 0, follow_up = 52,
                 analysis_time = c(26, 52), par_p = -log(0.9) / 52,
                 par_d = -log(0.7) / 52)
  pr <- list(bdribs = list(lambda_hp = -log(0.9) / 52, pt_hp = 2000))
  oc <- sim_oc(10, design, models = c("BDRIBS1", "BDRIBS2", "EXP_BERN"),
               priors = pr, t_star = 52, seed = 1, keep = TRUE)
  s <- oc$summary
  expect_identical(nrow(s), 6L)
  expect_identical(s$model, rep(c("BDRIBS1", "BDRIBS2", "EXP_BERN"), 2))
  expect_identical(s$look, rep(1:2, each = 3))
  expect_true(all(s$n_ok == 10L & s$n_failed == 0L))
  expect_true(all(s$prob_signal >= 0 & s$prob_signal <= 1))
  expect_equal(s$mean_n, rep(60, 6))
  expect_equal(oc$truth, c(placebo = 0.1, drug = 0.3, difference = 0.2))
  r <- oc$results
  expect_identical(nrow(r), 60L)
  # with the same events, the smaller exposure of BDRIBS2 gives a larger
  # posterior probability of r > 1
  expect_true(all(r$prob[r$model == "BDRIBS2"] >=
                    r$prob[r$model == "BDRIBS1"] - 1e-12))
  oc2 <- sim_oc(10, design, models = c("BDRIBS1", "BDRIBS2", "EXP_BERN"),
                priors = pr, t_star = 52, seed = 1)
  expect_identical(oc2$summary, oc$summary)
  expect_null(oc2$results)
})

test_that("sim_oc runs LL_BERN and checks its input", {
  design <- list(n = 40, k = 1, accrual = 0, follow_up = 104,
                 dist = "loglogistic", par_p = c(6, 0.8),
                 par_d = c(5.3, 0.8))
  pr <- list(llbern = list(mu_hp = 6, sigma_hp = 0.8, pt_hp = 200))
  oc <- sim_oc(2, design, models = "LL_BERN", priors = pr, t_star = 104,
               seed = 2)
  expect_identical(oc$summary$n_ok, 2L)
  expect_equal(unname(oc$truth["placebo"]), ll_cdf(104, 6, 0.8))
  expect_error(sim_oc(2, design, models = "LL_BERN", priors = pr), "t_star")
  expect_error(sim_oc(2, design, models = "BDRIBS1", priors = pr), "bdribs")
  expect_error(sim_oc(2, design, models = "XYZ", priors = pr), "subset")
  bad <- list(bdribs = list(lambda_hp = 0.001, pt_hp = 100, wrong = 1))
  expect_warning(sim_oc(2, design, models = "BDRIBS1", priors = bad,
                        seed = 1), "fits failed")
})
