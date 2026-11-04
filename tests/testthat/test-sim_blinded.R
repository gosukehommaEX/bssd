test_that("sim_blinded returns the expected structure", {
  set.seed(1)
  d <- sim_blinded(n = 30, k = 2, accrual = 10, follow_up = 52,
                   analysis_time = c(20, 62), par_p = 0.01, par_d = 0.02)
  expect_s3_class(d, "data.frame")
  expect_named(d, c("look", "analysis_time", "id", "group", "entry",
                    "follow_up", "event", "event_time", "time_at_risk"))
  expect_identical(sort(unique(d$look)), 1:2)
  l1 <- d[d$look == 1, ]
  l2 <- d[d$look == 2, ]
  expect_identical(nrow(l2), 30L)
  expect_identical(sum(l2$group == "drug"), 20L)
  expect_true(all(l1$entry < 20))
  expect_equal(l1$follow_up, 20 - l1$entry)
  expect_equal(l2$follow_up, rep(52, 30))
  ev <- d$event == 1
  expect_true(all(d$event_time[ev] <= d$follow_up[ev]))
  expect_true(all(is.na(d$event_time[!ev])))
  expect_equal(d$time_at_risk, ifelse(ev, d$event_time, d$follow_up))
  # an event at the first analysis is also an event at the second analysis
  expect_true(all(l2$event[match(l1$id[l1$event == 1], l2$id)] == 1))
})

test_that("sim_blinded generates the specified incidences", {
  set.seed(3)
  d <- sim_blinded(n = 20000, k = 1, accrual = 0, follow_up = 52,
                   par_p = -log(0.9) / 52, par_d = -log(0.7) / 52)
  ip <- tapply(d$event, d$group, mean)
  # Monte Carlo standard errors are about 0.003 and 0.005
  expect_equal(unname(ip["placebo"]), 0.1, tolerance = 0.1)
  expect_equal(unname(ip["drug"]), 0.3, tolerance = 0.05)
  dl <- sim_blinded(n = 20000, k = 1, accrual = 0, follow_up = 104,
                    dist = "loglogistic", par_p = c(7.3, 0.8),
                    par_d = c(5.85, 0.8))
  ip2 <- tapply(dl$event, dl$group, mean)
  expect_equal(unname(ip2["placebo"]), ll_cdf(104, 7.3, 0.8), tolerance = 0.2)
  expect_equal(unname(ip2["drug"]), ll_cdf(104, 5.85, 0.8), tolerance = 0.08)
})

test_that("sim_blinded handles dropout and checks its input", {
  set.seed(5)
  d <- sim_blinded(n = 2000, accrual = 0, follow_up = 52, par_p = 0.001,
                   par_d = 0.001, dropout = 0.01)
  # follow-up is min(52, time to dropout), whose mean is
  # (1 - exp(-0.01 * 52)) / 0.01, about 40.5
  expect_equal(mean(d$follow_up), (1 - exp(-0.01 * 52)) / 0.01,
               tolerance = 0.03)
  set.seed(2)
  a <- sim_blinded(n = 20, accrual = 5, follow_up = 10, par_p = 0.1,
                   par_d = 0.1)
  set.seed(2)
  b <- sim_blinded(n = 20, accrual = 5, follow_up = 10, par_p = 0.1,
                   par_d = 0.1)
  expect_identical(a, b)
  expect_error(sim_blinded(n = 10, accrual = 0, follow_up = 52,
                           dist = "loglogistic", par_p = 6,
                           par_d = c(6, 1)), "length 2")
  expect_error(sim_blinded(n = 10, accrual = 0, follow_up = 52,
                           dist = "loglogistic", par_p = c(6, -1),
                           par_d = c(6, 1)), "positive")
  expect_error(sim_blinded(n = 10, accrual = 0, follow_up = 52,
                           par_p = 0, par_d = 0.1), "out of range")
})
