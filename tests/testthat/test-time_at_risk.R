test_that("time_at_risk computes both exposure measures", {
  out <- time_at_risk(follow_up = c(52, 52, 30, 12), event = c(1, 0, 1, 0),
                      event_time = c(10, NA, 25, NA))
  expect_equal(out, c(exposure = 146, time_at_risk = 99, n_event = 2))
})

test_that("time_at_risk checks its input", {
  expect_error(time_at_risk(c(10, 20), c(1, 0), c(15, NA)), "not larger")
  expect_error(time_at_risk(c(10, 20), c(1, 0), c(NA, NA)), "positive")
  expect_error(time_at_risk(c(10, 20), c(1, 0), 5), "same length")
})
