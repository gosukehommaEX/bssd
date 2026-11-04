test_that("check_subject_data converts and validates the data", {
  d <- check_subject_data(c(TRUE, FALSE), c(1, 2))
  expect_identical(d$y, c(1L, 0L))
  expect_identical(d$t, c(1, 2))
  expect_error(check_subject_data(c(0, 2), c(1, 2)), "out of range")
  expect_error(check_subject_data(c(0, 1), c(1, 0)), "out of range")
  expect_error(check_subject_data(c(0, 1), 1), "length 2")
})
