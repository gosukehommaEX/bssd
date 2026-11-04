test_that("check_num accepts valid input and rejects invalid input", {
  expect_silent(check_num(1, "x"))
  expect_silent(check_num(c(0.5, 2), "x", len = 2L, lower = 0,
                          lower_open = TRUE))
  expect_error(check_num("a", "x"), "finite numeric")
  expect_error(check_num(NA_real_, "x"), "finite numeric")
  expect_error(check_num(c(1, 2), "x"), "length 1")
  expect_error(check_num(0, "x", lower = 0, lower_open = TRUE), "out of range")
  expect_error(check_num(1, "x", upper = 1, upper_open = TRUE), "out of range")
  expect_error(check_num(1.5, "x", integer = TRUE), "integer valued")
})
