test_that("merge_control merges settings and rejects unknown names", {
  m <- merge_control(list(b = 3), list(a = 1, b = 2))
  expect_identical(m, list(a = 1, b = 3))
  expect_error(merge_control(list(z = 1), list(a = 1)), "Unknown")
  expect_error(merge_control(1, list(a = 1)), "must be a list")
})
