test_that("bdribs_grid evaluates bdribs_post on all combinations", {
  g <- bdribs_grid(y = 6:7, exposure = c(80, 140), lambda_hp = 0.031,
                   p_prior = c(0.5, 0.5))
  expect_identical(nrow(g), 4L)
  expect_identical(names(g), c("y", "exposure", "lambda_hp", "prob"))
  p <- bdribs_post(y = 7, exposure = 140, k = 1, lambda_hp = 0.031,
                   pt_hp = Inf, p_prior = c(0.5, 0.5), details = FALSE)$prob
  expect_equal(g$prob[g$y == 7 & g$exposure == 140], p, tolerance = 1e-12)
  # more exposure with the same events lowers the probability
  expect_true(g$prob[g$y == 6 & g$exposure == 140] <
                g$prob[g$y == 6 & g$exposure == 80])
  expect_error(bdribs_grid(6, 80, 0.031, details = TRUE), "details")
})
