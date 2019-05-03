context("weighted_median_line")
library(robsurvey)

test_that("Function outputs correct result for example", {
  x <- c(1, 2, 4, 5)
  y <- c(3, 2, 7, 4)
  expect_equal(coef(weighted_median_line(y~x, type="prod")), c(3,0), tolerance=1e-05)
})
