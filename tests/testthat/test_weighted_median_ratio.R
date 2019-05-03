context("weighted_median_ratio")
library(robsurvey)

test_that("Function outputs correct result for example", {
  x <- c(1,2,4,5)
  y <- c(1,0,5,2)
  expect_equal(coef(weighted_median_ratio(y~x)), 0.4, tolerance=1e-05)
  expect_equal(fitted(weighted_median_ratio(y~x)), c(0.4, 0.8, 1.6, 2.0), tolerance=1e-05)
})
