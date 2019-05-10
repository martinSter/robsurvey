context("weighted_median_ratio")
library(robsurvey)

test_that("Function outputs correct result for example", {
  x <- c(1,2,4,5)
  y <- c(1,0,5,2)
  expect_equal(coef(weighted_median_ratio(y~x, w=rep(1, length(x)))), 0.7, tolerance=1e-05)
  expect_equal(fitted(weighted_median_ratio(y~x, w=rep(1, length(x)))), c(0.7, 1.4, 2.8, 3.5), tolerance=1e-05)
})
