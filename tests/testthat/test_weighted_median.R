context("weighted_median")
library(robsurvey)

test_that("Function outputs correct result for example", {
  x <- c(0.1, 0.35, 0.05, 0.1, 0.15, 0.05, 0.2)
  expect_equal(weighted_median(x, x), 0.2, tolerance=1e-05)
})

test_that("Function results in error if no weights are provided", {
  x <- c(0.1, 0.35, 0.05, 0.1, 0.15, 0.05, 0.2)
  expect_error(weighted_median(x))
})

test_that("Weighted median is equal to median if weights are all 1", {
  x <- c(0.1, 0.35, 0.05, 0.1, 0.15, 0.05, 0.2)
  w <- rep(1, length(x))
  expect_equal(weighted_median(x, w), median(x), tolerance=1e-05)
})

test_that("Weighted median and 50% weighted quantile are equal", {
  x <- c(0.1, 0.35, 0.05, 0.1, 0.15, 0.05, 0.2)
  expect_equal(weighted_median(x, x), weighted_quantile(x, x, 0.5), tolerance=1e-05)
})
