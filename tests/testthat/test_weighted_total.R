context("weighted_total")
library(robsurvey)

test_that("Function outputs correct result for example", {
  x <- c(0.1, 0.35, 0.05, 0.1, 0.15, 0.05, 0.2)
  expect_equal(weighted_total(x, x), sum(x*x), tolerance=1e-05)
})

test_that("Weighted total is equal to sum if weights are all 1", {
  x <- c(0.1, 0.35, 0.05, 0.1, 0.15, 0.05, 0.2)
  w <- c(1, 1, 1, 1, 1, 1, 1)
  expect_equal(weighted_total(x, w), sum(x), tolerance=1e-05)
})
