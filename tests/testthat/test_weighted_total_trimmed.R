context("weighted.total.trimmed")
library(robsurvey)

test_that("Function outputs correct result for example", {
  x <- c(0.1, 0.35, 0.05, 0.1, 0.15, 0.05, 0.2)
  expect_equal(weighted.total.trimmed(x, x, LB = 0.2), 0.1361111, tolerance=1e-05)
})

test_that("Weighted mean is equal to mean if weights are all 1", {
  x <- c(0.1, 0.35, 0.05, 0.1, 0.15, 0.05, 0.2)
  w <- c(1, 1, 1, 1, 1, 1, 1)
  expect_equal(weighted.total.trimmed(x, w, LB = 0.2), mean(x, trim = 0.2) * sum(w), tolerance=1e-05)
})
