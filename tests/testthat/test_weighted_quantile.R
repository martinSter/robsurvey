context("weighted_quantile")
library(robsurvey)

test_that("Function outputs correct result for example", {
  x <- c(0.1, 0.35, 0.05, 0.1, 0.15, 0.05, 0.2)
  expect_equal(weighted_quantile(x, x, probs = c(0.25, 0.5, 0.75)), c(0.05, 0.15, 0.20), tolerance=1e-05)
})

test_that("Weighted quantile is equal to quantile if weights are all 1", {
  x <- c(0.1, 0.35, 0.05, 0.1, 0.15, 0.05, 0.2)
  w <- rep(1, length(x))
  probs <- c(0.25, 0.5, 0.75)
  expect_equal(weighted_quantile(x, w, probs), as.vector(quantile(x, probs, type = 1)), tolerance=1e-05)
})
