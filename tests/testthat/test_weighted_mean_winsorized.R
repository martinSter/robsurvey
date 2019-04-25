context("weighted.mean.winsorized")
library(robsurvey)

test_that("Function outputs correct result for example", {
  x <- c(0.1, 0.35, 0.05, 0.1, 0.15, 0.05, 0.2)
  expect_equal(weighted.mean.winsorized(x, x, LB = 0.2), 0.09, tolerance=1e-05)
})
