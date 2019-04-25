context("weighted.total.huber")
library(robsurvey)

test_that("Function outputs correct result for example", {
  x <- c(0.1, 0.35, 0.05, 0.1, 0.15, 0.05, 0.2)
  expect_equal(weighted.total.huber(x, x, k = 1.34), 0.1881523, tolerance=1e-05)
})

