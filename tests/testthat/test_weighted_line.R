context("weighted.line")
library(robsurvey)

test_that("Function outputs correct result for example", {
  data(cars)
  out1 <- weighted.line(cars$speed, cars$dist)
  out2 <- weighted.line(cars$speed, cars$dist, w=rep(1:10, each=5))
  expect_equal(coef(out1), c(-29.333333, 4.666667), tolerance=1e-05)
  expect_equal(coef(out2), c(-26.666667, 4.222222), tolerance=1e-05)
})

test_that("weighted.line() results are identical to line() results", {
  data(cars)
  out1 <- weighted.line(cars$speed, cars$dist)
  out2 <- line(cars$speed, cars$dist)
  expect_equal(coef(out1), coef(out2), tolerance=1e-05)
})
