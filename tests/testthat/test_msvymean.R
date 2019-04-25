context("msvymean")
library(robsurvey)

test_that("Function outputs correct result for example", {

  suppressWarnings(library(survey))
  data(api)
  dstrat <- svydesign(id=~1, strata=~stype, weights=~pw, data=apistrat, fpc=~fpc)
  invisible(capture.output(out <- msvymean(~api00, dstrat, k = 2)))

  expect_equal(coef(out)[[1]], 662.9068, tolerance=1e-05)
  expect_equal(vcov(out)[[1]], 79.67986, tolerance=1e-05)
})

test_that("Results of barebone function and estimation methods are identical", {

  suppressWarnings(library(survey))
  data(api)
  dstrat <- svydesign(id=~1, strata=~stype, weights=~pw, data=apistrat, fpc=~fpc)
  w <- weights(dstrat)
  invisible(capture.output(out <- msvymean(~api00, dstrat, k = 2)))

  expect_equal(coef(out)[[1]], weighted.mean.huber(apistrat$api00, w, k = 2), tolerance=1e-05)
})
