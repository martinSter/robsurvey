context("wsvymean")
library(robsurvey)

test_that("Function outputs correct result for example", {

  suppressWarnings(library(survey))
  data(api)
  dstrat <- svydesign(id=~1, strata=~stype, weights=~pw, data=apistrat, fpc=~fpc)
  invisible(capture.output(out <- wsvymean(~api00, dstrat, LB = 0.05)))

  expect_equal(coef(out)[[1]], 640.5986, tolerance=1e-05)
  expect_equal(vcov(out)[[1]], 133.8212, tolerance=1e-05)
})

test_that("Results of barebone function and estimation methods are identical", {

  suppressWarnings(library(survey))
  data(api)
  dstrat <- svydesign(id=~1, strata=~stype, weights=~pw, data=apistrat, fpc=~fpc)
  w <- weights(dstrat)
  invisible(capture.output(out <- wsvymean(~api00, dstrat, LB = 0.05)))

  expect_equal(coef(out)[[1]], weighted.mean.winsorized(apistrat$api00, w, LB = 0.05), tolerance=1e-05)
})
