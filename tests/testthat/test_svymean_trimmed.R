context("svymean_trimmed")
library(robsurvey)

test_that("Function outputs correct result for example", {

  suppressWarnings(library(survey))
  data(api)
  dstrat <- svydesign(id=~1, strata=~stype, weights=~pw, data=apistrat, fpc=~fpc)
  invisible(capture.output(out <- svymean_trimmed(~api00, dstrat, LB = 0.05)))

  expect_equal(coef(out)[[1]], 655.362, tolerance=1e-05)
  expect_equal(vcov(out)[[1]], 133.8212, tolerance=1e-05)
})

test_that("Results of barebone function and estimation methods are identical", {

  suppressWarnings(library(survey))
  data(api)
  dstrat <- svydesign(id=~1, strata=~stype, weights=~pw, data=apistrat, fpc=~fpc)
  w <- weights(dstrat)
  invisible(capture.output(out <- svymean_trimmed(~api00, dstrat, LB = 0.05)))

  expect_equal(coef(out)[[1]], weighted_mean_trimmed(apistrat$api00, w, LB = 0.05), tolerance=1e-05)
})
