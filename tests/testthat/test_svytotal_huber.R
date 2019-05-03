context("svytotal_huber")
library(robsurvey)

test_that("Function outputs correct result for example", {

  suppressWarnings(library(survey))
  data(api)
  dstrat <- svydesign(id=~1, strata=~stype, weights=~pw, data=apistrat, fpc=~fpc)
  invisible(capture.output(out <- svytotal_huber(~api00, dstrat, k = 2)))

  expect_equal(coef(out)[[1]], 4106045, tolerance=1e-05)
  expect_equal(vcov(out)[[1]], 3056968590, tolerance=1e-05)
})
