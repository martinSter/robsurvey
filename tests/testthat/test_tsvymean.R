context("tsvymean")
library(robsurvey)

test_that("Function outputs correct result for example", {

  data(api)
  dstrat <- svydesign(id=~1, strata=~stype, weights=~pw, data=apistrat, fpc=~fpc)
  invisible(capture.output(tm1 <- tsvymean(~api00, dstrat, trim=c(0.01, 0.09), type="trim")))

  expect_equal(coef(tm1)[[1]], 641.5356, tolerance=1e-05)
  expect_equal(vcov(tm1)[[1]], 117.8851, tolerance=1e-05)
})
