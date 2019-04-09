context("msvymean")
library(robsurvey)

test_that("Function outputs correct result for example", {

  data(api)
  dstrat <- svydesign(id=~1, strata=~stype, weights=~pw, data=apistrat, fpc=~fpc)
  invisible(capture.output(rht1 <- msvymean(~api00, dstrat, type="rht", k=1.2)))

  expect_equal(coef(rht1)[[1]], 666.8978, tolerance=1e-05)
  expect_equal(vcov(rht1)[[1]], 88.28954, tolerance=1e-05)
})
