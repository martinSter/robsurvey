context("msvytotal")
library(robsurvey)

test_that("Function outputs correct result for example", {

  data(api)
  dstrat <- svydesign(id=~1, strata=~stype, weights=~pw, data=apistrat, fpc=~fpc)
  invisible(capture.output(rht1 <- msvytotal(~api00, dstrat, k=1.2)))

  expect_equal(coef(rht1)[[1]], 4050931, tolerance=1e-05)
  expect_equal(vcov(rht1)[[1]], 3525487144, tolerance=1e-05)
})
