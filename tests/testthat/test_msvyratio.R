context("msvyratio")
library(robsurvey)

test_that("Function outputs correct result for example", {

  data(api)
  dstrat <- svydesign(id=~1, strata=~stype, weights=~pw, data=apistrat, fpc=~fpc)
  invisible(capture.output(ratio1 <- msvyratio(~api00, ~api99, dstrat, k=1.2, na.rm=TRUE)))

  expect_equal(coef(ratio1)[[1]], 1.04975, tolerance=1e-05)
  expect_equal(vcov(ratio1)[[1]], 9.593178e-06, tolerance=1e-05)
})
