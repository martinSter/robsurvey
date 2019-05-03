context("svytotal_winsorized")
library(robsurvey)

test_that("Function outputs correct result for example", {

  suppressWarnings(library(survey))
  data(api)
  dstrat <- svydesign(id=~1, strata=~stype, weights=~pw, data=apistrat, fpc=~fpc)
  invisible(capture.output(out <- svytotal_winsorized(~api00, dstrat, LB = 0.05)))

  expect_equal(coef(out)[[1]], 3967868, tolerance=1e-05)
  expect_equal(vcov(out)[[1]], 5134136032, tolerance=1e-05)
})
