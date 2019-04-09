context("mer")
library(robsurvey)

test_that("Function outputs correct result for example", {

  data(api)
  dstrat <- svydesign(id=~1, strata=~stype, weights=~pw, data=apistrat, fpc=~fpc)
  invisible(capture.output(m1 <- msvymean(~api00, dstrat, type="rht", k=1.3)))
  invisible(capture.output(m1.mer <- mer(m1)))

  expect_equal(m1.mer[[1]], 771.8075)
})
