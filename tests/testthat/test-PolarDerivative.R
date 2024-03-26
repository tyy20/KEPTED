test_that("PolarDerivative works", {
  expect_equal(PolarDerivative(c(1/sqrt(2),1/sqrt(2))),
               matrix(c(1/sqrt(2),-1/sqrt(2)),nrow=1,ncol=2))
})
