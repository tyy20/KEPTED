test_that("PolarDerivatived works", {
  expect_equal(PolarDerivative(c(1/sqrt(2),1/2,1/2)),
               matrix(c(1/sqrt(2),-1/2,-1/2,
                        0,1,-1),
                      nrow=2,ncol=3,byrow=TRUE))
})
