test_that("kerGauss works", {
  expect_equal(kerGauss(1,c(1,2),c(2,1)),exp(-2))
})
