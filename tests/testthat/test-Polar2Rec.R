test_that("Polar2Rec works", {
  expect_equal(Polar2Rec(2,pi/4), list(X=c(sqrt(2),sqrt(2)),V=c(sqrt(2),sqrt(2))/2))
})
