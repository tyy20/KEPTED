test_that("Polar2Recd works", {
  expect_equal(Polar2Rec(4,c(pi/4,pi/4)),
               list(X=c(2*sqrt(2),2,2),V=c(2*sqrt(2),2,2)/4))
})
