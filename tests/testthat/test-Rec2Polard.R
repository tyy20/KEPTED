test_that("Rec2Polard works", {
  expect_equal(Rec2Polar(c(-sqrt(2),1,1)), list(R=2,Theta=c(-pi/4,pi/4)))
})
