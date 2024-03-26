test_that("Rec2Polar works", {
  expect_equal(Rec2Polar(c(-sqrt(2),sqrt(2))), list(R=2,Theta=-pi/4))
})
