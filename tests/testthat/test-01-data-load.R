context("Check data loading")

test_that("check rodents data", {
  hash_val <- digest(data(rodents))
  expect_equal(hash_val, "d1d3201c55cb74b44109513daf07bee4")
})