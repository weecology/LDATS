context("Check load of data against historic data load")

test_that("data(rodents) produces same data set", {
  hash_val <- digest::digest(data(rodents))
  expect_equal(hash_val, "d1d3201c55cb74b44109513daf07bee4")
})