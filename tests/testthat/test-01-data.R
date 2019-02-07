context("Check data loading")

test_that("check rodents data", {
  data(rodents)
  expect_equal(length(rodents), 2)
  expect_equal(names(rodents), c("document_term_table",
                                 "document_covariate_table"))
})

test_that("check jornada data", {
  data(jornada)
  expect_equal(length(rodents), 2)
  expect_equal(names(rodents), c("document_term_table",
                                 "document_covariate_table"))
})