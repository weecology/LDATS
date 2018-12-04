context("Check ptMCMC functions")

test_that("check prep_temp_sequence", {
  expect_equal(length(prep_temp_sequence()), 6)
  expect_equal(length(prep_temp_sequence(TS_controls_list(ntemps = 9))), 9)
  expect_equal(round(prep_temp_sequence()[3], 2), 8)
  expect_equal(round(prep_temp_sequence(TS_controls_list(q = 1))[3], 1), 2.8)
  expect_error(prep_temp_sequence(LDA_controls_list()))
})
