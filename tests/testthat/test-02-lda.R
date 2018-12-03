context("Check LDA functions")

data(rodents)
lda_data <- rodents$document_term_table
lda <- LDA_set(lda_data, c(2, 4), nseeds = 2, LDA_controls_list(quiet = TRUE))

test_that("check output from LDA_set", {
  expect_equal(length(lda), 4)
  expect_is(lda, "LDA_set")
  expect_is(lda[[1]], "LDA")
  expect_is(lda[[2]], "LDA")
  expect_is(lda[[3]], "LDA")
  expect_is(lda[[4]], "LDA")
})

test_that("check error catching of check_document_term_table", {
  dtt <- "a"
  expect_error(check_document_term_table(dtt))
  dtt <- matrix(1:100, 10, 10)
  expect_error(check_document_term_table(dtt, NA))
  dtt <- data.frame("dummy" = 1:100)
  expect_error(check_document_term_table(dtt, NA))
})

test_that("check error catching of check_topics", {
  expect_error(check_topics("a"))
  expect_error(check_topics(1.5))
  expect_error(check_topics(1))
  expect_error(check_topics(2), NA)
  expect_error(check_topics(c(2, 3, 4)), NA)
})

test_that("check output from prep_LDA_control", {
  expect_is(prep_LDA_control(1), "list")
  expect_equal(prep_LDA_control(1)$seed, 1)
  expect_equal(prep_LDA_control(1, LDA_controls_list(seed = 10))$seed, 1)
})

test_that("check selection via select_LDA", {
  expect_is(select_LDA(lda), "LDA_set")
  expect_equal(length(select_LDA(lda)), 1)
  expect_equal(select_LDA(lda)[1], lda[3])
})