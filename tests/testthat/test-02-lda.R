context("Check LDA functions")

data(rodents)
rem <- which(colnames(rodents) %in% c("newmoon", "date", "plots", "traps"))
lda_data <- rodents[ , -rem]
lda <- LDA_set(lda_data, c(2, 4), nseeds = 2)

test_that("check output from LDA_set", {
  expect_equal(length(lda), 4)
  expect_is(lda, "LDA_list")
  expect_is(lda[[1]], "LDA")
  expect_is(lda[[2]], "LDA")
  expect_is(lda[[3]], "LDA")
  expect_is(lda[[4]], "LDA")
})

test_that("check error catching of check_MV", {
  expect_error(check_MV("a"))
  expect_error(check_MV(matrix(1:100, 10, 10), NA))
  expect_error(check_MV(data.frame("dummy" = 1:100), NA))
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
  expect_equal(prep_LDA_control(1, list(seed = 10))$seed, 1)
})

test_that("check selection via select_LDA", {
  expect_is(select_LDA(lda), "LDA_list")
  expect_equal(length(select_LDA(lda)), 1)
  expect_equal(select_LDA(lda)[1], lda[3])
})