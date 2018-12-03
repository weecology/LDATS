context("Check utilities")

test_that("check modalvalue", {
  xx <- c(1, 2, 3, 4, 5, 4, 3, 2, 1, 2)
  expect_equal(modalvalue(xx), 2)
  expect_error(modalvalue("ok"))
})

test_that("check document_weights", {
  data(rodents)
  lda_data <- rodents$document_term_table
  doc_weights <- document_weights(lda_data)
  expect_equal(round(mean(doc_weights), 3), 0.282)
  expect_equal(length(doc_weights), 436)
  expect_error(document_weights("ok"))
})

test_that("check qprint", {
  expect_error(qprint())
  expect_error(qprint("ok"))
  expect_error(qprint("ok", ""))
  expect_output(qprint("ok", "", quiet = FALSE))
  expect_silent(qprint("ok", "", quiet = TRUE))
})

test_that("check mirror_vcov", {
  y <- 1:10
  x <- 101:110
  mod <- lm(y ~ x)
  vcv <- mirror_vcov(mod)  
  expect_equal(isSymmetric(vcv), TRUE)
  expect_error(mirror_vcov("ok"))
})

test_that("check normalize", {
  xx <- c(1, 2, 3, 4, 5, 4, 3, 2, 1, 2)
  expect_equal(mean(normalize(xx)), 0.425)
  expect_equal(normalize(xx)[1], 0)
  xx <- -1000:100
  expect_equal(mean(normalize(xx)), 0.5)
  expect_equal(round(sd(normalize(xx)), 3), 0.289)
  expect_error(normalize("ok"))
})

test_that("check memoise_fun", {
  expect_is(memoise_fun(sum, TRUE), "memoised")
  expect_is(memoise_fun(sum, FALSE), "function")  
  expect_error(memoise_fun(1, TRUE))
  expect_error(memoise_fun(sum, 1))
})

test_that("check check_control", {
  expect_silent(check_control(TS_controls_list(), "TS_controls"))
  expect_silent(check_control(LDA_TS_controls_list(), "LDA_TS_controls"))
  expect_silent(check_control(LDA_controls_list(), "LDA_controls"))
  expect_error(check_control(TS_controls_list(), "LDA_controls"))
  expect_error(check_control(LDA_controls_list(), "TS_controls"))
  expect_error(check_control(LDA_TS_controls_list(), "LDA_controls"))
})
