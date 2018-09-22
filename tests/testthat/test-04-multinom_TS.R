context("Check multinomial TS functions")

data(rodents)
rem <- which(colnames(rodents) %in% c("newmoon", "date", "plots", "traps"))
lda_data <- rodents[ , -rem]
lda <- LDA_set(lda_data, c(4), nseeds = 1)
dct <- data.frame(time = rodents[ , "newmoon"])
mts_data <- data.frame(dct)
mts_data$gamma <- lda[[1]]@gamma

test_that("check good output from multinom_TS", {
  mts <- multinom_TS(data = mts_data, formula_RHS = "1", 
           changepoints = c(20,50), weights = NULL, 
           control = TS_controls_list())
  expect_is(mts, "list")
  expect_equal(length(mts), 2)
  expect_equal(names(mts), c("chunk models", "logLik"))
  expect_equal(length(mts$"chunk models"), 3)
  expect_is(mts$logLik, "numeric")
})

test_that("check failed output from multinom_TS", {
  mts <- multinom_TS(data = mts_data, formula_RHS = "1", 
           changepoints = c(50,40), weights = NULL, 
           control = TS_controls_list())
  expect_is(mts, "list")
  expect_equal(length(mts), 2)
  expect_equal(names(mts), c("chunk models", "logLik"))
  expect_equal(mts$"chunk models", NA)
  expect_equal(mts$logLik, -Inf)
})

test_that("check bad changepoint catching of check_chunks", {
  expect_equal(check_chunks(mts_data, -1), FALSE)
  expect_equal(check_chunks(mts_data, 1e5), FALSE)
  expect_equal(check_chunks(mts_data, NULL), TRUE)
  expect_equal(check_chunks(mts_data, 100), TRUE)
  expect_equal(check_chunks(mts_data, c(10, 50, 100)), TRUE)
})

test_that("check memoization of multinom_TS_chunk", {
  expect_is(memoise_fun(multinom_TS_chunk, TRUE), "memoised")
})

test_that("check multinom_TS_chunk", {
  expect_is(multinom_TS_chunk(mts_data, "gamma ~ 1", 0, 40), "multinom")
})

test_that("check memoised multinom_TS_chunk", {
  multinom_TS_chunk_memo <- memoise_fun(multinom_TS_chunk, TRUE)
  expect_is(multinom_TS_chunk_memo(mts_data, "gamma ~ 1", 0, 40), "multinom")
})