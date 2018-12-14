context("Check multinomial TS functions")

data(rodents)
lda_data <- rodents$document_term_table
lda <- LDA_set(lda_data, c(4), nseeds = 1, LDA_controls_list(quiet = TRUE))
dct <- rodents$document_covariate_table
mts_data <- data.frame(dct)
mts_data$gamma <- lda[[1]]@gamma

test_that("check packaging of chunk fits", {
  control <- TS_controls_list()
  TS_chunk_memo <- memoise_fun(multinom_TS_chunk, control$memoise)
  chunks <- prep_chunks(data = mts_data, changepoints = c(20,50), 
                        timename = control$timename)
  nchunks <- nrow(chunks)
  fits <- vector("list", length = nchunks)
  for (i in 1:nchunks){
    fits[[i]] <- TS_chunk_memo(data = mts_data, formula = gamma ~ 1, 
                               chunk = chunks[i, ], weights = NULL)
  }
  packaged <- package_chunk_fits(chunks, fits)
  expect_is(packaged, "multinom_TS_fit")
  expect_equal(round(packaged$logLik, 2), -516.58)
})

test_that("check logLik for multinom_TS_fit", {
  mts <- multinom_TS(data = mts_data, formula = gamma~1, 
           changepoints = c(20,50), weights = NULL, 
           control = TS_controls_list())
  expect_is(logLik(mts), "logLik")
  expect_equal(round(as.numeric(logLik(mts))), -517)
})

test_that("check good output from multinom_TS", {
  mts <- multinom_TS(data = mts_data, formula = gamma~1, 
           changepoints = c(20,50), weights = NULL, 
           control = TS_controls_list())
  expect_is(mts, "list")
  expect_is(mts, "multinom_TS_fit")
  expect_equal(length(mts), 3)
  expect_equal(names(mts), c("chunk models", "logLik", "chunks"))
  expect_equal(length(mts$"chunk models"), 3)
  expect_is(mts$logLik, "numeric")
})

test_that("check check_changepoints", {
  expect_silent(check_changepoints())
  expect_silent(check_changepoints(1))
  expect_error(check_changepoints("ok"))
  expect_error(check_changepoints(0.3))
})

test_that("check failed output from multinom_TS", {
  mts <- multinom_TS(data = mts_data, formula = gamma~1, 
           changepoints = c(50,40), weights = NULL, 
           control = TS_controls_list())
  expect_is(mts, "list")
  expect_equal(length(mts), 3)
  expect_equal(names(mts), c("chunk models", "logLik", "chunks"))
  expect_equal(mts$"chunk models", NA)
  expect_equal(mts$logLik, -Inf)
})

test_that("check bad change point location catching of 
          verify_changepoint_locations", {
  expect_equal(verify_changepoint_locations(mts_data, -1), FALSE)
  expect_equal(verify_changepoint_locations(mts_data, 1e5), FALSE)
  expect_equal(verify_changepoint_locations(mts_data, NULL), TRUE)
  expect_equal(verify_changepoint_locations(mts_data, 100), TRUE)
  expect_equal(verify_changepoint_locations(mts_data, c(10, 50, 100)), TRUE)
})

test_that("check memoization of multinom_TS_chunk", {
  expect_is(memoise_fun(multinom_TS_chunk, TRUE), "memoised")
})

chunk <- data.frame(start = 0, end = 40)
test_that("check multinom_TS_chunk", {
  expect_is(multinom_TS_chunk(mts_data, "gamma ~ 1", chunk), "multinom")
})

test_that("check memoised multinom_TS_chunk", {
  multinom_TS_chunk_memo <- memoise_fun(multinom_TS_chunk, TRUE)
  expect_is(multinom_TS_chunk_memo(mts_data, "gamma ~ 1", chunk), "multinom")
})

test_that("check prepping of chunks", {
  expect_is(prep_chunks(mts_data, NULL, "newmoon"), "data.frame")
  expect_equal(prep_chunks(mts_data, NULL, "newmoon")$start, 1)
  expect_equal(prep_chunks(mts_data, NULL, "newmoon")$end, 467)
  expect_equal(prep_chunks(mts_data, c(10), "newmoon")$start, c(1, 11))
  expect_equal(prep_chunks(mts_data, c(10), "newmoon")$end, c(10, 467))
})
