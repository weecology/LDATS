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

test_that("check logLik for LDA_VEM", {
  expect_is(logLik(lda[[1]]), "logLik")
  expect_equal(round(as.numeric(logLik(lda[[1]]))), -47889)
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
  expect_silent(check_topics(5))
  expect_silent(check_topics(2:5))
})

test_that("check error catching of check_seeds", {
  expect_error(check_seeds("a"))
  expect_error(check_seeds(1.5))
  expect_error(check_seeds(2), NA)
  expect_error(check_seeds(c(2, 3, 4)), NA)
  expect_silent(check_seeds(5))
  expect_silent(check_seeds(1:5))
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

test_that("check check_LDA_set_inputs", {
  expect_silent(check_LDA_set_inputs(lda_data, 2, 1, LDA_controls_list()))
  expect_error(check_LDA_set_inputs(lda_data, 2, "ok", 2))
  expect_error(check_LDA_set_inputs(lda_data, "ok", 2, LDA_controls_list()))
  expect_error(check_LDA_set_inputs("ok", 2, 1, LDA_controls_list()))
})

test_that("check package_LDA_set", {
  document_term_table <- lda_data
  topics <- 2
  nseeds <- 1 
  control <- LDA_controls_list()
  check_LDA_set_inputs(document_term_table, topics, nseeds, control)
  mod_topics <- rep(topics, each = length(seq(2, nseeds * 2, 2)))
  mod_seeds <- rep(seq(2, nseeds * 2, 2), length(topics))
  nmods <- length(mod_topics)
  mods <- vector("list", length = nmods)
  for (i in 1:nmods){
    LDA_msg(mod_topics[i], mod_seeds[i], control)
    control_i <- prep_LDA_control(seed = mod_seeds[i], control = control)
    mods[[i]] <- LDA(document_term_table, k = mod_topics[i], 
                     control = control_i)
  }
  expect_is(package_LDA_set(mods, mod_topics, mod_seeds), "LDA_set")
  expect_error(package_LDA_set(mods, 0.2, mod_seeds))
  expect_error(package_LDA_set(mods, mod_topics, 0.2))
  expect_error(package_LDA_set("ok", mod_topics, mod_seeds))
})

test_that("check LDA_msg", {
  expect_output(LDA_msg(2, 1, LDA_controls_list()))
  expect_error(LDA_msg(2, 0.5, LDA_controls_list()))
})

test_that("Check LDA_controls_list", {
  expect_is(LDA_controls_list(), "LDA_controls")
  expect_equal(length(LDA_controls_list()), 3)
})
