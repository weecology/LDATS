context("Check LDA_TS functions")

data(rodents)
lda_data <- rodents$document_term_table
document_term_table <- rodents$document_term_table
document_covariate_table <- rodents$document_covariate_table

test_that("LDA_TS on 0 changepoints", {
  mod0 <- LDA_TS(document_term_table, document_covariate_table,
                 topics = 2, nseeds = 1, formulas = ~ 1, nchangepoints = 0,
                 weights = document_weights(document_term_table), 
                 control = LDA_TS_controls_list())
  expect_is(mod0, "LDA_TS")
  expect_equal(length(names(mod0)), 4)
  expect_is(mod0[[4]], "TS_fit")
})

test_that("LDA_TS on 1 changepoints", {
  mod1 <- LDA_TS(document_term_table, document_covariate_table,
                 topics = 2, nseeds = 1, formulas = ~ 1, nchangepoints = 1,
                 weights = document_weights(document_term_table), 
                 control = LDA_TS_controls_list(
                           TS_control = TS_controls_list(nit = 100)))
  expect_is(mod1, "LDA_TS")
  expect_equal(length(names(mod1)), 4)
  expect_is(mod1[[4]], "TS_fit")
})

test_that("check print on LDA_TS", {
  mod1 <- LDA_TS(document_term_table, document_covariate_table,
                 topics = 2, nseeds = 1, formulas = ~ 1, nchangepoints = 1,
                 weights = document_weights(document_term_table), 
                 control = LDA_TS_controls_list(
                           TS_control = TS_controls_list(nit = 100)))
  expect_output(print(mod1))
})

test_that("Check LDA_TS_controls_list", {
  expect_is(LDA_TS_controls_list(), "LDA_TS_controls")
  expect_equal(length(LDA_TS_controls_list()), 3)
})

test_that("Check package_LDA_TS", {
  topics <- 2
  nseeds <- 1
  formulas <- ~ 1
  nchangepoints <- 1
  weights <- document_weights(document_term_table)
  control <- LDA_TS_controls_list(TS_control = TS_controls_list(nit = 100))
  LDAs <- LDA_set(document_term_table, topics, nseeds, control$LDA_control)
  sel_LDA <- select_LDA(LDAs, control$LDA_control)
  TSs <- TS_on_LDA(sel_LDA, document_covariate_table, formulas, nchangepoints,
                   weights, control$TS_control)
  sel_TSs <- select_TS(TSs, control$TS_control)
  
  expect_is(package_LDA_TS(LDAs, sel_LDA, TSs, sel_TSs), "LDA_TS")
  expect_equal(length(package_LDA_TS(LDAs, sel_LDA, TSs, sel_TSs)), 4)
  expect_error(package_LDA_TS())
  expect_error(package_LDA_TS("ok"))
})