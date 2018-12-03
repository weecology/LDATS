context("Check TS functions")

data(rodents)
lda_data <- rodents$document_term_table
document_term_table <- rodents$document_term_table
document_covariate_table <- rodents$document_covariate_table
topics <- 2
nseeds <- 1
formulas <- ~ 1
nchangepoints <- 0
weights <- document_weights(document_term_table)
control <- LDA_TS_controls_list()
LDAs <- LDA_set(document_term_table, topics, nseeds, control$LDA_control)
LDA_models <- select_LDA(LDAs, control$LDA_control)
control <- control$TS_control
mods <- expand_TS(LDA_models, formulas, nchangepoints)
formula <- mods$formula[[1]]
nchangepoints <- mods$nchangepoints[1]
data <- prep_TS_data(document_covariate_table, LDA_models, mods, 1)
TSmod <- TS(data, formula, nchangepoints, weights, control)


test_that("Check check_formula", {
  expect_silent(check_formula(data, formula))
  expect_error(check_formula(data, ~1))
  expect_error(check_formula(document_covariate_table, formula))
})

test_that("Check TS_controls_list", {
  expect_is(TS_controls_list(), "TS_controls")
  expect_equal(length(TS_controls_list()), 17)
})