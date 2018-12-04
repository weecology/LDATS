context("Check ptMCMC functions")


data(rodents)
lda_data <- rodents$document_term_table
document_term_table <- rodents$document_term_table
document_covariate_table <- rodents$document_covariate_table
topics <- 2
nseeds <- 1
formulas <- ~ 1
nchangepoints <- 1
weights <- document_weights(document_term_table)
control <- LDA_TS_controls_list()
LDAs <- LDA_set(document_term_table, topics, nseeds, control$LDA_control)
LDA_models <- select_LDA(LDAs, control$LDA_control)
control <- TS_controls_list(nit = 1e3, seed = 1)
mods <- expand_TS(LDA_models, formulas, nchangepoints)
formula <- mods$formula[[1]]
nchangepoints <- mods$nchangepoints[1]
data <- prep_TS_data(document_covariate_table, LDA_models, mods, 1)
TSmod <- TS(data, formula, nchangepoints, weights, control)

rho_dist <- est_changepts(data, formula, nchangepoints, weights, control)
eta_dist <- est_regressors(rho_dist, data, formula, weights, control)
set.seed(1)



test_that("check prep_temp_sequence", {
  expect_equal(length(prep_temp_sequence()), 6)
  expect_equal(length(prep_temp_sequence(TS_controls_list(ntemps = 9))), 9)
  expect_equal(round(prep_temp_sequence()[3], 2), 8)
  expect_equal(round(prep_temp_sequence(TS_controls_list(q = 1))[3], 1), 2.8)
  expect_error(prep_temp_sequence(LDA_controls_list()))
})
