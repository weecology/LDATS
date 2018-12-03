# to do:
#   bug-catching through the pipeline at TS
#   tests 
#   tidy up documentation
#   full vignette
#

devtools::load_all()
data(rodents)
document_term_table <- rodents$document_term_table
document_covariate_table <- rodents$document_covariate_table


topics = 4
nseeds = 1
formulas = c(~ 1)
nchangepoints = 4
weights = document_weights(document_term_table)
control = LDA_TS_controls_list(TS_controls_list(nit = 1e4))
LDAs <- LDA_set(document_term_table, topics, nseeds, control$LDA_control)
LDA_models <- select_LDA(LDAs, control$LDA_control)
control <- control$TS_control

  mods <- expand_TS(LDA_models, formulas, nchangepoints)
  nmods <- nrow(mods)
  TSmods <- vector("list", nmods)
  i <- 1
  print_model_run_message(mods, i, LDA_models, control)
  formula_i <- mods$formula[[i]]
  nchangepoints_i <- mods$nchangepoints[i]
  data_i <- prep_TS_data(document_covariate_table, LDA_models, mods, i)
  TSmods[[i]] <- TS(data_i, formula_i, nchangepoints_i, weights, control)


