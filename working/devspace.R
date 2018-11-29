# TS is working now just great
# next to do: plot functions, tests, tidy up documentation
# after that: select_TS
# and then: LDA_TS
# and then: full vignette

devtools::load_all()
data(rodents)
document_term_table <- rodents$document_term_table
document_covariate_table <- rodents$document_covariate_table

mods <- LDA_TS(document_term_table, document_covariate_table,
               topics = 2:3, nseeds = 2, formulas = c(~ 1, ~newmoon), 
               nchangepoints = 0:1,
               weights = NULL, LDA_control = NULL, 
               TS_control = TS_controls_list())
