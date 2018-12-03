# to do:
#   tests 
#   tidy up documentation
#   full vignette
#


devtools::load_all()
data(rodents)
document_term_table <- rodents$document_term_table
document_covariate_table <- rodents$document_covariate_table


xx <- LDA_TS(document_term_table, document_covariate_table,
                   topics = 2, nseeds = 1, formulas = ~ 1, nchangepoints = 1,
                   weights = document_weights(document_term_table), 
                   control = LDA_TS_controls_list())