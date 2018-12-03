# to do:
#   tests 
#   tidy up documentation
#   full vignette
#
#
#  for tests, i'm just working through script-by-script
#  still to do:
#    TS
#    TS_on_LDA
#    TS_plots


devtools::load_all()
data(rodents)
document_term_table <- rodents$document_term_table
document_covariate_table <- rodents$document_covariate_table


xx <- LDA_TS(document_term_table, document_covariate_table,
                   topics = 2, nseeds = 1, formulas = ~ 1, nchangepoints = 1,
                   weights = document_weights(document_term_table), 
                   control = LDA_TS_controls_list())

xx <- LDA_TS(document_term_table, document_covariate_table,
                   topics = 4, nseeds = 10, formulas = ~ sin_year + cos_year, 
                   nchangepoints = 4,
                   weights = document_weights(document_term_table), 
                   control = LDA_TS_controls_list(
                             TS_control = TS_controls_list(nit = 1e3)))