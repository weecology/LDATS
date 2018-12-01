# currently bug checking the pipeline
#  in particular dealing with whatever weirdness is happening to do with
#  the normalization...i think some of it might have to do with needing
#  to standardize within chunks? idk, we might just need to put the 
#  standardizing on the "to do later" issues list
#
# to do:
#   plot functions
#   tests (add new functions to the list to do!)
#   tidy up documentation
#   full vignette
#
# also add better catching in the pipeline, so that everything doesn't fail 
# when one thing fails

devtools::load_all()
data(rodents)
document_term_table <- rodents$document_term_table
document_covariate_table <- rodents$document_covariate_table

topics = 2
nseeds = 1
formulas = c(~ 1)
nchangepoints = 1
weights = document_weights(document_term_table)
control = LDA_TS_controls_list()
LDAs <- LDA_set(document_term_table, topics, nseeds, control$LDA_control)
LDA_models <- select_LDA(LDAs, control$LDA_control)

control = TS_controls_list()

# into TS_on_LDA
#
#
# 2. work through TS ugggh ok lets go




mods <- LDA_TS(document_term_table, document_covariate_table,
               topics = 2, nseeds = 1, formulas = c(~ 1), nchangepoints = 1, 
               weights = document_weights(document_term_table), 
               control = LDA_TS_controls_list())
