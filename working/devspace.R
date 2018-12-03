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

topics = 4
nseeds = 1
formulas = c(~ sin_year + cos_year)
nchangepoints = 4
weights = document_weights(document_term_table)
control = LDA_TS_controls_list(TS_controls_list(nit = 1e4))
LDAs <- LDA_set(document_term_table, topics, nseeds, control$LDA_control)
LDA_models <- select_LDA(LDAs, control$LDA_control)

control = TS_controls_list()

x <- TS_on_LDA(LDA_models, document_covariate_table, formulas, nchangepoints,
               weights, control)


# into TS_on_LDA
#
#
# 2. work through TS ugggh ok lets go

save(x, file="ok2.RData")


mods <- LDA_TS(document_term_table, document_covariate_table,
               topics = 2, nseeds = 1, formulas = c(~ 1), nchangepoints = 1, 
               weights = document_weights(document_term_table), 
               control = LDA_TS_controls_list())

x<-TSmods[[1]]
class(x)
names(x)
cols = list(rho_cols = NULL, gamma_cols = NULL,
                                    rho_option = "D", gamma_option = "B", 
                                    rho_alpha = 0.4, gamma_alpha = 0.8)
TS_summary_plot(x, bin_width = 5, xlab = "New moon", "median", cols)

dim(x$rhos)
head(x$rhos)

range(data[ , control$timename])

