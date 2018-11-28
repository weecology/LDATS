devtools::load_all()
data(rodents)
rem <- which(colnames(rodents) %in% c("newmoon", "date", "plots", "traps"))
lda_data <- rodents[,-rem]
dct <- data.frame(newmoon = rodents[,"newmoon"])
r_lda <- LDA_set(lda_data, topics = 3, nseeds = 2)

LDA_models = select_LDA(r_lda)
document_covariate_table = dct
timename = "newmoon"
formulas = c(~ 1)
nchangepoints = 2
weights = NULL
control = TS_controls_list()


gamma = LDA_models[[1]]@gamma
data <- data.frame(document_covariate_table)
data$gamma <- gamma
formula<- gamma ~ newmoon
changepoints <- NULL


#
# working within the TS function
#   everything is set and working great within est_changepts
#   now need to include the unconditional estimates of the regressors
#   also want to include tests on the all existing functions and probably
#   want to break up the TS function script a bit soon...


t1 <- system.time(
  rho_dist <- est_changepts(data, formula, nchangepoints, weights, control)
)

