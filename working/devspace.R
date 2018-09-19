devtools::load_all()
data(rodents)
rem <- which(colnames(rodents) %in% c("newmoon", "date", "plots", "traps"))
lda_data <- rodents[,-rem]
dct <- data.frame(newmoon = rodents[,"newmoon"])

r_lda <- LDA_set(lda_data, topics = 2, nseeds = 2)
smod <- select_LDA(r_lda)

slotNames(x1)

yy <- smod[[1]]@gamma
XX(yy~1)


TS_set_on_LDA(r_lda, dct, "newmoon")

LDA_models <- smod
document_covariate_table <- dct
timename <- "newmoon"
formula = ~ 1
changepoints = 0
weights = NULL
ptMCMC_controls = ptMCMC_controls_list()

TS(gamma_i, formula_i, nchangepoints_i, ptMCMC_controls)






