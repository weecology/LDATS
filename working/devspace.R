devtools::load_all()
data(rodents)
rem <- which(colnames(rodents) %in% c("newmoon", "date", "plots", "traps"))
lda_data <- rodents[,-rem]
dct <- data.frame(newmoon = rodents[,"newmoon"])
r_lda <- LDA_set(lda_data, topics = 2, nseeds = 2)

LDA_models = select_LDA(r_lda)
document_covariate_table = dct
timename = "newmoon"
formula = ~ 1
changepoints = 0
weights = NULL
control = TS_controls_list()


i <- 1
gamma = LDA_models[[mods$LDA[i]]]@gamma



