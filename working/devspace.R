devtools::load_all()
data(rodents)
rem <- which(colnames(rodents) %in% c("newmoon", "date", "plots", "traps"))
lda_data <- rodents[,-rem]
dct <- data.frame(newmoon = rodents[,"newmoon"])
r_lda <- LDA_set(lda_data, topics = 3, nseeds = 2)

LDA_models <- select_LDA(r_lda)
data <- data.frame(dct)
data$gamma <- LDA_models[[1]]@gamma



mod0 <- TS(data = data, formula = gamma ~1, nchangepoints= 0, weights = NULL, 
               control = TS_controls_list())
mod1 <- TS(data = data, formula = gamma ~1, nchangepoints= 1, weights = NULL, 
               control = TS_controls_list())
mod2 <- TS(data = data, formula = gamma ~1, nchangepoints= 2, weights = NULL, 
               control = TS_controls_list())

