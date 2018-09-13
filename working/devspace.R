devtools::load_all()
data(rodents)
rem <- which(colnames(rodents) %in% c("newmoon", "date", "plots", "traps"))
lda_data <- rodents[,-rem]
dct <- data.frame(newmoon = rodents[,"newmoon"])

r_lda <- LDA_set(MV = lda_data, topics = 2, nseeds = 2)
smod <- select_LDA(r_lda)

slotNames(x1)

yy <- smod[[1]]@gamma
XX(yy~1)


TS_set_on_LDA(r_lda, dct, "newmoon")