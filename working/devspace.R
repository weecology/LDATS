devtools::load_all()
data(rodents)
rem <- which(colnames(rodents) %in% c("newmoon", "date", "plots", "traps"))
lda_data <- rodents[,-rem]

r_lda <- LDA_set(MV = lda_data, topics = 2, nseeds = 2)
smod <- select_LDA(r_lda)

slotNames(x1)

yy <- smod[[1]]@gamma
XX(yy~1)