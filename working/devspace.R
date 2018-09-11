rem <- which(colnames(rodents) %in% c("newmoon", "date", "plots", "traps"))
lda_data <- rodents[,-rem]

x <- LDA_set(lda_data, c(2,5))
x1 <- x[[1]]

slotNames(x1)


test LDA functions.
fill out any more function doc info. 