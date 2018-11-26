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
nchangepoints = 0
weights = NULL
control = TS_controls_list()


gamma = LDA_models[[1]]@gamma
data <- data.frame(document_covariate_table)
data$gamma <- gamma
formula<- gamma ~ newmoon
changepoints <- NULL

mts <- multinom_TS(data, formula, changepoints = c(5, 100))

# multinom_TS and its set of underlying functions are all good to go at this 
# point. so now let's jump back to TS and get going there!
#  presently working within TS, have all of the prep work done, leading up to
#  the main for loop

xx <- multinom_TS(data, formula, changepoints = c(5, 100))[[1]][[1]]
vcov(xx)


