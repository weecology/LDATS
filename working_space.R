# an example based on the Portal rodent data

"%>%" <- magrittr::"%>%"
data(rodents)
lda_data <- rodents %>%
            dplyr::select(-c(newmoonnumber, newmoondate, nplots, ntraps))
r_LDA <- LDATS::LDA(data = lda_data, ntopics = 2:5, nseeds = 2, ncores = 4)
ldamodel <- r_LDA[[1]]
plot(r_LDA)
plot(ldamodel)

ts_data <- data.frame(rodents[ , "newmoonnumber"])
colnames(ts_data) <- "time"
ts_data$gamma <- ldamodel@gamma
preds <- "time"
sample_sizes <- apply(lda_data, 1, sum)
wts <- round(sample_sizes/max(sample_sizes), 3)

LDATS::MTS(formula = preds, data = ts_data, nchangepoints = 1, weights = wts, 
  nit = 100)
