# an example based on the Portal rodent data

"%>%" <- magrittr::"%>%"
data(rodents)
lda_data <- rodents %>%
            dplyr::select(-c(newmoonnumber, newmoondate, nplots, ntraps))
ts_data <- data.frame(rodents[ , "newmoonnumber"])
colnames(ts_data) <- "time"

r_LDA <- LDATS::LDA(data = lda_data, ntopics = 2:5, nseeds = 2, ncores = 4)
lda_mods <- 
LDATS::LDA_TS(lda_data, ts_data, ntopics = 2:5, nseeds = 2, 
                          ncores = 4, nit = 100)



wts <- doc_weights(lda_data)

ldamodel <- r_LDA[[1]]
plot(r_LDA)
plot(ldamodel)

LDATS::MTS(formula = "1", data = lda_mods[[1]], nchangepoints = 1, weights = wts, 
  nit = 100)


LDATS::MTS_set(prepped_data = lda_mods, formula = "1", nchangepoints = 1, weights = wts, 
  nit = 100)