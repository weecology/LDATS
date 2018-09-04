data(rodents)
lda_data <- dplyr::select(rodents, -c(newmoon, date, plots, traps))
ts_data <- data.frame("time" = data.frame(rodents)[ , "newmoon"])

r_LDA <- LDATS::LDA(data = lda_data, ntopics = 2:5, nseeds = 2, ncores = 4)
r_LDATS <- LDATS::LDA_TS(lda_data, ts_data, formula = c("1", "time"),
                         ntopics = 2:5, nseeds = 2, ncores = 4, nit = 100)






