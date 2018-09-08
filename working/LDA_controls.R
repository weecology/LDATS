  data(rodents)
  lda_data <- dplyr::select(rodents, -c(newmoon, date, plots, traps))
  r_LDA <- parLDA(data = lda_data, ntopics = 2:7, nseeds = 10, ncores = 4)


AICv <- as.numeric(lapply(r_LDA, AIC))
ntops <- rep(2:5, each = 10)
plot(ntops, AICv)