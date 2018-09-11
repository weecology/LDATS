  data(rodents)
  lda_data <- dplyr::select(rodents, -c(newmoon, date, plots, traps))
  r_LDA <- parLDA(data = lda_data, ntopics = 2:7, nseeds = 10, ncores = 6,
                   )


AICv <- as.numeric(lapply(r_LDA, AIC))
ntops <- rep(2:7, each = 10)
postalph <- NULL
for(i in 1:length(r_LDA))
postalph[i]<-r_LDA[[i]]@alpha


par(mfrow=c(2,1))
plot(ntops, AICv)
plot(ntops, postalph)



  max_cores <- detectCores(logical = TRUE)
  seed_in <- rep(seq(2, nseeds * 2, 2), length(ntopics))
  k_in <- rep(ntopics, each = length(seq(2, nseeds * 2, 2)))
  nruns <- length(seed_in)

  if (ncores > max_cores){
    capped_cores <- floor(0.75 * parallel::detectCores(logical = TRUE))
    msg <- paste(ncores, " is larger than max cores available (", max_cores,
             "); capped at ", capped_cores, " cores.", sep = "")
    warning(msg)
    ncores <- capped_cores
  }
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  i <- 1
  mods <- foreach(i = 1:nruns, .packages = "topicmodels",
                  .errorhandling = "pass") %dopar% {

    topicmodels::LDA(data, k = k_in[i], control = list(seed = seed_in[i]))
  }

  stopCluster(cl)
  names(mods) <- paste("k: ", k_in, ", seed: ", seed_in, sep = "")
  class(mods) <- c("LDA_list", "list")



AICv <- as.numeric(lapply(mods, AIC))
ntops <- rep(2:7, each = 10)
postalph <- NULL
for(i in 1:length(mods))
postalph[i]<-mods[[i]]@alpha


par(mfrow=c(2,1))
plot(ntops, AICv)
plot(ntops, postalph)


