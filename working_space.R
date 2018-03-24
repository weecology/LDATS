# working through an example based on the Portal rodent data

"%>%" <- dplyr::"%>%"
"period" <- lubridate::"period"
"%dopar%" <- foreach::"%dopar%"
"%dorng%" <- doRNG::"%dorng%"

rodent_data <- LDATS::create_rodent_table()
lda_data <- rodent_data %>%
            dplyr::select(-c(period, censusdate, nplots, ntraps))

# next up: formatting and such with the batch_LDA function
# def need to add some progress bar or something, yeah?

rodent_LDA <- LDATS::batch_LDA(data = lda_data, ntopics = 2:6, nseeds = 100, 
                method = "VEM", sort = TRUE, sortby = "aicc",
                parallel = TRUE, ncores = 8)

  model_summary <- rodent_LDA$ModelSummaries
  min_AICc <- min(model_summary[ , "aicc"])
  delta_AICc <- model_summary[ , "aicc"] - min_AICc

  exp_half_d_aic <- exp(-0.5 * delta_AICc)
  sum_exp_half_d_aic <- sum(exp_half_d_aic)
  model_wt <- exp_half_d_aic / sum_exp_half_d_aic 
  plot(model_wt)
 
  length(which(model_wt > 0.001)) / length(model_wt)

  plot(model_summary[ , 2], model_wt)


  rodent_cp <- cp_models(data = dat, dates = dates, ntopics = 3, SEED = 68, 
                 weights = rep(1, nrow(rodent_data)), maxit = 1e4, maxcps = 4)





