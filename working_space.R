# for code dev: working through an example based on the rodent data

  # prep the data

    # grab the control plot rodent data

      dat <- create_rodent_table(period_first = 1, period_last = 436,
                                 selected_plots = c(2, 4, 8, 11, 12, 14, 
                                                    17, 22),
                                 selected_species = c('BA', 'DM', 'DO', 'DS', 
                                                      'NA', 'OL', 'OT', 'PB', 
                                                      'PE', 'PF', 'PH', 'PI', 
                                                      'PL', 'PM', 'PP', 
                                                      'RF', 'RM', 'RO', 
                                                      'SF', 'SH', 'SO')) 

    # grab the dates to go with count data

      moondat <- read.csv(text = RCurl::getURL(paste(
                            "https://raw.githubusercontent.com/",
                             "weecology/PortalData/master/Rodents/",
                             "moon_dates.csv", sep = "")),
                          stringsAsFactors = F)

      moondat$date <- as.Date(moondat$censusdate)
      period_dates <- dplyr::filter(moondat, period %in% rownames(dat)) %>% 
                      dplyr::select(period, date)
      dates <- period_dates$date

    # combine the counts and the dates to create the full data set
    
      data_full <- data.frame(dates, dat)

  # run the LDA 

    rodent_LDA <- batch_LDA(data = dat, ntopics = 2:6, nseeds = 100, 
                            method = "VEM", sort = TRUE, sortby = "aicc",
                            parallel = TRUE, ncores = 8)


  # looking at model weights for considering how best to progress

    MS <- rodent_LDA$ModelSummaries
    minAICc <- min(MS[, "aicc"])
    deltAICc <- MS[, "aicc"] - minAICc


    expnhdaic <- exp(-0.5 * deltAICc)
    sumexpnhdaic <- sum(expnhdaic)

    mw <- expnhdaic / sumexpnhdaic 
    plot(mw)
 
    length(which(mw > 0.001)) / length(mw)

    plot(MS[,2], mw)

  
  # run a the cp model

    rodent_cp <- cp_models(data = dat, dates = dates, ntopics = 3, SEED = 68, 
                           weights = rep(1, nrow(dat)), 
                           maxit = 1e4, maxcps = 4)





