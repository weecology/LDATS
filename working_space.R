# for code dev

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

  # run an ldats



  # run the LDA model

      # set the seeds to use

        seeds <- 2 * seq(200)

      # define topic range

        mint <- 2
        maxt <- 8

        full_n_topics <- repeat_VEM(dat, seeds,
                                    topic_min = mint, topic_max = maxt)

  # use the best LDA model to run the change point model

        full_cp <- cp_models(data = data_full, ntopics = 4, SEED = 206, 
                                     weights = rep(1, nrow(data_full)), 
                                     maxit = 1e4, maxcps = 5)




data = data_full
ntopics = 4
SEED = 206
weights = rep(1, nrow(data_full))
maxit = 1e4
maxcps = 5




