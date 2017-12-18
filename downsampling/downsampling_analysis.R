

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


  # downsample the data

    # quarterly (frequency = 4)

      # first sample of the quarter

        quarterly_first <- downsample(data = data_full, frequency = 4, 
                                      span = c("1978-01-01", "2015-01-01"), 
                                      selection_random = FALSE, 
                                      selection_order = 1)

      # second sample of the quarter

        quarterly_second  <- downsample(data = data_full, frequency = 4, 
                                        span = c("1978-01-01", "2015-01-01"), 
                                        selection_random = FALSE, 
                                        selection_order = 2)

        # some quarters only had one sample, so we get NAs
        #  drop the NAs

          quarterly_second <- na.omit(quarterly_second)

      # random samples from the quarter

        set.seed(123)
        quarterly_random_1  <- downsample(data = data_full, frequency = 4, 
                                          span = c("1978-01-01", 
                                                   "2015-01-01"), 
                                          selection_random = TRUE)

        set.seed(321)
        quarterly_random_2  <- downsample(data = data_full, frequency = 4, 
                                          span = c("1978-01-01", 
                                                   "2015-01-01"), 
                                          selection_random = TRUE)

        set.seed(456)
        quarterly_random_3  <- downsample(data = data_full, frequency = 4, 
                                          span = c("1978-01-01", 
                                                   "2015-01-01"), 
                                          selection_random = TRUE)

    # semi-yearly (frequency = 2)

      # first sample of the half year

        semi_first <- downsample(data = data_full, frequency = 2, 
                                 span = c("1978-01-01", "2015-01-01"), 
                                 selection_random = FALSE, 
                                 selection_order = 1)

      # third sample of the half year

        semi_third <- downsample(data = data_full, frequency = 2, 
                                 span = c("1978-01-01", "2015-01-01"), 
                                 selection_random = FALSE, 
                                 selection_order = 3)

      # random samples from the half year

        set.seed(123)
        semi_random_1 <- downsample(data = data_full, frequency = 2, 
                                      span = c("1978-01-01", "2015-01-01"), 
                                      selection_random = TRUE)

        set.seed(321)
        semi_random_2 <- downsample(data = data_full, frequency = 2, 
                                      span = c("1978-01-01", "2015-01-01"), 
                                      selection_random = TRUE)

        set.seed(456)
        semi_random_3 <- downsample(data = data_full, frequency = 2, 
                                      span = c("1978-01-01", "2015-01-01"), 
                                      selection_random = TRUE)


    # yearly (frequency = 1)

      # first sample of the year

        yearly_first <- downsample(data = data_full, frequency = 1, 
                                   span = c("1978-01-01", "2015-01-01"), 
                                   selection_random = FALSE, 
                                   selection_order = 1)

      # sixth sample of the year

        yearly_sixth <- downsample(data = data_full, frequency = 1, 
                                   span = c("1978-01-01", "2015-01-01"), 
                                   selection_random = FALSE, 
                                   selection_order = 6)

      # random samples from the year

        set.seed(123)
        yearly_random_1 <- downsample(data = data_full, frequency = 1, 
                                      span = c("1978-01-01", "2015-01-01"), 
                                      selection_random = TRUE)

        set.seed(321)
        yearly_random_2 <- downsample(data = data_full, frequency = 1, 
                                      span = c("1978-01-01", "2015-01-01"), 
                                      selection_random = TRUE)

        set.seed(456)
        yearly_random_3 <- downsample(data = data_full, frequency = 1, 
                                      span = c("1978-01-01", "2015-01-01"), 
                                      selection_random = TRUE)


  # Determine the number of topics

    # prep work

      # set the seeds to use

        seeds <- 2 * seq(200)

      # define topic range

        mint <- 2
        maxt <- 8


    # full

        full_n_topics <- repeat_VEM(data_full[ , -1], seeds,
                                    topic_min = mint, topic_max = maxt)

    # quarterly (frequency = 4)

      # first sample of the quarter

        quarterly_first_n_topics <- repeat_VEM(quarterly_first[ , -1],
                                               seeds, topic_min = mint,
                                               topic_max = maxt)

      # second sample of the quarter

        quarterly_second_n_topics <- repeat_VEM(quarterly_second[ , -1], 
                                                seeds, topic_min = mint,
                                                topic_max = maxt)

      # random samples from the quarter

        quarterly_r1_n_topics <- repeat_VEM(quarterly_random_1[ , -1], seeds,
                                           topic_min = mint, topic_max = maxt)

        quarterly_r2_n_topics <- repeat_VEM(quarterly_random_2[ , -1], seeds,
                                           topic_min = mint, topic_max = maxt)

        quarterly_r3_n_topics <- repeat_VEM(quarterly_random_3[ , -1], seeds,
                                           topic_min = mint, topic_max = maxt)

    # semi-yearly (frequency = 2)

      # first sample of the half year

        semi_first_n_topics <- repeat_VEM(semi_first[ , -1], seeds,
                                          topic_min = mint, topic_max = maxt)

      # third sample of the half year

        semi_third_n_topics <- repeat_VEM(semi_third[ , -1], seeds,
                                          topic_min = mint, topic_max = maxt)

      # random samples from the half year

        semi_r1_n_topics <- repeat_VEM(semi_random_1[ , -1], seeds,
                                       topic_min = mint, topic_max = maxt)

        semi_r2_n_topics <- repeat_VEM(semi_random_2[ , -1],seeds,
                                       topic_min = mint, topic_max = maxt)

        semi_r3_n_topics <- repeat_VEM(semi_random_3[ , -1], seeds,
                                       topic_min = mint, topic_max = maxt)

    # yearly (frequency = 1)

      # first sample of the year

        yearly_first_n_topics <- repeat_VEM(yearly_first[ , -1], seeds,
                                           topic_min = mint, topic_max = maxt)

      # sixth sample of the year

        yearly_sixth_n_topics <- repeat_VEM(yearly_sixth[ , -1], seeds,
                                           topic_min = mint, topic_max = maxt)

      # random samples from the year

        yearly_r1_n_topics <- repeat_VEM(yearly_random_1[ , -1], seeds,
                                         topic_min = mint, topic_max = maxt)

        yearly_r2_n_topics <- repeat_VEM(yearly_random_2[ , -1], seeds,
                                         topic_min = mint, topic_max = maxt)

        yearly_r3_n_topics <- repeat_VEM(yearly_random_3[ , -1], seeds,
                                         topic_min = mint, topic_max = maxt)


  # save the data and LDA ntopic models

    topackage <- c("data_full", "quarterly_first", "quarterly_second",   
                   "quarterly_random_1", "quarterly_random_2", 
                   "quarterly_random_3", "semi_first", "semi_third", 
                   "semi_random_1", "semi_random_2", "semi_random_3", 
                   "yearly_first", "yearly_sixth", "yearly_random_1", 
                   "yearly_random_2", "yearly_random_3",                 
                   "full_n_topics", "quarterly_first_n_topics",
                   "quarterly_second_n_topics", "quarterly_r1_n_topics",
                   "quarterly_r2_n_topics", "quarterly_r3_n_topics",
                   "semi_first_n_topics", "semi_third_n_topics",
                   "semi_r1_n_topics", "semi_r2_n_topics", "semi_r3_n_topics",
                   "yearly_first_n_topics", "yearly_sixth_n_topics",
                   "yearly_r1_n_topics", "yearly_r2_n_topics",
                   "yearly_r3_n_topics")

    save(list = topackage, file = "data_and_LDA_ntopic_models.RData")


  # output figures

    # ntopics

      ntopic_hist_file(dd = full_n_topics, nn = "Full_Data_Set")

      ntopic_hist_file(dd = quarterly_first_n_topics, nn = "Quarterly_1st")
      ntopic_hist_file(dd = quarterly_second_n_topics, nn = "Quarterly_2nd")
      ntopic_hist_file(dd = quarterly_r1_n_topics, nn = "Quarterly_R1")
      ntopic_hist_file(dd = quarterly_r2_n_topics, nn = "Quarterly_R2")
      ntopic_hist_file(dd = quarterly_r3_n_topics, nn = "Quarterly_R3")

      ntopic_hist_file(dd = semi_first_n_topics, nn = "Semi_1st")
      ntopic_hist_file(dd = semi_third_n_topics, nn = "Semi_3rd")
      ntopic_hist_file(dd = semi_r1_n_topics, nn = "Semi_R1")
      ntopic_hist_file(dd = semi_r2_n_topics, nn = "Semi_R2")
      ntopic_hist_file(dd = semi_r3_n_topics, nn = "Semi_R3")

      ntopic_hist_file(dd = yearly_first_n_topics, nn = "Yearly_1st")
      ntopic_hist_file(dd = yearly_sixth_n_topics, nn = "Yearly_6th")
      ntopic_hist_file(dd = yearly_r1_n_topics, nn = "Yearly_R1")
      ntopic_hist_file(dd = yearly_r2_n_topics, nn = "Yearly_R2")
      ntopic_hist_file(dd = yearly_r3_n_topics, nn = "Yearly_R3")




  # Run the change point models
  # 
  # Current issues with NAs generated by the change point model function
  #  these haven't been successfully run yet

    # full

        full_cp <- cp_models(data = data_full, ntopics = 4, SEED = 206, 
                                     weights = rep(1, nrow(data_full)), 
                                     maxit = 1e4, maxcps = 5)

    # quarterly (frequency = 4)

      # first sample of the quarter

        quarterly_first_cp <- cp_models(data = quarterly_first, ntopics = 4,  
                                     weights = rep(1, nrow(quarterly_first)), 
                                     SEED = 318, maxit = 1e4, maxcps = 5)

      # second sample of the quarter

        quarterly_second_cp <- cp_models(data = quarterly_second, ntopics = 4,  
                                     weights = rep(1, nrow(quarterly_second)), 
                                     SEED = 276, maxit = 1e4, maxcps = 5)

      # random samples from the quarter

        quarterly_r1_cp <- cp_models(data = quarterly_random_1, ntopics = 4,  
                                     weights = rep(1, nrow(quarterly_first)), 
                                     SEED = 352, maxit = 1e4, maxcps = 5)

        quarterly_r2_cp <- cp_models(data = quarterly_random_2, ntopics = 4,  
                                     weights = rep(1, nrow(quarterly_first)), 
                                     SEED = 28, maxit = 1e4, maxcps = 5)

        quarterly_r3_cp <- cp_models(data = quarterly_random_3, ntopics = 4,  
                                     weights = rep(1, nrow(quarterly_first)), 
                                     SEED = 214, maxit = 1e4, maxcps = 5)

    # semi-yearly (frequency = 2)

      # first sample of the half year

        semi_first_cp <- cp_models(data = semi_first, ntopics = 4, 
                                   weights = rep(1, nrow(semi_first)), 
                                   SEED = 176, maxit = 1e4, maxcps = 5)

      # third sample of the half year

        semi_third_cp <- cp_models(data = semi_third, ntopics = 4, 
                                   weights = rep(1, nrow(semi_third)), 
                                   SEED = 314, maxit = 1e4, maxcps = 5)

      # random samples from the half year

        semi_r1_cp <- cp_models(data = semi_random_1, ntopics = 4, 
                                weights = rep(1, nrow(semi_random_1)), 
                                SEED = 62, maxit = 1e4, maxcps = 5)

        semi_r2_cp <- cp_models(data = semi_random_2, ntopics = 3, 
                                weights = rep(1, nrow(semi_random_2)), 
                                SEED = 186, maxit = 1e4, maxcps = 5)

        semi_r3_cp <- cp_models(data = semi_random_3, ntopics = 4, 
                                weights = rep(1, nrow(semi_random_3)), 
                                SEED = 324, maxit = 1e4, maxcps = 5)

    # yearly (frequency = 1)

      # first sample of the year

        yearly_first_cp <- cp_models(data = yearly_first, ntopics = 3, 
                                     weights = rep(1, nrow(yearly_first)), 
                                     SEED = 82, maxit = 1e4, maxcps = 5)

      # sixth sample of the year

        yearly_sixth_cp <- cp_models(data = yearly_sixth, ntopics = 3, 
                                     weights = rep(1, nrow(yearly_sixth)), 
                                     SEED = 6, maxit = 1e4, maxcps = 5)

      # random samples from the year

        yearly_r1_cp <- cp_models(data = yearly_random_1, ntopics = 3, 
                                  weights = rep(1, nrow(yearly_random_1)), 
                                  SEED = 330, maxit = 1e4, maxcps = 5)

        yearly_r2_cp <- cp_models(data = yearly_random_2, ntopics = 3, 
                                  weights = rep(1, nrow(yearly_random_2)), 
                                  SEED = 18, maxit = 1e4, maxcps = 5)

        yearly_r3_cp <- cp_models(data = yearly_random_3, ntopics = 4, 
                                  weights = rep(1, nrow(yearly_random_3)), 
                                  SEED = 92, maxit = 1e4, maxcps = 5)



  # save the data and LDA ntopic and change point models

    topackage2 <- c("full_cp", "quarterly_first_cp", "quarterly_second_cp", 
                    "quarterly_r1_cp", "quarterly_r2_cp", "quarterly_r3_cp")


    save(list = topackage2, 
         file = "changepoint_models_full_and_quarterly.RData")


    topackage3 <- c("semi_first_cp", "semi_third_cp", "semi_r1_cp",  
                   "semi_r2_cp","semi_r3_cp", "yearly_first_cp",  
                   "yearly_sixth_cp", "yearly_r1_cp", "yearly_r2_cp", 
                   "yearly_r3_cp")


    save(list = topackage3, 
         file = "changepoint_models_semi_and_yearly.RData")




  # output figures

    # ntopics

      ntopic_hist_file(dd = full_n_topics, nn = "Full_Data_Set")

      ntopic_hist_file(dd = quarterly_first_n_topics, nn = "Quarterly_1st")
      ntopic_hist_file(dd = quarterly_second_n_topics, nn = "Quarterly_2nd")
      ntopic_hist_file(dd = quarterly_r1_n_topics, nn = "Quarterly_R1")
      ntopic_hist_file(dd = quarterly_r2_n_topics, nn = "Quarterly_R2")
      ntopic_hist_file(dd = quarterly_r3_n_topics, nn = "Quarterly_R3")

      ntopic_hist_file(dd = semi_first_n_topics, nn = "Semi_1st")
      ntopic_hist_file(dd = semi_third_n_topics, nn = "Semi_3rd")
      ntopic_hist_file(dd = semi_r1_n_topics, nn = "Semi_R1")
      ntopic_hist_file(dd = semi_r2_n_topics, nn = "Semi_R2")
      ntopic_hist_file(dd = semi_r3_n_topics, nn = "Semi_R3")

      ntopic_hist_file(dd = yearly_first_n_topics, nn = "Yearly_1st")
      ntopic_hist_file(dd = yearly_sixth_n_topics, nn = "Yearly_6th")
      ntopic_hist_file(dd = yearly_r1_n_topics, nn = "Yearly_R1")
      ntopic_hist_file(dd = yearly_r2_n_topics, nn = "Yearly_R2")
      ntopic_hist_file(dd = yearly_r3_n_topics, nn = "Yearly_R3")

    # change points

      cp_hists_file(dd = data_full, cpms = full_cp, 
                    ntopics = 4, nn = "Full_Data_Set")


