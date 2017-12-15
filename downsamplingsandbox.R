#
# downsampling sandbox
#

  # data prepping

    dat <- create_rodent_table(period_first = 1,
                          period_last = 436,
                          selected_plots = c(2, 4, 8, 11, 12, 14, 17, 22),
                          selected_species = c('BA', 'DM', 'DO', 'DS', 'NA',
                                               'OL', 'OT', 'PB', 'PE', 'PF',
                                               'PH', 'PI', 'PL', 'PM', 'PP', 
                                               'RF', 'RM', 'RO', 'SF', 'SH',
                                               'SO'))

    # dates to go with count data

      moondat <- read.csv(text = 
                   RCurl::getURL(paste("https://raw.githubusercontent.com/",
                                     "weecology/PortalData/master/Rodents/",
                                     "moon_dates.csv", sep = "")),
                   stringsAsFactors = F)
      moondat$date <- as.Date(moondat$censusdate)
      period_dates <- dplyr::filter(moondat, period %in% rownames(dat)) %>% 
                      dplyr::select(period, date)
      dates <- period_dates$date

    # combine full data set
    
      data_full <- data.frame(dates, dat)

    # downsample

      nspecies <- ncol(dat)
      yr <- format(dates, "%Y")
      mo <- format(dates, "%m")
      years <- 1977:2015

      # year   


        nyears <- length(years)
        dat_year <- matrix(NA, nrow = nyears, ncol = nspecies)
        specific_date <- as.Date(rep(NA, nyears))

        for(i in 1:nyears){
 
          dates_within <- dates[which(yr == years[i])]
          specific_date[i] <- dates_within[which.min(dates_within)]
          specific_data <- which(dates == specific_date[i])
          dat_year[i, ] <- as.numeric(dat[specific_data, ])
        }

        colnames(dat_year) <- colnames(dat)
        data_year <- data.frame(date = specific_date, dat_year)

      # semiyear

        dat_semiyear <- matrix(NA, nrow = nyears * 2, ncol = nspecies)
        specific_date <- as.Date(rep(NA, nyears * 2))

        for(i in 1:nyears){
 
          
          dates_within <- dates[which(yr == years[i])]
          mo_within <- as.numeric(format(dates_within, "%m"))

          sem1 <- which(mo_within >= 1 & mo_within < 6) 
          sem2 <- which(mo_within >= 6) 


          if(length(sem1) > 0){

            e1 <- (i - 1) * 2 + 1
            specific_date[e1] <- (dates_within[sem1])[which.min(dates_within[sem1])]
            specific_data <- which(dates == specific_date[e1])
            dat_semiyear[e1, ] <- as.numeric(dat[specific_data, ])

          }

          if(length(sem2) > 0){

            e2 <- (i - 1) * 2 + 2
            specific_date[e2] <- (dates_within[sem2])[which.min(dates_within[sem2])]
            specific_data <- which(dates == specific_date[e2])
            dat_semiyear[e2, ] <- as.numeric(dat[specific_data, ])

          }

        }

        colnames(dat_semiyear) <- colnames(dat)
        data_semiyear <- data.frame(date = specific_date, dat_semiyear)
        data_semiyear <- na.omit(data_semiyear)

      # quarterly

        dat_quarter <- matrix(NA, nrow = nyears * 4, ncol = nspecies)
        specific_date <- as.Date(rep(NA, nyears * 4))

        for(i in 1:nyears){
 
          
          dates_within <- dates[which(yr == years[i])]
          mo_within <- as.numeric(format(dates_within, "%m"))

          q1 <- which(mo_within >= 1 & mo_within < 4) 
          q2 <- which(mo_within >= 4 & mo_within < 7) 
          q3 <- which(mo_within >= 7 & mo_within < 10) 
          q4 <- which(mo_within >= 10) 


          if(length(q1) > 0){

            e1 <- (i - 1) * 4 + 1
            specific_date[e1] <- (dates_within[q1])[which.min(dates_within[q1])]
            specific_data <- which(dates == specific_date[e1])
            dat_quarter[e1, ] <- as.numeric(dat[specific_data, ])

          }

          if(length(q2) > 0){

            e2 <- (i - 1) * 4 + 2
            specific_date[e2] <- (dates_within[q2])[which.min(dates_within[q2])]
            specific_data <- which(dates == specific_date[e2])
            dat_quarter[e2, ] <- as.numeric(dat[specific_data, ])

          }


          if(length(q3) > 0){

            e3 <- (i - 1) * 4 + 3
            specific_date[e3] <- (dates_within[q3])[which.min(dates_within[q3])]
            specific_data <- which(dates == specific_date[e3])
            dat_quarter[e3, ] <- as.numeric(dat[specific_data, ])

          }


          if(length(q4) > 0){

            e4 <- (i - 1) * 4 + 4
            specific_date[e4] <- (dates_within[q4])[which.min(dates_within[q4])]
            specific_data <- which(dates == specific_date[e4])
            dat_quarter[e4, ] <- as.numeric(dat[specific_data, ])

          }
        }

        colnames(dat_quarter) <- colnames(dat)
        data_quarter <- data.frame(date = specific_date, dat_quarter)
        data_quarter <- na.omit(data_quarter)



  # select the number of topics

    # set seeds

      seeds <- 2 * seq(20)

    # define topic range

      mint <- 2
      maxt <- 8

    # full data

      dat_full <- data_full[, -1]
      full_n_topics <- repeat_VEM(dat_full,
                                  seeds,
                                  topic_min = mint,
                                  topic_max = maxt)

    # quartely data

      dat_quarter <- data_quarter[, -1]
      quarter_n_topics <- repeat_VEM(dat_quarter,
                                  seeds,
                                  topic_min = mint,
                                  topic_max = maxt)

    # semi annual data

      dat_semiyear <- data_semiyear[, -1]
      semi_n_topics <- repeat_VEM(dat_semiyear,
                                  seeds,
                                  topic_min = mint,
                                  topic_max = maxt)

    # annual data 

      dat_year <- data_year[, -1]
      year_n_topics <- repeat_VEM(dat_year,
                                  seeds,
                                  topic_min = mint,
                                  topic_max = maxt)


    par(mfrow=c(4,1))
    hist(full_n_topics$k, breaks= seq(0.5, 9.5, 1),
         main = "Month", xlab = "")
    hist(quarter_n_topics$k, breaks= seq(0.5, 9.5, 1),
         main = "Quarter", xlab = "")
    hist(semi_n_topics$k, breaks= seq(0.5, 9.5, 1),
         main = "Semi-Year", xlab = "")
    hist(year_n_topics$k, breaks= seq(0.5, 9.5, 1),
         main = "Year", xlab = "Topics")

