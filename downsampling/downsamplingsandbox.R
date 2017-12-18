#
# downsampling sandbox
#  deleting once the relevant code is included in the "downsampling_analysis.R" script 


    tiff("ntopics.tif", height = 8, width = 6, unit = "in", res = 200)
    par(mfrow=c(4,1))
    hist(full_n_topics$k, breaks= seq(0.5, 9.5, 1),
         main = "Month", xlab = "")
    hist(quarter_n_topics$k, breaks= seq(0.5, 9.5, 1),
         main = "Quarter", xlab = "")
    hist(semi_n_topics$k, breaks= seq(0.5, 9.5, 1),
         main = "Semi-Year", xlab = "")
    hist(year_n_topics$k, breaks= seq(0.5, 9.5, 1),
         main = "Year", xlab = "Topics")
    dev.off()

    # it's 4 for full, quarter, and semi-year, but 2 (barely) for year


  # changepoint model

    # full

      # Baseline LDA model

        # inputs [seed: 206]

          ntopics_full <- 4
          SEED_full <- full_n_topics$SEED[which.min(full_n_topics$aic)]

        # run model 

          full_bl_lda <- topicmodels::LDA(dat_full, ntopics_full, 
                                          control = list(seed = SEED_full), 
                                          method = 'VEM')

      # Change point model

        # set up time for model

          full_year_continuous <- 1970 + 
                                   as.integer(julian(data_full[,1])) / 365.25
          full_x <- data.frame(year_continuous = full_year_continuous,
                               sin_year = sin(full_year_continuous * 2 * pi),
                               cos_year = cos(full_year_continuous * 2 * pi))


        # run models with 1, 2, 3, 4, 5 changepoints

          full_cp_1 <- changepoint_model(full_bl_lda, full_x, 1, maxit = 1e4,
                              weights = rep(1, length(full_year_continuous)))
          full_cp_2 <- changepoint_model(full_bl_lda, full_x, 2, maxit = 1e4,
                              weights = rep(1, length(full_year_continuous)))
          full_cp_3 <- changepoint_model(full_bl_lda, full_x, 3, maxit = 1e4,
                              weights = rep(1, length(full_year_continuous)))
          full_cp_4 <- changepoint_model(full_bl_lda, full_x, 4, maxit = 1e4,
                              weights = rep(1, length(full_year_continuous)))
          full_cp_5 <- changepoint_model(full_bl_lda, full_x, 5, maxit = 1e4,
                              weights = rep(1, length(full_year_continuous)))

      # output
      
        tiff("cp_full.tif", height = 10, width = 6, unit = "in", res = 200)
        par(mfrow = c(5,1))
        annual_hist(full_cp_1, full_year_continuous)
        text(1973, 1.1e4, cex = 1.25, adj = -1, xpd = T,
         paste("AIC = ", 
         round(mean(full_cp_1$saved_lls * -2) + 2 * (3 * (ntopics - 1) * (1 + 1) + (1)), 2),
                   sep = ""))
        mtext(side = 3, "Monthly Censuses", line = 1)
        annual_hist(full_cp_2, full_year_continuous)
        text(1973, 1.1e4, cex = 1.25, adj = -1, xpd = T,
         paste("AIC = ", 
         round(mean(full_cp_2$saved_lls * -2) + 2 * (3 * (ntopics - 1) * (2 + 1) + (2)), 2),
                   sep = ""))
        annual_hist(full_cp_3, full_year_continuous)
        text(1973, 1.1e4, cex = 1.25, adj = -1, xpd = T,
         paste("AIC = ", 
         round(mean(full_cp_3$saved_lls * -2) + 2 * (3 * (ntopics - 1) * (3 + 1) + (3)), 2),
                   sep = ""))
        annual_hist(full_cp_4, full_year_continuous)
        text(1973, 1.1e4, cex = 1.25, adj = -1, xpd = T,
         paste("AIC = ", 
         round(mean(full_cp_4$saved_lls * -2) + 2 * (3 * (ntopics - 1) * (4 + 1) + (4)), 2),
                   sep = ""))
        annual_hist(full_cp_5, full_year_continuous)
        text(1973, 1.1e4, cex = 1.25, adj = -1, xpd = T,
         paste("AIC = ", 
         round(mean(full_cp_5$saved_lls * -2) + 2 * (3 * (ntopics - 1) * (5 + 1) + (5)), 2),
                   sep = ""))
        dev.off()

        ntopics <- 4
        mean(full_cp_1$saved_lls * -2) + 2 * (3 * (ntopics - 1) * (1 + 1) + (1))
        mean(full_cp_2$saved_lls * -2) + 2 * (3 * (ntopics - 1) * (2 + 1) + (2))
        mean(full_cp_3$saved_lls * -2) + 2 * (3 * (ntopics - 1) * (3 + 1) + (3))
        mean(full_cp_4$saved_lls * -2) + 2 * (3 * (ntopics - 1) * (4 + 1) + (4))
        mean(full_cp_5$saved_lls * -2) + 2 * (3  *(ntopics - 1) * (5 + 1) + (5))

    # quarter

      # Baseline LDA model

        # inputs [seed 318]

          ntopics_quarter <- 4
          SEED_quarter <- quarter_n_topics$SEED[which.min(quarter_n_topics$aic)]

        # run model 

          quarter_bl_lda <- topicmodels::LDA(dat_quarter, ntopics_quarter, 
                             control = list(seed = SEED_quarter), 
                             method = 'VEM')

      # Change point model

        # set up time for model

          quarter_year_continuous <- 1970 + 
                                 as.integer(julian(data_quarter[,1])) / 365.25
          quarter_x <- data.frame(
                    year_continuous = quarter_year_continuous,
                    sin_year = sin(quarter_year_continuous * 2 * pi),
                    cos_year = cos(quarter_year_continuous * 2 * pi)
                    )

        # run models with 1, 2, 3, 4, 5 changepoints

          quarter_cp_1 <- changepoint_model(quarter_bl_lda, quarter_x, 1, 
                                            maxit = 1e4,
                                            weights = rep(1, 
                                             length(quarter_year_continuous)))
          quarter_cp_2 <- changepoint_model(quarter_bl_lda, quarter_x, 2,
                                            maxit = 1e4,
                                            weights = rep(1, 
                                             length(quarter_year_continuous)))
          quarter_cp_3 <- changepoint_model(quarter_bl_lda, quarter_x, 3,
                                            maxit = 1e4,
                                            weights = rep(1, 
                                             length(quarter_year_continuous)))
          quarter_cp_4 <- changepoint_model(quarter_bl_lda, quarter_x, 4,
                                            maxit = 1e4,
                                            weights = rep(1, 
                                             length(quarter_year_continuous)))
          quarter_cp_5 <- changepoint_model(quarter_bl_lda, quarter_x, 5,
                                            maxit = 1e4,
                                            weights = rep(1, 
                                             length(quarter_year_continuous)))

      # output

        ntopics <- 4

        tiff("cp_quarter.tif", height = 10, width = 6, unit = "in", res = 200)
        par(mfrow = c(5,1))
        annual_hist(quarter_cp_1, quarter_year_continuous)
        text(1973, 1e4, cex = 1.25, adj = -1, xpd = T,
         paste("AIC = ", 
         round(mean(quarter_cp_1$saved_lls * -2) + 2 * (3 * (ntopics - 1) * (1 + 1) + (1)), 2),
                   sep = ""))
        mtext(side = 3, "Quarterly Censuses", line = 1)
        annual_hist(quarter_cp_2, quarter_year_continuous)
        text(1973, 1e4, cex = 1.25, adj = -1, xpd = T,
         paste("AIC = ", 
         round(mean(quarter_cp_2$saved_lls * -2) + 2 * (3 * (ntopics - 1) * (2 + 1) + (2)), 2),
                   sep = ""))
        annual_hist(quarter_cp_3, quarter_year_continuous)
        text(1973, 1e4, cex = 1.25, adj = -1, xpd = T,
         paste("AIC = ", 
         round(mean(quarter_cp_3$saved_lls * -2) + 2 * (3 * (ntopics - 1) * (3 + 1) + (3)), 2),
                   sep = ""))
        annual_hist(quarter_cp_4, quarter_year_continuous)
        text(1973, 1e4, cex = 1.25, adj = -1, xpd = T,
         paste("AIC = ", 
         round(mean(quarter_cp_4$saved_lls * -2) + 2 * (3 * (ntopics - 1) * (4 + 1) + (4)), 2),
                   sep = ""))
        annual_hist(quarter_cp_5, quarter_year_continuous)
        text(1973, 1e4, cex = 1.25, adj = -1, xpd = T,
         paste("AIC = ", 
         round(mean(quarter_cp_5$saved_lls * -2) + 2 * (3 * (ntopics - 1) * (5 + 1) + (5)), 2),
                   sep = ""))
        dev.off()

        mean(quarter_cp_1$saved_lls * -2) + 2 * (3 * (ntopics - 1) * (1 + 1) + (1))
        mean(quarter_cp_2$saved_lls * -2) + 2 * (3 * (ntopics - 1) * (2 + 1) + (2))
        mean(quarter_cp_3$saved_lls * -2) + 2 * (3 * (ntopics - 1) * (3 + 1) + (3))
        mean(quarter_cp_4$saved_lls * -2) + 2 * (3 * (ntopics - 1) * (4 + 1) + (4))
        mean(quarter_cp_5$saved_lls * -2) + 2 * (3  *(ntopics - 1) * (5 + 1) + (5))

    # semi

      # Baseline LDA model

        # inputs [seed 306]

          ntopics_semi <- 4
          SEED_semi <- semi_n_topics$SEED[which.min(semi_n_topics$aic)]

        # run model 

          semi_bl_lda <- topicmodels::LDA(dat_semiyear, ntopics_semi, 
                             control = list(seed = SEED_semi), 
                             method = 'VEM')

      # Change point model

        # set up time for model

          semi_year_continuous <- 1970 + as.integer(julian(data_semiyear[,1])) / 365.25
          semi_x <- data.frame(
                    year_continuous = semi_year_continuous,
                    sin_year = sin(semi_year_continuous * 2 * pi),
                    cos_year = cos(semi_year_continuous * 2 * pi)
                    )

        # run models with 1, 2, 3, 4, 5 changepoints

          semi_cp_1 <- changepoint_model(semi_bl_lda, semi_x, 1, maxit = 1e4,
                                   weights = rep(1, length(semi_year_continuous)))
          semi_cp_2 <- changepoint_model(semi_bl_lda, semi_x, 2, maxit = 1e4,
                                   weights = rep(1, length(semi_year_continuous)))
          semi_cp_3 <- changepoint_model(semi_bl_lda, semi_x, 3, maxit = 1e4,
                                   weights = rep(1, length(semi_year_continuous)))
          semi_cp_4 <- changepoint_model(semi_bl_lda, semi_x, 4, maxit = 1e4,
                                   weights = rep(1, length(semi_year_continuous)))
          semi_cp_5 <- changepoint_model(semi_bl_lda, semi_x, 5, maxit = 1e4,
                                   weights = rep(1, length(semi_year_continuous)))

      # output

        ntopics <- 4
        tiff("cp_semi.tif", height = 10, width = 6, unit = "in", res = 200)
        par(mfrow = c(5,1))
        annual_hist(semi_cp_1, semi_year_continuous)
        text(1973, 1e4, cex = 1.25, adj = -1, xpd = T,
          paste("AIC = ", 
           round(mean(semi_cp_1$saved_lls * -2) + 2 * (3 * (ntopics - 1) * (1 + 1) + (1)), 2),
                   sep = ""))
        mtext(side = 3, "Semi-Yearly Censuses", line = 1)
        annual_hist(semi_cp_2, semi_year_continuous)
        text(1973, 1e4, cex = 1.25, adj = -1, xpd = T,
          paste("AIC = ", 
           round(mean(semi_cp_2$saved_lls * -2) + 2 * (3 * (ntopics - 1) * (2 + 1) + (2)), 2),
                   sep = ""))
        annual_hist(semi_cp_3, semi_year_continuous)
        text(1973, 1e4, cex = 1.25, adj = -1, xpd = T,
          paste("AIC = ", 
           round(mean(semi_cp_3$saved_lls * -2) + 2 * (3 * (ntopics - 1) * (3 + 1) + (3)), 2),
                   sep = ""))
        annual_hist(semi_cp_4, semi_year_continuous)
        text(1973, 1e4, cex = 1.25, adj = -1, xpd = T,
          paste("AIC = ", 
           round(mean(semi_cp_4$saved_lls * -2) + 2 * (3 * (ntopics - 1) * (4 + 1) + (4)), 2),
                   sep = ""))
        annual_hist(semi_cp_5, semi_year_continuous)
        text(1973, 1e4, cex = 1.25, adj = -1, xpd = T,
          paste("AIC = ", 
           round(mean(semi_cp_5$saved_lls * -2) + 2 * (3 * (ntopics - 1) * (5 + 1) + (5)), 2),
                   sep = ""))
        dev.off()


        mean(semi_cp_1$saved_lls * -2) + 2 * (3 * (ntopics - 1) * (1 + 1) + (1))
        mean(semi_cp_2$saved_lls * -2) + 2 * (3 * (ntopics - 1) * (2 + 1) + (2))
        mean(semi_cp_3$saved_lls * -2) + 2 * (3 * (ntopics - 1) * (3 + 1) + (3))
        mean(semi_cp_4$saved_lls * -2) + 2 * (3 * (ntopics - 1) * (4 + 1) + (4))
        mean(semi_cp_5$saved_lls * -2) + 2 * (3  *(ntopics - 1) * (5 + 1) + (5))


    # year

      # Baseline LDA model

        # inputs [seed 82]

          ntopics_year <- 2
          SEED_year <- year_n_topics$SEED[which.min(year_n_topics$aic)]

        # run model 

          year_bl_lda <- topicmodels::LDA(dat_year, ntopics_year, 
                             control = list(seed = SEED_year), 
                             method = 'VEM')

      # Change point model

        # set up time for model

          year_year_continuous <- 1970 + as.integer(julian(data_year[,1])) / 365.25
          year_x <- data.frame(
                    year_continuous = year_year_continuous,
                    sin_year = sin(year_year_continuous * 2 * pi),
                    cos_year = cos(year_year_continuous * 2 * pi)
                    )

        # run models with 1, 2, 3, 4, 5 changepoints

          year_cp_1 <- changepoint_model(year_bl_lda, year_x, 1, maxit = 1e4,
                                   weights = rep(1, length(year_year_continuous)))
          year_cp_2 <- changepoint_model(year_bl_lda, year_x, 2, maxit = 1e4,
                                   weights = rep(1, length(year_year_continuous)))
          year_cp_3 <- changepoint_model(year_bl_lda, year_x, 3, maxit = 1e4,
                                   weights = rep(1, length(year_year_continuous)))
          year_cp_4 <- changepoint_model(year_bl_lda, year_x, 4, maxit = 1e4,
                                   weights = rep(1, length(year_year_continuous)))
          year_cp_5 <- changepoint_model(year_bl_lda, year_x, 5, maxit = 1e4,
                                   weights = rep(1, length(year_year_continuous)))

      # output

        ntopics <- 2
        tiff("cp_year.tif", height = 10, width = 6, unit = "in", res = 200)
        par(mfrow = c(5,1))
        annual_hist(year_cp_1, year_year_continuous)
        text(1973, 1e4, cex = 1.25, adj = -1, xpd = T,
          paste("AIC = ", 
           round(mean(year_cp_1$saved_lls * -2) + 2 * (3 * (ntopics - 1) * (1 + 1) + (1)), 2),
                   sep = ""))
        mtext(side = 3, "Yearly Censuses", line = 1)
        annual_hist(year_cp_2, year_year_continuous)
        text(1973, 1e4, cex = 1.25, adj = -1, xpd = T,
          paste("AIC = ", 
           round(mean(year_cp_2$saved_lls * -2) + 2 * (3 * (ntopics - 1) * (2 + 1) + (2)), 2),
                   sep = ""))
        annual_hist(year_cp_3, year_year_continuous)
        text(1973, 1e4, cex = 1.25, adj = -1, xpd = T,
          paste("AIC = ", 
           round(mean(year_cp_3$saved_lls * -2) + 2 * (3 * (ntopics - 1) * (3 + 1) + (3)), 2),
                   sep = ""))
        annual_hist(year_cp_4, year_year_continuous)
        text(1973, 1e4, cex = 1.25, adj = -1, xpd = T,
          paste("AIC = ", 
           round(mean(year_cp_4$saved_lls * -2) + 2 * (3 * (ntopics - 1) * (4 + 1) + (4)), 2),
                   sep = ""))
        annual_hist(year_cp_5, year_year_continuous)
        text(1973, 1e4, cex = 1.25, adj = -1, xpd = T,
          paste("AIC = ", 
           round(mean(year_cp_5$saved_lls * -2) + 2 * (3 * (ntopics - 1) * (5 + 1) + (5)), 2),
                   sep = ""))
        dev.off()


        mean(year_cp_1$saved_lls * -2) + 2 * (3 * (ntopics - 1) * (1 + 1) + (1))
        mean(year_cp_2$saved_lls * -2) + 2 * (3 * (ntopics - 1) * (2 + 1) + (2))
        mean(year_cp_3$saved_lls * -2) + 2 * (3 * (ntopics - 1) * (3 + 1) + (3))
        mean(year_cp_4$saved_lls * -2) + 2 * (3 * (ntopics - 1) * (4 + 1) + (4))
        mean(year_cp_5$saved_lls * -2) + 2 * (3  *(ntopics - 1) * (5 + 1) + (5))


