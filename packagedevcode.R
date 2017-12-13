##############################################################################
#
# package development code for LDA.pointbreak
#
# Latent Dirichlet Allocation with Break Point Analysis
#
# version 0.0.1 December 2017
#
# held under MIT
#
# #RequisiteJohnnyUtahReference
#
# Erica Christensen, David Harris, Hao Ye, and Juniper Simonis
#
##############################################################################

# This script is a work in progress transitioning from Erica and Dave's 
#  original LDA + break point analysis to a more general package of R code 
#  that works with and builds upon existing R topic model packages.
#
# The goal is to translate this script (originally a direct copy of the R
#  script for Erica and Dave's analysis) into the components of the package
#  and vignettes/examples based on the rodent data.
#  
# 

##############################################################################
#
# Preliminary Development Code
#
##############################################################################

  # load devtools

    library(devtools)

  # set the package location

    pkgloc <- getwd()

  # add dependencies to the description and load them here 

    pkgdpns <- c("topicmodels", "RCurl", "multipanelfigure", "reshape2",  
                 "dplyr", "memoise", "lubridate", "progress", "ggplot2",  
                 "viridis", "nnet", "RColorBrewer", "Rcpp",  
                 "tidyverse", "gridExtra", "topicmodels")

    for(i in 1:length(pkgdpns)){
      devtools::use_package(pkgdpns[i], "Imports", pkgloc)
    }

  # load the package

    devtools::load_all(devtools::as.package(pkgloc))

  # populate the documentation

    devtools::document(devtools::as.package(pkgloc))

  # define the pipe operator 

    `%>%` <- dplyr::`%>%`
    `period` <- lubridate::`period`

##############################################################################
# 
# Rodent Data Example
# 
##############################################################################

 # Prepare Data

  # species counts

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



 # Select number of topics
 #  Fit a bunch of LDA models with different seeds

  # set seeds 
  #  Only use even numbers for seeds because consecutive seeds give 
  #  identical results

    seeds = 2 * seq(200)

  # repeat LDA model fit and AIC calculation with a bunch of different 
  #  seeds to test robustness of the analysis

    best_ntopic <- repeat_VEM(dat,
                              seeds,
                              topic_min = 2,
                              topic_max = 8)

  # histogram of how many seeds chose how many topics

    hist(best_ntopic$k, breaks= seq(0.5, 9.5, 1), xlab = 'best # of topics', 
         main = '')

 # variability among seeds in composition with 4 topics

  # get the best 100 seeds where 4 topics was the best LDA model

    seeds_4topics <- best_ntopic %>%
                      dplyr::filter(k == 4) %>%
                      dplyr::arrange(aic) %>%
                      head(100) %>%
                      dplyr::pull(SEED)

  # choose seed with highest log likelihood for all following analyses
  #  (also produces plot of community composition for 'best' run compared 
  #  to 'worst')

    best_seed <- calculate_LDA_distance(dat, seeds_4topics)
    mean_dist <- unlist(best_seed)[2]
    max_dist <- unlist(best_seed)[3]

 # Baseline LDA model

  # inputs

    ntopics <- 4
    SEED <- unlist(best_seed)[1]

  # run model 

    ldamodel <- topicmodels::LDA(dat, ntopics, control = list(seed = SEED), 
                                 method = 'VEM')

 # Change point model

  # set up time for model

    year_continuous <- 1970 + as.integer(julian(dates)) / 365.25
    x <- data.frame(
                    year_continuous = year_continuous,
                    sin_year = sin(year_continuous * 2 * pi),
                    cos_year = cos(year_continuous * 2 * pi)
                    )

  # run models with 1, 2, 3, 4, 5 changepoints

    cp_results_rodent1 = changepoint_model(ldamodel, x, 1, 
                             weights = rep(1, length(year_continuous)))
    cp_results_rodent2 = changepoint_model(ldamodel, x, 2, 
                             weights = rep(1, length(year_continuous)))
    cp_results_rodent3 = changepoint_model(ldamodel, x, 3, 
                             weights = rep(1, length(year_continuous)))
    cp_results_rodent4 = changepoint_model(ldamodel, x, 4, 
                             weights = rep(1, length(year_continuous)))
    cp_results_rodent5 = changepoint_model(ldamodel, x, 5, 
                             weights = rep(1, length(year_continuous)))

  # some quick histograms of changepoint model results

    hist(year_continuous[cp_results_rodent4$saved[,1,]], 
         breaks = seq(1977,2016,.25), xlab = '', main = 'Changepoint Estimate')
    annual_hist(cp_results_rodent4, year_continuous)

  # turn changepoint results into data frame

    df_4 <- as.data.frame(t(cp_results_rodent1$saved[,1,])) %>% reshape2::melt()
    df_4$value <- year_continuous[df_4$value]

  # find 95% confidence intervals on each change point:

    quantile(df_4[df_4$variable == 'V1', 'value'], probs=c(.025, .975)) %>% 
        date_decimal() %>% format('%d-%m-%Y')
    quantile(df_4[df_4$variable == 'V2', 'value'], probs=c(.025, .975)) %>% 
        date_decimal() %>% format('%d-%m-%Y')
    quantile(df_4[df_4$variable == 'V3', 'value'], probs=c(.025, .975)) %>% 
        date_decimal() %>% format('%d-%m-%Y')
    quantile(df_4[df_4$variable == 'V4', 'value'], probs=c(.025, .975)) %>% 
        date_decimal() %>% format('%d-%m-%Y')

  # change point model selection
  #   mean deviance ( -2 * log likelihood) + 2*(#parameters)

    mean(cp_results_rodent1$saved_lls * -2) + 
       2 * (3 * (ntopics - 1) * (1 + 1) + (1))
    mean(cp_results_rodent2$saved_lls * -2) + 
       2 * (3 * (ntopics - 1) * (2 + 1) + (2))
    mean(cp_results_rodent3$saved_lls * -2) + 
       2 * (3 * (ntopics - 1) * (3 + 1) + (3))
    mean(cp_results_rodent4$saved_lls * -2) + 
       2 * (3 * (ntopics - 1) * (4 + 1) + (4))
    mean(cp_results_rodent5$saved_lls * -2) + 
       2 * (3 * (ntopics - 1) * (5 + 1) + (5))

# =================================================================
# 5. figures
# =================================================================

# plot community compositions
beta1 = community_composition(ldamodel)
# put columns in order of largest species to smallest
composition = beta1[,c('NA','DS','SH','SF','SO','DO','DM','PB','PH','OL','OT','PL','PM','PE','PP','PI','RF','RM','RO','BA','PF')]
plot_community_composition(composition,c(3,4,1,2))


# community composition with grassland communities highlighted
P = plot_community_composition_gg(composition,c(3,4,1,2),ylim=c(0,.8))

(figure_spcomp <- multi_panel_figure(
  width = c(70,70,70,70),
  height = c(70,10),
  panel_label_type = "none",
  column_spacing = 0))
figure_spcomp %<>% fill_panel(
  P[[1]],
  row = 1, column = 1)
figure_spcomp %<>% fill_panel(
  P[[2]],
  row = 1, column = 2)
figure_spcomp %<>% fill_panel(
  P[[3]],
  row = 1, column = 3)
figure_spcomp %<>% fill_panel(
  P[[4]],
  row = 1, column = 4)
figure_spcomp


# plot of component communities over time
cc = plot_component_communities(ldamodel,ntopics,dates)


# changepoint histogram w 4 cpts
H_4 = ggplot(data = df_4, aes(x=value)) +
  geom_histogram(data=subset(df_4,variable=='V1'),aes(y=..count../sum(..count..)),binwidth = .5,fill='black',alpha=.2) +
  geom_histogram(data=subset(df_4,variable=='V2'),aes(y=..count../sum(..count..)),binwidth = .5,fill='black',alpha=.4) +
  geom_histogram(data=subset(df_4,variable=='V3'),aes(y=..count../sum(..count..)),binwidth = .5,fill='black',alpha=.6) +
  geom_histogram(data=subset(df_4,variable=='V4'),aes(y=..count../sum(..count..)),binwidth = .5,fill='black',alpha=1) +
  labs(x='',y='') +
  xlim(range(year_continuous)) +
  scale_y_continuous(labels=c('0.00','0.20','0.40','0.60','0.80'),breaks = c(0,.2,.4,.6,.8)) +
  theme(axis.text=element_text(size=12),
        panel.border=element_rect(colour='black',fill=NA),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour='grey90'),
        panel.grid.minor = element_line(colour='grey90'))
  #theme_bw()
H_4


# changepoint model plot
cpts = find_changepoint_location(cp_results_rodent4)
cpt_plot = get_ll_non_memoized_plot(ldamodel,x,cpts,make_plot=T,weights=rep(1,length(year_continuous)))


# Figure 3 -- community composition, LDA model, changepoint histogram, changepoint timeseries
(figure <- multi_panel_figure(
  width = c(70,70,70,70),
  height = c(60,60,60,60),
  column_spacing = 0))
figure %<>% fill_panel(
  figure_spcomp,
  row = 1, column = 1:4)
figure %<>% fill_panel(
  cc,
  row = 2, column = 1:4)
figure %<>% fill_panel(
  H_4,
  row = 3, column = 1:4)
figure %<>% fill_panel(
  cpt_plot,
  row = 4, column = 1:4)
figure


