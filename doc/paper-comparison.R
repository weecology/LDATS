## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)


## ---- eval=FALSE---------------------------------------------------------
#  install.packages("devtools")
#  devtools::install_github("weecology/LDATS")

## ---- eval = T-----------------------------------------------------------
library(LDATS)
library(dplyr)
load(here::here('vignettes', 'christensenetal-comparison-files',  'data_files.Rdata'))
load(here::here('vignettes', 'christensenetal-comparison-files', 'model-cache',  'ldas.RData'))
load(here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'ldats_changepoints.RData'))
load(here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'paper_changepoints.RData'))

## ----download Christensen, eval = F--------------------------------------
#  
#  paper_filepath <- here::here('vignettes', 'christensenetal-comparison-files', 'christensenetal-2018-files')
#  test_filepath <- here::here('vignettes', 'christensenetal-comparison-files', 'christensenetal-2018-files', 'rodent_LDA_analysis.r')
#  
#  if(!file.exists(test_filepath)) {
#  
#  files_to_download <- c('rodent_LDA_analysis.r', 'rodent_data_for_LDA.r', 'AIC_model_selection.R', 'changepointmodel.r', 'LDA-distance.R', 'Rodent_table_dat.csv', 'LDA_figure_scripts.R')
#  
#  for(i in 1:length(files_to_download)) {
#    download.file(url = paste0("https://raw.githubusercontent.com/emchristensen/Extreme-events-LDA/master/", files_to_download[i]),
#                  destfile = paste0(paper_filepath, '/', files_to_download[i]))
#  }
#  
#  rm(files_to_download)
#  rm(i)
#  
#  }
#  
#  rm(test_filepath)
#  rm(paper_filepath)

## ----LDATS data----------------------------------------------------------
data(rodents)

head(rodents[[1]])

## ----Paper data, eval = F------------------------------------------------
#  source(here::here('vignettes', 'christensenetal-comparison-files', 'christensenetal-2018-files', 'rodent_data_for_LDA.r'))
#  
#  dat = create_rodent_table(period_first = 1,
#                            period_last = 436,
#                            selected_plots = c(2,4,8,11,12,14,17,22),
#                            selected_species = c('BA','DM','DO','DS','NA','OL','OT','PB','PE','PF','PH','PI','PL','PM','PP','RF','RM','RO','SF','SH','SO'))
#  
#  # dates to go with count data
#  moondat = read.csv(text=getURL("https://raw.githubusercontent.com/weecology/PortalData/master/Rodents/moon_dates.csv"),stringsAsFactors = F)
#  moondat$date = as.Date(moondat$censusdate)
#  
#  period_dates = filter(moondat,period %in% rownames(dat)) %>% select(period,date)
#  dates = period_dates$date
#  
#  paper_dat <- dat
#  
#  paper_dates <- dates
#  
#  rm(list = c('dat', 'dates', 'period_dates',
#              'create_rodent_table'))
#  

## ----rodent data comparison----------------------------------------------

compare <- rodents[[1]] == paper_dat

length(which(rowSums(compare) < ncol(compare)))

## ----adjust LDATS data after Christensen et al, eval = F-----------------
#  
#  trap_table = read.csv('https://raw.githubusercontent.com/weecology/PortalData/master/Rodents/Portal_rodent_trapping.csv')
#  

## ----continue adjust LDATS after paper-----------------------------------
 trap_table_controls =dplyr::filter(trap_table, plot %in% c(2,4,8,11,12,14,17,22))
  nplots_controls = aggregate(trap_table_controls$sampled,by=list(period = trap_table_controls$period),FUN=sum)
   # adjust species counts by number of plots trapped that month
  ldats_rodents_adjusted = as.data.frame.matrix(rodents[[1]])
  for (n in 1:436) {
    #divide by number of control plots actually trapped (should be 8) and multiply by 8 to estimate captures as if all plots were trapped
    ldats_rodents_adjusted[n,] = round(ldats_rodents_adjusted[n,]/nplots_controls$x[n]*8)
  }


compare_raw <- rodents[[1]] == ldats_rodents_adjusted
length(which(rowSums(compare_raw) < ncol(compare_raw)))

compare_adjusted <- ldats_rodents_adjusted == paper_dat
length(which(rowSums(compare_adjusted) < ncol(compare_adjusted)))


## ----data adjustment cleanup, include = FALSE----------------------------
rm(list = c('compare', 'compare_adjusted', 'compare_raw', 'n', 'ldats_rodents_adjusted', 'nplots_controls', 'trap_table', 'trap_table_controls'))

## ----switch to adjusted rodents------------------------------------------
rodents[[1]] <- paper_dat

## ----add dates to covariate table----------------------------------------

head(rodents$document_covariate_table)

moondat <- dplyr::select(moondat, newmoonnumber, censusdate)
colnames(moondat) <- c('newmoon', 'censusdate')

new_cov_table <- dplyr::left_join(rodents$document_covariate_table, moondat, by = 'newmoon')

rodents$document_covariate_table <- new_cov_table


## ----cleanup dates, include = F------------------------------------------
rm(moondat)
rm(new_cov_table)

## ----LDATS LDAs, eval = F------------------------------------------------
#  
#  ldats_ldas <- LDATS::LDA_set(document_term_table = rodents$document_term_table, topics = c(2:6), nseeds = 100)
#  ldats_lda_selected <- LDATS::select_LDA(LDA_models = ldats_ldas)
#  
#  save(ldats_ldas, file = here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'ldats_ldas.Rds'))
#  save(ldats_lda_selected, file = here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'ldats_lda_selected.Rds'))
#  

## ----create dat for paper lda, include = F-------------------------------
dat = paper_dat

## ----paper LDAs, eval = F------------------------------------------------
#  source(here::here('vignettes', 'christensenetal-comparison-files', 'christensenetal-2018-files', 'AIC_model_selection.R'))
#  source(here::here('vignettes', 'christensenetal-comparison-files', 'christensenetal-2018-files', 'LDA-distance.R'))
#  
#  # Fit a bunch of LDA models with different seeds
#  # Only use even numbers for seeds because consecutive seeds give identical results
#  seeds = 2*seq(200)
#  
#  # repeat LDA model fit and AIC calculation with a bunch of different seeds to test robustness of the analysis
#  best_ntopic = repeat_VEM(paper_dat,
#                           seeds,
#                           topic_min=2,
#                           topic_max=6)
#  hist(best_ntopic$k,breaks=c(0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5),xlab='best # of topics', main='')
#  
#  # 2b. how different is species composition of 4 community-types when LDA is run with different seeds?
#  # ==================================================================
#  # get the best 100 seeds where 4 topics was the best LDA model
#  seeds_4topics = best_ntopic %>%
#    filter(k == 4) %>%
#    arrange(aic) %>%
#    head(100) %>%
#    pull(SEED)
#  
#  # choose seed with highest log likelihood for all following analyses
#  #    (also produces plot of community composition for 'best' run compared to 'worst')
#  best_seed = calculate_LDA_distance(paper_dat,seeds_4topics)
#  mean_dist = unlist(best_seed)[2]
#  max_dist = unlist(best_seed)[3]
#  
#  # ==================================================================
#  # 3. run LDA model
#  # ==================================================================
#  ntopics = 4
#  SEED = unlist(best_seed)[1]  # For the paper, use seed 206
#  ldamodel = LDA(paper_dat,ntopics, control = list(seed = SEED),method='VEM')
#  
#  save(ldamodel,file = here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'paper_ldamodel.Rds'))

## ----plot paper LDA, fig.width=7, fig.height=6---------------------------
# Paper
plot(ldamodel, cols = NULL, option = "D")


## ----plot LDATS LDA, fig.width=7, fig.height=6---------------------------
# LDATS
plot(ldats_lda_selected[[1]], cols = NULL, option = "D")

## ----run LDATS LDA and paper cpt, eval = F-------------------------------
#  #### Run changepoint ####
#  source(here::here('vignettes', 'christensenetal-comparison-files', 'christensenetal-2018-files', 'changepointmodel.r'))
#  
#  # set up parameters for model
#  year_continuous = 1970 + as.integer(julian(paper_dates)) / 365.25
#  x = data.frame(
#    year_continuous = year_continuous,
#    sin_year = sin(year_continuous * 2 * pi),
#    cos_year = cos(year_continuous * 2 * pi)
#  )
#  
#  # run models with 1, 2, 3, 4, 5, 6 changepoints
#  ldats_paper_cp_results_rodent = changepoint_model(ldats_lda_selected[[1]], x, 1, weights = rep(1,length(year_continuous)))
#  save(ldats_paper_cp_results_rodent, file = here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'ldats_paper_1.Rds'))
#  rm(ldats_paper_cp_results_rodent)
#  
#  ldats_paper_cp_results_rodent2 = changepoint_model(ldats_lda_selected[[1]], x, 2, weights = rep(1,length(year_continuous)))
#  save(ldats_paper_cp_results_rodent2, file = here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'ldats_paper_2.Rds'))
#  rm(ldats_paper_cp_results_rodent2)
#  
#  ldats_paper_cp_results_rodent3 = changepoint_model(ldats_lda_selected[[1]], x, 3, weights = rep(1,length(year_continuous)))
#  save(ldats_paper_cp_results_rodent3, file = here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'ldats_paper_3.Rds'))
#  rm(ldats_paper_cp_results_rodent3)
#  
#  ldats_paper_cp_results_rodent4 = changepoint_model(ldats_lda_selected[[1]], x, 4, weights = rep(1,length(year_continuous)))
#  save(ldats_paper_cp_results_rodent4, file = here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'ldats_paper_4.Rds'))
#  rm(ldats_paper_cp_results_rodent4)
#  
#  ldats_paper_cp_results_rodent5 = changepoint_model(ldats_lda_selected[[1]], x, 5, weights = rep(1,length(year_continuous)))
#  save(ldats_paper_cp_results_rodent5, file = here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'ldats_paper_5.Rds'))
#  rm(ldats_paper_cp_results_rodent5)
#  
#  ldats_paper_cp_results_rodent6 = changepoint_model(ldats_lda_selected[[1]], x, 6, weights = rep(1,length(year_continuous)))
#  save(ldats_paper_cp_results_rodent6, file = here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'ldats_paper_6.Rds'))
#  rm(ldats_paper_cp_results_rodent6)

## ----select LDATS LDA and paper cpt, eval = F----------------------------
#  #### Changepoint model selection ####
#  load(here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'ldats_paper_1.Rds'))
#  load(here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'ldats_paper_2.Rds'))
#  load(here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'ldats_paper_3.Rds'))
#  load(here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'ldats_paper_4.Rds'))
#  load(here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'ldats_paper_5.Rds'))
#  load(here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'ldats_paper_6.Rds'))
#  

## ----select LDATS-paper, eval = T----------------------------------------
ntopics = ldats_lda_selected[[1]]@k

# change point model selection
# mean deviance ( -2 * log likelihood) + 2*(#parameters)
mean(ldats_paper_cp_results_rodent$saved_lls * -2) + 2*(3*(ntopics-1)*(1+1)+(1))
mean(ldats_paper_cp_results_rodent2$saved_lls * -2)+ 2*(3*(ntopics-1)*(2+1)+(2))
mean(ldats_paper_cp_results_rodent3$saved_lls * -2)+ 2*(3*(ntopics-1)*(3+1)+(3))
mean(ldats_paper_cp_results_rodent4$saved_lls * -2)+ 2*(3*(ntopics-1)*(4+1)+(4))
mean(ldats_paper_cp_results_rodent5$saved_lls * -2)+ 2*(3*(ntopics-1)*(5+1)+(5))
mean(ldats_paper_cp_results_rodent6$saved_lls * -2)+ 2*(3*(ntopics-1)*(6+1)+(6))

# The lowest deviance is for 4 changepoints.

ldats_paper_cpt_selected = ldats_paper_cp_results_rodent4

df_4 = as.data.frame(t(ldats_paper_cpt_selected$saved[,1,])) %>% reshape::melt()
year_continuous = 1970 + as.integer(julian(paper_dates)) / 365.25
df_4$value = year_continuous[df_4$value]


ldats_paper_cpt_dates <- vector(length = 4)
ldats_paper_cpt_dates[1] <- mean(df_4[df_4$variable=='V1','value']) %>% lubridate::date_decimal() %>% format('%Y-%m-%d')
ldats_paper_cpt_dates[2] <-mean(df_4[df_4$variable=='V2','value']) %>% lubridate::date_decimal() %>% format('%Y-%m-%d')
ldats_paper_cpt_dates[3] <-mean(df_4[df_4$variable=='V3','value']) %>% lubridate::date_decimal() %>% format('%Y-%m-%d')
ldats_paper_cpt_dates[4] <-mean(df_4[df_4$variable=='V4','value']) %>% lubridate::date_decimal() %>% format('%Y-%m-%d')

ldats_paper_cpt_dates


## ----load paper LDA for paper cpt, eval = F------------------------------
#  
#  #### Load LDA ####
#  load(here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'paper_ldamodel.Rds'))
#  

## ----run paper LDA and paper cpt, eval = F-------------------------------
#  #### Run changepoint ####
#  source(here::here('vignettes', 'christensenetal-comparison-files', 'christensenetal-2018-files', 'changepointmodel.r'))
#  
#  # set up parameters for model
#  year_continuous = 1970 + as.integer(julian(paper_dates)) / 365.25
#  x = data.frame(
#    year_continuous = year_continuous,
#    sin_year = sin(year_continuous * 2 * pi),
#    cos_year = cos(year_continuous * 2 * pi)
#  )
#  
#  # run models with 1, 2, 3, 4, 5, 6 changepoints
#  paper_paper_cp_results_rodent = changepoint_model(ldamodel, x, 1, weights = rep(1,length(year_continuous)))
#  save(paper_paper_cp_results_rodent, file = here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'paper_paper_1.Rds'))
#  rm(paper_paper_cp_results_rodent)
#  
#  paper_paper_cp_results_rodent2 = changepoint_model(ldamodel, x, 2, weights = rep(1,length(year_continuous)))
#  save(paper_paper_cp_results_rodent2, file = here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'paper_paper_2.Rds'))
#  rm(paper_paper_cp_results_rodent2)
#  
#  paper_paper_cp_results_rodent3 = changepoint_model(ldamodel, x, 3, weights = rep(1,length(year_continuous)))
#  save(paper_paper_cp_results_rodent3, file = here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'paper_paper_3.Rds'))
#  rm(paper_paper_cp_results_rodent3)
#  
#  paper_paper_cp_results_rodent4 = changepoint_model(ldamodel, x, 4, weights = rep(1,length(year_continuous)))
#  save(paper_paper_cp_results_rodent4, file = here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'paper_paper_4.Rds'))
#  rm(paper_paper_cp_results_rodent4)
#  
#  paper_paper_cp_results_rodent5 = changepoint_model(ldamodel, x, 5, weights = rep(1,length(year_continuous)))
#  save(paper_paper_cp_results_rodent5, file = here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'paper_paper_5.Rds'))
#  rm(paper_paper_cp_results_rodent5)
#  
#  paper_paper_cp_results_rodent6 = changepoint_model(ldamodel, x, 6, weights = rep(1,length(year_continuous)))
#  save(paper_paper_cp_results_rodent6, file = here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'paper_paper_6.Rds'))
#  rm(paper_paper_cp_results_rodent6)

## ----select paper LDA and paper cpt, eval = F----------------------------
#  #### Changepoint model selection ####
#  load(here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'paper_paper_1.Rds'))
#  load(here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'paper_paper_2.Rds'))
#  load(here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'paper_paper_3.Rds'))
#  load(here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'paper_paper_4.Rds'))
#  load(here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'paper_paper_5.Rds'))
#  load(here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'paper_paper_6.Rds'))

## ----select paper-paper, eval = T----------------------------------------
ntopics = ldamodel@k

# change point model selection
# mean deviance ( -2 * log likelihood) + 2*(#parameters)
mean(paper_paper_cp_results_rodent$saved_lls * -2) + 2*(3*(ntopics-1)*(1+1)+(1))
mean(paper_paper_cp_results_rodent2$saved_lls * -2)+ 2*(3*(ntopics-1)*(2+1)+(2))
mean(paper_paper_cp_results_rodent3$saved_lls * -2)+ 2*(3*(ntopics-1)*(3+1)+(3))
mean(paper_paper_cp_results_rodent4$saved_lls * -2)+ 2*(3*(ntopics-1)*(4+1)+(4))
mean(paper_paper_cp_results_rodent5$saved_lls * -2)+ 2*(3*(ntopics-1)*(5+1)+(5))
mean(paper_paper_cp_results_rodent6$saved_lls * -2)+ 2*(3*(ntopics-1)*(6+1)+(6))

# The lowest deviance is for 4 changepoints.

paper_paper_cpt_selected = paper_paper_cp_results_rodent4

df_4 = as.data.frame(t(paper_paper_cpt_selected$saved[,1,])) %>%  reshape::melt()
df_4$value = year_continuous[df_4$value]


paper_paper_cpt_dates <- vector(length = 4)
paper_paper_cpt_dates[1] <- mean(df_4[df_4$variable=='V1','value']) %>% lubridate::date_decimal() %>% format('%Y-%m-%d')
paper_paper_cpt_dates[2] <-mean(df_4[df_4$variable=='V2','value']) %>% lubridate::date_decimal() %>% format('%Y-%m-%d')
paper_paper_cpt_dates[3] <-mean(df_4[df_4$variable=='V3','value']) %>% lubridate::date_decimal() %>% format('%Y-%m-%d')
paper_paper_cpt_dates[4] <-mean(df_4[df_4$variable=='V4','value']) %>% lubridate::date_decimal() %>% format('%Y-%m-%d')

paper_paper_cpt_dates



## ----run LDATS LDA and LDATS changepoint, eval = F-----------------------
#  #### Run LDATS changepoint ####
#  
#  ldats_ldats_cpt <- TS_on_LDA(LDA_models = ldats_lda_selected,
#                               document_covariate_table = rodents$document_covariate_table,
#                               formulas = ~ sin_year + cos_year,
#                               nchangepoints = 1:6,
#                               weights = NULL,
#                               timename = 'newmoon')
#  
#  
#  save(ldats_ldats_cpt, file = here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'ldats_ldats.Rds'))
#  rm(ldats_ldats_cpt)
#  rm(ldats_lda_selected)

## ----select ldats lda + ldats cpt----------------------------------------

ldats_ldats_cpt_selected <- select_TS(ldats_ldats_cpt)

ldats_ldats_cpt_selected$nchangepoints


## ----get ldats lda + ldats cpt dates-------------------------------------
ldats_ldats_cpt_dates <- ldats_ldats_cpt_selected$rho_summary$Mean
ldats_ldats_cpt_dates <- round(ldats_ldats_cpt_dates)

ldats_ldats_cpt_dates <- as.data.frame(ldats_ldats_cpt_dates)
colnames(ldats_ldats_cpt_dates) <- 'newmoon'
ldats_ldats_cpt_dates <- dplyr::left_join(ldats_ldats_cpt_dates, rodents$document_covariate_table, by = 'newmoon')
ldats_ldats_cpt_dates <- ldats_ldats_cpt_dates$censusdate


## ----load paper LDA for LDATS cpt, eval = F------------------------------
#  
#  #### Load LDA ####
#  load(here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'paper_ldamodel.Rds'))
#  

## ----run paper LDA and LDATS changepoint, eval = F-----------------------
#  #### Run LDATS changepoint ####
#  
#  paper_ldats_cpt <- TS_on_LDA(LDA_models = ldamodel,
#                               document_covariate_table = rodents$document_covariate_table,
#                               formulas = ~ sin_year + cos_year,
#                               nchangepoints = 1:6,
#                               weights = NULL,
#                               timename = 'newmoon')
#  
#  
#  save(paper_ldats_cpt, file = here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'paper_ldats.Rds'))
#  rm(paper_ldats_cpt)
#  rm(ldamodel)

## ----select paper lda + ldats cpt----------------------------------------

paper_ldats_cpt_selected <- select_TS(paper_ldats_cpt)

paper_ldats_cpt_selected$nchangepoints


## ----get paper lda + ldats cpt dates-------------------------------------
paper_ldats_cpt_dates <- paper_ldats_cpt_selected$rho_summary$Mean
paper_ldats_cpt_dates <- round(paper_ldats_cpt_dates)

paper_ldats_cpt_dates <- as.data.frame(paper_ldats_cpt_dates)
colnames(paper_ldats_cpt_dates) <- 'newmoon'
paper_ldats_cpt_dates <- dplyr::left_join(paper_ldats_cpt_dates, rodents$document_covariate_table, by = 'newmoon')

halfinterval = as.integer((as.Date(rodents$document_covariate_table$censusdate[ which(rodents$document_covariate_table$newmoon == paper_ldats_cpt_dates$newmoon[1] + 1)]) - as.Date(rodents$document_covariate_table$censusdate[ which(rodents$document_covariate_table$newmoon == paper_ldats_cpt_dates$newmoon[1] - 1)]))/2)

newdate <- lubridate::ymd(rodents$document_covariate_table$censusdate[ which(rodents$document_covariate_table$newmoon == paper_ldats_cpt_dates$newmoon[1] -1 )]) + lubridate::days(halfinterval)

paper_ldats_cpt_dates$censusdate[1] <- as.character(newdate)



halfinterval =as.integer((as.Date(rodents$document_covariate_table$censusdate[ which(rodents$document_covariate_table$newmoon == paper_ldats_cpt_dates$newmoon[3] + 2)]) - as.Date(rodents$document_covariate_table$censusdate[ which(rodents$document_covariate_table$newmoon == paper_ldats_cpt_dates$newmoon[3] - 1)]))/3)

newdate <- lubridate::ymd(rodents$document_covariate_table$censusdate[ which(rodents$document_covariate_table$newmoon == paper_ldats_cpt_dates$newmoon[3] -2 )]) + lubridate::days(halfinterval)

paper_ldats_cpt_dates$censusdate[3] <- as.character(newdate)

paper_ldats_cpt_dates <- paper_ldats_cpt_dates$censusdate

## ----load cpts to plot, include = F--------------------------------------

year_continuous = 1970 + as.integer(julian(paper_dates)) / 365.25
x = data.frame(
  year_continuous = year_continuous,
  sin_year = sin(year_continuous * 2 * pi),
  cos_year = cos(year_continuous * 2 * pi)
)
ntopics = 4


## ----plot paper LDA and LDATS cpts, fig.width=7, fig.height=6------------

plot(paper_ldats_cpt_selected)


## ----plot ldats LDA and LDATS cpt, fig.width=7, fig.height=6-------------
plot(ldats_ldats_cpt_selected)

## ----plot paper LDA and paper cpt, fig.width=7, fig.height=6-------------

source(here::here('vignettes', 'christensenetal-comparison-files', 'christensenetal-2018-files', 'LDA_figure_scripts.R'))
source(here::here('vignettes', 'christensenetal-comparison-files', 'christensenetal-2018-files', 'changepointmodel.r'))
paper_cpts = find_changepoint_location(paper_paper_cpt_selected)
paper_cpt_plot = get_ll_non_memoized_plot(ldamodel,x,paper_cpts,make_plot=T,weights=rep(1,length(year_continuous)))
annual_hist(paper_paper_cpt_selected, year_continuous)
paper_cpt_plot


## ----plot LDATS lda and paper cpt, fig.width=7, fig.height=6-------------
ntopics = 6
ldats_cpts = find_changepoint_location(ldats_paper_cpt_selected)
ldats_cpt_plot = get_ll_non_memoized_plot(ldats_lda_selected[[1]],x,ldats_cpts,make_plot=T,weights=rep(1,length(year_continuous)))
annual_hist(ldats_paper_cpt_selected, year_continuous)
ldats_cpt_plot

## ----report cpt dates, include = F---------------------------------------

cpt_dates <- as.data.frame(paper_paper_cpt_dates)
cpt_dates <- cbind(cpt_dates, ldats_paper_cpt_dates)
cpt_dates <- cbind(cpt_dates, ldats_ldats_cpt_dates)
cpt_dates <- cbind(cpt_dates, paper_ldats_cpt_dates)

colnames(cpt_dates) <- c('paperLDA_papercpt', 'ldatsLDA_papercpt', 'ldatsLDA_ldatscpt', 'paperLDA_ldatscpt')


## ----print cpt dates-----------------------------------------------------
cpt_dates

