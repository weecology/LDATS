## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- eval=FALSE---------------------------------------------------------
## install.packages("devtools")
## devtools::install_github("weecology/LDATS")

## ---- eval = T-----------------------------------------------------------
library(LDATS)


## ----download Christensen------------------------------------------------

paper_filepath <- here::here('vignettes', 'christensenetal-comparison-files', 'christensenetal-2018-files')
test_filepath <- here::here('vignettes', 'christensenetal-ristensenetal-comparison-files', 'christensenetal-2018-files', 'rodent_LDA_analysis.r')

if(!file.exists(test_filepath)) {

files_to_download <- c('rodent_LDA_analysis.r', 'rodent_data_for_LDA.r', 'AIC_model_selection.R', 'changepointmodel.r', 'LDA-distance.R', 'Rodent_table_dat.csv', 'LDA_figure_scripts.R')

for(i in 1:length(files_to_download)) {
  download.file(url = paste0("https://raw.githubusercontent.com/emchristensen/Extreme-events-LDA/master/", files_to_download[i]),
                destfile = paste0(paper_filepath, '/', files_to_download[i]))
}

rm(files_to_download)
rm(i)

}

rm(test_filepath)
rm(paper_filepath)

## ----LDATS data----------------------------------------------------------

data(rodents)

head(rodents[[1]])


## ----Paper data----------------------------------------------------------
source(here::here('vignettes', 'christensenetal-comparison-files', 'christensenetal-2018-files', 'rodent_data_for_LDA.r'))

dat = create_rodent_table(period_first = 1,
                          period_last = 436,
                          selected_plots = c(2,4,8,11,12,14,17,22),
                          selected_species = c('BA','DM','DO','DS','NA','OL','OT','PB','PE','PF','PH','PI','PL','PM','PP','RF','RM','RO','SF','SH','SO'))

# dates to go with count data
moondat = read.csv(text=getURL("https://raw.githubusercontent.com/weecology/PortalData/master/Rodents/moon_dates.csv"),stringsAsFactors = F)
moondat$date = as.Date(moondat$censusdate)

period_dates = filter(moondat,period %in% rownames(dat)) %>% select(period,date)
dates = period_dates$date

paper_dat <- dat

paper_dates <- dates

rm(list = c('dat', 'dates', 'period_dates',
            'create_rodent_table'))


## ----rodent data comparison----------------------------------------------

compare <- rodents[[1]] == paper_dat

length(which(rowSums(compare) < ncol(compare)))

## ----adjust LDATS data after Christensen et al---------------------------

trap_table = read.csv('https://raw.githubusercontent.com/weecology/PortalData/master/Rodents/Portal_rodent_trapping.csv')
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
## 
## ldats_ldas <- LDATS::LDA_set(document_term_table = rodents$document_term_table, topics = c(2:6), nseeds = 100)
## ldats_lda_selected <- LDATS::select_LDA(LDA_models = ldats_ldas)
## 
## save(ldats_ldas, file = here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'ldats_ldas.Rds'))
## save(ldats_lda_selected, file = here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'ldats_lda_selected.Rds'))
## 

## ----rm LDATS LDAS, include = F, eval = F--------------------------------
## rm(list = c('ldats_ldas', 'ldats_lda_selected'))

## ----create dat for paper lda, include = F-------------------------------
dat = paper_dat

## ----paper LDAs, eval = F------------------------------------------------
## source(here::here('vignettes', 'christensenetal-comparison-files', 'christensenetal-2018-files', 'AIC_model_selection.R'))
## source(here::here('vignettes', 'christensenetal-comparison-files', 'christensenetal-2018-files', 'LDA-distance.R'))
## 
## # Fit a bunch of LDA models with different seeds
## # Only use even numbers for seeds because consecutive seeds give identical results
## seeds = 2*seq(200)
## 
## # repeat LDA model fit and AIC calculation with a bunch of different seeds to test robustness of the analysis
## best_ntopic = repeat_VEM(paper_dat,
##                          seeds,
##                          topic_min=2,
##                          topic_max=6)
## hist(best_ntopic$k,breaks=c(0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5),xlab='best # of topics', main='')
## 
## # 2b. how different is species composition of 4 community-types when LDA is run with different seeds?
## # ==================================================================
## # get the best 100 seeds where 4 topics was the best LDA model
## seeds_4topics = best_ntopic %>%
##   filter(k == 4) %>%
##   arrange(aic) %>%
##   head(100) %>%
##   pull(SEED)
## 
## # choose seed with highest log likelihood for all following analyses
## #    (also produces plot of community composition for 'best' run compared to 'worst')
## best_seed = calculate_LDA_distance(paper_dat,seeds_4topics)
## mean_dist = unlist(best_seed)[2]
## max_dist = unlist(best_seed)[3]
## 
## # ==================================================================
## # 3. run LDA model
## # ==================================================================
## ntopics = 4
## SEED = unlist(best_seed)[1]  # For the paper, use seed 206
## ldamodel = LDA(paper_dat,ntopics, control = list(seed = SEED),method='VEM')
## 
## save(ldamodel,file = here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'paper_ldamodel.Rds'))

## ----cleanup paper LDAS, include = F, eval = F---------------------------
## rm(list = c('dat', 'ldamodel', 'max_dist', 'mean_dist',
##             'ntopics', 'SEED', 'seeds', 'seeds_4topics', 'aic_model', 'aic_model_gibbs', 'calculate_LDA_distance', 'Hellinger', 'min_H', 'repeat_VEM', 'best_ntopic', 'best_seed'))

## ----reload LDAS, include = FALSE----------------------------------------

load(here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'ldats_lda_selected.Rds'))
load(here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'paper_ldamodel.Rds'))


## ----plot paper LDA, fig.width=6, fig.height=6---------------------------
# Paper
plot(ldamodel, cols = NULL, option = "D")


## ----plot LDATS LDA, fig.width=6, fig.height=6---------------------------
# LDATS
plot(ldats_lda_selected[[1]], cols = NULL, option = "D")

## ----remove LDAS after plotting, include = F-----------------------------
rm(ldamodel)
rm(ldats_lda_selected)

## ----load LDATS LDA for paper cpt----------------------------------------

#### Load LDA ####
load(here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'ldats_lda_selected.Rds'))


## ----run LDATS LDA and paper cpt, eval = F-------------------------------
## #### Run changepoint ####
## source(here::here('vignettes', 'christensenetal-comparison-files', 'christensenetal-2018-files', 'changepointmodel.r'))
## 
## # set up parameters for model
## year_continuous = 1970 + as.integer(julian(paper_dates)) / 365.25
## x = data.frame(
##   year_continuous = year_continuous,
##   sin_year = sin(year_continuous * 2 * pi),
##   cos_year = cos(year_continuous * 2 * pi)
## )
## 
## # run models with 1, 2, 3, 4, 5, 6 changepoints
## cp_results_rodent = changepoint_model(ldats_lda_selected[[1]], x, 1, weights = rep(1,length(year_continuous)))
## save(cp_results_rodent, file = here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'ldats_paper_1.Rds'))
## rm(cp_results_rodent)
## 
## cp_results_rodent2 = changepoint_model(ldats_lda_selected[[1]], x, 2, weights = rep(1,length(year_continuous)))
## save(cp_results_rodent2, file = here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'ldats_paper_2.Rds'))
## rm(cp_results_rodent2)
## 
## cp_results_rodent3 = changepoint_model(ldats_lda_selected[[1]], x, 3, weights = rep(1,length(year_continuous)))
## save(cp_results_rodent3, file = here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'ldats_paper_3.Rds'))
## rm(cp_results_rodent3)
## 
## cp_results_rodent4 = changepoint_model(ldats_lda_selected[[1]], x, 4, weights = rep(1,length(year_continuous)))
## save(cp_results_rodent4, file = here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'ldats_paper_4.Rds'))
## rm(cp_results_rodent4)
## 
## cp_results_rodent5 = changepoint_model(ldats_lda_selected[[1]], x, 5, weights = rep(1,length(year_continuous)))
## save(cp_results_rodent5, file = here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'ldats_paper_5.Rds'))
## rm(cp_results_rodent5)
## 
## cp_results_rodent6 = changepoint_model(ldats_lda_selected[[1]], x, 6, weights = rep(1,length(year_continuous)))
## save(cp_results_rodent6, file = here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'ldats_paper_6.Rds'))
## rm(cp_results_rodent6)

## ----select LDATS LDA and paper cpt--------------------------------------
#### Changepoint model selection ####
load(here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'ldats_paper_1.Rds'))
load(here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'ldats_paper_2.Rds'))
load(here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'ldats_paper_3.Rds'))
load(here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'ldats_paper_4.Rds'))
load(here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'ldats_paper_5.Rds'))
load(here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'ldats_paper_6.Rds'))

ntopics = ldats_lda_selected[[1]]@k

# change point model selection
# mean deviance ( -2 * log likelihood) + 2*(#parameters)
mean(cp_results_rodent$saved_lls * -2) + 2*(3*(ntopics-1)*(1+1)+(1))
mean(cp_results_rodent2$saved_lls * -2)+ 2*(3*(ntopics-1)*(2+1)+(2))
mean(cp_results_rodent3$saved_lls * -2)+ 2*(3*(ntopics-1)*(3+1)+(3))
mean(cp_results_rodent4$saved_lls * -2)+ 2*(3*(ntopics-1)*(4+1)+(4))
mean(cp_results_rodent5$saved_lls * -2)+ 2*(3*(ntopics-1)*(5+1)+(5))
mean(cp_results_rodent6$saved_lls * -2)+ 2*(3*(ntopics-1)*(6+1)+(6))

# The lowest deviance is for 4 changepoints.


df_4 = as.data.frame(t(cp_results_rodent4$saved[,1,])) %>% reshape::melt()
year_continuous = 1970 + as.integer(julian(paper_dates)) / 365.25
df_4$value = year_continuous[df_4$value]


ldats_paper_cpt_dates <- vector(length = 4)
ldats_paper_cpt_dates[1] <- mean(df_4[df_4$variable=='V1','value']) %>% lubridate::date_decimal() %>% format('%Y-%m-%d')
ldats_paper_cpt_dates[2] <-mean(df_4[df_4$variable=='V2','value']) %>% lubridate::date_decimal() %>% format('%Y-%m-%d')
ldats_paper_cpt_dates[3] <-mean(df_4[df_4$variable=='V3','value']) %>% lubridate::date_decimal() %>% format('%Y-%m-%d')
ldats_paper_cpt_dates[4] <-mean(df_4[df_4$variable=='V4','value']) %>% lubridate::date_decimal() %>% format('%Y-%m-%d')

ldats_paper_cpt_dates


## ----clean up LDATS LDA and paper changepoint, include = F---------------
rm(list = c('cp_results_rodent', 'cp_results_rodent2', 'cp_results_rodent3', 'cp_results_rodent4', 'cp_results_rodent5', 'cp_results_rodent6', 'ntopics', 'ldats_lda_selected', 'df_4'))

## ----load paper LDA for paper cpt----------------------------------------

#### Load LDA ####
load(here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'paper_ldamodel.Rds'))


## ----run paper LDA and paper cpt, eval = F-------------------------------
## #### Run changepoint ####
## source(here::here('vignettes', 'christensenetal-comparison-files', 'christensenetal-2018-files', 'changepointmodel.r'))
## 
## # set up parameters for model
## year_continuous = 1970 + as.integer(julian(paper_dates)) / 365.25
## x = data.frame(
##   year_continuous = year_continuous,
##   sin_year = sin(year_continuous * 2 * pi),
##   cos_year = cos(year_continuous * 2 * pi)
## )
## 
## # run models with 1, 2, 3, 4, 5, 6 changepoints
## cp_results_rodent = changepoint_model(ldamodel, x, 1, weights = rep(1,length(year_continuous)))
## save(cp_results_rodent, file = here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'paper_paper1.Rds'))
## rm(cp_results_rodent)
## 
## cp_results_rodent2 = changepoint_model(ldamodel, x, 2, weights = rep(1,length(year_continuous)))
## save(cp_results_rodent2, file = here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'paper_paper2.Rds'))
## rm(cp_results_rodent2)
## 
## cp_results_rodent3 = changepoint_model(ldamodel, x, 3, weights = rep(1,length(year_continuous)))
## save(cp_results_rodent3, file = here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'paper_paper3.Rds'))
## rm(cp_results_rodent3)
## 
## cp_results_rodent4 = changepoint_model(ldamodel, x, 4, weights = rep(1,length(year_continuous)))
## save(cp_results_rodent4, file = here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'paper_paper4.Rds'))
## rm(cp_results_rodent4)
## 
## cp_results_rodent5 = changepoint_model(ldamodel, x, 5, weights = rep(1,length(year_continuous)))
## save(cp_results_rodent5, file = here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'paper_paper5.Rds'))
## rm(cp_results_rodent5)
## 
## cp_results_rodent6 = changepoint_model(ldamodel, x, 6, weights = rep(1,length(year_continuous)))
## save(cp_results_rodent6, file = here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'paper_paper6.Rds'))
## rm(cp_results_rodent6)

## ----select paper LDA and paper cpt--------------------------------------
#### Changepoint model selection ####
load(here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'paper_paper_1.Rds'))
load(here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'paper_paper_2.Rds'))
load(here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'paper_paper_3.Rds'))
load(here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'paper_paper_4.Rds'))
load(here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'paper_paper_5.Rds'))
load(here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'paper_paper_6.Rds'))

ntopics = ldamodel@k

# change point model selection
# mean deviance ( -2 * log likelihood) + 2*(#parameters)
mean(cp_results_rodent$saved_lls * -2) + 2*(3*(ntopics-1)*(1+1)+(1))
mean(cp_results_rodent2$saved_lls * -2)+ 2*(3*(ntopics-1)*(2+1)+(2))
mean(cp_results_rodent3$saved_lls * -2)+ 2*(3*(ntopics-1)*(3+1)+(3))
mean(cp_results_rodent4$saved_lls * -2)+ 2*(3*(ntopics-1)*(4+1)+(4))
mean(cp_results_rodent5$saved_lls * -2)+ 2*(3*(ntopics-1)*(5+1)+(5))
mean(cp_results_rodent6$saved_lls * -2)+ 2*(3*(ntopics-1)*(6+1)+(6))

# The lowest deviance is for 4 changepoints.

df_4 = as.data.frame(t(cp_results_rodent4$saved[,1,])) %>%  reshape::melt()
df_4$value = year_continuous[df_4$value]


paper_paper_cpt_dates <- vector(length = 4)
paper_paper_cpt_dates[1] <- mean(df_4[df_4$variable=='V1','value']) %>% lubridate::date_decimal() %>% format('%Y-%m-%d')
paper_paper_cpt_dates[2] <-mean(df_4[df_4$variable=='V2','value']) %>% lubridate::date_decimal() %>% format('%Y-%m-%d')
paper_paper_cpt_dates[3] <-mean(df_4[df_4$variable=='V3','value']) %>% lubridate::date_decimal() %>% format('%Y-%m-%d')
paper_paper_cpt_dates[4] <-mean(df_4[df_4$variable=='V4','value']) %>% lubridate::date_decimal() %>% format('%Y-%m-%d')

paper_paper_cpt_dates


## ----clean up paper LDA and paper changepoint, include = F---------------
rm(list = c('cp_results_rodent', 'cp_results_rodent2', 'cp_results_rodent3', 'cp_results_rodent4', 'cp_results_rodent5', 'cp_results_rodent6', 'ntopics', 'ldamodel', 'df_4'))

## ----load LDATS LDA for LDATS cpt----------------------------------------

#### Load LDA ####
load(here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'ldats_lda_selected.Rds'))


## ----run LDATS LDA and LDATS changepoint, eval = F-----------------------
## #### Run LDATS changepoint ####
## 
## ldats_ldats_cpt <- TS_on_LDA(LDA_models = ldats_lda_selected,
##                              document_covariate_table = rodents$document_covariate_table,
##                              formulas = ~ sin_year + cos_year,
##                              nchangepoints = 1:6,
##                              weights = NULL)
## 
## 
## save(ldats_ldats_cpt, file = here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'ldats_ldats.Rds'))
## rm(ldats_ldats_cpt)
## rm(ldats_lda_selected)

## ----load LDATS LDA and LDATS changepoint, include = F-------------------
load(here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'ldats_ldats.Rds'))

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


## ----cleanup LDATS LDA and LDATS changepoint, include = F----------------

save(ldats_ldats_cpt_selected, file = here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'ldats_ldats_selected.Rds'))

rm(ldats_ldats_cpt_selected)
rm(ldats_ldats_cpt)

## ----load paper LDA for LDATS cpt----------------------------------------

#### Load LDA ####
load(here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'paper_ldamodel.Rds'))


## ----run paper LDA and LDATS changepoint, eval = F-----------------------
## #### Run LDATS changepoint ####
## 
## paper_ldats_cpt <- TS_on_LDA(LDA_models = ldamodel,
##                              document_covariate_table = rodents$document_covariate_table,
##                              formulas = ~ sin_year + cos_year,
##                              nchangepoints = 1:6,
##                              weights = NULL)
## 
## 
## save(paper_ldats_cpt, file = here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'paper_ldats.Rds'))
## rm(paper_ldats_cpt)
## rm(ldamodel)

## ----load paper LDA and LDATS changepoint, include = F-------------------
load(here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'paper_ldats.Rds'))

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

## ----cleanup paper LDA and LDATS changepoint, include = F----------------

save(paper_ldats_cpt_selected, file = here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'paper_ldats_selected.Rds'))

rm(paper_ldats_cpt_selected)
rm(paper_ldats_cpt)
rm(list = c('newdate', 'halfinterval'))

## ----load cpts to plot, include = F--------------------------------------
load(here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'paper_ldats_selected.Rds'))
load(here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'ldats_ldats_selected.Rds'))
load(here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'paper_paper_4.Rds'))
paper_paper_cpt_selected <- cp_results_rodent4
load(here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'ldats_paper_4.Rds'))
ldats_paper_cpt_selected <- cp_results_rodent4

year_continuous = 1970 + as.integer(julian(paper_dates)) / 365.25
x = data.frame(
  year_continuous = year_continuous,
  sin_year = sin(year_continuous * 2 * pi),
  cos_year = cos(year_continuous * 2 * pi)
)
ntopics = 4


## ----plot paper LDA and LDATS cpts, fig.width=6, fig.height=6------------

plot(paper_ldats_cpt_selected)


## ----plot ldats LDA and LDATS cpt, fig.width=6, fig.height=6-------------
plot(ldats_ldats_cpt_selected)

## ----plot paper LDA and paper cpt, fig.width=6, fig.height=6-------------

source(here::here('vignettes', 'christensenetal-comparison-files', 'christensenetal-2018-files', 'LDA_figure_scripts.R'))
source(here::here('vignettes', 'christensenetal-comparison-files', 'christensenetal-2018-files', 'changepointmodel.r'))
paper_cpts = find_changepoint_location(paper_paper_cpt_selected)
paper_cpt_plot = get_ll_non_memoized_plot(ldamodel,x,paper_cpts,make_plot=T,weights=rep(1,length(year_continuous)))
annual_hist(paper_paper_cpt_selected, year_continuous)
paper_cpt_plot


## ----plot LDATS lda and paper cpt, fig.width=6, fig.height=6-------------
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

## ----cleanup cpt_dates, include = F--------------------------------------
rm(list= 'ldats_ldats_cpt_dates', 'ldats_paper_cpt_dates',  'paper_ldats_cpt_dates', 'paper_paper_cpt_dates', 'cpt_dates')

