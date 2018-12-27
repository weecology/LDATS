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

## ----download Christensen, eval = FALSE----------------------------------
## 
## files_to_download <- c('rodent_LDA_analysis.r', 'rodent_data_for_LDA.r', 'AIC_model_selection.R', 'changepointmodel.r', 'LDA-distance.R', 'Rodent_table_dat.csv', 'LDA_figure_scripts.R')
## 
## for(i in 1:length(files_to_download)) {
##   download.file(url = paste0("https://raw.githubusercontent.com/emchristensen/Extreme-events-LDA/master/", files_to_download[i]),
##                 destfile = paste0('christensen-ecology-files/', files_to_download[i]))
## }
## 
## rm(files_to_download)
## rm(i)
## 

## ----LDATS data----------------------------------------------------------

data(rodents)

ldats_dat <- rodents[[1]]
ldats_dates <- rodents[[2]]
head(ldats_dat)


## ----Ch data-------------------------------------------------------------
source('christensen-ecology-files/rodent_data_for_LDA.r')
dat = create_rodent_table(period_first = 1,
                          period_last = 436,
                          selected_plots = c(2,4,8,11,12,14,17,22),
                          selected_species = c('BA','DM','DO','DS','NA','OL','OT','PB','PE','PF','PH','PI','PL','PM','PP','RF','RM','RO','SF','SH','SO'), diagnose = F)

# dates to go with count data
moondat = read.csv(text=getURL("https://raw.githubusercontent.com/weecology/PortalData/master/Rodents/moon_dates.csv"),stringsAsFactors = F)
moondat$date = as.Date(moondat$censusdate)

period_dates = filter(moondat,period %in% rownames(dat)) %>% select(period,date)
dates = period_dates$date

ch_dat <- dat

ch_dates <- dates

rm(list = c('dat', 'dates', 'period_dates', 'moondat',
            'create_rodent_table'))


## ----rodent data comparison, eval = F------------------------------------
## 
## compare <- ldats_dat == ch_dat
## 
## compare_rows <- vector(length = nrow(compare))
## 
## for(i in 1:nrow(compare)) {
##   if(sum(compare[i, ]) == 21) {
##     compare_rows[i] <- TRUE
##   } else {
##     compare_rows[i] <- FALSE
##   }
## }
## 
## length(which(compare_rows == F))
## ldats_dat[ which(compare_rows == F), ]
## ch_dat[ which(compare_rows == F), ]

## ----Ch rodent adjustment, eval = F--------------------------------------
## 
## source('christensen-ecology-files/rodent_data_for_LDA.r')
## diagnose_dat = create_rodent_table(period_first = 1,
##                           period_last = 436,
##                           selected_plots = c(2,4,8,11,12,14,17,22),
##                           selected_species = c('BA','DM','DO','DS','NA','OL','OT','PB','PE','PF','PH','PI','PL','PM','PP','RF','RM','RO','SF','SH','SO'), diagnose = T)
## 
## adjusted_table <- diagnose_dat[[1]]
## raw_table <- diagnose_dat[[2]]
## nplots <- diagnose_dat[[3]]
## nplots <- nplots[ which(nplots$period <= 436), ]
## 
## compare_raw <- raw_table == ldats_dat
## length(which(compare_raw == F))
## 
## which(compare_rows == F)
## which(nplots$x < 8)
## 

## ----data adjustment cleanup, eval = F, include = FALSE------------------
## rm(list = c('diagnose_dat', 'nplots', 'compare', 'compare_raw', 'adjusted_table', 'i', 'compare_rows', 'raw_table', 'create_rodent_table'))

## ----switch to adjusted rodent data--------------------------------------

ldats_dat <- as.matrix(ch_dat)
ldats_dat <- apply(ldats_dat, c(1,2), FUN = as.integer)
rodents$document_term_table <- ldats_dat
rm(ldats_dat)
rm(ldats_dates)
rm(ch_dat)

## ----LDATS LDAs, eval = F------------------------------------------------
## 
## ldats_ldas <- LDATS::LDA_set(document_term_table = rodents$document_term_table, topics = c(2:6), nseeds = 100)
## ldats_lda_selected <- LDATS::select_LDA(LDA_models = ldats_ldas)
## 

## ----paper LDAs, eval = F------------------------------------------------
## source('christensen-ecology-files/AIC_model_selection.R')
## source('christensen-ecology-files/LDA-distance.R')
## 
## # Fit a bunch of LDA models with different seeds
## # Only use even numbers for seeds because consecutive seeds give identical results
## seeds = 2*seq(200)
## 
## # repeat LDA model fit and AIC calculation with a bunch of different seeds to test robustness of the analysis
## best_ntopic = repeat_VEM(rodents[[1]],
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
## best_seed = calculate_LDA_distance(rodents[[1]],seeds_4topics)
## mean_dist = unlist(best_seed)[2]
## max_dist = unlist(best_seed)[3]
## 
## # ==================================================================
## # 3. run LDA model
## # ==================================================================
## ntopics = 4
## SEED = unlist(best_seed)[1]  # For the paper, I use seed 206
## ldamodel = LDA(rodents[[1]],ntopics, control = list(seed = SEED),method='VEM')
## 
## 

## ----save LDAs, eval = F, include = F------------------------------------
## save(ldamodel, file = 'model-stash/paper_lda.Rds')
## save(ldats_lda_selected, file = 'model-stash/ldats_lda.Rds')

## ----reload LDAS, include = FALSE----------------------------------------

load('model-stash/ldats_lda.Rds')
load('paper-comparison/model-stash/paper_lda.Rds')


## ----plot LDAs-----------------------------------------------------------
# Paper
plot(ldamodel, cols = NULL, option = "D")

# LDATS
plot(ldats_lda_selected[[1]], cols = NULL, option = "D")

## ---- include = F--------------------------------------------------------
# remove paper LDA (for memory)
rm(ldamodel)

## ----LDATS LDA + LDATS cpt, eval = F-------------------------------------
## # Run changepoint
## 
## ldats_ldats <- LDATS::TS_on_LDA(LDA_models = ldats_lda_selected,
##                                 document_covariate_table = rodents[[2]],
##                                 formulas =c(~ sin_year + cos_year),
##                                 nchangepoints = c(2, 3),
##                                 weights = NULL)
## 
## save(ldats_ldats, file = 'model-stash/ldats_ldats.Rds')
## 
## # Select best changepoint
## 
## ldats_ldats_selected <- LDATS::select_TS(TS_models = ldats_ldats)
## 
## save(ldats_ldats_selected, file = '/Users/renatadiaz/Documents/model-stash/ldats_ldats_selected.Rds')
## 

## ---- eval = F-----------------------------------------------------------
## rm(ldats_ldats)
## rm(ldats_ldats_selected)

## ----ldats paper, eval = F-----------------------------------------------
source('christensen-ecology-files/changepointmodel.r')
# set up parameters for model
year_continuous = 1970 + as.integer(julian(ch_dates)) / 365.25
x = data.frame(
  year_continuous = year_continuous,
  sin_year = sin(year_continuous * 2 * pi),
  cos_year = cos(year_continuous * 2 * pi)
)

# run models with 1, 2, 3, 4, 5 changepoints
cp_results_rodent = changepoint_model(ldats_lda_selected[[1]], x, 1, weights = rep(1,length(year_continuous)))

save(cp_results_rodent, file = '/Users/renatadiaz/Documents/model-stash/paper_cpt_ldats_lda1.Rds')
rm(cp_results_rodent)


cp_results_rodent2 = changepoint_model(ldats_lda_selected[[1]], x, 2, weights = rep(1,length(year_continuous)))


save(cp_results_rodent2, file = '/Users/renatadiaz/Documents/model-stash/paper_cpt_ldats_lda2.Rds')
rm(cp_results_rodent2)


cp_results_rodent3 = changepoint_model(ldats_lda_selected[[1]], x, 3, weights = rep(1,length(year_continuous)))
save(cp_results_rodent3, file = '/Users/renatadiaz/Documents/model-stash/paper_cpt_ldats_lda3.Rds')
rm(cp_results_rodent3)

cp_results_rodent4 = changepoint_model(ldats_lda_selected[[1]], x, 4, weights = rep(1,length(year_continuous)))
save(cp_results_rodent4, file = '/Users/renatadiaz/Documents/model-stash/paper_cpt_ldats_lda4.Rds')
rm(cp_results_rodent4)

cp_results_rodent5 = changepoint_model(ldats_lda_selected[[1]], x, 5, weights = rep(1,length(year_continuous)))
save(cp_results_rodent5, file = '/Users/renatadiaz/Documents/model-stash/paper_cpt_ldats_lda5.Rds')
rm(cp_results_rodent5)

cp_results_rodent6 = changepoint_model(ldats_lda_selected[[1]], x, 6, weights = rep(1,length(year_continuous)))
save(cp_results_rodent6, file = '/Users/renatadiaz/Documents/model-stash/paper_cpt_ldats_lda6.Rds')
rm(cp_results_rodent6)

rm(ldats_lda_selected)

## ----paper paper, eval = F-----------------------------------------------
## source('christensen-ecology-files/changepointmodel.r')
## # set up parameters for model
## year_continuous = 1970 + as.integer(julian(ch_dates)) / 365.25
## x = data.frame(
## year_continuous = year_continuous,
## sin_year = sin(year_continuous * 2 * pi),
## cos_year = cos(year_continuous * 2 * pi)
## )
## # run models with 1, 2, 3, 4, 5 changepoints
## cp_results_rodent = changepoint_model(ldamodel, x, 1, weights = rep(1,length(year_continuous)))
## save(cp_results_rodent, file = '/Users/renatadiaz/Documents/model-stash/paper_cpt_paper_lda1.Rds')
## rm(cp_results_rodent)
## 
## 
## cp_results_rodent2 = changepoint_model(ldamodel, x, 2, weights = rep(1,length(year_continuous)))
## save(cp_results_rodent2, file = '/Users/renatadiaz/Documents/model-stash/paper_cpt_paper_lda2.Rds')
## rm(cp_results_rodent2)
## 
## cp_results_rodent3 = changepoint_model(ldamodel, x, 3, weights = rep(1,length(year_continuous)))
## save(cp_results_rodent3, file = '/Users/renatadiaz/Documents/model-stash/paper_cpt_paper_lda3.Rds')
## rm(cp_results_rodent3)
## 
## cp_results_rodent4 = changepoint_model(ldamodel, x, 4, weights = rep(1,length(year_continuous)))
## save(cp_results_rodent4, file = '/Users/renatadiaz/Documents/model-stash/paper_cpt_paper_lda4.Rds')
## rm(cp_results_rodent4)
## 
## cp_results_rodent5 = changepoint_model(ldamodel, x, 5, weights = rep(1,length(year_continuous)))
## save(cp_results_rodent5, file = '/Users/renatadiaz/Documents/model-stash/paper_cpt_paper_lda5.Rds')
## rm(cp_results_rodent5)
## 
## cp_results_rodent6 = changepoint_model(ldamodel, x, 6, weights = rep(1,length(year_continuous)))
## save(cp_results_rodent6, file = '/Users/renatadiaz/Documents/model-stash/paper_cpt_paper_lda6.Rds')
## rm(cp_results_rodent6)

## ----paper paper selection-----------------------------------------------
ntopics = ldamodel@k
# Load models and compare deviance to select the best one
load('/Users/renatadiaz/Documents/model-stash/paper_cpt_paper_lda1.Rds')
load('/Users/renatadiaz/Documents/model-stash/paper_cpt_paper_lda2.Rds')
load('/Users/renatadiaz/Documents/model-stash/paper_cpt_paper_lda3.Rds')
load('/Users/renatadiaz/Documents/model-stash/paper_cpt_paper_lda4.Rds')
load('/Users/renatadiaz/Documents/model-stash/paper_cpt_paper_lda5.Rds')
load('/Users/renatadiaz/Documents/model-stash/paper_cpt_paper_lda6.Rds')

# change point model selection
# mean deviance ( -2 * log likelihood) + 2*(#parameters)
mean(cp_results_rodent$saved_lls * -2) + 2*(3*(ntopics-1)*(1+1)+(1))
mean(cp_results_rodent2$saved_lls * -2)+ 2*(3*(ntopics-1)*(2+1)+(2))
mean(cp_results_rodent3$saved_lls * -2)+ 2*(3*(ntopics-1)*(3+1)+(3))
mean(cp_results_rodent4$saved_lls * -2)+ 2*(3*(ntopics-1)*(4+1)+(4))
mean(cp_results_rodent5$saved_lls * -2)+ 2*(3*(ntopics-1)*(5+1)+(5))
mean(cp_results_rodent6$saved_lls * -2)+ 2*(3*(ntopics-1)*(5+1)+(5))


## ----paper ldats, eval = F-----------------------------------------------
## # Run changepoint
## 
## paper_ldats <- LDATS::TS_on_LDA(LDA_models = ldamodel,
##                                 document_covariate_table = rodents[[2]],
##                                 formulas =c(~ sin_year + cos_year),
##                                 nchangepoints = c(2, 3, 4, 5, 6),
##                                 weights = NULL)
## 
## save(paper_ldats, file = 'model-stash/paper_ldats.Rds')
## 
## # Select best changepoint
## 
## paper_ldats_selected <- LDATS::select_TS(TS_models = paper_ldats)
## 
## save(paper_ldats_selected, file = '/Users/renatadiaz/Documents/model-stash/paper_ldats_selected.Rds')
## 

