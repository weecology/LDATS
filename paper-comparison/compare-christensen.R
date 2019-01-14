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

## ----adjusted rodent data------------------------------------------------

ldats_dat_adj <- as.matrix(ch_dat)
ldats_dat_adj <- apply(ldats_dat, c(1,2), FUN = as.integer)
rodents_adj <- rodents
rodents_adj$document_term_table <- ldats_dat_adj
rm(ldats_dat_adj)
rm(ch_dat)

## ----adj vs non adj LDAS-------------------------------------------------

# paper LDA w adjusted data
source('christensen-ecology-files/AIC_model_selection.R')
source('christensen-ecology-files/LDA-distance.R')

# Fit a bunch of LDA models with different seeds
# Only use even numbers for seeds because consecutive seeds give identical results
seeds = 2*seq(200)

# repeat LDA model fit and AIC calculation with a bunch of different seeds to test robustness of the analysis
best_ntopic = repeat_VEM(rodents_adj[[1]],
                         seeds,
                         topic_min=2,
                         topic_max=6)
hist(best_ntopic$k,breaks=c(0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5),xlab='best # of topics', main='')

# 2b. how different is species composition of 4 community-types when LDA is run with different seeds?
# ==================================================================
# get the best 100 seeds where 4 topics was the best LDA model
seeds_4topics = best_ntopic %>%
  filter(k == 4) %>%
  arrange(aic) %>%
  head(100) %>%
  pull(SEED)

# choose seed with highest log likelihood for all following analyses
#    (also produces plot of community composition for 'best' run compared to 'worst')
best_seed = calculate_LDA_distance(rodents[[1]],seeds_4topics)
mean_dist = unlist(best_seed)[2]
max_dist = unlist(best_seed)[3]

# ==================================================================
# 3. run LDA model
# ==================================================================
ntopics = 4
SEED = unlist(best_seed)[1]  # For the paper, I use seed 206
ldamodel_adj = LDA(rodents[[1]],ntopics, control = list(seed = SEED),method='VEM')

plot(ldamodel_adj)

# paper LDA w nonadjusted data

# Fit a bunch of LDA models with different seeds
# Only use even numbers for seeds because consecutive seeds give identical results
seeds = 2*seq(200)

# repeat LDA model fit and AIC calculation with a bunch of different seeds to test robustness of the analysis
best_ntopic = repeat_VEM(rodents[[1]],
                         seeds,
                         topic_min=2,
                         topic_max=6)
hist(best_ntopic$k,breaks=c(0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5),xlab='best # of topics', main='')

# 2b. how different is species composition of 4 community-types when LDA is run with different seeds?
# ==================================================================
# get the best 100 seeds where 4 topics was the best LDA model
seeds_4topics = best_ntopic %>%
  filter(k == 4) %>%
  arrange(aic) %>%
  head(100) %>%
  pull(SEED)

# choose seed with highest log likelihood for all following analyses
#    (also produces plot of community composition for 'best' run compared to 'worst')
best_seed = calculate_LDA_distance(rodents[[1]],seeds_4topics)
mean_dist = unlist(best_seed)[2]
max_dist = unlist(best_seed)[3]

# ==================================================================
# 3. run LDA model
# ==================================================================
ntopics = 4
SEED = unlist(best_seed)[1]  # For the paper, I use seed 206
ldamodel_nonadj = LDA(rodents[[1]],ntopics, control = list(seed = SEED),method='VEM')

plot(ldamodel_nonadj)

# ldats LDA w adjusted data

ldats_adj_ldas <- LDATS::LDA_set(document_term_table = rodents_adj$document_term_table, topics = c(2:6), nseeds = 100)
ldats_adj_lda_selected <- LDATS::select_LDA(LDA_models = ldats_adj_ldas)

plot(ldats_adj_lda_selected)

# ldats LDA w nonadjusted data

ldats_nonadj_ldas <- LDATS::LDA_set(document_term_table = rodents$document_term_table, topics = c(2:6), nseeds = 100)
ldats_nonadj_lda_selected <- LDATS::select_LDA(LDA_models = ldats_ldas)

plot(ldats_nonadj_lda_selected)


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
## save(ldamodel, file = '~/Dropbox/ldats-models/paper_lda.Rds')
## save(ldats_lda_selected, file = '~/Dropbox/ldats-models/ldats_lda.Rds')

## ----reload LDAS, include = FALSE----------------------------------------

load('~/Dropbox/ldats-models/ldats_lda.Rds')
load('~/Dropbox/ldats-models/paper_lda.Rds')


## ----plot LDAs-----------------------------------------------------------
# Paper
plot(ldamodel, cols = NULL, option = "D")

# LDATS
plot(ldats_lda_selected[[1]], cols = NULL, option = "D")

## ----read model results--------------------------------------------------
ldats_ldats <- read.csv('ldats_ldats_dates.csv')
ldats_paper <- read.csv('ldats_paper_dates.csv')
paper_ldats <- read.csv('paper_ldats_dates.csv')
paper_paper <- read.csv('paper_paper_dates.csv')


estimates <- c(1:4)
estimates <- as.data.frame(estimates)

# Naming convention: ldamethod_cptmethod
# so ldats_paper is the ldats lda and paper cpt, etc
estimates$ldats_ldats <- ldats_ldats$newmoondate
estimates$ldats_paper <- ldats_paper$cpt_dates
estimates$paper_ldats <- paper_ldats$newmoondate
estimates$paper_paper <- paper_paper$cpt_dates

estimates
