library(LDATS)

setwd('paper-comparison')

data(rodents)

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


ldats_dat_adj <- as.matrix(ch_dat)
ldats_dat_adj <- apply(ldats_dat_adj, c(1,2), FUN = as.integer)
rodents_adj <- rodents
rodents_adj$document_term_table <- ldats_dat_adj
rm(ldats_dat_adj)
rm(ch_dat)


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
best_seed = calculate_LDA_distance(rodents_adj[[1]],seeds_4topics)
mean_dist = unlist(best_seed)[2]
max_dist = unlist(best_seed)[3]

# ==================================================================
# 3. run LDA model
# ==================================================================
ntopics = 4
SEED = unlist(best_seed)[1]  # For the paper, I use seed 206
ldamodel_adj = LDA(rodents_adj[[1]],ntopics, control = list(seed = SEED),method='VEM')

save(ldamodel_adj, file = '~/Dropbox/ldats-models/paper_lda_adj_data.Rds')

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

save(ldamodel_nonadj, file = '~/Dropbox/ldats-models/paper_lda_nonadj_data.Rds')
plot(ldamodel_nonadj)
ldamodel_nonadj@k
head(ldamodel_nonadj@gamma)



# ldats LDA w adjusted data

ldats_adj_ldas <- LDATS::LDA_set(document_term_table = rodents_adj$document_term_table, topics = c(2:6), nseeds = 100)
ldats_adj_lda_selected <- LDATS::select_LDA(LDA_models = ldats_adj_ldas)

save(ldats_adj_lda_selected, file = '~/Dropbox/ldats-models/ldats_lda_adj_data.Rds')

plot(ldats_adj_lda_selected)

# ldats LDA w nonadjusted data

ldats_nonadj_ldas <- LDATS::LDA_set(document_term_table = rodents$document_term_table, topics = c(2:6), nseeds = 100)
ldats_nonadj_lda_selected <- LDATS::select_LDA(LDA_models = ldats_nonadj_ldas)

save(ldats_nonadj_lda_selected, file = '~/Dropbox/ldats-models/ldats_lda_nonadj_data.Rds')


plot(ldats_nonadj_lda_selected)
