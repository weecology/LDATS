# LDA and changepoint analysis pipeline on rodent data -- uses VEM method
#
#  1. prepare data
#      - specific script for each data set
#  2a. run LDA with different number of topics and use AIC to select best model
#  2b. run LDA with best number of topics (determined by 2a)
#  3. run changepoint model
#  4. produce figures


library(topicmodels)
library(RCurl)


source('rodent_data_for_LDA.R')
source('AIC_model_selection.R')
source('LDA_figure_scripts.R')


# ===================================================================
# prepare rodent data
dat = create_rodent_table(period_first = 1,
                          period_last = 436,
                          selected_plots = c(2,4,8,11,12,14,17,22),
                          selected_species = c('BA','DM','DO','DS','NA','OL','OT','PB','PE','PF','PH','PI','PL','PM','PP','RF','RM','RO','SF','SH','SO'))

# dates to go with count data
moondat = read.csv(text=getURL("https://raw.githubusercontent.com/weecology/PortalData/master/Rodents/moon_dates.csv"),stringsAsFactors = F)
moondat$date = as.Date(moondat$CensusDate)

period_dates = filter(moondat,Period %in% rownames(dat)) %>% select(Period,date)
dates = period_dates$date



# ===================================
# model parameters:
SEED = 2010
topic_min = 2
topic_max = 9

# ==================================================================
# select number of topics

# choose number of topics -- model selection using AIC
aic_values = aic_model(dat,SEED=2010,topic_min,topic_max)
ntopics = filter(aic_values,aic==min(aic)) %>% select(k) %>% as.numeric()

# repeat aic calculation with a bunch of different seeds to test robustness of the analysis
best_ntopic = repeat_VEM(dat,
                         ntimes=20,
                         topic_min=2,
                         topic_max=9)

# ==================================================================
# run LDA model

ldamodel = LDA(dat,ntopics, control = list(seed = SEED),method='VEM')

# ==================================================================
# change point model -- Not working today.  Not sure why.

# year_continuous = 1970 + as.integer(julian(dates)) / 365.25
# x = data.frame(
#   year_continuous = year_continuous,
#   sin_year = sin(year_continuous * 2 * pi),
#   cos_year = cos(year_continuous * 2 * pi)
# )
# cp_results = changepoint_model(ldamodel, x, n_chpoints, maxit)

# =================================================================
# figures

beta1=matrix(apply(ldamodel$beta,2,mean),ntopics,nspp)
plot_community_composition(beta1)

plot_component_communities_gibbs_credible(ldamodel,ntopics,dates)
