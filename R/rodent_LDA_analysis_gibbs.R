# LDA and changepoint analysis pipeline on rodent data -- uses Denis Valle's Gibbs method
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
dat = create_rodent_table(1,436,c(2,4,8,11,12,14,17,22),
                          selected_species = c('BA','DM','DO','DS','NA','OL','OT','PB','PE','PF','PH','PI','PL','PM','PP','RF','RM','RO','SF','SH','SO'))

# dates to go with count data
moondat = read.csv(text=getURL("https://raw.githubusercontent.com/weecology/PortalData/master/Rodents/moon_dates.csv"),stringsAsFactors = F)
moondat$date = as.Date(moondat$CensusDate)

period_dates = filter(moondat,Period %in% rownames(dat)) %>% select(Period,date)
dates = period_dates$date



# ===================================
# model parameters -- could be inputs for future wrapper function:
ngibbs=500
test_topics = c(2,8) # only comparing models with 2 or 3 topics for now -- faster
n_chpoints = 2
maxit = 1000
nspp = 21
nplots = 436

# ==================================================================
# select number of topics

# choose number of topics -- model selection using AIC
aic_values = aic_model_gibbs(dat,ngibbs,test_topics[1],test_topics[2],F)

# ==================================================================
# run LDA model
ntopics = filter(aic_values,aic==min(aic)) %>% select(k) %>% as.numeric()
ldamodel = gibbs.samp(dat.agg=dat,ngibbs=ngibbs,ncommun=ntopics,a.betas=1,a.theta=1)

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
