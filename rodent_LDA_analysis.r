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
SEED = 113052032
topic_min = 2
topic_max = 9
nspp=length(dat)
# ==================================================================
# select number of topics

# choose number of topics -- model selection using AIC
aic_values = aic_model(dat,SEED,topic_min,topic_max)
ntopics = filter(aic_values,aic==min(aic)) %>% select(k) %>% as.numeric()

# repeat aic calculation with a bunch of different seeds to test robustness of the analysis
best_ntopic = repeat_VEM(dat,
                         ntimes=200,
                         topic_min=2,
                         topic_max=9)
#hist(best_ntopic[,1],breaks=c(0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5),xlab='best # of topics')

# ==================================================================
# run LDA model

ldamodel = LDA_analysis_VEM(dat,SEED,c(topic_min,topic_max))

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

beta1 = community_composition(ldamodel)
# put columns in order of largest species to smallest
composition = beta1[,c('NA','DS','SH','SF','SO','DO','DM','PB','PH','OL','OT','PL','PM','PE','PP','PI','RF','RM','RO','BA','PF')]
plot_community_composition(composition)

plot_component_communities(ldamodel,ntopics,dates)
